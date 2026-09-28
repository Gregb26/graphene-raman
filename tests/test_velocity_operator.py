"""
Tests of the M0 building blocks of the velocity operator (optics/velocity_operator.py).

Everything is checked on the analytic graphene tight-binding model, whose answers are known:
    - the model itself: Hermiticity in R space, bond lengths, invariance of the bonds under the
      relabelling shift_B, number of R vectors;
    - the reciprocal lattice;
    - the Fourier transform: closed-form H_AB(k), Dirac point, Hermiticity of H(k), dH(k) and A(k),
      and dH(k) against centred finite differences;
    - dagger / hermitize.

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture (see `tb`): every property that does not depend on how the B orbital is labelled
must hold in both.

Run with:  .venv/bin/python -m pytest tests/test_velocity_operator.py -v
"""

import pytest
import numpy as np
from electron_defect_interaction.optics.velocity_operator import *
from types import SimpleNamespace

N = 100    # the k grid has N x N = 1e4 points: fast, yet covers the whole Brillouin zone
DK = 1e-5  # finite-difference step (1/Angstrom); truncation ~DK^2 and round-off ~eps/DK balance here

@pytest.fixture
def tb(request):
    """
    Graphene tight-binding model (WannierTB).

    Built with shift_B = (0,0,0) by default. A test can request other gauges with an indirect
    parametrization,

        @pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)

    pytest then passes each value to THIS fixture as request.param (instead of passing the raw
    tuple to the test), and runs the test once per value (ids tb0, tb1). For tests without that
    decorator request.param does not exist, hence getattr with a default.
    """
    shift = getattr(request, "param", (0,0,0))
    return make_graphene_tb(shift_B=shift)

@pytest.fixture
def grid(tb):
    """
    Quantities derived from `tb` that the tests reuse, bundled in a SimpleNamespace.

    It depends on the `tb` fixture, so pytest rebuilds it for every gauge of a parametrized test
    (the shifted model has 7 R vectors instead of 5).

    Fields:
        B      : (3, 3) reciprocal vectors in columns, 1/Angstrom
        k_red  : (N^2, 3) reduced k points
        k_cart : (N^2, 3) Cartesian k points, 1/Angstrom
        K      : (3,) Dirac point, 1/Angstrom
        R_cart : (nR, 3) Cartesian R vectors in rows, Angstrom
        index  : dict, integer triplet tuple(R) -> row of R in tb.R_int
        minus  : (nR,) int, minus[iR] = row of -R_int[iR], so that X_R[minus] lists X(-R)
    """
    B = reciprocal(tb.lattice)
    k_red, k_cart, K = k_grid(B,N)
    R_cart = tb.R_int @ tb.lattice.T # rows of integers -> rows of cartesian vectors
    index = {tuple(R): iR for iR, R in enumerate(tb.R_int)}
    minus = np.array([index[tuple(-R)] for R in tb.R_int]) # KeyError if some -R is missing

    return SimpleNamespace(B=B, k_red=k_red, k_cart=k_cart, K=K, R_cart=R_cart, index=index, minus=minus)


def make_graphene_tb_analytic(tb, grid):
    """
    Closed-form off-diagonal element of the graphene Bloch Hamiltonian, reference for `fourier`.

        H_AB(k) = -t (1 + e^{-ik.a1} + e^{-ik.a2}) = -t sum_{R in S} e^{ik.R},
        S = {0, -a1, -a2}  (cells of the three B neighbours of A when shift_B = 0).

    Valid only for shift_B = (0,0,0): with a shift L, H_AB(k) picks up the gauge phase e^{-ik.L}
    and only |H_AB| would still match.

    Inputs:
        tb : WannierTB
            Provides t (eV) and the lattice (a_i in columns, Angstrom).
        grid : SimpleNamespace from the `grid` fixture
            Provides k_cart, (Nk, 3), 1/Angstrom.
    Returns:
        H_AB : (Nk,) complex, eV
    """

    t = tb.t
    A = tb.lattice
    k = grid.k_cart

    return -t*(1+np.exp(-1j*k @ A[:,0].T) + np.exp(-1j*k @ A[:,1].T))

def test_fourier_analytic(tb, grid):
    """
    H_AB(k) from `fourier` equals the closed form on the whole grid (shift_B = 0 only, see
    make_graphene_tb_analytic). Checks at once the sign of the phase, the ndegen division, the
    conversion R_cart = R_int @ lattice.T and the orbital ordering (A = 0, B = 1).
    """

    H_tb = make_graphene_tb_analytic(tb, grid)
    H_k = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.k_cart)

    assert np.allclose(H_tb, H_k[:,0,1], atol=1e-12), 'tight binding model not correct'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_dirac_point(tb, grid):
    """
    The whole H(K) vanishes at the Dirac point K = (2 b1 + b2)/3: the three neighbour phases cancel,
    1 + e^{-i4pi/3} + e^{-i2pi/3} = 0, and the on-site energies are zero. Holds in both gauges,
    since the shift only multiplies H_AB by a phase. K is passed as a (1, 3) array (K[None, :])
    because `fourier` expects a list of k points.
    """

    H_K = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.K[None, :])

    assert np.allclose(H_K, 0, atol=1e-12), 'Hamiltonian not zero at Dirac point'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_hermitian(tb, grid):
    """
    Hermiticity in real space, H(-R) = H(R)^dagger, i.e. H_ij(R) = H_ji(-R)^* for every R.

    This is what makes H(k) Hermitian at every k. H_R[minus] lists H(-R) and is compared with the
    conjugate transpose of H_R. The dagger goes on ONE side only: on both sides the test would check
    H(-R) = H(R), an inversion symmetry that graphene's H(R) does not have (H_AB(-a1) = -t but
    H_AB(+a1) = 0).

    Also checked:
        - minus is an involution (minus[minus] = identity): -(-R) = R, and every R has its -R;
        - the same relation for every Cartesian component of r(R) (swapaxes(-1, -2) only acts on
          the two orbital axes, so the Cartesian axis 1 is untouched);
        - ndegen(-R) = ndegen(R): otherwise H(k) would not be Hermitian even with Hermitian H(R).
    """

    assert np.allclose(grid.minus[grid.minus], np.arange(len(tb.R_int))), 'R list not closed under R -> -R'
    assert np.allclose(tb.H_R, tb.H_R[grid.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'TB Hamiltonian not hermitian'
    assert np.allclose(tb.r_R, tb.r_R[grid.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'positian operator not hermitian'
    assert np.allclose(tb.ndegen, tb.ndegen[grid.minus], atol=1e-12), 'R and -R dont have the same weight'

def test_reciprocal(tb, grid):
    """
    a_i . b_j = 2 pi delta_ij, checked in the equivalent completeness form
    sum_i a_i b_i^T = lattice @ B.T = 2 pi I (for square matrices A^T B = 2 pi I <=> A B^T = 2 pi I).
    """

    assert np.allclose(tb.lattice @ grid.B.T, 2*np.pi*np.eye(3), atol=1e-12)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fourier_hermitian(tb, grid):
    """
    H(k), A(k) and dH(k)/dk are Hermitian at every k, in both gauges.

    dH is Hermitian because it is the derivative of a Hermitian matrix. In M0, A(k) is exactly
    Hermitian (and k-independent): r(R) only has the real Wannier centres on the diagonal of r(0).
    """

    H_k = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.k_cart)
    A_k = fourier(tb.r_R, grid.R_cart, tb.ndegen, grid.k_cart)
    dH_k = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.k_cart, deriv=True)

    assert np.allclose(H_k, dagger(H_k), atol=1e-12), 'H_k not hermitian'
    assert np.allclose(A_k, dagger(A_k), atol=1e-12), 'A_k not hermitian'
    assert np.allclose(dH_k, dagger(dH_k), atol=1e-12), 'dH_k not hermitian'

def test_hermitize(tb, grid):
    """
    `hermitize` really symmetrizes, is idempotent, and leaves a Hermitian input unchanged.

    A random complex (hence non-Hermitian) array is needed: on an already Hermitian input, wrong
    implementations such as the identity or the plain dagger would pass. It has the shape of A(k),
    (Nk, 3, nW, nW), to check that only the last two axes are touched. The generator is seeded so
    that a failure is reproducible. The Hermitian input is A(k) as returned by `fourier`, not a
    matrix symmetrized with `hermitize` itself (the test would then check itself).
    """

    rng = np.random.default_rng(0)
    A_k = fourier(tb.r_R, grid.R_cart, tb.ndegen, grid.k_cart)

    A = rng.random(A_k.shape) + 1j*rng.random(A_k.shape)

    assert np.allclose(hermitize(A), dagger(hermitize(A)), atol=1e-12), 'hermitize function does not hermitize the object'
    assert np.allclose(hermitize(A), hermitize(hermitize(A)), atol=1e-12), 'hermitize function not idempotent'
    assert np.allclose(A_k, hermitize(A_k), atol=1e-12), 'hermitize function breaks hermicity of hermitian object'

def finite_difference_hamiltonian(tb, grid):
    """
    Centred finite-difference gradient of H(k), reference for fourier(..., deriv=True).

        dH/dk_mu (k) ~ [H(k + DK e_mu) - H(k - DK e_mu)] / (2 DK),    mu = x, y, z

    The displaced points are built for all (k, mu) pairs at once by broadcasting,
    (Nk, 1, 3) + (1, 3, 3) -> (Nk, 3, 3), then flattened to (3 Nk, 3) because `fourier` takes a
    list of k points. Flattening in C order keeps the pairs in the order (k0, x), (k0, y), (k0, z),
    (k1, x), ...

    Returns:
        (3 Nk, nW, nW) complex, eV*Angstrom, still flattened: reshape it to (Nk, 3, nW, nW) to
        compare with dH(k). The z component is exactly 0 since all R lie in the plane.
    """

    dk = DK*np.eye(3) # row mu = DK e_mu
    k_plus = grid.k_cart[:, None, :] + dk[None, ...]
    k_minus = grid.k_cart[:, None, :] - dk[None, ...]

    Nk = grid.k_cart.shape[0]

    H_plus = fourier(tb.H_R, grid.R_cart, tb.ndegen, k_plus.reshape(3*Nk, 3))
    H_minus = fourier(tb.H_R, grid.R_cart, tb.ndegen, k_minus.reshape(3*Nk, 3))

    dH_k_fd = (H_plus - H_minus) / (2*DK)

    return dH_k_fd

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fourier_finite_difference(tb, grid):
    """
    Analytic derivative (i R_mu H(R) inside `fourier`) against centred finite differences, in both
    gauges. Tolerance 1e-6 eV*Angstrom: the observed error is ~1e-9 to 1e-8 (larger R vectors in
    the shifted gauge), against |dH| ~ 10 eV*Angstrom. A sign or 2 pi error would show up at O(1).
    """

    dH_k_df = finite_difference_hamiltonian(tb, grid)
    dH_k = fourier(tb.H_R, grid.R_cart  , tb.ndegen, grid.k_cart, deriv=True)

    assert np.allclose(dH_k_df.reshape(dH_k.shape), dH_k, atol=1e-6), 'analytical and numerical dH dont match'


def bond_vectors(tb):
    """
    Nearest-neighbour bond vectors reconstructed from the model output only.

    A hop A(0) -> B(R) links the A centre r_AA(0) to the B centre of cell R, R + r_BB(0), so

        delta = R_cart + r_BB(0) - r_AA(0)

    for every R where H_AB(R) != 0. Nothing is read from inside make_graphene_tb: the function
    checks the geometry that the WannierTB actually encodes. With shift_B = L, the hop label
    R = s - L and the stored B centre tau_B + L change together and L cancels, so the bonds must
    come out identical in every gauge. The same reconstruction applies to Wannier90 data (R list
    and Wannier centres from `_tb.dat`).

    R_cart and the row of R = 0 are recomputed from `tb` (not taken from the `grid` fixture), so the
    function can be applied to two different models in the same test.

    Returns:
        deltas : (n_hops, 3) float, Angstrom, one bond vector per row (3 for graphene)
    """
    R_cart = tb.R_int @ tb.lattice.T
    index = {tuple(R): iR for iR, R in enumerate(tb.R_int)}

    iR0 = index[(0,0,0)]
    r_0AA = tb.r_R[iR0,:,0,0].real; r_0BB = tb.r_R[iR0,:,1,1].real # centres are real vectors
    # A -> B hops: rows where the AB element is non-zero (abs + threshold rather than != 0, since
    # the values are complex and real data are never exactly zero)
    mask = np.abs(tb.H_R[:,0,1]) > 1e-12

    deltas = R_cart[mask] + r_0BB - r_0AA

    return deltas

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_bond_length(tb):
    """
    Exactly three A -> B hops, each of length a_cc, in both gauges. The norm is taken along axis=1
    (the x, y, z components of each row), giving one length per bond.
    """
    deltas = bond_vectors(tb)

    assert len(deltas) == 3, 'must have three non-zero hopping for first neighbhour graphene tb'
    assert np.allclose(np.linalg.norm(deltas, axis=1), tb.a_cc, atol=1e-12), 'carbon-carbon bond length not correct'

def sort_by_angle(d):
    """
    Reorder a list of in-plane vectors (rows) by polar angle, to compare two lists that contain the
    same vectors in different orders.

    np.argsort on a single key followed by row indexing keeps every vector intact, whereas
    np.sort(d, axis=0) would sort the x and y columns independently and mix vectors. The angle is
    arctan2(y, x) (not arctan, whose second positional argument is the output array). For graphene
    the angles are 30, 150 and -90 degrees, far from each other and from the +-180 degree branch
    cut, so the order is robust to round-off.
    """

    angles = np.arctan2(d[:,1], d[:,0])
    order = np.argsort(angles)
    return d[order]

def test_tb_shift_same_bonds():
    """
    The relabelling shift_B changes the R list but not the physics.

    Without shift the R list is {0, +-a1, +-a2} (nR = 5); with shift_B = (1,0,0) the hops move to
    S - L and L - S and R = 0 only carries the centres (nR = 7). The three bond vectors must be the
    same in both models, up to their order. This test needs both models at once, so it builds them
    itself instead of using the parametrized `tb` fixture.
    """
    tb = make_graphene_tb()
    tb_shifted = make_graphene_tb(shift_B=(1,0,0))

    assert len(tb.R_int) == 5
    assert len(tb_shifted.R_int) == 7

    deltas = bond_vectors(tb)
    deltas_shifted = bond_vectors(tb_shifted)

    assert np.allclose(sort_by_angle(deltas), sort_by_angle(deltas_shifted), atol=1e-12)
