"""
Tests of the M0 building blocks of the velocity operator (optics/velocity_operator.py).

Everything is checked on the analytic graphene tight-binding model, whose answers are known:
    - the model itself: Hermiticity in R space, bond lengths, invariance of the bonds under the
      relabelling shift_B, number of R vectors;
    - the reciprocal lattice;
    - the Fourier transform: closed-form H_AB(k), Dirac point, Hermiticity of H(k), dH(k) and A(k),
      and dH(k) against centred finite differences;
    - dagger / hermitize;
    - the velocity operator (F5): diagonalization, Hermiticity of hbar v, inputs left untouched,
      diagonal = gradient of the bands, Fermi velocity 3 t a_cc / 2 at the Dirac cone, and gauge
      invariance of |hbar v_mn|^2 under shift_B with the Berry term (broken without it).

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture (see `tb`): every property that does not depend on how the B orbital is labelled
must hold in both.

Run with:  .venv/bin/python -m pytest tests/test_velocity_operator.py -v
"""

import pytest
import numpy as np
from electron_defect_interaction.optics.velocity_operator import *


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
    k grid of `tb` (GridTB: B, k_red, k_cart, K) with N x N points.

    It depends on the `tb` fixture, so pytest rebuilds it for every gauge of a parametrized test.
    """
    return make_grid_tb(tb, N)

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
        grid : dataclass from the `grid` fixture
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
    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart)

    assert np.allclose(H_tb, H_k[:,0,1], atol=1e-12), 'tight binding model not correct'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_dirac_point(tb, grid):
    """
    The whole H(K) vanishes at the Dirac point K = (2 b1 + b2)/3: the three neighbour phases cancel,
    1 + e^{-i4pi/3} + e^{-i2pi/3} = 0, and the on-site energies are zero. Holds in both gauges,
    since the shift only multiplies H_AB by a phase. K is passed as a (1, 3) array (K[None, :])
    because `fourier` expects a list of k points.
    """

    H_K = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.K[None, :])

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

    assert np.allclose(tb.minus[tb.minus], np.arange(len(tb.R_int))), 'R list not closed under R -> -R'
    assert np.allclose(tb.H_R, tb.H_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'TB Hamiltonian not hermitian'
    assert np.allclose(tb.r_R, tb.r_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'positian operator not hermitian'
    assert np.allclose(tb.ndegen, tb.ndegen[tb.minus], atol=1e-12), 'R and -R dont have the same weight'

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

    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart)
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, grid.k_cart)
    dH_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart, deriv=True)

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
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, grid.k_cart)

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

    H_plus = fourier(tb.H_R, tb.R_cart, tb.ndegen, k_plus.reshape(3*Nk, 3))
    H_minus = fourier(tb.H_R, tb.R_cart, tb.ndegen, k_minus.reshape(3*Nk, 3))

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
    dH_k = fourier(tb.H_R, tb.R_cart  , tb.ndegen, grid.k_cart, deriv=True)

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

# ---------------------------------------------------------------------------------------------
# Velocity operator (F5)
# ---------------------------------------------------------------------------------------------

def velocity_from_tb(tb, k, berry=True):
    """
    Whole chain from the tight-binding model to hbar v at the k points `k`.

    Three Fourier transforms (H, dH/dk, A) followed by `velocity`, with or without the Berry term.
    This is what the driver will do on every block of k points.

    Inputs:
        tb : WannierTB
        k : (Nk, 3) float, 1/Angstrom, Cartesian k points in rows (any list, not only the grid)
        berry : bool
            True: hbar v = V^dag (dH + i[H, A]) V. False: the Berry term is dropped.
    Returns:
        H_k : (Nk, nW, nW) complex, eV, H(k) in the Wannier gauge (needed by test_velocity_diagonalizes)
        eps, V, hv : as returned by `velocity`
    """

    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, k) # (nk, nW, nW)
    dH_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, k, deriv=True) # (nk, 3, nW, nW)
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, k) # (nk, 3, nW, nW)
    A_k = hermitize(A_k)

    if berry:
        eps, V, hv = velocity(H_k, dH_k, A_k)
    else:
        eps, V, hv = velocity(H_k, dH_k)

    return H_k, eps, V, hv

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_velocity_diagonalizes(tb, grid):
    """
    The eigh step inside `velocity`, in both gauges: V is unitary, V^dag H V = diag(eps), and the
    bands are sorted in ascending energy.

    These are the checks that would cost a diagonalization per block if they sat inside `velocity`.
    diag(eps) is built as a stack of diagonal matrices by broadcasting, (Nk, nW, 1) * (nW, nW) ->
    [k, m, n] = eps_m(k) delta_mn (np.diag only handles a single matrix). The sorting is what the
    other tests rely on: in M0, v = band 0 and c = band 1 around mu = 0, and the gap is
    eps[:, 1] - eps[:, 0].
    """
    H_k, eps, V, hv = velocity_from_tb(tb, grid.k_cart)
    _, nW = eps.shape

    assert np.allclose(dagger(V) @ V, np.eye(nW), atol=1e-12), 'V not unitary'
    assert np.allclose(dagger(V) @ H_k @ V, eps[..., None] * np.eye(nW), atol=1e-12), 'some other problem with V'
    assert np.all(np.diff(eps, axis=1) >= 0), 'eigenvalues must be increasing (why?)'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("berry", [True, False])
def test_velocity_hermitian(tb, grid, berry):
    """
    hbar v is Hermitian, for the 4 combinations gauge x (with, without Berry).

    Two stacked parametrize decorators give the Cartesian product of their values: `tb` is indirect
    (the shift goes through the fixture), `berry` is direct (True/False arrives as is in the
    argument). This is the test that catches a missing factor i in the Berry term: [H, A] alone is
    anti-Hermitian, and the diagonal test cannot see it since the Berry term vanishes on the
    diagonal. Observed defect ~4e-15 eV*Angstrom.
    """

    hv = velocity_from_tb(tb, grid.k_cart, berry)[-1] # last returned value; the others are not needed

    assert np.allclose(dagger(hv), hv, atol=1e-12), 'velocity operator not hermitian'

def test_velocity_inputs_unchanged(tb, grid):
    """
    `velocity` does not modify its input arrays (regression test).

    With W = dH_k instead of dH_k.copy(), W and dH_k are two names for the same array, and
    W += i[H, A] silently overwrites the caller's dH_k: any later use of it (the "without Berry"
    comparison, a finite difference) would then be wrong. The references are taken with .copy()
    before the call (a plain `=` would alias them too and the test could never fail), and compared
    with np.array_equal: nothing may change, not even at 1e-16. Only the Berry branch can overwrite
    dH_k, so the test calls `velocity` with A_k, in the default gauge.
    """

    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart) # (nk, nW, nW)
    dH_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart, deriv=True) # (nk, 3, nW, nW)
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, grid.k_cart) # (nk, 3, nW, nW)

    H_k_ref = H_k.copy(); dH_k_ref = dH_k.copy(); A_k_ref = A_k.copy()
    _ = velocity(H_k, dH_k, A_k)

    assert np.array_equal(H_k, H_k_ref)
    assert np.array_equal(dH_k, dH_k_ref)
    assert np.array_equal(A_k, A_k_ref)

def finite_difference_eps(tb, grid):
    """
    Centred finite-difference gradient of the band energies, reference for the diagonal of hbar v.

        d eps_n/dk_mu (k) ~ [eps_n(k + DK e_mu) - eps_n(k - DK e_mu)] / (2 DK)

    Same displaced k points as finite_difference_hamiltonian, (Nk, 3, 3) flattened to (3 Nk, 3),
    but only the eigenvalues are needed, hence eigvalsh (cheaper than eigh). Both eigvalsh calls
    return the bands in ascending order, so band n is subtracted from band n: fine as long as the
    two bands do not cross, i.e. away from K. The result is reshaped back to (Nk, 3, nW), the shape
    of the diagonal of hbar v.

    Returns:
        (Nk, 3, nW) float, eV*Angstrom
    """
    nk, d = grid.k_cart.shape
    dk = DK*np.eye(d) # (3,3)

    k_plus  = grid.k_cart[:, None, :] + dk[None, ...] # (nk, 3, 3)
    k_minus = grid.k_cart[:, None, :] - dk[None, ...]  # (nk, 3, 3)

    H_k_plus = fourier(tb.H_R, tb.R_cart, tb.ndegen, k_plus.reshape(d*nk, d))
    H_k_minus = fourier(tb.H_R, tb.R_cart, tb.ndegen, k_minus.reshape(d*nk, d))

    eps_plus = np.linalg.eigvalsh(H_k_plus)
    eps_minus = np.linalg.eigvalsh(H_k_minus)

    eps_fd = (eps_plus - eps_minus) / (2*DK)

    return eps_fd.reshape(nk, d, -1) # (nk, 3, nW)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("berry", [True, False])
def test_velocity_diagonal_gradient(tb, grid, berry):
    """
    Diagonal of hbar v = group velocity d eps_n/dk (Hellmann-Feynman), with and without Berry, in
    both gauges.

    Passing without Berry too is expected: the Berry term i(eps_m - eps_n) Abar_mn vanishes for
    m = n, so this test checks dH/dk and the rotation V^dag ... V, not the Berry term (see
    test_velocity_hermitian and test_velocity_gauge for that). The diagonal of a Hermitian matrix is
    real, which is checked first.

    Tolerance: atol = 1e-7 eV*Angstrom with rtol = 0 (np.allclose tests |a - b| <= atol + rtol |b|,
    and the default rtol = 1e-5 would silently loosen it to ~6e-5 here). Observed error 1.7e-8 on
    the masked points. Near K the cone eps = +-hbar v_F |q| curves like 1/q and the truncation error
    of the centred difference grows like DK^2 v_F / q^2: 1.5e-6 at the two grid points closest to K
    and K' (gap 0.098 eV, q ~ 8.5e-3 1/Angstrom). The mask gap > 1 eV removes them, in both valleys
    at once, without any geometry.
    """

    eps_fd = finite_difference_eps(tb, grid) # (nk, 3, nW)
    _, eps, _, hv = velocity_from_tb(tb, grid.k_cart, berry) # (nk, nW), (nk, 3, nW, nW)

    # "far from K": gap above 1 eV. Boolean (nk,) mask, applied on the k axis of (nk, 3, nW) arrays
    mask = eps[:,1] - eps[:,0] > 1
    hv_diag = np.diagonal(hv, axis1=-2, axis2=-1) # (nk, 3, nW): the diagonal axis goes last

    assert np.allclose(np.imag(hv_diag), 0, atol=1e-12), 'intraband velocity operator eigenvalues must be real !'
    assert np.allclose(np.real(hv_diag[mask]), eps_fd[mask], atol=1e-7, rtol=0), ' diagonal of velocity operator does not match finite differenced eigenvalues'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fermi_velocity(tb, grid):
    """
    Fermi velocity of the Dirac cone, hbar v_F = 3 t a_cc / 2 = 5.751 eV*Angstrom, in both gauges.

    On a circle of radius q = 1e-3 1/Angstrom around K (12 equally spaced angles, never K itself,
    where c and v are degenerate), two independent quantities must equal hbar v_F:
        - interband: sqrt(|hbar v^x_cv|^2 + |hbar v^y_cv|^2), element [c, v] = [1, 0]. For the
          Dirac Hamiltonian hbar v_F (sigma . q) this is exactly hbar v_F at every angle, while each
          component alone varies like sin or cos of theta (the node of test_node, F6);
        - intraband: |grad eps_c| = |(Re hbar v^x_cc, Re hbar v^y_cc)|, the group velocity.
    With Berry, both are gauge invariant, hence the two gauges.

    Tolerances are relative (rtol, atol = 0) because the deviation is proportional to the quantity:
        - point by point, rtol = 2e-3: trigonal warping gives a correction ~ q a_cc cos(3 theta),
          observed 7.1e-4 at worst;
        - averaged over theta, rtol = 1e-5: the cos(3 theta) term cancels exactly on 12 equally
          spaced angles and only the O(q^2) term remains, observed 3.8e-7. This is the real test of
          hbar v_F; the point-by-point one catches an aberrant angle.
    """

    ntheta = 12
    q = 1e-3 # 1/Angstrom

    # k = K + q (cos theta, sin theta, 0): column_stack puts the three (ntheta,) arrays in columns,
    # one k point per row; the z component must be an array of zeros, not the scalar 0
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False) # (ntheta,)
    k = grid.K[None, :] + q*np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))) # (ntheta, 3)

    hv = velocity_from_tb(tb, k, berry=True)[-1] # (ntheta, 3, nW, nW)
    # in M0 bands are sorted and mu = 0 lies between them: v = 0, c = 1
    interband_norm = np.sqrt(np.abs(hv[:, 0, 1, 0])**2 + np.abs(hv[:, 1, 1, 0])**2) # (ntheta,)
    group_velocity = np.sqrt(np.abs(np.real(hv[:, 0, 1, 1]))**2 + np.abs(np.real(hv[:, 1, 1, 1]))**2) # (ntheta, )

    vF = 3 *tb.t *tb.a_cc / 2 # hbar v_F, eV*Angstrom, read from the model rather than hard-coded

    assert np.allclose(interband_norm, vF, rtol=2e-3, atol=0)
    assert np.allclose(group_velocity, vF, rtol=2e-3, atol=0)
    assert np.isclose(interband_norm.mean(), vF, rtol=1e-5, atol=0)
    assert np.isclose(group_velocity.mean(), vF, rtol=1e-5, atol=0)

def test_velocity_gauge(grid):
    """
    Gauge invariance of the velocity under the relabelling shift_B = L = a1: exact with the Berry
    term, broken without it. It is the only test that checks the VALUE of i[H, A], not just its
    form (the F5 counterpart of test_gauge_shift, without the Kubo sum).

    Why it is exact: moving B by L gives H'(k) = U H(k) U^dag with U = diag(1, e^{ik.L}). The
    derivative of U adds i L_mu [P, H] to dH' (P = projector on B), and the shifted centre
    tau_B + L adds exactly the opposite i L_mu [H, P] to i[H', A']. So dH' + i[H', A'] =
    U (dH + i[H, A]) U^dag, V' = U V, and hbar v' = hbar v up to the arbitrary phase eigh gives each
    eigenvector: hbar v_mn may change phase, |hbar v_mn|^2 may not, hence the comparison of squared
    moduli. Observed difference 3.7e-13 eV^2*Angstrom^2 for |hbar v|^2 up to ~59.

    Without Berry the extra term i L_mu [P, H] stays: the x component changes (by up to ~800
    eV^2*Angstrom^2), the y component does not, because L = a1 has no y component. The `not
    allclose` makes sure the shift really reaches the velocity (a shift_B ignored somewhere would
    pass the Berry part and fail here).

    The test needs both models at once, so it builds them itself. The `grid` fixture is built from
    the default `tb` fixture (shift 0), unrelated to the local `tb`; that is fine because the grid
    only depends on the lattice, identical in both models.
    """
    tb = make_graphene_tb()
    tb_shifted = make_graphene_tb(shift_B=(1,0,0))

    k = grid.k_cart

    # with Berry: same bands, same |hbar v_mn|^2 for every (k, mu, m, n)
    _, eps, _, hv = velocity_from_tb(tb, k, berry=True)
    _, eps_shifted, _, hv_shifted = velocity_from_tb(tb_shifted, k, berry=True)

    assert np.allclose(eps, eps_shifted, atol=1e-10)
    assert np.allclose(np.abs(hv)**2, np.abs(hv_shifted)**2, atol=1e-10)

    # without Berry: x (mu = 0) must differ, y (mu = 1) must not
    _, eps, _, hv = velocity_from_tb(tb, k, berry=False)
    _, eps_shifted, _, hv_shifted = velocity_from_tb(tb_shifted, k, berry=False)

    assert not np.allclose(np.abs(hv[:, 0, ...])**2, np.abs(hv_shifted[:, 0, ...])**2, atol=1e-10)
    assert np.allclose(np.abs(hv[:, 1, ...])**2, np.abs(hv_shifted[:, 1, ...])**2, atol=1e-10)



    




