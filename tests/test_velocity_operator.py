"""
Tests of the M0 building blocks of optics/velocity_operator.py, on the analytic graphene model:
the model (R-space Hermiticity, bonds, shift_B), the reciprocal lattice, `fourier` (closed form,
Dirac point, Hermiticity, finite differences), dagger/hermitize, `velocity` (F5) and `ring` (F6).

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture: gauge-independent properties must hold in both. Reference values and the
derivations behind the tolerances: memoire/EM/EM.md (M0, implementation notes).

Run with:  .venv/bin/python -m pytest tests/test_velocity_operator.py -v
"""

import pytest
import numpy as np
from electron_defect_interaction.optics.velocity_operator import *


N = 100    # the k grid has N x N = 1e4 points: fast, yet covers the whole Brillouin zone
DK = 1e-5  # finite-difference step (1/Angstrom); truncation ~DK^2 and round-off ~eps/DK balance here
HW = 2.33  # eV, 532 nm laser

@pytest.fixture
def tb(request):
    """
    Graphene model, shift_B = (0,0,0) by default. Other gauges through
    @pytest.mark.parametrize("tb", [...], indirect=True), which puts each value in request.param.
    """
    shift = getattr(request, "param", (0,0,0))
    return make_graphene_tb(shift_B=shift)

@pytest.fixture
def grid(tb):
    """
    N x N GridTB of `tb`, rebuilt for every gauge of a parametrized test.
    """
    return make_grid_tb(tb, N)

def make_graphene_tb_analytic(tb, grid):
    """
    Closed form H_AB(k) = -t (1 + e^{-ik.a1} + e^{-ik.a2}), reference for `fourier`. Valid for
    shift_B = 0 only (a shift L adds the phase e^{-ik.L}). Returns (Nk,) complex, eV.
    """

    t = tb.t
    A = tb.lattice
    k = grid.k_cart

    return -t*(1+np.exp(-1j*k @ A[:,0].T) + np.exp(-1j*k @ A[:,1].T))

def test_fourier_analytic(tb, grid):
    """
    `fourier` gives the closed-form H_AB(k) on the whole grid: checks the phase sign, ndegen,
    R_cart and the orbital order at once.
    """

    H_tb = make_graphene_tb_analytic(tb, grid)
    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart)

    assert np.allclose(H_tb, H_k[:,0,1], atol=1e-12), 'tight binding model not correct'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_dirac_point(tb, grid):
    """
    H(K) = 0 at K = (2 b1 + b2)/3, in both gauges.
    """

    H_K = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.K[None, :])

    assert np.allclose(H_K, 0, atol=1e-12), 'Hamiltonian not zero at Dirac point'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_hermitian(tb, grid):
    """
    R-space Hermiticity X_ij(R) = X_ji(-R)* for H and r (dagger on one side only), R list closed
    under R -> -R, and ndegen(-R) = ndegen(R): together they make H(k) and A(k) Hermitian.
    """

    assert np.allclose(tb.minus[tb.minus], np.arange(len(tb.R_int))), 'R list not closed under R -> -R'
    assert np.allclose(tb.H_R, tb.H_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'TB Hamiltonian not hermitian'
    assert np.allclose(tb.r_R, tb.r_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'positian operator not hermitian'
    assert np.allclose(tb.ndegen, tb.ndegen[tb.minus], atol=1e-12), 'R and -R dont have the same weight'

def test_reciprocal(tb, grid):
    """
    a_i . b_j = 2 pi delta_ij, written as lattice @ B.T = 2 pi I.
    """

    assert np.allclose(tb.lattice @ grid.B.T, 2*np.pi*np.eye(3), atol=1e-12)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fourier_hermitian(tb, grid):
    """
    H(k), A(k) and dH/dk are Hermitian at every k, in both gauges.
    """

    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart)
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, grid.k_cart)
    dH_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, grid.k_cart, deriv=True)

    assert np.allclose(H_k, dagger(H_k), atol=1e-12), 'H_k not hermitian'
    assert np.allclose(A_k, dagger(A_k), atol=1e-12), 'A_k not hermitian'
    assert np.allclose(dH_k, dagger(dH_k), atol=1e-12), 'dH_k not hermitian'

def test_hermitize(tb, grid):
    """
    `hermitize` symmetrizes a random (non-Hermitian) array of the shape of A(k), is idempotent, and
    leaves the Hermitian A(k) of M0 unchanged.
    """

    rng = np.random.default_rng(0)
    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, grid.k_cart)

    A = rng.random(A_k.shape) + 1j*rng.random(A_k.shape)

    assert np.allclose(hermitize(A), dagger(hermitize(A)), atol=1e-12), 'hermitize function does not hermitize the object'
    assert np.allclose(hermitize(A), hermitize(hermitize(A)), atol=1e-12), 'hermitize function not idempotent'
    assert np.allclose(A_k, hermitize(A_k), atol=1e-12), 'hermitize function breaks hermicity of hermitian object'

def finite_difference_hamiltonian(tb, grid):
    """
    Centred finite difference [H(k + DK e_mu) - H(k - DK e_mu)] / (2 DK) for all (k, mu) at once.
    Returns (3 Nk, nW, nW) complex, eV*Angstrom, ordered (k0, x), (k0, y), (k0, z), (k1, x), ...:
    reshape to (Nk, 3, nW, nW).
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
    Analytic dH/dk of `fourier` against finite differences, in both gauges. Observed 1e-9 to 1e-8
    eV*Angstrom for |dH| ~ 10; a sign or 2 pi error would show at O(1).
    """

    dH_k_df = finite_difference_hamiltonian(tb, grid)
    dH_k = fourier(tb.H_R, tb.R_cart  , tb.ndegen, grid.k_cart, deriv=True)

    assert np.allclose(dH_k_df.reshape(dH_k.shape), dH_k, atol=1e-6), 'analytical and numerical dH dont match'


def bond_vectors(tb):
    """
    Bond vectors delta = R_cart + r_BB(0) - r_AA(0) for every R with H_AB(R) != 0, rebuilt from the
    WannierTB alone (a shift L cancels between R and the B centre). Returns (n_hops, 3), Angstrom.
    """
    R_cart = tb.R_int @ tb.lattice.T
    index = {tuple(R): iR for iR, R in enumerate(tb.R_int)}

    iR0 = index[(0,0,0)]
    r_0AA = tb.r_R[iR0,:,0,0].real; r_0BB = tb.r_R[iR0,:,1,1].real # centres are real vectors
    # A -> B hops: rows where the AB element is non-zero
    mask = np.abs(tb.H_R[:,0,1]) > 1e-12

    deltas = R_cart[mask] + r_0BB - r_0AA

    return deltas

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_bond_length(tb):
    """
    Three A -> B hops of length a_cc, in both gauges.
    """
    deltas = bond_vectors(tb)

    assert len(deltas) == 3, 'must have three non-zero hopping for first neighbhour graphene tb'
    assert np.allclose(np.linalg.norm(deltas, axis=1), tb.a_cc, atol=1e-12), 'carbon-carbon bond length not correct'

def sort_by_angle(d):
    """
    Rows of d sorted by polar angle arctan2(y, x), to compare lists of vectors given in any order.
    """

    angles = np.arctan2(d[:,1], d[:,0])
    order = np.argsort(angles)
    return d[order]

def test_tb_shift_same_bonds():
    """
    shift_B = (1,0,0) changes the R list (nR = 5 -> 7) but not the bonds. Builds both models itself.
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
    fourier (H, dH, A) -> hermitize (A) -> velocity at the k points `k`, with or without the Berry
    term: what the driver will do on each block. Returns H_k and the (eps, V, hv) of `velocity`.
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
    V unitary, V^dag H V = diag(eps) and bands in ascending order, in both gauges (the checks kept
    out of `velocity`).
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
    hbar v Hermitian for each gauge, with and without Berry. Catches a missing i in the Berry term,
    which the diagonal test cannot see. Observed 4e-15.
    """

    hv = velocity_from_tb(tb, grid.k_cart, berry)[-1]

    assert np.allclose(dagger(hv), hv, atol=1e-12), 'velocity operator not hermitian'

def test_velocity_inputs_unchanged(tb, grid):
    """
    `velocity` leaves its inputs untouched (regression: W = dH_k without copy overwrote dH_k).
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
    Centred finite difference of the band energies (eigvalsh at k +/- DK e_mu), reference for the
    diagonal of hbar v. Valid away from K, where the sorted bands do not cross. Returns
    (Nk, 3, nW) float, eV*Angstrom.
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
    Diagonal of hbar v = d eps_n/dk, with and without Berry (the Berry term vanishes on the
    diagonal), in both gauges. Points with a gap below 1 eV are masked: near K the finite-difference
    error grows like DK^2 v_F / q^2 (1.5e-6 next to K, 1.7e-8 elsewhere). rtol = 0, since the default
    rtol would loosen atol to ~6e-5.
    """

    eps_fd = finite_difference_eps(tb, grid) # (nk, 3, nW)
    _, eps, _, hv = velocity_from_tb(tb, grid.k_cart, berry) # (nk, nW), (nk, 3, nW, nW)

    mask = eps[:,1] - eps[:,0] > 1 # far from K and K'
    hv_diag = np.diagonal(hv, axis1=-2, axis2=-1) # (nk, 3, nW)

    assert np.allclose(np.imag(hv_diag), 0, atol=1e-12), 'intraband velocity operator eigenvalues must be real !'
    assert np.allclose(np.real(hv_diag[mask]), eps_fd[mask], atol=1e-7, rtol=0), ' diagonal of velocity operator does not match finite differenced eigenvalues'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fermi_velocity(tb, grid):
    """
    hbar v_F = 3 t a_cc / 2 at q = 1e-3 1/Angstrom around K (12 angles, both gauges), from the
    interband |hbar v_cv| and from the group velocity |grad eps_c|. Point by point, trigonal warping
    gives up to 7e-4 (rtol 2e-3); it cancels in the mean over theta (3.8e-7, rtol 1e-5).
    """

    ntheta = 12
    q = 1e-3 # 1/Angstrom

    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False) # (ntheta,)
    k = grid.K[None, :] + q*np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))) # (ntheta, 3)

    hv = velocity_from_tb(tb, k, berry=True)[-1] # (ntheta, 3, nW, nW)
    # v = band 0, c = band 1 in M0
    interband_norm = np.sqrt(np.abs(hv[:, 0, 1, 0])**2 + np.abs(hv[:, 1, 1, 0])**2) # (ntheta,)
    group_velocity = np.sqrt(np.abs(np.real(hv[:, 0, 1, 1]))**2 + np.abs(np.real(hv[:, 1, 1, 1]))**2) # (ntheta, )

    vF = 3 *tb.t *tb.a_cc / 2 # hbar v_F, eV*Angstrom

    assert np.allclose(interband_norm, vF, rtol=2e-3, atol=0)
    assert np.allclose(group_velocity, vF, rtol=2e-3, atol=0)
    assert np.isclose(interband_norm.mean(), vF, rtol=1e-5, atol=0)
    assert np.isclose(group_velocity.mean(), vF, rtol=1e-5, atol=0)

def test_velocity_gauge(grid):
    """
    Gauge invariance under shift_B = a1. With Berry, eps and |hbar v_mn|^2 are unchanged (3.7e-13;
    moduli because the eigenvector phases are arbitrary); without Berry, the x component changes (up
    to ~800 eV^2 Angstrom^2) and y does not, since L has no y component. Proof in EM.md. The grid
    comes from the default fixture: it only depends on the lattice, identical in both models.
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

# ---------------------------------------------------------------------------------------------
# resonant ring (F6)
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_on_shell(tb, grid, hw):
    """
    `ring` lies on the resonance: |k - K| = q, and eps_c - eps_v = hw recomputed through
    `velocity_from_tb`, not `_gap` (observed 5e-12 eV).
    """

    k, theta, q, q0 = ring(tb, grid.K, hw)

    assert np.allclose(q, np.linalg.norm(k-grid.K[None, :], axis=1), atol=1e-12)

    eps = velocity_from_tb(tb, k)[1]
    gap = eps[:, 1] - eps[:, 0]

    assert np.allclose(gap, hw, atol=1e-10, rtol=0)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_ring_dirac_limit(tb, grid):
    """
    At hw = 0.1 eV the ring is the Dirac circle q = q0 (observed q/q0 in [0.997, 1.003]).
    """

    k, theta, q, q0 = ring(tb, grid.K, hw=0.1)

    assert np.allclose(q/q0, np.ones(len(q)), rtol=5e-3, atol=0)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_symmetry(tb, grid, hw):
    """
    q(theta) has the C3 symmetry and the mirror symmetry theta -> -theta of the bands around K.
    """

    # test C3 symmetry q(theta + 120) = q(theta)
    k, theta, q, q0 = ring(tb, grid.K, hw)
    ntheta = len(theta)

    assert np.allclose(np.roll(q, ntheta//3), q, atol=1e-12)

    # test mirror symmetry q(-theta) = q(theta): index of -theta_j is (n - j) mod n
    assert np.allclose(np.roll(q[::-1], 1), q, atol=1e-12)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_node(tb, grid, hw):
    """
    With Berry, |hbar v^x_cv|^2 vanishes at theta = 0 and pi (q // x): the node of Gruneis 2003.
    """
    k = ring(tb, grid.K, hw)[0]
    ntheta = k.shape[0]
    hv = velocity_from_tb(tb, k)[-1]
    assert np.allclose(np.abs(hv[0, 0, 1, 0])**2, 0, atol=1e-12) and np.allclose(np.abs(hv[ntheta // 2, 0, 1, 0])**2, 0, atol=1e-12)

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_isotropy(tb, grid, hw):
    """
    With Berry, <|v^x_cv|^2> = <|v^y_cv|^2> and <Re(v^x_cv* v^y_cv)> = 0 on the ring (C3), which gives
    sigma_xx = sigma_yy and sigma_xy = 0 in F7. |.|^2 is taken before the mean over theta.
    """
    k = ring(tb, grid.K, hw)[0]
    hv_cv = velocity_from_tb(tb, k)[-1][:, :, 1, 0] # (ntheta, 3): element [c, v] = [1, 0]

    vx2 = np.mean(np.abs(hv_cv[:, 0])**2) # eV^2 Angstrom^2
    vy2 = np.mean(np.abs(hv_cv[:, 1])**2)
    vxy = np.mean(np.real(np.conj(hv_cv[:, 0]) * hv_cv[:, 1]))

    assert np.isclose(vx2, vy2, rtol=1e-10, atol=0), '<|v_x|^2> != <|v_y|^2> on the ring'
    assert abs(vxy) < 1e-10 * vx2, '<Re(v_x* v_y)> != 0 on the ring'

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw, ratio", [(0.1, 1.0), (2.33, 0.98909)])
def test_ring_average(tb, grid, hw, ratio):
    """
    <|hbar v^x_cv|^2> / ((hbar v_F)^2 / 2) on the ring against independent references: 1 in the
    Dirac limit (0.1 eV), 0.989 at 2.33 eV (EM.md, trigonal warping). No reference at 1.0 or 2.54 eV.
    """
    k = ring(tb, grid.K, hw)[0]
    hv = velocity_from_tb(tb, k)[-1] # (ntheta, 3, nW, nW)
    hv_F = 3*tb.t*tb.a_cc/2 # eV Angstrom

    vx2 = np.mean(np.abs(hv[:, 0, 1, 0])**2) # |.|^2 first, then the mean over theta

    assert np.isclose(vx2/(hv_F**2/2), ratio, rtol=1e-4, atol=0), 'ring average of |v_x|^2 off its reference'

def test_ring_without_berry(tb, grid):
    """
    Without Berry, at 2.33 eV with B in cell 0 (EM.md): no node at theta = 0, pi (|v^x_cv|^2 = 2.05);
    two nodes at 193.6 and 343.4 degrees (samples 193.5 and 343.5); <|v^x|^2> without / with 1.125;
    |v_cv| without / with in [0.664, 1.271]. Gauge-dependent numbers, hence the default `tb` only.
    """
    k, theta = ring(tb, grid.K, 2.33)[:2]
    ntheta = len(theta)
    hv = velocity_from_tb(tb, k, berry=True)[-1][:, :2, 1, 0] # (ntheta, 2): x, y of [c, v]
    hv_nb = velocity_from_tb(tb, k, berry=False)[-1][:, :2, 1, 0]

    vx2_nb = np.abs(hv_nb[:, 0])**2 # (ntheta,)

    assert vx2_nb[0] > 1 and vx2_nb[ntheta // 2] > 1, 'node still at theta = 0 or pi without Berry'

    # local minima: lower than both neighbours (np.roll wraps around theta = 2 pi)
    minima = (vx2_nb < np.roll(vx2_nb, 1)) & (vx2_nb < np.roll(vx2_nb, -1))
    nodes = np.degrees(theta[minima])
    assert len(nodes) == 2, f'expected 2 nodes without Berry, found {nodes} degrees'
    assert np.allclose(nodes, [193.5, 343.5], atol=0.5), f'nodes without Berry at {nodes} degrees'

    ratio_x = np.mean(vx2_nb) / np.mean(np.abs(hv[:, 0])**2)
    ratio_norm = np.linalg.norm(hv_nb, axis=1) / np.linalg.norm(hv, axis=1) # |v_cv| in the plane, per theta

    assert np.isclose(ratio_x, 1.125, rtol=1e-3, atol=0), '<|v_x|^2> without / with Berry off 1.125'
    assert np.isclose(ratio_norm.min(), 0.664, atol=1e-3), '|v_cv| without / with Berry: minimum off 0.664'
    assert np.isclose(ratio_norm.max(), 1.271, atol=1e-3), '|v_cv| without / with Berry: maximum off 1.271'
