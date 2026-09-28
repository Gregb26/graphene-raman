"""
Tests of velocity_operator: `fourier` (closed form, Dirac point, Hermiticity, finite differences),
dagger/hermitize, and `velocity` (F5).
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import dagger, fourier, hermitize, make_graphene_tb, velocity
from conftest import velocity_from_tb  # temporary: will be replaced by the src chain function

DK = 1e-5  # finite-difference step (1/Angstrom); truncation ~DK^2 and round-off ~eps/DK balance here

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
