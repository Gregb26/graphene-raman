"""
Tests of the tight-binding model (tb_model): R-space Hermiticity, bonds, relabelling shift_B of the
M0 graphene model; construction of the real 27 x 27 model by make_wannier_tb and its centres-only
version centres_only (M1, fixtures tb_w90 and eig_w90 of conftest.py).
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import (centres_only, compute_velocity, dagger, fourier,
                                                         hermitize, make_graphene_tb, make_grid_tb,
                                                         reciprocal, ring)

E_D = -4.238895  # eV, Dirac point of the 27 x 27 data at K
FROZ_MAX = -1.74 # eV, top of the frozen window of the 27 x 27 wannierisation
HW = 2.33        # eV, 532 nm laser

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


def test_wannier_tb_fields(tb_w90):
    """
    make_wannier_tb fills the derived fields consistently: R_cart = R_int @ lattice.T, index and
    minus invert the R list, and H is Hermitian in R space (r is not, see the Berry test).
    """
    tb = tb_w90

    assert tb.H_R.shape == (741, 5, 5) and tb.r_R.shape == (741, 3, 5, 5)
    assert np.allclose(tb.R_cart, tb.R_int @ tb.lattice.T, rtol=0, atol=1e-12)
    assert all(tb.index[tuple(R)] == iR for iR, R in enumerate(tb.R_int))
    assert np.array_equal(tb.R_int[tb.minus], -tb.R_int)
    assert np.allclose(tb.H_R, tb.H_R[tb.minus].swapaxes(-1, -2).conj(), rtol=0, atol=1e-12)


def test_wannier_tb_bands_match_eig(tb_w90, eig_w90):
    """
    Whole Cartesian chain on the real model: k_cart = k_red @ B.T, then compute_velocity, reproduces
    the .eig in the frozen window at the 729 coarse points (1.3e-5 eV, as the transform in reduced
    coordinates of test_wannier_io). Checks k_cart . R_cart = 2 pi k_red . R_int on the real lattice.
    """
    k_red, E_dft = eig_w90
    B = reciprocal(tb_w90.lattice)

    _, eps, _, _ = compute_velocity(tb_w90, k_red @ B.T, mode='berry')

    nfrozen = np.sum(E_dft < FROZ_MAX, axis=1)
    err = max(np.abs(eps[k, :n] - E_dft[k, :n]).max() for k, n in enumerate(nfrozen))
    assert err < 3e-5


def test_wannier_tb_dirac_point(tb_w90):
    """
    The K of GridTB, (2 b1 + b2)/3, is the Dirac point of the real lattice too (a1, a2 at 60 degrees,
    oriented differently from M0): the p_z bands 4 and 5 are degenerate at E_D.
    """
    K = make_grid_tb(tb_w90, 3).K

    _, eps, _, _ = compute_velocity(tb_w90, K[None, :], mode='berry')

    assert abs(eps[0, 4] - eps[0, 3]) < 1e-6
    assert np.allclose(eps[0, 3:], E_D, rtol=0, atol=1e-5)


def test_wannier_tb_berry_connection_raw(tb_w90):
    """
    r_R is kept raw: A(K) = fourier(r_R) is not Hermitian (2.3e-3 Angstrom in x, 4.9e-3 in y, EM.md
    section 2), and hermitize removes the defect.
    """
    K = make_grid_tb(tb_w90, 3).K

    A = fourier(tb_w90.r_R, tb_w90.R_cart, tb_w90.ndegen, K[None, :]) # (1, 3, nW, nW)
    defect = np.abs(A - dagger(A)).max(axis=(0, 2, 3)) # (3,), x y z
    assert np.allclose(defect[:2], [2.31e-3, 4.86e-3], rtol=0.02, atol=0)

    A_h = hermitize(A)
    assert np.abs(A_h - dagger(A_h)).max() < 1e-15


def test_wannier_tb_ring(tb_w90):
    """
    The M0 defaults t = 2.7 eV and a_cc = 1.42 Angstrom suffice for `ring` on the real data: at 2.33 eV
    and mu = E_D the bracket [0, 2 q0] holds (q/q0 = 0.95-1.29; real hbar v_F = 5.47 eV Angstrom
    against 5.75 for the defaults), and every point is on shell between the p_z bands 4 and 5.
    """
    K = make_grid_tb(tb_w90, 3).K

    k, _, q, q0 = ring(tb_w90, K, HW, mu=E_D)
    _, eps, _, _ = compute_velocity(tb_w90, k, mode='berry')

    assert np.all(np.sum(eps < E_D, axis=1) == 4) # c = band 4, v = band 3 (0-based) everywhere
    assert np.allclose(eps[:, 4] - eps[:, 3], HW, rtol=0, atol=1e-9)
    assert 0.9 < q.min()/q0 and q.max()/q0 < 1.35


def random_k(tb, n, seed=0):
    """
    n Cartesian k points drawn uniformly in the in-plane Brillouin zone of tb, (n, 3), 1/Angstrom.
    """
    k_red = np.random.default_rng(seed).random((n, 3)) * [1, 1, 0]
    return k_red @ reciprocal(tb.lattice).T


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_centres_only_graphene(tb):
    """
    The M0 model has only centres: centres_only changes nothing, in both gauges.
    """
    assert np.array_equal(centres_only(tb).r_R, tb.r_R)


def test_centres_only_structure(tb_w90):
    """
    Only the diagonal of r(0) survives, unchanged; tb itself is not modified (tb_w90 is shared by
    the whole session) and the other fields are the same objects.
    """
    r_before = tb_w90.r_R.copy()
    tb_c = centres_only(tb_w90)
    i0 = tb_w90.index[(0, 0, 0)]
    nW = tb_w90.r_R.shape[-1]

    assert tb_c.r_R.shape == tb_w90.r_R.shape
    assert np.count_nonzero(tb_c.r_R) == 3 * nW
    diag = np.diagonal(tb_c.r_R[i0], axis1=-2, axis2=-1)
    assert np.array_equal(diag, np.diagonal(tb_w90.r_R[i0], axis1=-2, axis2=-1))

    assert np.array_equal(tb_w90.r_R, r_before), 'centres_only modified its input'
    assert tb_c.H_R is tb_w90.H_R and tb_c.R_cart is tb_w90.R_cart and tb_c.ndegen is tb_w90.ndegen


def test_centres_only_berry_connection(tb_w90):
    """
    A(k) = diag(tau_n) at every k: constant, and Hermitian without hermitize.
    """
    tb_c = centres_only(tb_w90)
    i0 = tb_w90.index[(0, 0, 0)]
    k = random_k(tb_w90, 200)

    A = fourier(tb_c.r_R, tb_c.R_cart, tb_c.ndegen, k) # (Nk, 3, nW, nW)

    assert np.allclose(A, tb_c.r_R[i0][None], rtol=0, atol=1e-14)
    assert np.abs(A - dagger(A)).max() < 1e-14


def test_centres_only_diagonal_unchanged(tb_w90):
    """
    The Berry term vanishes on the band diagonal: eps and hbar v_nn = d eps_n / dk are the same with
    the full r and with the centres only.
    """
    k = random_k(tb_w90, 200)

    _, eps, _, hv = compute_velocity(tb_w90, k, mode='berry')
    _, eps_c, _, hv_c = compute_velocity(centres_only(tb_w90), k, mode='berry')

    assert np.allclose(eps_c, eps, rtol=0, atol=1e-12)
    diag = lambda X: np.diagonal(X, axis1=-2, axis2=-1)
    assert np.allclose(diag(hv_c), diag(hv), rtol=0, atol=1e-12)


def velocity_atomic_gauge(tb, k):
    """
    hbar v from dH_at/dk alone, with H in the atomic gauge,
    H_at,mn(k) = sum_R e^{ik.(R + tau_n - tau_m)} H_mn(R) / ndegen(R), tau from the diagonal of r(0):
    no Berry term, no `fourier`, no commutator. Returns hv (Nk, 3, nW, nW), eV*Angstrom.
    """
    tau = np.diagonal(tb.r_R[tb.index[(0, 0, 0)]].real, axis1=-2, axis2=-1).T # (nW, 3)
    d = tb.R_cart[:, None, None, :] + tau[None, None, :, :] - tau[None, :, None, :] # (nR, m, n, 3)
    phase = np.exp(1j * np.einsum("kc,rmnc->krmn", k, d)) / tb.ndegen[None, :, None, None]

    H_at = np.einsum("krmn,rmn->kmn", phase, tb.H_R)
    dH_at = np.einsum("krmn,rmnc,rmn->kcmn", phase, 1j*d, tb.H_R)
    _, V = np.linalg.eigh(H_at)
    return dagger(V)[:, None] @ dH_at @ V[:, None]


def test_centres_only_atomic_gauge(tb_w90):
    """
    With the centres only, the Berry term turns dH/dk into its atomic-gauge form (EM.md,
    implementation notes): |hbar v_mn|^2 equals that of velocity_atomic_gauge (measured 2e-12
    eV^2 Angstrom^2), while the full r differs by up to ~30 eV^2 Angstrom^2. Moduli only: the phases
    of the two sets of eigenvectors are unrelated.
    """
    k = random_k(tb_w90, 200)

    _, _, _, hv_c = compute_velocity(centres_only(tb_w90), k, mode='berry')
    _, _, _, hv = compute_velocity(tb_w90, k, mode='berry')
    hv_at = velocity_atomic_gauge(tb_w90, k)

    assert np.allclose(np.abs(hv_c)**2, np.abs(hv_at)**2, rtol=0, atol=1e-9)
    assert np.abs(np.abs(hv)**2 - np.abs(hv_at)**2).max() > 1
