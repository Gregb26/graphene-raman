"""
Tests of the tight-binding model (tb_model): R-space Hermiticity, bonds, relabelling shift_B of the
M0 graphene model; construction of the real 27 x 27 model by make_wannier_tb, its centres-only
version centres_only (M1, fixtures tb_w90, eig_w90 and w90_ref of conftest.py), extract_block and
its p_z block pz_block (M3).
"""

import pytest
import numpy as np
from graphene_raman.electron_photon import (centres_only, compute_velocity, dagger,
                                                         fermi_velocity, fourier, hermitize,
                                                         make_graphene_tb, make_grid_tb, pz_block,
                                                         reciprocal, ring)
from graphene_raman.electron_photon.tb_model import extract_block

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


def test_wannier_tb_fields(tb_w90, w90_ref):
    """R_cart, index and minus consistent with R_int; H Hermitian in R space."""
    tb = tb_w90

    nR, nW = w90_ref.nR, w90_ref.nW
    assert tb.H_R.shape == (nR, nW, nW) and tb.r_R.shape == (nR, 3, nW, nW)
    assert np.allclose(tb.R_cart, tb.R_int @ tb.lattice.T, rtol=0, atol=1e-12)
    assert all(tb.index[tuple(R)] == iR for iR, R in enumerate(tb.R_int))
    assert np.array_equal(tb.R_int[tb.minus], -tb.R_int)
    assert np.allclose(tb.H_R, tb.H_R[tb.minus].swapaxes(-1, -2).conj(), rtol=0, atol=1e-12)


def test_wannier_tb_bands_match_eig(tb_w90, eig_w90, w90_ref):
    """Cartesian chain (k_red @ B.T, compute_velocity) reproduces the .eig in the frozen window."""
    k_red, E_dft = eig_w90
    B = reciprocal(tb_w90.lattice)

    _, eps, _, _ = compute_velocity(tb_w90, k_red @ B.T, mode='berry')

    nfrozen = np.sum(E_dft < w90_ref.froz_max, axis=1)
    err = max(np.abs(eps[k, :n] - E_dft[k, :n]).max() for k, n in enumerate(nfrozen))
    assert err < w90_ref.eig_tol


def test_wannier_tb_dirac_point(tb_w90, w90_ref):
    """The K of GridTB is the Dirac point of the real lattice: pi and pi* degenerate at E_D."""
    K = make_grid_tb(tb_w90, 3).K
    v, c = w90_ref.bands_pi

    _, eps, _, _ = compute_velocity(tb_w90, K[None, :], mode='berry')

    assert abs(eps[0, c] - eps[0, v]) < 1e-6
    assert np.allclose(eps[0, [v, c]], w90_ref.E_D, rtol=0, atol=1e-5)


def test_wannier_tb_berry_connection_raw(tb_w90, w90_ref):
    """A(K) from the raw r is not Hermitian (EM.md section 2); hermitize fixes it."""
    K = make_grid_tb(tb_w90, 3).K

    A = fourier(tb_w90.r_R, tb_w90.R_cart, tb_w90.ndegen, K[None, :]) # (1, 3, nW, nW)
    defect = np.abs(A - dagger(A)).max(axis=(0, 2, 3)) # (3,), x y z
    assert np.allclose(defect[:2], w90_ref.r_defect_K, rtol=0.02, atol=0)

    A_h = hermitize(A)
    assert np.abs(A_h - dagger(A_h)).max() < 1e-15


def test_wannier_tb_ring(tb_w90, w90_ref):
    """The default t, a_cc bracket the 2.33 eV ring on the real data; points on shell."""
    K = make_grid_tb(tb_w90, 3).K
    v, c = w90_ref.bands_pi
    lo, hi = w90_ref.ring_q_over_q0

    k, _, q, q0 = ring(tb_w90, K, HW, mu=w90_ref.E_D)
    _, eps, _, _ = compute_velocity(tb_w90, k, mode='berry')

    assert np.all(np.sum(eps < w90_ref.E_D, axis=1) == c) # occupied up to pi everywhere
    assert np.allclose(eps[:, c] - eps[:, v], HW, rtol=0, atol=1e-9)
    assert lo < q.min()/q0 and q.max()/q0 < hi


def random_k(tb, n, seed=0):
    """n random Cartesian k points in the in-plane Brillouin zone of tb, (n, 3), 1/Angstrom."""
    k_red = np.random.default_rng(seed).random((n, 3)) * [1, 1, 0]
    return k_red @ reciprocal(tb.lattice).T


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_centres_only_graphene(tb):
    """The M0 model has only centres: centres_only changes nothing, in both gauges."""
    assert np.array_equal(centres_only(tb).r_R, tb.r_R)


def test_centres_only_structure(tb_w90):
    """Only the diagonal of r(0) survives; tb is not modified (session fixture); other fields shared."""
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
    """A(k) = diag(tau_n): constant in k, Hermitian without hermitize."""
    tb_c = centres_only(tb_w90)
    i0 = tb_w90.index[(0, 0, 0)]
    k = random_k(tb_w90, 200)

    A = fourier(tb_c.r_R, tb_c.R_cart, tb_c.ndegen, k) # (Nk, 3, nW, nW)

    assert np.allclose(A, tb_c.r_R[i0][None], rtol=0, atol=1e-14)
    assert np.abs(A - dagger(A)).max() < 1e-14


def test_centres_only_diagonal_unchanged(tb_w90):
    """The Berry term vanishes on the band diagonal: eps and hbar v_nn unchanged."""
    k = random_k(tb_w90, 200)

    _, eps, _, hv = compute_velocity(tb_w90, k, mode='berry')
    _, eps_c, _, hv_c = compute_velocity(centres_only(tb_w90), k, mode='berry')

    assert np.allclose(eps_c, eps, rtol=0, atol=1e-12)
    diag = lambda X: np.diagonal(X, axis1=-2, axis2=-1)
    assert np.allclose(diag(hv_c), diag(hv), rtol=0, atol=1e-12)


def velocity_atomic_gauge(tb, k):
    """
    hbar v = V^dagger dH_at/dk V with H_at,mn(k) = sum_R e^{ik.(R + tau_n - tau_m)} H_mn(R)/ndegen(R):
    no Berry term, no `fourier`. Returns hv (Nk, 3, nW, nW), eV*Angstrom.
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
    Centres only = atomic gauge: |hbar v_mn|^2 agree to rounding, while the full r clearly differs.
    Moduli only, the eigenvector phases are unrelated.
    """
    k = random_k(tb_w90, 200)

    _, _, _, hv_c = compute_velocity(centres_only(tb_w90), k, mode='berry')
    _, _, _, hv = compute_velocity(tb_w90, k, mode='berry')
    hv_at = velocity_atomic_gauge(tb_w90, k)

    assert np.allclose(np.abs(hv_c)**2, np.abs(hv_at)**2, rtol=0, atol=1e-9)
    assert np.abs(np.abs(hv)**2 - np.abs(hv_at)**2).max() > 1


def test_extract_block():
    """rows x columns block of the last two axes, leading axes kept (reference: np.ix_)."""
    X = np.random.default_rng(1).random((2, 3, 5, 5))
    rows, cols = [0, 2], [1, 3, 4]

    block = extract_block(X, rows, cols)

    assert block.shape == (2, 3, 2, 3)
    assert np.array_equal(block, X[(slice(None), slice(None)) + np.ix_(rows, cols)])


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_pz_block_graphene(tb):
    """The M0 model has only its two p_z: pz_block keeps it whole, in both gauges."""
    tb2 = pz_block(tb, [0, 1])

    assert np.array_equal(tb2.H_R, tb.H_R) and np.array_equal(tb2.r_R, tb.r_R)


def test_pz_block_structure(tb_w90, w90_ref):
    """pz x pz blocks of H and r (centres included), as copies; tb intact; R list and ndegen shared."""
    pz = w90_ref.pz
    H_before, r_before = tb_w90.H_R.copy(), tb_w90.r_R.copy()
    ix = np.ix_(pz, pz)

    tb2 = pz_block(tb_w90, pz)

    assert tb2.H_R.shape == (w90_ref.nR, 2, 2) and tb2.r_R.shape == (w90_ref.nR, 3, 2, 2)
    assert np.array_equal(tb2.H_R, tb_w90.H_R[(slice(None),) + ix])
    assert np.array_equal(tb2.r_R, tb_w90.r_R[(slice(None), slice(None)) + ix])
    assert not np.shares_memory(tb2.H_R, tb_w90.H_R) and not np.shares_memory(tb2.r_R, tb_w90.r_R)

    assert np.array_equal(tb_w90.H_R, H_before) and np.array_equal(tb_w90.r_R, r_before), 'pz_block modified its input'
    assert tb2.R_int is tb_w90.R_int and tb2.ndegen is tb_w90.ndegen and tb2.minus is tb_w90.minus


def test_pz_block_near_K(tb_w90, w90_ref):
    """
    Circle q = 0.1 around K: eps and |hbar v^{x,y}_cv|^2 of pi, pi* are those of the 5 x 5
    (measured 5e-15 eV and 1.5e-13 eV^2 Angstrom^2). In-plane light does not mix sigma and pi.
    """
    K = make_grid_tb(tb_w90, 3).K
    v, c = w90_ref.bands_pi
    theta = np.linspace(0, 2*np.pi, 360, endpoint=False)
    k = K + 0.1 * np.stack([np.cos(theta), np.sin(theta), np.zeros_like(theta)], axis=1)

    _, eps5, _, hv5 = compute_velocity(tb_w90, k, mode='berry')
    _, eps2, _, hv2 = compute_velocity(pz_block(tb_w90, w90_ref.pz), k, mode='berry')

    assert np.allclose(eps2, eps5[:, [v, c]], rtol=0, atol=1e-12)
    assert np.allclose(np.abs(hv2[:, :2, 1, 0])**2, np.abs(hv5[:, :2, c, v])**2, rtol=0, atol=1e-10)


def test_pz_block_bands_in_5x5(tb_w90, w90_ref):
    """
    Whole zone: each band of the 2 x 2 is a band of the 5 x 5, matched by energy. Pairing by index
    fails (pi crosses sigma bands far from K), hence the second assert.
    """
    v, c = w90_ref.bands_pi
    k = random_k(tb_w90, 500)

    _, eps5, _, _ = compute_velocity(tb_w90, k, mode='no_berry')
    _, eps2, _, _ = compute_velocity(pz_block(tb_w90, w90_ref.pz), k, mode='no_berry')

    nearest = np.abs(eps2[:, :, None] - eps5[:, None, :]).min(axis=2) # (Nk, 2), eV
    assert nearest.max() < 1e-12
    assert np.abs(eps2 - eps5[:, [v, c]]).max() > 1


def test_pz_block_fermi_velocity_and_ring(tb_w90, w90_ref):
    """Same hbar v_F as the 5 x 5, and same <|hbar v_cv|^2> on the 2.33 eV ring (31.443 eV^2 Angstrom^2)."""
    K = make_grid_tb(tb_w90, 3).K
    tb2 = pz_block(tb_w90, w90_ref.pz)

    assert np.isclose(fermi_velocity(tb2, K, w90_ref.E_D)['inter_avg'], w90_ref.hv_F, rtol=1e-5, atol=0)

    means = []
    for model, (v, c) in [(tb_w90, w90_ref.bands_pi), (tb2, (0, 1))]:
        k, _, _, _ = ring(model, K, HW, mu=w90_ref.E_D)
        _, _, _, hv = compute_velocity(model, k, mode='berry')
        means.append(np.mean(np.abs(hv[:, 0, c, v])**2 + np.abs(hv[:, 1, c, v])**2))
    assert np.isclose(means[1], means[0], rtol=1e-10, atol=0)
