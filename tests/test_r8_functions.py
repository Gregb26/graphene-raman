"""Minimal tests of the R8 functions (2026-09-28): disorder_average.tbar_reduce (Q1), tbar_k (Q2), green_k (Q3), dos_average (Q4),
spectral_path / spectral_maxima (Q5) and utils.lattice.build_k_path (Q6).

Reference = direct resolvent of a finite periodic lattice: nearest-neighbour p_z graphene (electron_photon.make_graphene_tb) on 6 x 6 cells,
a random Hermitian defect potential on a 5-cell cluster, T = V [1 - G0 V]^-1 by inversion of the 72 x 72 real-space matrices. On the
6 x 6 grid, local_green_batch gives the exact cluster block of G0 of that finite lattice, so every identity holds to round-off.
Run: PYTHONPATH=src pytest tests/test_r8_functions.py"""
import numpy as np
import pytest

from graphene_raman.defects.many_body import cluster_tmatrix as ct
from graphene_raman.defects.many_body import disorder_average as da
from graphene_raman.defects.many_body.pole_criterion import cluster_t_cache
from graphene_raman.electron_photon import make_graphene_tb
from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk
from graphene_raman.utils.lattice import build_k_path

N = 6                                                        # finite lattice N x N cells (periodic)
ETA = 0.1                                                    # eV, broadening of g0 (and of G_k where the identity needs it)
EGRID = np.array([-2.3, -0.7, -0.05, 0.4, 1.9])
R_LOC = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])


def _pz_model():
    tb = make_graphene_tb(t=2.7)
    return tb.H_R, tb.R_int, tb.ndegen                       # p_z(A), p_z(B)


def _cell(R):
    return (int(R[0]) % N) * N + (int(R[1]) % N)


def _real_space_H(Hpi, Rw, nd):
    nw = Hpi.shape[1]
    H = np.zeros((N * N * nw, N * N * nw), complex)
    for i in range(N):
        for j in range(N):
            L = _cell((i, j))
            for iR, R in enumerate(Rw):
                M = _cell((i + R[0], j + R[1]))
                H[L * nw:(L + 1) * nw, M * nw:(M + 1) * nw] += Hpi[iR] / nd[iR]     # H(R)[m, n] = <m 0|H|n R>
    return H


def _random_hermitian(n, rng, scale=1.0):
    X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return scale * 0.5 * (X + X.conj().T)


@pytest.fixture(scope="module")
def setup():
    rng = np.random.default_rng(8)
    Hpi, Rw, nd = _pz_model()
    nw = 2
    V = _random_hermitian(len(R_LOC) * nw, rng, scale=1.5)                  # defect potential on the cluster, index L*nw + w
    Hfull = _real_space_H(Hpi, Rw, nd)
    idx = np.array([_cell(R) * nw + w for R in R_LOC for w in range(nw)])
    Vfull = np.zeros_like(Hfull); Vfull[np.ix_(idx, idx)] = V
    n = len(Hfull)
    Tfull = []
    for e in EGRID:
        G0 = np.linalg.inv((e + 1j * ETA) * np.eye(n) - Hfull)
        Tfull.append(Vfull @ np.linalg.inv(np.eye(n) - G0 @ Vfull))
    k_int = ct.mp_grid(N, N, 1)
    Hk = Hwr_to_Hwk(Hpi, Rw, k_int, ndegen=nd)[0]
    g0 = ct.cluster_green_batch(Hk, k_int, R_LOC, EGRID, ETA)
    t = cluster_t_cache(V, g0)
    return dict(Hpi=Hpi, Rw=Rw, nd=nd, V=V, idx=idx, Tfull=np.array(Tfull), k_int=k_int, Hk=Hk, g0=g0, t=t, nw=nw)


def _tbar_direct(S, k):
    """Tbar_k[w, w'] = sum_{L, L'} e^{-2 pi i k.R_L} T_{(L,w),(L',w')} e^{+2 pi i k.R_L'} from the real-space T (unwrapped labels)."""
    nw = S["nw"]
    T = S["Tfull"][:, S["idx"]][:, :, S["idx"]].reshape(len(EGRID), len(R_LOC), nw, len(R_LOC), nw)
    ph = np.exp(2j * np.pi * (np.asarray(k) @ R_LOC.T))                     # (nk, nL)
    return np.einsum("kL,eLaMb,kM->ekab", ph.conj(), T, ph)


def test_tbar_matches_direct_resolvent(setup):
    S = setup
    assert np.allclose(S["t"], S["Tfull"][:, S["idx"]][:, :, S["idx"]], rtol=0, atol=1e-12)       # local t = cluster block of T
    rng = np.random.default_rng(1)
    k = np.vstack([rng.random((5, 3)) * [1, 1, 0], [[2 / 3, 1 / 3, 0]], S["k_int"][:7]])
    Du, tau = da.tbar_reduce(S["t"], R_LOC, S["nw"])
    Tb = da.tbar_k(tau, Du, k, k_chunk=4)
    ref = _tbar_direct(S, k)
    assert Tb.shape == (len(EGRID), len(k), 2, 2)
    assert np.abs(Tb - ref).max() < 1e-12 * np.abs(ref).max()
    # single energy, and band basis U^dagger Tbar U
    Du1, tau1 = da.tbar_reduce(S["t"][2], R_LOC, S["nw"])
    assert np.allclose(da.tbar_k(tau1, Du1, k), Tb[2], rtol=0, atol=1e-12)
    _, _, U = Hwr_to_Hwk(S["Hpi"], S["Rw"], k, ndegen=S["nd"])
    Tn = da.tbar_k(tau, Du, k, U=U)
    assert np.abs(Tn - np.einsum("kwn,ekwv,kvm->eknm", U.conj(), ref, U)).max() < 1e-12 * np.abs(ref).max()


def test_tbar_reduce_keeps_selected_wfs(setup):
    """5-WF layout of the production (sigma 0-2 decoupled at -100 eV, p_z = 3, 4), V on the p_z only: wfs=(3, 4) gives the pi-only result."""
    S = setup
    Hwr = np.zeros((len(S["Rw"]), 5, 5), complex)
    Hwr[:, 3:5, 3:5] = S["Hpi"]
    Hwr[int(np.flatnonzero((S["Rw"] == 0).all(axis=1))[0]), [0, 1, 2], [0, 1, 2]] = -100.0
    Rw, nd = S["Rw"], S["nd"]
    nL = len(R_LOC)
    V5 = np.zeros((nL * 5, nL * 5), complex)
    ip = np.array([L * 5 + w for L in range(nL) for w in (3, 4)])
    V5[np.ix_(ip, ip)] = S["V"]
    Hk5 = Hwr_to_Hwk(Hwr, Rw, S["k_int"], ndegen=nd)[0]
    t5 = cluster_t_cache(V5, ct.cluster_green_batch(Hk5, S["k_int"], R_LOC, EGRID, ETA))
    Du, tau5 = da.tbar_reduce(t5, R_LOC, 5, wfs=(3, 4))
    Du2, tau2 = da.tbar_reduce(S["t"], R_LOC, 2)
    assert np.array_equal(Du, Du2) and np.abs(tau5 - tau2).max() < 1e-12 * np.abs(tau2).max()


def test_green_and_dos_match_direct_inversion(setup):
    S = setup
    k, Hk = S["k_int"], S["Hk"]
    Du, tau = da.tbar_reduce(S["t"], R_LOC, S["nw"])
    cs, eta_G = np.array([0.0, 0.01, 0.05]), 0.07
    out = da.dos_average(Hk, tau, Du, k, cs, EGRID, eta_G, e_chunk=2, k_chunk=7)
    Tref = _tbar_direct(S, k)
    ref = np.array([[-np.trace(np.linalg.inv((e + 1j * eta_G) * np.eye(2) - Hk - c * Tref[j]), axis1=1, axis2=2).imag.sum() / (np.pi * len(k))
                     for j, e in enumerate(EGRID)] for c in cs])
    assert out["rho"].shape == (3, len(EGRID))
    assert np.abs(out["rho"] - ref).max() < 1e-12 * np.abs(ref).max()
    assert np.allclose(out["rho"][0], out["rho0"], rtol=0, atol=1e-14)
    # green_k itself, and a scalar concentration
    G = da.green_k(Hk, 0.01 * Tref[1], EGRID[1], eta_G)
    assert np.allclose(G, np.linalg.inv((EGRID[1] + 1j * eta_G) * np.eye(2) - Hk - 0.01 * Tref[1]), rtol=0, atol=1e-12)
    assert np.allclose(da.dos_average(Hk, tau, Du, k, 0.01, EGRID, eta_G)["rho"], out["rho"][1], rtol=0, atol=1e-14)


def test_linear_term_equals_lloyd(setup):
    """k = internal grid of g0 and eta_G = eta_t: -(1/(pi N_k)) sum_k Im Tr[g0_k Tbar_k g0_k] = (1/pi) Im Tr[t dg0/deps]."""
    S = setup
    Du, tau = da.tbar_reduce(S["t"], R_LOC, S["nw"])
    lin = da.dos_average(S["Hk"], tau, Du, S["k_int"], 0.0, EGRID, ETA, linear=True)["drho_lin"]
    g0p = ct.cluster_green_batch(S["Hk"], S["k_int"], R_LOC, EGRID, ETA, deriv=True)
    lloyd = np.array([np.trace(S["t"][j] @ g0p[j]).imag / np.pi for j in range(len(EGRID))])
    assert np.abs(lin - lloyd).max() < 1e-12 * np.abs(lloyd).max()


def test_spectral_function_and_maxima(setup):
    S = setup
    B = 2 * np.pi * np.linalg.inv(np.array([[np.sqrt(3) / 2, np.sqrt(3) / 2, 0], [-0.5, 0.5, 0], [0, 0, 5.0]])).T
    kp, _, _, _ = build_k_path([("G", [0, 0, 0]), ("K", [2 / 3, 1 / 3, 0]), ("M", [0.5, 0, 0])], 13, B)
    Du, tau = da.tbar_reduce(S["t"], R_LOC, S["nw"])
    Hp = Hwr_to_Hwk(S["Hpi"], S["Rw"], kp, ndegen=S["nd"])[0]
    A = da.spectral_path(Hp, da.tbar_k(tau, Du, kp), 0.02, EGRID, 0.05)
    Tref = _tbar_direct(S, kp)
    ref = np.array([-2 * np.trace(np.linalg.inv((e + 0.05j) * np.eye(2) - Hp - 0.02 * Tref[j]), axis1=1, axis2=2).imag
                    for j, e in enumerate(EGRID)])
    assert A.shape == (len(EGRID), len(kp)) and np.abs(A - ref).max() < 1e-12 * np.abs(ref).max()
    # maxima: exact vertex for parabolic peaks, grid maximum without refinement
    eg = np.linspace(-1.0, 1.0, 201)
    a = np.column_stack([np.maximum(0.0, 1.0 - (eg - 0.1234) ** 2 / 0.01), np.maximum(0.0, 2.0 - (eg + 0.4321) ** 2 / 0.04)])
    pk = da.spectral_maxima(a, eg)
    assert abs(pk[0][0, 0] - 0.1234) < 1e-12 and abs(pk[0][0, 1] - 1.0) < 1e-12
    assert abs(pk[1][0, 0] + 0.4321) < 1e-12 and abs(pk[1][0, 1] - 2.0) < 1e-12
    raw = da.spectral_maxima(a, eg, refine=False)
    assert abs(raw[0][0, 0] - 0.12) < 1e-12


def test_build_k_path_proportional():
    a = 2.4659
    A = a * np.array([[np.sqrt(3) / 2, np.sqrt(3) / 2, 0], [-0.5, 0.5, 0], [0, 0, 8.0 / a]])   # columns a1, a2 (60 degrees), a3
    B = 2 * np.pi * np.linalg.inv(A).T
    corners = [("G", [0, 0, 0]), ("K", [2 / 3, 1 / 3, 0]), ("M", [0.5, 0, 0])]
    k, labels, idx, s = build_k_path(corners, 601, B)
    assert len(k) == 601 and labels == ["G", "K", "M"] and idx == [0, 400, 600]                # |GK| : |KM| = 2 : 1
    for (_, p), j in zip(corners, idx):
        assert np.array_equal(k[j], np.asarray(p, float))
    assert np.all(np.diff(s) > 0) and abs(s[400] - 4 * np.pi / (3 * a)) < 1e-12                 # |GK| = 4 pi / (3 a)
    step = np.diff(s)
    assert np.ptp(step[:400]) < 1e-12 and np.ptp(step[400:]) < 1e-12 and abs(step[0] - step[-1]) < 1e-12
    k2, _, idx2, _ = build_k_path(corners, 7, B)
    assert len(k2) == 7 and idx2 == [0, 4, 6]
    with pytest.raises(ValueError):
        build_k_path(corners, 2, B)


def test_sigma_eff_is_the_schur_complement():
    """Q7: 1/(z - eps_n - Sigma^eff_n) equals the diagonal of the direct 2 x 2 inverse."""
    rng = np.random.default_rng(7)
    e = np.linspace(-1.0, 1.0, 9); eps_k = np.array([-0.3, 0.4]); eta = 0.02
    S = rng.standard_normal((9, 2, 2)) * 0.1 + 1j * rng.standard_normal((9, 2, 2)) * 0.1
    Se = da.sigma_eff(eps_k, S, e, eta)
    G = np.linalg.inv((e + 1j * eta)[:, None, None] * np.eye(2) - np.diag(eps_k)[None] - S)
    for n in (0, 1):
        assert np.allclose(1.0 / ((e + 1j * eta) - eps_k[n] - Se[:, n]), G[:, n, n], rtol=0, atol=1e-12)


def test_dirac_model():
    """Q8: G0bar = Hilbert transform of the Dirac site DOS (quadrature, 1e-9) ; pole and cutoff are inverse of each other (1e-12)."""
    from scipy.integrate import quad
    hv, Ac, Lam, gv = 5.459, 5.266, 7.0, 2
    C = Ac * gv / (4 * np.pi * hv ** 2)
    for x in (-3.0, -0.2, 0.15, 2.5):
        pv = quad(lambda y: C * abs(y), -Lam, Lam, weight="cauchy", wvar=x, epsabs=1e-14, epsrel=1e-13, limit=200)[0]
        g = da.dirac_g0bar(x, Lam, hv, Ac, gv)
        assert abs(g.real - (-pv)) < 1e-9 * max(1.0, abs(pv))            # quad cauchy integrates f(y)/(y - x)
        assert abs(g.imag + np.pi * C * abs(x)) < 1e-14
    for V0, L in ((28.8, 1e4), (29.7, 1e3), (-10.0, 1e4)):
        ep = da.dirac_pole(V0, L, hv, Ac)
        assert np.sign(ep) == -np.sign(V0)
        assert abs(da.dirac_g0bar(ep, L, hv, Ac).real - 1.0 / V0) < 1e-12 / abs(V0) * 1e2
        assert abs(da.dirac_lambda_for_pole(ep, V0, hv, Ac) / L - 1.0) < 1e-9
        t0 = da.dirac_t0(ep * np.array([0.5, 1.5]), V0, L, hv, Ac)
        assert np.all(np.isfinite(t0))
