"""
disorder_average.py
    Disorder-averaged Green's function of a crystal with a dilute concentration of identical point defects, in the T-matrix
    approximation of Kaasbjerg, PRB 101, 045433 (2020), Eqs. (22)-(24), (28), (34) (R8, 2026-09-28).

    The single-defect local t-matrix t(eps) on the cluster R_local (local_tmatrix.local_t / pole_criterion.local_t_cache) is
    brought to k space, where the average over random defect positions restores translation symmetry:

        Tbar_k(eps)[w, w'] = sum_{L, L'} exp(-2 pi i k.R_L) t_{(L,w),(L',w')}(eps) exp(+2 pi i k.R_L')         (Wannier gauge)
                           = sum_D exp(-2 pi i k.D) tau_{ww'}(D; eps),   tau(D) = sum_{R_L - R_L' = D} t_{(L,.),(L',.)}
        Sigma_k(eps) = c_cell Tbar_k(eps)            (c_cell = defects per unit cell; Tbar = N T is intensive, Eq. (24))
        G_k(eps) = [(eps + i eta) 1 - H_k - Sigma_k(eps)]^{-1}                                                  (Eq. (22))
        rho(eps) = -(1/(pi N_k)) sum_k Im Tr G_k(eps)        states / eV / unit cell / spin                      (Eq. (34))
        A_k(eps) = -2 Im Tr G_k(eps)                                                                            (Eq. (28))

    Tbar_k is the Bloch-sum matrix element <<wk|T|w'k>> with Bloch sums normalized to one unit cell, i.e. the same normalization
    as the unit-cell M (bloch_norm = 'unit_cell') and as scattering_rate: <nk|t|nk> = U^dagger Tbar_k U. H_k is the Wannier-gauge
    Hamiltonian (Hwr_to_Hwk) with the same "R only" phase convention, so Tr G_k is gauge invariant.

    Nothing here recomputes g0, t or the cluster Green's function: those come from local_tmatrix (local_green_batch, local_t,
    cluster_ldos) and pole_criterion (local_t_cache).
"""
import numpy as np
from scipy.signal import find_peaks

from electron_defect_interaction.defects.many_body.local_tmatrix import _diff_table


def tbar_reduce(t, R_local, nw, wfs=None):
    """
    Q1. Reduce the cluster t-matrix to the distinct cell differences D = R_L - R_L':

        tau_{ab}(D; eps) = sum_{L, L' : R_L - R_L' = D} t_{(L, wfs[a]), (L', wfs[b])}(eps)

    so that the Wannier-gauge k-diagonal element is a plain lattice Fourier sum, Tbar_k = sum_D exp(-2 pi i k.D) tau(D) (tbar_k).
    Exact rearrangement of the double sum over the cluster (no approximation).

    Inputs:
        t: (nE, nL*nw, nL*nw) or (nL*nw, nL*nw) complex, t-matrix on the cluster, flat index L*nw + w (local_t / local_t_cache).
        R_local: (nL, 3) ints, cluster cells (reduced coordinates), same order as the flat index of t.
        nw: int, number of Wannier functions per cell in the layout of t.
        wfs: sequence of Wannier-function indices to keep (e.g. (3, 4) for the pi block of the 5-WF layout); None keeps all nw.
    Returns:
        Du: (nD, 3) ints, distinct differences D (local_tmatrix._diff_table order).
        tau: (nE, nD, nb, nb) complex (nb = len(wfs)); (nD, nb, nb) if t was a single matrix.
    """
    t = np.asarray(t)
    single = t.ndim == 2
    if single:
        t = t[None]
    R_local = np.asarray(R_local, int)
    nL = len(R_local)
    nE = t.shape[0]
    assert t.shape[1:] == (nL * nw, nL * nw), f"t {t.shape[1:]} vs nL*nw = {nL * nw}"
    wfs = np.arange(nw) if wfs is None else np.asarray(wfs, int)
    nb = len(wfs)
    Du, inv = _diff_table(R_local)                                     # inv[L, L'] = index of R_L - R_L' in Du
    nD = len(Du)
    # t[e, L, a, L', b] on the kept orbitals -> (L, L') pairs as rows
    t5 = t.reshape(nE, nL, nw, nL, nw)[:, :, wfs][:, :, :, :, wfs]      # (nE, nL, nb, nL, nb)
    pairs = np.transpose(t5, (1, 3, 0, 2, 4)).reshape(nL * nL, nE * nb * nb)
    # one-hot (nD, nL^2) selection matrix: tau = S @ pairs sums the pairs sharing the same difference D
    S = np.zeros((nD, nL * nL))
    S[inv.ravel(), np.arange(nL * nL)] = 1.0
    tau = (S @ pairs).reshape(nD, nE, nb, nb).transpose(1, 0, 2, 3)
    tau = np.ascontiguousarray(tau)
    return (Du, tau[0]) if single else (Du, tau)


def tbar_k(tau, Du, k, U=None, k_chunk=65536):
    """
    Q2. k-diagonal matrix element of the single-defect T matrix, intensive (unit-cell normalized Bloch sums):

        Tbar^{(W)}_k(eps) = sum_D exp(-2 pi i k.D) tau(D; eps)                       (Wannier gauge, sublattice basis for p_z)
        Tbar^{nm}_k(eps) = [U(k)^dagger Tbar^{(W)}_k(eps) U(k)]_{nm}                  (band basis, off-diagonal n != m included)

    With U the eigenvectors of H_W(k) (Hwr_to_Hwk), Tbar^{nn}_k(eps_nk) = <nk|t(eps_nk)|nk> of scattering_rate and the on-shell
    rate is Gamma_nk = -2 Im Tbar^{nn}_k(eps_nk).

    Inputs:
        tau: (nE, nD, nb, nb) or (nD, nb, nb) complex, from tbar_reduce.
        Du: (nD, 3) ints, differences D matching tau.
        k: (nk, 3) floats, k-points in reduced coordinates (any list: grid or path).
        U: optional (nk, nb, nb) complex, columns = band eigenvectors in the Wannier gauge of the same nb orbitals.
        k_chunk: k-points per BLAS call (memory: k_chunk x nD phases).
    Returns:
        Tbar: (nE, nk, nb, nb) complex, or (nk, nb, nb) if tau was a single energy.
    """
    tau = np.asarray(tau)
    single = tau.ndim == 3
    if single:
        tau = tau[None]
    nE, nD, nb, _ = tau.shape
    k = np.asarray(k, float)
    Du = np.asarray(Du, float)
    out = np.empty((nE, len(k), nb, nb), dtype=complex)
    flat = tau.transpose(1, 0, 2, 3).reshape(nD, nE * nb * nb)          # (nD, nE*nb*nb)
    for s in range(0, len(k), k_chunk):
        ph = np.exp(-2j * np.pi * (k[s:s + k_chunk] @ Du.T))            # (nk_c, nD)
        out[:, s:s + k_chunk] = (ph @ flat).reshape(-1, nE, nb, nb).transpose(1, 0, 2, 3)
    if U is not None:
        U = np.asarray(U)
        out = np.einsum("kwn,ekwv,kvm->eknm", U.conj(), out, U, optimize=True)
    return out[0] if single else out


def green_k(Hk, Sigma, e, eta):
    """
    Q3. Dyson Green's function at one energy, general batched inversion (np.linalg.inv, no closed form):

        G_k = [(e + i eta) 1 - H_k - Sigma_k]^{-1}

    Inputs:
        Hk: (nk, nb, nb) complex, Hamiltonian at each k (Wannier gauge or band basis, same basis as Sigma).
        Sigma: (nk, nb, nb) complex or None (pristine: G_k = g0_k).
        e: float, energy (eV); eta: float, broadening (eV).
    Returns:
        G: (nk, nb, nb) complex.
    """
    nb = Hk.shape[-1]
    A = (e + 1j * eta) * np.eye(nb)[None] - Hk
    if Sigma is not None:
        A = A - Sigma
    return np.linalg.inv(A)


def dos_average(Hk, tau, Du, k, c_cell, egrid, eta, linear=False, e_chunk=16, k_chunk=32768):
    """
    Q4. Disorder-averaged DOS in the T-matrix approximation (Kaasbjerg Eqs. (22), (24), (34)):

        Sigma_k(eps) = c_cell Tbar_k(eps),   G_k = [(eps + i eta) - H_k - Sigma_k]^{-1},
        rho(eps)  = -(1/(pi N_k)) sum_k Im Tr G_k(eps),     rho0 = same with Sigma = 0,
    and, if linear, the first-order term in c_cell:
        drho_lin(eps) = -(1/(pi N_k)) sum_k Im Tr[g0_k Tbar_k g0_k],   g0_k = G_k(Sigma = 0),
    which equals (1/pi) Im Tr[t d g0/d eps] (Lloyd, local_green_batch(deriv=True)) when k is the internal grid of g0 and eta is
    the same. Units: states / eV / unit cell / spin, for the nb orbitals kept in tau.

    Inputs:
        Hk: (nk, nb, nb) complex, Wannier-gauge Hamiltonian on the k list (same orbitals and order as tau).
        tau, Du: from tbar_reduce, tau (nE, nD, nb, nb) on egrid.
        k: (nk, 3) floats, output k-grid (reduced coordinates), weights 1/N_k.
        c_cell: float or 1-D array of concentrations (defects per unit cell).
        egrid: (nE,) energies of tau (eV); eta: broadening of G_k (eV).
        linear: also return drho_lin.
        e_chunk, k_chunk: energies and k-points per block (memory ~ e_chunk x k_chunk x nb^2 complex).
    Returns:
        dict(rho=(nE,) or (nc, nE), rho0=(nE,), [drho_lin=(nE,)], c_cell=array, eta=eta, nk=N_k).
    """
    egrid = np.asarray(egrid, float)
    tau = np.asarray(tau)
    nE = len(egrid)
    assert tau.shape[0] == nE, f"tau has {tau.shape[0]} energies, egrid {nE}"
    cs = np.atleast_1d(np.asarray(c_cell, float))
    k = np.asarray(k, float)
    nk = len(k)
    rho = np.zeros((len(cs), nE))
    rho0 = np.zeros(nE)
    lin = np.zeros(nE) if linear else None
    for ks in range(0, nk, k_chunk):
        kc, Hc = k[ks:ks + k_chunk], Hk[ks:ks + k_chunk]
        for es in range(0, nE, e_chunk):
            T = tbar_k(tau[es:es + e_chunk], Du, kc, k_chunk=k_chunk)   # (ec, nkc, nb, nb)
            for j in range(T.shape[0]):
                e = egrid[es + j]
                g0 = green_k(Hc, None, e, eta)
                rho0[es + j] += np.trace(g0, axis1=1, axis2=2).imag.sum()
                for ic, c in enumerate(cs):
                    G = green_k(Hc, c * T[j], e, eta)
                    rho[ic, es + j] += np.trace(G, axis1=1, axis2=2).imag.sum()
                if linear:
                    lin[es + j] += np.einsum("kab,kbc,kca->", g0, T[j], g0).imag
    f = -1.0 / (np.pi * nk)
    out = dict(rho=f * (rho[0] if np.ndim(c_cell) == 0 else rho), rho0=f * rho0, c_cell=cs, eta=eta, nk=nk)
    if linear:
        out["drho_lin"] = f * lin
    return out


def spectral_path(Hk, Tbar, c_cell, egrid, eta):
    """
    Q5a. Spectral function along a k list (Kaasbjerg Eq. (28), summed over bands):

        A_k(eps) = -2 Im Tr G_k(eps),   G_k = [(eps + i eta) - H_k - c_cell Tbar_k(eps)]^{-1},
    normalized as in the article: integral of A_k d eps / (2 pi) = number of bands nb.

    Inputs:
        Hk: (nk, nb, nb) complex (same basis as Tbar); Tbar: (nE, nk, nb, nb) from tbar_k on egrid; c_cell: float;
        egrid: (nE,) eV; eta: eV.
    Returns:
        A: (nE, nk) floats, 1/eV.
    """
    egrid = np.asarray(egrid, float)
    A = np.empty((len(egrid), Hk.shape[0]))
    for j, e in enumerate(egrid):
        G = green_k(Hk, c_cell * Tbar[j], e, eta)
        A[j] = -2.0 * np.trace(G, axis1=1, axis2=2).imag
    return A


def spectral_maxima(A, egrid, prominence=0.0, refine=True):
    """
    Q5b. Maxima of A_k(eps) along the energy axis at fixed k (the "white dots" of Kaasbjerg Fig. 14).

    Local maxima by scipy.signal.find_peaks (with the given prominence, same units as A); if refine, the vertex of the parabola
    through the three grid points around an interior maximum (exact for a parabolic peak), height = parabola value there.

    Inputs:
        A: (nE, nk) floats; egrid: (nE,) uniform or not; prominence: float; refine: bool.
    Returns:
        list of nk arrays of shape (npeaks, 2): columns (eps_peak, A_peak), sorted by energy.
    """
    egrid = np.asarray(egrid, float)
    out = []
    for ik in range(A.shape[1]):
        a = A[:, ik]
        idx, _ = find_peaks(a, prominence=prominence)
        pk = []
        for i in idx:
            x, y = egrid[i], a[i]
            if refine and 0 < i < len(a) - 1:
                x0, x1, x2 = egrid[i - 1:i + 2]
                y0, y1, y2 = a[i - 1:i + 2]
                # parabola y = p x^2 + q x + r through the three points; vertex at -q / (2 p)
                p, q, r = np.polyfit([x0, x1, x2], [y0, y1, y2], 2)
                if p < 0:
                    x = -q / (2 * p)
                    y = r - q * q / (4 * p)
            pk.append((x, y))
        out.append(np.array(pk, float).reshape(-1, 2))
    return out


def sigma_eff(eps_k, Sigma, e, eta=0.0):
    """
    Q7. Effective self-energy of each band of a 2 x 2 subspace (Kaasbjerg Eqs. (44)-(45)): eliminating the other band n' from
    G_k = [(e + i eta) - diag(eps_k) - Sigma]^{-1} (Schur complement) gives

        G^{nn}_k(e) = 1 / (e + i eta - eps_nk - Sigma^eff_nk(e)),
        Sigma^eff_nk(e) = Sigma_nn(e) + Sigma_nn'(e) Sigma_n'n(e) / (e + i eta - eps_n'k - Sigma_n'n'(e)),   n' != n.

    The second term is the defect-induced coupling between the two bands (Kaasbjerg: responsible for the band-gap opening at K).
    Quasiparticle solutions: e - eps_nk - Re Sigma^eff_nk(e) = 0 (pole_criterion.sign_changes on a grid).

    Inputs:
        eps_k: (2,) band energies at k (band basis of Sigma).
        Sigma: (nE, 2, 2) complex self-energy in the band basis (e.g. c_cell * tbar_k(..., U=U)[:, ik]).
        e: (nE,) energies; eta: broadening added to e (0 in Eq. (45)).
    Returns:
        Seff: (nE, 2) complex, column n = Sigma^eff of band n.
    """
    e = np.asarray(e, float); S = np.asarray(Sigma); eps_k = np.asarray(eps_k, float)
    z = e + 1j * eta
    out = np.empty((len(e), 2), dtype=complex)
    for n, m in ((0, 1), (1, 0)):
        out[:, n] = S[:, n, n] + S[:, n, m] * S[:, m, n] / (z - eps_k[m] - S[:, m, m])
    return out


def dirac_g0bar(eps, Lam, hbar_vF, A_cell, g_v=2):
    """
    Q8a. k-summed Green's function of the Dirac model per sublattice site (Kaasbjerg Eq. (48)):

        G0bar(eps) = A_cell (rho0bar / 2) [ eps ln| eps^2 / (eps^2 - Lam^2) | - i pi |eps| theta(Lam - |eps|) ],
        rho0bar = g_v / (2 pi (hbar vF)^2),

    i.e. the Hilbert transform of the site DOS D(eps) = A_cell (rho0bar/2) |eps| on [-Lam, Lam].
    Units: eps, Lam (eV), hbar_vF (eV Angstrom), A_cell (Angstrom^2) -> 1/eV.
    """
    eps = np.asarray(eps, float)
    C = A_cell * g_v / (2 * np.pi * hbar_vF ** 2) / 2
    with np.errstate(divide="ignore"):
        re = C * eps * np.log(np.abs(eps ** 2 / (eps ** 2 - Lam ** 2)))
    re = np.where(eps == 0.0, 0.0, re)
    im = -np.pi * C * np.abs(eps) * (np.abs(eps) < Lam)
    return re + 1j * im


def dirac_t0(eps, V0, Lam, hbar_vF, A_cell, g_v=2):
    """Q8b. T0(eps) = V0 / (1 - V0 G0bar(eps)) (Kaasbjerg Eq. (47))."""
    return V0 / (1.0 - V0 * dirac_g0bar(eps, Lam, hbar_vF, A_cell, g_v))


def dirac_pole(V0, Lam, hbar_vF, A_cell, g_v=2):
    """
    Q8c. Pole of T0 closest to the Dirac point: root of 1/V0 = Re G0bar(eps) with the smallest |eps|, on the side sign(eps) = -sign(V0)
    (below the Dirac point for a repulsive V0 > 0). Re G0bar vanishes at 0 and at |eps| = Lam/sqrt(2) and has one extremum in between,
    so the root is bracketed between 0 and that extremum (log-spaced search, then brentq). Returns NaN if |1/V0| exceeds the extremum.
    """
    from scipy.optimize import brentq
    s = -np.sign(V0)
    f = lambda x: dirac_g0bar(s * x, Lam, hbar_vF, A_cell, g_v).real - 1.0 / V0
    xs = np.geomspace(1e-12 * Lam, Lam / np.sqrt(2) * (1 - 1e-12), 4000)
    fx = f(xs)
    i = np.where(np.sign(fx[:-1]) != np.sign(fx[1:]))[0]
    if len(i) == 0:
        return float("nan")
    return float(s * brentq(f, xs[i[0]], xs[i[0] + 1], xtol=1e-15, rtol=1e-14))


def dirac_lambda_for_pole(eps_p, V0, hbar_vF, A_cell, g_v=2):
    """
    Q8d. Cutoff Lam for which the Dirac-model pole sits at eps_p: from 1/V0 = C eps_p ln(eps_p^2 / (Lam^2 - eps_p^2)) (|eps_p| < Lam),
        Lam = |eps_p| sqrt(1 + exp(-1 / (V0 C eps_p))),   C = A_cell g_v / (4 pi (hbar vF)^2).
    """
    C = A_cell * g_v / (2 * np.pi * hbar_vF ** 2) / 2
    return float(abs(eps_p) * np.sqrt(1.0 + np.exp(-1.0 / (V0 * C * eps_p))))
