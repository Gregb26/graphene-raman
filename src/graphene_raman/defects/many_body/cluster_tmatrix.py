"""
cluster_tmatrix.py
    Single-defect T-matrix on a finite cluster of cells around the defect (the set C of retained cells), in the Wannier
    basis, and the on-shell scattering rate it gives on any output k-grid without forming the dense Bloch T-matrix.

        M_cluster = M_W(R, R') for R, R' in the cluster, flat index L*nw + w            (cluster_potential)
        g0(eps)   = (1/N) sum_k e^{2 pi i k.(R - R')} [eps + i eta - H_W(k)]^{-1}       (cluster_green, cluster_green_batch)
        t(eps)    = M_cluster [1 - g0 M_cluster]^{-1}, exact if the cluster covers the support of M      (cluster_t)
        Gamma_nk  = -2 Im phi_nk^dag t(eps_nk) phi_nk,  phi_nk[(L, w)] = U_wn(k) e^{2 pi i k.R_L}  (per defect)

    R_cut is a norm on the reduced cell labels, not a Cartesian distance. The internal grid of g0 is decoupled from the
    output grid. Energies in eV.
"""
import numpy as np

from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk
from graphene_raman.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr


def mp_grid(n1, n2=None, n3=1):
    """Unshifted Gamma-centered MP grid in [0,1), reduced coords, C-order (i,j,l)."""
    n2 = n1 if n2 is None else n2
    i, j, l = np.meshgrid(np.arange(n1), np.arange(n2), np.arange(n3), indexing="ij")
    return np.stack([i.ravel() / n1, j.ravel() / n2, l.ravel() / n3], axis=1).astype(float)


def _phase(k, R):
    """exp(2 pi i k.R): (nk,3),(nL,3) -> (nk, nL)."""
    return np.exp(2j * np.pi * (k @ np.asarray(R, float).T))


def mwr_locality(Mwr, R_mwr, R0=(0, 0, 0)):
    """
    Block norm ||M_W(R, R0)|| against |R - R0| (reduced labels), sorted by distance. Guardrail a: raises if the largest
    block is not at R0 (defect not at the origin).
    Inputs:
        Mwr:   (nw, nR, nw, nR) complex, M_W(R, R').
        R_mwr: (nR, 3) ints, labels of Mwr, recentred on the defect.
        R0:    (3,) ints, defect cell.
    Returns:
        dist:   (nR,) floats, |R - R0| in reduced labels, sorted.
        weight: (nR,) floats, ||M_W(R, R0)||_F in the same order.
    """
    R_mwr = np.asarray(R_mwr, int)
    i0 = int(np.argmin(np.abs(R_mwr - np.asarray(R0)).sum(axis=1)))
    # Mwr shape (nw, nr, nw, nr); block norm between R and R0
    w = np.array([np.linalg.norm(Mwr[:, i, :, i0]) for i in range(len(R_mwr))])
    dist = np.linalg.norm(R_mwr - R_mwr[i0], axis=1)
    if w.argmax() != i0:
        raise AssertionError(f"guardrail a: max |Mwr(R,R0)| not at R0 (at R={R_mwr[w.argmax()]}); "
                             "defect not centered at origin / stray recentering phase.")
    order = np.argsort(dist)
    return dist[order], w[order]

def cluster_cells(Rn, R_cut):
    """
    Cells of the cluster: labels whose reduced-coordinate norm is <= R_cut (not a Cartesian distance), in the order of Rn.
    Inputs:
        Rn:    (nR, 3) ints, cell labels recentred on the defect (defect_mwr).
        R_cut: float, cutoff in cells (tolerance 1e-9).
    Returns:
        R_cluster: (nL, 3) ints, the retained labels.
    """
    mask = np.linalg.norm(Rn, axis=1) <= R_cut + 1e-9 # (nR,) bool, True for the cluster cells

    R_cluster = Rn[mask] # (nL, 3)
    return R_cluster
    
def cluster_potential(Mwr, R_mwr, R_cluster, herm_atol=1e-10):
    """
    Block of M_W on the cluster, flat index L*nw + w. Guardrail b: if the Hermiticity residual exceeds herm_atol, it is
    printed and the block is symmetrized.
    Inputs:
        Mwr:       (nw, nR, nw, nR) complex, M_W(R, R') (eV).
        R_mwr:     (nR, 3) ints, labels of Mwr.
        R_cluster: (nL, 3) ints, cluster cells (each must be in R_mwr).
        herm_atol: float, tolerance on max|M - M^dag|.
    Returns:
        M_cluster: (nL*nw, nL*nw) complex, Hermitian (eV).
        herm_res:  float, residual before symmetrization.
    """
    R_mwr = np.asarray(R_mwr, int)
    R_cluster = np.asarray(R_cluster, int)
    nw = Mwr.shape[0]
    idx = []
    for R in R_cluster:
        hit = np.where((R_mwr == R).all(axis=1))[0]
        if len(hit) == 0:
            raise ValueError(f"R_cluster vector {R} not present in Mwr R grid")
        idx.append(int(hit[0]))
    sub = Mwr[:, idx, :, :][:, :, :, idx]               # (nw, nL, nw, nL)
    nL = len(idx)
    M_cluster = np.transpose(sub, (1, 0, 3, 2)).reshape(nL * nw, nL * nw)   # index L*nw+w
    res = float(np.max(np.abs(M_cluster - M_cluster.conj().T)))
    if res > herm_atol:
        print(f"[guardrail b] M_cluster Hermiticity residual {res:.2e} > {herm_atol:.0e}; symmetrizing.", flush=True)
        M_cluster = 0.5 * (M_cluster + M_cluster.conj().T)
    return M_cluster, res


def recenter_mwr(Mwr, R_mwr, MP_grid):
    """
    Labels recentred on the defect: R_d is the cell with the largest on-site block ||M_W(R, R)||, and R - R_d is folded
    back into the MP-dual box. Only the labels change (Gamma is translation invariant).
    Inputs:
        Mwr:     (nw, nR, nw, nR) complex, M_W(R, R').
        R_mwr:   (nR, 3) ints, raw labels (Mwk_to_Mwr).
        MP_grid: (3,) ints, MP grid (D1, D2, D3).
    Returns:
        R_new: (nR, 3) ints, recentred labels, same order as R_mwr.
        R_d:   (3,) ints, defect cell in the raw labels.
    """
    R_mwr = np.asarray(R_mwr, int)
    w = np.array([np.linalg.norm(Mwr[:, i, :, i]) for i in range(len(R_mwr))])
    i_d = int(np.argmax(w))
    R_d = R_mwr[i_d].copy()
    Nv = np.asarray(MP_grid, int)
    R_new = ((R_mwr - R_d + Nv // 2) % Nv) - Nv // 2
    return R_new, R_d


def cluster_green(Hwk, k_int, R_cluster, eps, eta):
    """
    Cluster block of the pristine lattice Green's function at one energy, by direct inversion (reference of
    cluster_green_batch).
    Inputs:
        Hwk:       (nk, nw, nw) complex, Wannier Hamiltonian on k_int (eV).
        k_int:     (nk, 3) floats, internal k-grid (reduced coordinates).
        R_cluster: (nL, 3) ints, cluster cells.
        eps:       float, energy (eV).
        eta:       float, broadening (eV).
    Returns:
        g0: (nL*nw, nL*nw) complex, flat index L*nw + w (1/eV).
    """
    nki, nw, _ = Hwk.shape
    nL = len(R_cluster)
    A = (eps + 1j * eta) * np.eye(nw)[None] - Hwk               # (nki, nw, nw)
    G0k = np.linalg.inv(A)                                       # (nki, nw, nw)
    ph = _phase(k_int, R_cluster)                                 # (nki, nL)
    # g0[L,w,L',w'] = (1/nki) sum_k ph[k,L] G0k[k,w,w'] conj(ph[k,L'])
    g0 = np.einsum("kL,kwv,kM->LwMv", ph, G0k, np.conj(ph), optimize=True) / nki
    return g0.reshape(nL * nw, nL * nw)


def _diff_table(R_cluster):
    """Distinct lattice differences D = R_L - R_M of a cluster and the (L,M) -> D index map."""
    R = np.asarray(R_cluster, int)
    nL = len(R)
    D = (R[:, None, :] - R[None, :, :]).reshape(-1, 3)
    Du, inv = np.unique(D, axis=0, return_inverse=True)
    return Du, np.asarray(inv).reshape(nL, nL)


def cluster_green_batch(Hwk, k_int, R_cluster, egrid, eta, k_chunk=8192, e_chunk=512, deriv=False):
    """
    cluster_green on a whole energy grid, exactly restructured: g0 depends on R_L - R_L' only (sum over the distinct
    differences D), and G0(k) is expanded in the eigenbasis of H(k), so the k-sum is one matrix product per chunk.
    Inputs:
        Hwk:       (nk, nw, nw) complex, Wannier Hamiltonian on k_int (eV).
        k_int:     (nk, 3) floats, internal k-grid (reduced coordinates).
        R_cluster: (nL, 3) ints, cluster cells.
        egrid:     (nE,) floats, energies (eV).
        eta:       float, broadening (eV).
        k_chunk:   int, k-points per chunk (memory only).
        e_chunk:   int, energies per chunk (memory only).
        deriv:     bool, if True returns dg0/d(eps) (Lloyd: delta_rho = (1/pi) Im Tr[t dg0/d(eps)]).
    Returns:
        g0: (nE, nL*nw, nL*nw) complex, flat index L*nw + w.
    """
    nki, nw, _ = Hwk.shape
    R_cluster = np.asarray(R_cluster, int)
    nL = len(R_cluster)
    Du, inv = _diff_table(R_cluster)
    nD = len(Du)
    egrid = np.asarray(egrid, float)
    nE = len(egrid)
    eps, U = np.linalg.eigh(Hwk)                                   # (nki, nw), (nki, nw, nw)
    gD = np.zeros((nE, nD * nw * nw), dtype=complex)
    for s in range(0, nki, k_chunk):
        k = k_int[s:s + k_chunk]; e = eps[s:s + k_chunk]; u = U[s:s + k_chunk]
        ph = _phase(k, Du)                                         # (nk, nD) = ph[k,L] conj(ph[k,M])
        # W[(k,n), (D,w,v)] = e^{2 pi i k.D} U[k,w,n] conj(U[k,v,n])   (spectral weight of G0(k) at eps_kn)
        W = np.einsum("kD,kwn,kvn->knDwv", ph, u, u.conj(), optimize=True).reshape(-1, nD * nw * nw)
        ek = e.reshape(-1)
        for t in range(0, nE, e_chunk):
            den = 1.0 / (egrid[t:t + e_chunk, None] + 1j * eta - ek[None, :])   # (ne, nk*nw)
            if deriv:
                den = -den * den
            gD[t:t + e_chunk] += den @ W
    gD = gD.reshape(nE, nD, nw, nw) / nki
    g0 = gD[:, inv]                                                # (nE, nL, nL, nw, nw)
    return np.ascontiguousarray(g0.transpose(0, 1, 3, 2, 4)).reshape(nE, nL * nw, nL * nw)


def cluster_t(M_cluster, g0):
    """
    Cluster T-matrix t = M_cluster [1 - g0 M_cluster]^{-1}, i.e. t = M + M g0 t.
    Inputs:
        M_cluster: (n, n) complex, defect potential on the cluster (eV).
        g0:        (n, n) complex, cluster block of the Green's function at one energy (1/eV).
    Returns:
        t: (n, n) complex (eV).
    """
    n = M_cluster.shape[0]
    return M_cluster @ np.linalg.solve(np.eye(n) - g0 @ M_cluster, np.eye(n))


def cluster_ldos(g0, M_cluster, idx):
    """
    Local DOS of the defective lattice on cluster orbitals, rho_i = -(1/pi) Im G_ii with G = g0 + g0 t g0 = [1 - g0 M]^{-1} g0
    (one linear solve per energy, exact on the cluster), and of the pristine lattice, rho0_i = -(1/pi) Im g0_ii.
    Inputs:
        g0:        (nE, n, n) or (n, n) complex, cluster block of the Green's function (cluster_green_batch).
        M_cluster: (n, n) complex, defect potential on the cluster (eV).
        idx:       sequence of ints, cluster indices L*nw + w.
    Returns:
        rho:  (nE, len(idx)) floats, or (len(idx),) for a single energy; states / eV / orbital / spin.
        rho0: same shape, pristine lattice.
    """
    g0 = np.asarray(g0)
    single = g0.ndim == 2
    if single:
        g0 = g0[None]
    idx = np.asarray(idx, int)
    n = M_cluster.shape[0]
    A = np.eye(n)[None] - g0 @ M_cluster                                   # (nE, n, n): 1 - g0 M
    G = np.linalg.solve(A, g0)                                         # (nE, n, n): [1 - g0 M]^{-1} g0 = g0 + g0 t g0
    rho = -np.diagonal(G, axis1=1, axis2=2)[:, idx].imag / np.pi
    rho0 = -np.diagonal(g0, axis1=1, axis2=2)[:, idx].imag / np.pi
    return (rho[0], rho0[0]) if single else (rho, rho0)


def scattering_rate(Hwr, Rw, ndegen, M_cluster, R_cluster, k_out, eta,
                    k_int=None, positivity_atol=1e-8):
    """
    On-shell rate per defect, Gamma_nk = -2 Im phi_nk^dag t(eps_nk) phi_nk, with g0 computed exactly at each eps_nk
    (reference of scattering_rate_fast; slow). Raises if min Gamma < -positivity_atol (sign or gauge error).
    Inputs:
        Hwr:             (nRw, nw, nw) complex, H(R) of the Wannier model (eV).
        Rw:              (nRw, 3) ints, its R vectors.
        ndegen:          (nRw,) ints, Wigner-Seitz degeneracies.
        M_cluster:       (nL*nw, nL*nw) complex, defect potential on the cluster (eV).
        R_cluster:       (nL, 3) ints, cluster cells.
        k_out:           (nko, 3) floats, output k-points (reduced coordinates).
        eta:             float, broadening (eV).
        k_int:           (nki, 3) floats, internal grid of g0 (default 300 x 300).
        positivity_atol: float or None, tolerance of the positivity check.
    Returns:
        gamma: (nw, nko) floats, eV.
    """
    if k_int is None:
        k_int = mp_grid(300, 300, 1)
    R_cluster = np.asarray(R_cluster, int)
    nL = len(R_cluster)
    nw = Hwr.shape[1]
    assert M_cluster.shape == (nL * nw, nL * nw), f"M_cluster {M_cluster.shape} vs (nL*nw)={nL*nw}"

    Hwk_int, _, _ = Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=ndegen)              # (nki, nw, nw)
    _, E_out, U_out = Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)            # (nko, nw), (nko, nw, nw)
    ph_out = _phase(k_out, R_cluster)                                        # (nko, nL)
    # phi[k, n, (L,w)] = U_out[k,w,n] * ph_out[k,L]
    phi = np.einsum("kL,kwn->knLw", ph_out, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)

    nko = len(k_out)
    gamma = np.zeros((nw, nko))
    for ik in range(nko):
        for n in range(nw):
            eps = float(E_out[ik, n])
            g0 = cluster_green(Hwk_int, k_int, R_cluster, eps, eta)
            t = cluster_t(M_cluster, g0)
            v = phi[ik, n]
            sigma = v.conj() @ t @ v                                       # <nk|t|nk>
            gamma[n, ik] = -2.0 * sigma.imag
    if positivity_atol is not None and np.min(gamma) < -positivity_atol:
        raise AssertionError(f"positivity: min Gamma = {np.min(gamma):.2e} < 0 (sign/gauge bug)")
    return gamma


def scattering_rate_fast(Hwr, Rw, ndegen, M_cluster, R_cluster, k_out, eta, k_int=None,
                         e_window=None, ne_per_eta=8, positivity_atol=1e-8):
    """
    Same rate as scattering_rate, with g0 computed once on an energy grid of step eta / ne_per_eta and each state taking
    the nearest grid energy (controlled approximation, check C8). States outside e_window get NaN.
    Inputs:
        e_window:   (2,) floats or None, absolute energy window (eV).
        ne_per_eta: int, grid points per eta.
        others:     as in scattering_rate.
    Returns:
        gamma: (nw, nko) floats, eV; NaN outside e_window.
    """
    if k_int is None:
        k_int = mp_grid(300, 300, 1)
    R_cluster = np.asarray(R_cluster, int)
    nL, nw = len(R_cluster), Hwr.shape[1]
    assert M_cluster.shape == (nL * nw, nL * nw)
    Hwk_int, _, _ = Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=ndegen)
    _, E_out, U_out = Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)
    ph_out = _phase(k_out, R_cluster)
    phi = np.einsum("kL,kwn->knLw", ph_out, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)

    sel = np.ones_like(E_out, dtype=bool)
    if e_window is not None:
        sel = (E_out >= e_window[0]) & (E_out <= e_window[1])
    gamma = np.full((nw, len(k_out)), np.nan)
    if not sel.any():
        return gamma
    E_sel = E_out[sel]
    de = eta / ne_per_eta
    egrid = np.arange(E_sel.min() - eta, E_sel.max() + eta + de, de)
    g0_all = cluster_green_batch(Hwk_int, k_int, R_cluster, egrid, eta)   # exact, batched (see cluster_green_batch)
    t_cache = [cluster_t(M_cluster, g0_all[j]) for j in range(len(egrid))]
    del g0_all
    for ik in range(len(k_out)):
        for n in range(nw):
            if not sel[ik, n]:
                continue
            j = int(np.argmin(np.abs(egrid - E_out[ik, n])))
            v = phi[ik, n]
            gamma[n, ik] = -2.0 * (v.conj() @ t_cache[j] @ v).imag
    g = gamma[np.isfinite(gamma)]
    if positivity_atol is not None and g.size and g.min() < -positivity_atol:
        raise AssertionError(f"positivity: min Gamma = {g.min():.2e} < 0 (sign/gauge bug)")
    return gamma


def defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N=0.0):
    """
    M_W(R, R') of the defect: rotation V^dag M V, double Fourier transform, Kumagai-Oba alignment M_W(R, R) - C_N for every
    Wannier function and every cell of the N x N supercell box (approximation (i)), then labels recentred on the defect.
    Inputs:
        Mbk:    (nb, nk, nb, nk) complex, M in Bloch gauge on the MP grid k (eV).
        U:      (nk, nb, nw) or (nk, nw, nw) complex, Wannier gauge matrix (see Mbk_to_Mwk).
        U_dis:  (nk, nb, nw) complex or None, disentanglement matrix.
        k:      (nk, 3) floats, full unshifted MP grid in reduced coordinates.
        MP:     (3,) ints, MP grid (D1, D2, 1): D = N (coarse) or pN (dense).
        n_box:  int, supercell size N; box = cells with (R mod D) in [0, N)^2, on the raw labels.
        C_N:    float, alignment constant Delta V_PA^(N) in the unit of Mbk (config "alignment"); 0 = unaligned.
    Returns: dict with
        Mwr:    (nw, nR, nw, nR) complex, aligned M_W(R, R').
        R:      (nR, 3) ints, raw labels of Mwk_to_Mwr.
        Rn:     (nR, 3) ints, labels recentred on the defect (recenter_mwr).
        R_d:    (3,) ints, defect cell in the raw labels.
        in_box: (nR,) bool, cells of the supercell box.
    """

    Mwk = Mbk_to_Mwk(Mbk, U, U_dis) # (nW, nk, nW, nk)
    Mwr, R = Mwk_to_Mwr(Mwk, k, MP) # (nW, nR, nW, nR), (nR, 3)

    nW = Mwk.shape[0]

    c = R.copy()
    c[:, :2] %= MP[:2] # (nR, 3)

    in_box = np.all(c[:,:2] < n_box, axis=1)

    for w in range(nW):
        Mwr[w, in_box, w, in_box] -= C_N

    Rn, R_d = recenter_mwr(Mwr, R, MP)

    return {"Mwr": Mwr, "R": R, "Rn": Rn, "R_d": R_d, "in_box": in_box}

