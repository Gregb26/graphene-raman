"""
local_tmatrix.py
    Single-defect scattering rate via the EXACT local t-matrix, evaluated on an arbitrarily fine
    output k-grid without ever forming the (nw*nkf)^2 dense T-matrix.

    For a defect localized on a finite set of Wannier sites R_local (centered on R0 = 0):

        V_loc = P Mwr P                              (local block of the Wannier-real-space M)
        g0(eps) = local block of the FULL lattice Green's function,
                  g0[wR, w'R'] = (1/Nk_int) sum_k e^{2pi i k.(R-R')} [(eps+i eta) - H_W(k)]^{-1}
        t(eps)  = V_loc [1 - g0(eps) V_loc]^{-1}     (exact if R_local covers supp(V))
        Gamma_perdef_nk = -2 Im < nk | t(eps_nk) | nk >,  with
                  <wR|nk> = U(k)_{wn} e^{2pi i k.R}  (phi_nk),  a gauge-invariant diagonal element.

    The internal grid Nk_int (for g0) is DECOUPLED from the output grid k_out and should be dense.
    Gamma_perdef is the per-defect (concentration-normalized) rate; multiply by 1/Nk for a given
    array concentration.
"""
import numpy as np

from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk


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
    Frobenius weight ||Mwr(R,R0)|| vs |R-R0| (documents the defect localization). Returns
    (dist, weight) sorted by distance, and asserts the weight peaks at R=R0 (guardrail a: the
    defect must sit at the origin, no stray recentering phase).
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


def extract_V_loc(Mwr, R_mwr, R_local, herm_atol=1e-10):
    """
    V_loc = P Mwr P for R,R' in R_local, flattened as index=L*nw+w. Guardrail b: check Hermiticity
    of V_loc; if beyond herm_atol, symmetrize and log the residual (interpolation can break it).
    Returns (V_loc, herm_residual).
    """
    R_mwr = np.asarray(R_mwr, int)
    R_local = np.asarray(R_local, int)
    nw = Mwr.shape[0]
    idx = []
    for R in R_local:
        hit = np.where((R_mwr == R).all(axis=1))[0]
        if len(hit) == 0:
            raise ValueError(f"R_local vector {R} not present in Mwr R grid")
        idx.append(int(hit[0]))
    sub = Mwr[:, idx, :, :][:, :, :, idx]               # (nw, nL, nw, nL)
    nL = len(idx)
    V = np.transpose(sub, (1, 0, 3, 2)).reshape(nL * nw, nL * nw)   # index L*nw+w
    res = float(np.max(np.abs(V - V.conj().T)))
    if res > herm_atol:
        print(f"[guardrail b] V_loc Hermiticity residual {res:.2e} > {herm_atol:.0e}; symmetrizing.", flush=True)
        V = 0.5 * (V + V.conj().T)
    return V, res


def local_green(Hwk, k_int, R_local, eps, eta):
    """
    Local block g0[(L,w),(L',w')] of the lattice Green's function at energy eps.
    Hwk: (nki, nw, nw) Wannier Hamiltonian on the internal grid k_int (nki,3).
    Returns (nL*nw, nL*nw) complex, flattened as index = L*nw + w.
    """
    nki, nw, _ = Hwk.shape
    nL = len(R_local)
    A = (eps + 1j * eta) * np.eye(nw)[None] - Hwk               # (nki, nw, nw)
    G0k = np.linalg.inv(A)                                       # (nki, nw, nw)
    ph = _phase(k_int, R_local)                                 # (nki, nL)
    # g0[L,w,L',w'] = (1/nki) sum_k ph[k,L] G0k[k,w,w'] conj(ph[k,L'])
    g0 = np.einsum("kL,kwv,kM->LwMv", ph, G0k, np.conj(ph), optimize=True) / nki
    return g0.reshape(nL * nw, nL * nw)


def local_t(V_loc, g0):
    """t = V_loc [1 - g0 V_loc]^{-1}."""
    n = V_loc.shape[0]
    return V_loc @ np.linalg.solve(np.eye(n) - g0 @ V_loc, np.eye(n))


def scattering_rate(Hwr, Rw, ndegen, V_loc, R_local, k_out, eta,
                    k_int=None, perdef=True, positivity_atol=1e-8):
    """
    On-shell per-defect scattering rate Gamma[n,k] on the output grid k_out.
    V_loc: (nL*nw, nL*nw) Hermitian, flattened as L*nw+w. Returns (nw, nk_out).
    """
    if k_int is None:
        k_int = mp_grid(300, 300, 1)
    R_local = np.asarray(R_local, int)
    nL = len(R_local)
    nw = Hwr.shape[1]
    assert V_loc.shape == (nL * nw, nL * nw), f"V_loc {V_loc.shape} vs (nL*nw)={nL*nw}"

    Hwk_int, _, _ = Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=ndegen)              # (nki, nw, nw)
    _, E_out, U_out = Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)            # (nko, nw), (nko, nw, nw)
    ph_out = _phase(k_out, R_local)                                        # (nko, nL)
    # phi[k, n, (L,w)] = U_out[k,w,n] * ph_out[k,L]
    phi = np.einsum("kL,kwn->knLw", ph_out, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)

    nko = len(k_out)
    gamma = np.zeros((nw, nko))
    for ik in range(nko):
        for n in range(nw):
            eps = float(E_out[ik, n])
            g0 = local_green(Hwk_int, k_int, R_local, eps, eta)
            t = local_t(V_loc, g0)
            v = phi[ik, n]
            sigma = v.conj() @ t @ v                                       # <nk|t|nk>
            gamma[n, ik] = -2.0 * sigma.imag
    if positivity_atol is not None and np.min(gamma) < -positivity_atol:
        raise AssertionError(f"positivity: min Gamma = {np.min(gamma):.2e} < 0 (sign/gauge bug)")
    return gamma


def scattering_rate_from_wannier(M_coarse, k_coarse, U, U_dis, Hwr, Rw, ndegen,
                                 R_local, k_out, eta, k_int=None, R0=(0, 0, 0)):
    """
    Production path: interpolate the coarse Bloch-gauge M to Wannier real space, extract the local
    defect block, and evaluate the per-defect on-shell rate on k_out. Enforces (RAISE):
      - nw consistency between U and H(R) (same wannierization);
      - k-grid consistency between U and the coarse M grid;
      - guardrail a: Mwr weight centered at R0 (mwr_locality);
      - guardrail b: V_loc Hermitian (extract_V_loc symmetrizes + logs).
    The caller must have already passed the hard gauge/provenance gate (wannier_provenance.load_
    wannier_checked) before reading U and H(R).
    """
    from electron_defect_interaction.wannier.wannier_interpolation import (
        Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order)

    nw_H = Hwr.shape[1]
    nw_U = U.shape[-1]
    if nw_H != nw_U:
        raise ValueError(f"gauge check: nw mismatch U({nw_U}) vs H(R)({nw_H}) -- different wannierizations.")
    if U.shape[0] != len(k_coarse):
        raise ValueError(f"gauge check: U has {U.shape[0]} k, coarse M has {len(k_coarse)} -- grid mismatch.")

    Mwk = Mbk_to_Mwk(M_coarse, U, U_dis)                        # (nw, nk, nw, nk)
    Mwr, R_mwr = Mwk_to_Mwr(Mwk, k_coarse, _infer_mp_grid(k_coarse))
    mwr_locality(Mwr, R_mwr, R0=R0)                             # guardrail a (hard)
    V_loc, _ = extract_V_loc(Mwr, R_mwr, R_local)              # guardrail b
    return scattering_rate(Hwr, Rw, ndegen, V_loc, R_local, k_out, eta, k_int=k_int)
