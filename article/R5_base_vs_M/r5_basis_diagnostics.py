"""
r5_basis_diagnostics.py  (R5, 2026-09-25 ; module de campagne, répertoire de travail, hors dépôt)
    Diagnostics of a truncated folded basis against reference (QE) supercell states.

    basis_diagnostics -- for H (hermitian, layout n*nk + k) and overlaps c: delta = 1 - c^dag c, Rayleigh quotient, residual
                         ||(H - e_ref) c|| / ||c||, weights on the eigenvectors of H (largest first)
    interlacing       -- sorted spectra of a reference and a model: counts below a bound, side-by-side lists in a window,
                         min_k [lambda_k(model) - lambda_{k-s}(ref)] for shifts s
    band_subblock     -- restrict (M, eps) to the first n bands
"""
import numpy as np


def basis_diagnostics(H, c, e_ref, evals=None, evecs=None, top=5):
    """
    H (n, n) hermitian (any energy scale), c (n,) overlaps of the reference state with the basis (reference state of unit norm),
    e_ref the reference energy on H's scale. Returns dict(norm_in = c^dag c, delta = 1 - c^dag c, rayleigh = c^dag H c / c^dag c,
    residual = ||(H - e_ref) c|| / ||c||, top = [(e_j, |<phi_j|c>|^2 / c^dag c)] for the `top` largest, sum_top).
    """
    c = np.asarray(c, complex)
    cc = float(np.vdot(c, c).real)
    if evals is None or evecs is None:
        evals, evecs = np.linalg.eigh(0.5 * (H + H.conj().T))
    Hc = H @ c
    R = float(np.vdot(c, Hc).real / cc)
    resid = float(np.linalg.norm(Hc - e_ref * c) / np.sqrt(cc))
    proj = np.abs(evecs.conj().T @ c) ** 2 / cc
    order = np.argsort(-proj)[:top]
    return dict(norm_in=cc, delta=1.0 - cc, rayleigh=R, residual=resid,
                top=[(float(evals[j]), float(proj[j])) for j in order], sum_top=float(proj[order].sum()))


def interlacing(eps_ref, eps_model, lo, hi, shifts=(0, 1, 2)):
    """
    Sorted spectra (same energy scale). Counts below lo, window lists [lo, hi], and for each shift s:
    d_k = lambda_k(model) - lambda_{k-s}(ref) over the common index range; the minimum over all k and over the k whose
    model value lies in [lo, hi] (with the index and the two values at the minimum).
    """
    r = np.sort(np.asarray(eps_ref, float)); m = np.sort(np.asarray(eps_model, float))
    out = dict(n_ref=int(len(r)), n_model=int(len(m)), n_ref_below_lo=int((r < lo).sum()), n_model_below_lo=int((m < lo).sum()),
               ref_window=[float(x) for x in r[(r >= lo) & (r <= hi)]], model_window=[float(x) for x in m[(m >= lo) & (m <= hi)]], shifts={})
    for s in shifts:
        kmax = min(len(m), len(r) + s)
        if kmax <= s:
            continue
        d = m[s:kmax] - r[:kmax - s]; ks = np.arange(s, kmax)
        j = int(d.argmin())
        rec = dict(min_all=float(d[j]), k_all=int(ks[j]), model_all=float(m[ks[j]]), ref_all=float(r[ks[j] - s]))
        inw = (m[s:kmax] >= lo) & (m[s:kmax] <= hi)
        if inw.any():
            jw = int(np.where(inw)[0][d[inw].argmin()])
            rec.update(min_window=float(d[jw]), k_window=int(ks[jw]), model_window=float(m[ks[jw]]), ref_window=float(r[ks[jw] - s]))
        out["shifts"][int(s)] = rec
    return out


def band_subblock(M, eps, n):
    """First n bands of M (nb, nk, nb, nk) and eps (nb, nk)."""
    return M[:n, :, :n, :], eps[:n]
