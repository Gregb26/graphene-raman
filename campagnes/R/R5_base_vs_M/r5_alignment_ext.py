"""
r5_alignment_ext.py  (R5, 2026-09-25 ; module de campagne, répertoire de travail, hors dépôt)
    Additions planned for defects/alignment.py and wannier/supercell_fold.py (kept here, outside the repository):

    check_k_reference   -- refuse a k list that differs from the k of the .save providing the coefficients by a non-zero integer
                           vector (the R4 D4 error of 0.3: plane-wave index built with k81 in [-1/2, 1/2) for coefficients referred to
                           k27 in [0, 1)); returns the integer offsets
    dirac_quadruplet    -- E_D of a perfect N x N supercell (N = 3m): mean of the n states closest to E_F, their spread
    gamma_shift         -- rigid shift between two Gamma spectra on their nlow lowest states (median, max residual)
    lowest_state_shift  -- shift between the lowest states of two spectra
    size_fits           -- eps(N) = eps_inf + a / N^p by least squares, per family
"""
import numpy as np


def check_k_reference(k_used, k_save, tol=1e-6):
    """k_used and k_save (nk, 3): same modulo 1 is required; if they differ by a non-zero integer vector -> ValueError."""
    d = np.asarray(k_used, float) - np.asarray(k_save, float)
    frac = np.abs(d - np.rint(d)).max()
    if frac > tol:
        raise ValueError(f"check_k_reference: k lists differ by a non-integer vector (max frac dev {frac:.2e})")
    off = np.rint(d).astype(int)
    nbad = int((np.abs(off).sum(axis=1) > 0).sum())
    if nbad:
        raise ValueError(f"check_k_reference: {nbad} of {len(d)} k differ from the .save k by a non-zero integer vector "
                         f"(offsets {np.unique(off, axis=0).tolist()}); the plane-wave index must use the .save k")
    return off


def dirac_quadruplet(eps, ef, n=4):
    """E_D = mean of the n eigenvalues closest to ef; returns (E_D, sorted indices, spread = max - min)."""
    eps = np.asarray(eps, float)
    q = np.argsort(np.abs(eps - ef))[:n]
    return float(eps[q].mean()), sorted(int(i) for i in q), float(eps[q].max() - eps[q].min())


def gamma_shift(eps_a, eps_b, nlow=8):
    """Rigid shift a - b on the nlow lowest states (sorted): (median, max |residual|, nlow)."""
    a = np.sort(np.asarray(eps_a, float))[:nlow]; b = np.sort(np.asarray(eps_b, float))[:nlow]
    d = a - b; s = float(np.median(d))
    return s, float(np.abs(d - s).max()), int(len(d))


def lowest_state_shift(eps_a, eps_b):
    """min(eps_a) - min(eps_b)."""
    return float(np.min(eps_a) - np.min(eps_b))


def size_fits(N, eps, powers=(1, 2)):
    """
    Least-squares eps(N) = eps_inf + a / N^p for each p in powers. N, eps: 1-d arrays (same length, >= 2 points).
    Returns {p: dict(eps_inf, a, rms, max_resid, n)}.
    """
    N = np.asarray(N, float); e = np.asarray(eps, float); out = {}
    for p in powers:
        X = np.stack([np.ones_like(N), 1.0 / N ** p], axis=1)
        coef, *_ = np.linalg.lstsq(X, e, rcond=None)
        r = e - X @ coef
        out[int(p)] = dict(eps_inf=float(coef[0]), a=float(coef[1]), rms=float(np.sqrt(np.mean(r ** 2))), max_resid=float(np.abs(r).max()), n=int(len(N)))
    return out
