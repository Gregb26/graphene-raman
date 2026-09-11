"""EPW phonselfen post-processing (read-only).
EPW's gamma___ / 'Phonon linewidth (meV)' in linewidth.phself.<T>K is documented in selfen.f90 (selfen_phon_q) as Im Pi_qnu, a
HALF width. P7 (scripts/epw_d2_extract.py, 2026-09-11) showed EMPIRICALLY that for graphene this number is 2.0 x the Dirac-cone Im Pi
built from EPW's own g and v_F at Gamma (brute-force checked, sin^2/cos^2 vertex structure confirmed on the q-zoom), and that it
coincides with the literature FWHM (10.8 cm^-1 at Gamma, 21.3 at K) while <D^2> matches Piscanec 2004 to 1-5 %. Convention adopted:
  gamma_epw  = raw EPW value (kept),  gamma_fwhm = gamma_epw (numerically the FWHM),  gamma_hwhm = gamma_epw / 2.
The origin of the factor inside EPW was not located in the source; the calibration is empirical (see NOTES_EPW.md 1d).
read_phself(run_dir) -> dict(T (nT,), q (nq,3) crystal, s (nq,) path coordinate (0..1, Cartesian lengths), omega (nq,6) meV,
                            gamma_hwhm (nT,nq,6) meV, lam (nT,nq,6)), q coordinates parsed from epw.out.
"""
import os, re, glob
import numpy as np

MEV2CM = 8.06554
# reciprocal basis of the 60-degree cell (b1, b2 at 120 deg), Cartesian in units 2pi/a: b1 = (1/sqrt3, -1)... only lengths matter
_B = np.array([[1 / np.sqrt(3.0), -1.0], [1 / np.sqrt(3.0), 1.0]]).T   # columns b1, b2 (2pi/a units, up to a common factor)


def path_coordinate(qc):
    kc = (_B @ qc[:, :2].T).T
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(kc, axis=0), axis=1))])
    return s / s[-1] if s[-1] > 0 else s


def read_phself(run_dir):
    txt = open(os.path.join(run_dir, "epw.out")).read()
    qs = {}
    for m in re.finditer(r"iq =\s*(\d+) coord\.:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", txt):
        qs[int(m.group(1))] = tuple(float(m.group(i)) for i in (2, 3, 4))
    nq = max(qs); q = np.array([qs[i] for i in range(1, nq + 1)])
    files = sorted(glob.glob(os.path.join(run_dir, "linewidth.phself.*K")), key=lambda f: float(re.search(r"phself\.([\d.]+)K", f).group(1)))
    T, om, gam = [], None, []
    for f in files:
        T.append(float(re.search(r"phself\.([\d.]+)K", f).group(1)))
        d = np.loadtxt(f, comments="#"); iq = d[:, 0].astype(int); im = d[:, 1].astype(int)
        w = np.full((nq, 6), np.nan); g = np.full((nq, 6), np.nan); w[iq - 1, im - 1] = d[:, 2]; g[iq - 1, im - 1] = d[:, 3]
        om = w if om is None else om; gam.append(g)
    lam = []
    for f in sorted(glob.glob(os.path.join(run_dir, "lambda.phself.*K")), key=lambda f: float(re.search(r"phself\.([\d.]+)K", f).group(1))):
        d = np.loadtxt(f, comments="#"); lam.append(d[:, 1:7] if d.ndim == 2 and d.shape[1] >= 7 else np.full((nq, 6), np.nan))
    gam = np.array(gam)
    return dict(T=np.array(T), q=q, s=path_coordinate(q), omega=om, gamma_epw=gam, gamma_hwhm=gam / 2.0, gamma_fwhm=gam, lam=np.array(lam) if lam else None,
                nkf=int(re.search(r"Using uniform k-mesh:\s*(\d+)", txt).group(1)), degaussw=float(re.search(r"Gaussian Broadening:\s*([\d.]+)", txt).group(1)),
                E_F=float(re.search(r"read from the input file: Ef =\s*([-\d.]+)", txt).group(1)))


def path_position(q, corners=((0, 0), (2 / 3, 1 / 3), (0.5, 0), (0, 0)), tol=2e-3):
    """position s in [0, 1] along the polyline through `corners` (crystal coords) for each q lying on it (nan otherwise).
    A point matching several segments (e.g. Γ at both ends) takes the candidate closest to the previous point's s, in list order."""
    C = np.array([(_B @ np.array(c)) for c in corners]); L = np.linalg.norm(np.diff(C, axis=0), axis=1); L0 = np.concatenate([[0], np.cumsum(L)])
    kc = (_B @ np.asarray(q)[:, :2].T).T; s = np.full(len(kc), np.nan); prev = 0.0
    for i, k in enumerate(kc):
        cands = []
        for j in range(len(L)):
            d = C[j + 1] - C[j]; t = np.dot(k - C[j], d) / L[j] ** 2; perp = np.linalg.norm(k - C[j] - t * d)
            if -1e-6 <= t <= 1 + 1e-6 and perp < tol: cands.append((L0[j] + t * L[j]) / L0[-1])
        if cands: s[i] = min(cands, key=lambda c: abs(c - prev)); prev = s[i]
    return s


def merge_runs(runs):
    """concatenate q lists of several runs (same T, nkf, degaussw), ordered along the Γ–K–M–Γ polyline; duplicates (same s) keep the first."""
    T = runs[0]["T"]
    for r in runs[1:]: assert np.allclose(r["T"], T) and r["nkf"] == runs[0]["nkf"] and r["degaussw"] == runs[0]["degaussw"], "incompatible runs"
    q = np.concatenate([r["q"] for r in runs]); s = np.concatenate([path_position(r["q"]) for r in runs])
    om = np.concatenate([r["omega"] for r in runs]); gam = np.concatenate([r["gamma_epw"] for r in runs], axis=1)
    lam = np.concatenate([r["lam"] for r in runs], axis=1) if all(r["lam"] is not None for r in runs) else None
    keep = ~np.isnan(s); order = np.argsort(s[keep], kind="stable"); idx = np.where(keep)[0][order]
    idx = idx[np.concatenate([[True], np.diff(s[idx]) > 1e-9])]
    return dict(T=T, q=q[idx], s=s[idx], omega=om[idx], gamma_epw=gam[:, idx], gamma_hwhm=gam[:, idx] / 2.0, gamma_fwhm=gam[:, idx], lam=lam[:, idx] if lam is not None else None, nkf=runs[0]["nkf"], degaussw=runs[0]["degaussw"], E_F=runs[0]["E_F"])


def special_points(q, tol=1e-4):
    """indices of Gamma, K (2/3,1/3 mod G, either valley) and M (1/2,0 and equivalents) in a crystal q list."""
    def eq(a, b): d = np.mod(np.asarray(a)[:2] - np.asarray(b)[:2] + 0.5, 1.0) - 0.5; return np.linalg.norm(d) < tol
    out = {}
    for i, qq in enumerate(q):
        for name, cands in (("G", [(0, 0)]), ("K", [(2 / 3, 1 / 3), (1 / 3, 2 / 3), (-1 / 3, 1 / 3), (1 / 3, -1 / 3)]), ("M", [(0.5, 0), (0, 0.5), (0.5, 0.5), (-0.5, 0.5)])):
            if any(eq(qq, c) for c in cands): out.setdefault(name, []).append(i)
    return out
