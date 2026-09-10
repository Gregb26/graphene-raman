"""EPW elecselfen post-processing (read-only).
read_selfen(epw_out)      -> list of dict(k (nk,3) crystal, E (nstate,) eV ABSOLUTE (epw.out prints E - E_F; E_F added back once),
                            ReS, ImS (nstate,) eV, ik, ibnd, T (K), degaussw, E_F) parsed from epw.out
                            ('E( ib )= ... Re[Sigma]= ... Im[Sigma]= ...'; linewidth.elself.<T>K holds the same Im Sigma in meV);
                            Gamma^ep = 2 Im Sigma is formed here, once.
mp_weights(k, n)          -> weights of an irreducible k list of an n x n MP mesh (hexagonal cell, 12 in-plane ops + TR)
gamma_of_energy(...)      -> Lorentzian-weighted average (eta) of the on-shell Gamma_nk on an energy grid (chapter-4 recipe)
window_median(eg, G, E_D, half) -> median of Gamma(eps) over |eps - E_D| <= half (same statistic as resonance_metrics.py)
"""
import re
import numpy as np

RY2EV = 13.605693122994


def _sym_ops_2d():
    """Integer 2x2 matrices acting on k in RECIPROCAL crystal coords (b1, b2 at 120 degrees for the 60-degree real cell):
    S^T G* S = G*, 12 ops (C6v incl. mirrors). These are the inverse-transposes of the real-space ops."""
    G = np.array([[1.0, -0.5], [-0.5, 1.0]])         # reciprocal metric: b1.b2 = -|b|^2/2
    ops = []
    for a, b, c, d in np.ndindex(3, 3, 3, 3):
        S = np.array([[a - 1, b - 1], [c - 1, d - 1]])
        if abs(round(np.linalg.det(S))) == 1 and np.allclose(S.T @ G @ S, G):
            ops.append(S)
    assert len(ops) == 12, len(ops)
    return np.array(ops)


def _canon(kk, n):
    """canonical integer label of k (crystal, in units of 1/n) under the point group + time reversal."""
    ops = _sym_ops_2d(); ops = np.concatenate([ops, -ops])
    ki = np.rint(kk[:, :2] * n).astype(int)                       # (N, 2)
    imgs = np.einsum("sij,nj->nsi", ops, ki) % n                    # (N, nops, 2)
    lab = imgs[..., 0] * n + imgs[..., 1]
    return lab.min(axis=1)


def mp_weights(k, n):
    """weights (sum = 1) of the irreducible points k (crystal) of the n x n Gamma-centred MP mesh, whatever subgroup EPW used:
    every full-mesh point is shared equally among the listed k of its full-group (12 ops + time reversal) orbit; raises if a
    full-mesh orbit has no listed representative (incomplete list)."""
    k = np.asarray(k, float)
    if not np.allclose(k[:, :2] * n, np.rint(k[:, :2] * n), atol=1e-4):
        raise ValueError("k list is not commensurate with the mesh")
    full = np.array([(i / n, j / n, 0.0) for i in range(n) for j in range(n)])
    cf, ci = _canon(full, n), _canon(k, n)
    members = {}
    for i, c in enumerate(ci): members.setdefault(c, []).append(i)
    w = np.zeros(len(k))
    for c in cf:
        if c not in members:
            raise ValueError("full-mesh point without irreducible representative")
        for i in members[c]: w[i] += 1.0 / len(members[c])
    return w / w.sum()


def read_selfen(epw_out):
    """Parse the elecselfen tables of epw.out. Returns per-temperature dict list."""
    txt = open(epw_out).read()
    m = re.search(r"Fermi energy is read from the input file: Ef =\s*([-\d.]+)\s*eV", txt)
    EF = float(m.group(1)) if m else None
    m = re.search(r"Broadening:\s*([\d.]+)\s*eV", txt); dg = float(m.group(1)) if m else None
    nkf = re.search(r"nkf1\s*=\s*(\d+)", open(epw_out.replace("epw.out", "epw.in")).read()); n_mesh = int(nkf.group(1)) if nkf else None
    res = []
    parts = re.split(r"Temperature:\s*([\d.]+)K", txt)                 # [pre, T1, blk1, T2, blk2, ...]
    for T_K, blk in zip(parts[1::2], parts[2::2]):
        T_K = float(T_K); ks, rows = [], []
        for ik_m in re.finditer(r"ik =\s*(\d+) coord\.:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)(.*?)(?=\n\s*ik =|\Z)", blk, re.S):
            ik = int(ik_m.group(1)); kk = tuple(float(ik_m.group(i)) for i in (2, 3, 4)); ks.append((ik, kk))
            for r in re.finditer(r"E\(\s*(\d+)\s*\)=\s*([-\d.]+)\s*eV\s+Re\[Sigma\]=\s*([-\d.Ee+]+)\s*meV\s+Im\[Sigma\]=\s*([-\d.Ee+]+)\s*meV", ik_m.group(5)):
                # epw.out prints E relative to E_F (fermi_energy input): stored absolute (eV), the single reference shift
                rows.append((ik, int(r.group(1)), float(r.group(2)) + EF, float(r.group(3)) * 1e-3, float(r.group(4)) * 1e-3))
        if not rows: continue
        rows = np.array(rows); kd = dict(ks); kidx = sorted(kd)
        k = np.array([kd[i] for i in kidx]); kpos = {i: j for j, i in enumerate(kidx)}
        res.append(dict(T=T_K, degaussw=dg, E_F=EF, n_mesh=n_mesh, k=k, ik=np.array([kpos[int(i)] for i in rows[:, 0]]), ibnd=rows[:, 1].astype(int),
                        E=rows[:, 2], ReS=rows[:, 3], ImS=rows[:, 4], Gamma=2.0 * rows[:, 4]))
    if not res: raise RuntimeError(f"no elecselfen table in {epw_out}")
    return res


def gamma_of_energy(E, G, w, eg, eta):
    """Lorentzian-weighted average of on-shell widths: Gamma(eps) = sum_s w_s G_s L(eps - E_s) / sum_s w_s L(eps - E_s)."""
    L = (eta / np.pi) / ((eg[:, None] - E[None, :]) ** 2 + eta ** 2) * w[None, :]
    return (L * G[None, :]).sum(1) / L.sum(1)


def window_median(eg, Ge, E_D, half):
    m = np.abs(eg - E_D) <= half
    return float(np.nanmedian(Ge[m])), int(m.sum())
