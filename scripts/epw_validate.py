#!/usr/bin/env python
"""
epw_validate.py -- chapter-5 validations of the EPW pipeline (read-only post-processing), raw numbers + npz.
  bands   : bands.x (DFT) vs EPW band.dat on the same path, constant reference offset removed, per band; window +-3 eV of E_D
  phonons : matdyn (cm^-1) vs EPW freq.dat (meV) on the same path, per mode
  decay   : decay.{H,dynmat,epmate,epmatp} -> max|.|(R_max)/max|.|(0)
  g       : |g| DFPT (ph.x electron_phonon='prt', k fixed, q on a path) vs EPW (prtgkk, same k, q on a fine path), with
            band-ORDER-INDEPENDENT matching: at each common q the states at k and k+q are grouped into degenerate subspaces
            by energy (tol_deg), DFPT and EPW subspaces are paired by energy (tol_match), phonon modes are grouped by
            frequency (tol_w), and the gauge-invariant sums  G^2 = sum_{i in A, j in B, nu in M} |g_ij^nu|^2  are compared.
Energies: eV; |g| meV; E_D explicit (min pi/pi* gap on the path). Output: results/epw/validation_<tag>.npz
Usage: python scripts/epw_validate.py --root <grid dir> --tag 24k24q [--bands bands/graphene.qe.bands] ...
"""
import argparse, os, re, numpy as np
HA2EV = 27.211386245988; CM2MEV = 1.0 / 8.06554

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--bands", default="bands/graphene.qe.bands"); ap.add_argument("--epw-band", default="band_freq_interp/band.eig")
ap.add_argument("--epw-kpt", default="band_freq_interp/graphene_band.kpt"); ap.add_argument("--freq", default="phonons/graphene.freq"); ap.add_argument("--epw-freq", default="band_freq_interp/phband.freq")
ap.add_argument("--decay-dir", default="epw1"); ap.add_argument("--decay-prefix", default="decay")
ap.add_argument("--g-dirs", default="G:dfpt_g_G:epw_g_G,K:dfpt_g_K:epw_g_K"); ap.add_argument("--g-prefix", default="epw_g")
ap.add_argument("--tol-deg", type=float, default=2e-3); ap.add_argument("--tol-match", type=float, default=0.03); ap.add_argument("--tol-w", type=float, default=0.3)
ap.add_argument("--out", default=None); ap.add_argument("--fsthick", type=float, default=3.5, help="window |E-E_D| (eV) for the |g| statistics relevant downstream"); a = ap.parse_args()
R = a.root; out = {}
def P(x): return os.path.join(R, x)

# ---------------- cell (for path coordinates)
a1 = np.array([4.0354919061, -2.3298923383]); a2 = np.array([4.0354919061, 2.3298923383]); alat = np.linalg.norm(a1)
A = np.array([a1, a2]).T; B = 2 * np.pi * np.linalg.inv(A).T

def read_bandsx(f, nb=None):
    txt = open(f).read().split("\n"); hdr = re.search(r"nbnd=\s*(\d+),\s*nks=\s*(\d+)", txt[0]); nb = int(hdr.group(1)); ks, vals = [], []; i = 1
    while i < len(txt):
        p = txt[i].split()
        if len(p) == 3 and i + 1 < len(txt) and len(txt[i + 1].split()) > 3:
            ks.append([float(x) for x in p]); e = []; j = i + 1
            while len(e) < nb: e += [float(x) for x in txt[j].split()]; j += 1
            vals.append(e); i = j
        else: i += 1
    ks = np.array(ks); s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(ks, axis=0), axis=1))]); return s / s[-1], np.array(vals)
def read_epw_plot(f, nkf):
    """EPW band_plot output: either the QE '&plot nbnd=, nks=' format (band.eig / phband.freq: per k a coordinate line then
    the values) or the 2-column (s, value) gnuplot format; returns (nkf, nb)."""
    txt = open(f).read().split("\n")
    if txt[0].strip().startswith("&plot"):
        hdr = re.search(r"nbnd=\s*(\d+),\s*nks=\s*(\d+)", txt[0]); nb, nks = int(hdr.group(1)), int(hdr.group(2)); assert nks == nkf, (nks, nkf)
        vals = []; i = 1
        while i < len(txt) and len(vals) < nks:
            if len(txt[i].split()) == 3:
                e = []; j = i + 1
                while len(e) < nb: e += [float(x) for x in txt[j].split()]; j += 1
                vals.append(e); i = j
            else: i += 1
        return np.array(vals)
    d = np.loadtxt(f); nb = len(d) // nkf; return d[:, 1].reshape(nb, nkf).T
def path_coord_epw(kptfile):
    k = np.loadtxt(kptfile, skiprows=1)[:, :2]; kc = (B @ k.T).T; s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(kc, axis=0), axis=1))]); return s / s[-1], len(k)

# ---------------- bands
if os.path.exists(P(a.bands)) and os.path.exists(P(a.epw_band)):
    s_d, E_d = read_bandsx(P(a.bands)); s_e, nkf = path_coord_epw(P(a.epw_kpt)); E_e = read_epw_plot(P(a.epw_band), nkf); nw = E_e.shape[1]
    gap = E_d[:, 4] - E_d[:, 3]; iD = int(np.argmin(gap)); ED_d = 0.5 * (E_d[iD, 3] + E_d[iD, 4])
    Epi = np.interp(s_d, s_e, E_e[:, 3]); Eps = np.interp(s_d, s_e, E_e[:, 4]); iDe = int(np.argmin(Eps - Epi)); ED_e = 0.5 * (Epi[iDe] + Eps[iDe]); shift = ED_e - ED_d
    print(f"[bands] E_D bands.x {ED_d:.4f} eV, EPW {ED_e:.4f} eV, reference offset {shift:+.4f} eV removed; DFT nk {len(s_d)}, EPW nkf {nkf}, nw {nw}")
    rows = []
    for n in range(nw):
        Ee = np.interp(s_d, s_e, E_e[:, n]) - shift; err = np.min(np.abs(E_d - Ee[:, None]), axis=1); w3 = np.abs(Ee - ED_d) <= 3.0
        rows.append((err.max(), err[w3].max() if w3.any() else np.nan, np.sqrt((err[w3] ** 2).mean()) if w3.any() else np.nan))
        print(f"   WF band {n+1}: max|dE| path {rows[-1][0]*1e3:7.2f} meV; |E-E_D|<=3 eV: max {rows[-1][1]*1e3:7.2f} meV rms {rows[-1][2]*1e3:6.2f} meV")
    out.update(bands_ED_dft=ED_d, bands_ED_epw=ED_e, bands_offset=shift, bands_err=np.array(rows), bands_s=s_d, bands_E_dft=E_d, bands_E_epw_on_dft=np.stack([np.interp(s_d, s_e, E_e[:, n]) - shift for n in range(nw)], 1))

# ---------------- phonons
if os.path.exists(P(a.freq)) and os.path.exists(P(a.epw_freq)):
    txt = open(P(a.freq)).read().split("\n"); qs, fr = [], []; i = 1
    while i < len(txt):
        p = txt[i].split()
        if len(p) in (3, 4) and i + 1 < len(txt) and len(txt[i + 1].split()) >= 6: qs.append([float(x) for x in p[:3]]); fr.append([float(x) for x in txt[i + 1].split()[:6]]); i += 2
        else: i += 1
    qs = np.array(qs); F_md = np.sort(np.array(fr), axis=1); s_md = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(qs, axis=0), axis=1))]); s_md /= s_md[-1]
    s_e, nkf = path_coord_epw(P(a.epw_kpt)); F_e = np.sort(read_epw_plot(P(a.epw_freq), nkf), axis=1) / CM2MEV       # meV -> cm^-1
    rows = []
    for m in range(6):
        Fe = np.interp(s_md, s_e, F_e[:, m]); e = np.abs(Fe - F_md[:, m]); rows.append((e.max(), np.sqrt((e ** 2).mean()))); print(f"[phonons] mode {m+1}: max|dw| {e.max():6.2f} cm^-1 ({e.max()*CM2MEV:5.2f} meV), rms {rows[-1][1]:5.2f} cm^-1; matdyn {F_md[:,m].min():.0f}..{F_md[:,m].max():.0f} cm^-1")
    out.update(ph_err=np.array(rows), ph_s=s_md, ph_F_matdyn=F_md, ph_F_epw_on_matdyn=np.stack([np.interp(s_md, s_e, F_e[:, m]) for m in range(6)], 1))

# ---------------- decay
for name in ("H", "dynmat", "epmate", "epmatp"):
    f = P(f"{a.decay_dir}/{a.decay_prefix}.{name}")
    if not os.path.exists(f): continue
    d = np.loadtxt(f, comments="#"); r, v = d[:, 0], d[:, 1]; v0 = v[np.isclose(r, r.min())].max(); vf = v[r >= 0.9 * r.max()].max()
    print(f"[decay] {name:7s}: max|.|(R=0) {v0:.3e} Ry, max|.|(R>=0.9 R_max={r.max():.1f} A) {vf:.3e} Ry, ratio {vf/v0:.2e}"); out.update(**{f"decay_{name}_r": r, f"decay_{name}_v": v, f"decay_{name}_ratio": vf / v0})

# ---------------- |g| DFPT vs EPW, band-order independent
def read_dfpt_prt(path):
    """{q_cart: {'E': {ib: e_k}, 'Eq': {jb: e_kq}, 'w': {nu: omega}, 'g': {(ib,jb,nu): |g|}}} in the order of appearance."""
    blocks = {}; order = []; q = None
    for line in open(path):
        if "q coord.:" in line:
            q = tuple(round(float(x), 6) for x in line.split(":")[1].split()); blocks.setdefault(q, {"E": {}, "Eq": {}, "w": {}, "g": {}}); order.append(q) if q not in order else None
        p = line.split()
        if q is not None and len(p) == 10 and p[0].isdigit() and p[1].isdigit() and p[2].isdigit():
            ib, jb, nu = int(p[0]), int(p[1]), int(p[2]); b = blocks[q]; b["E"][ib] = float(p[3]); b["Eq"][jb] = float(p[4]); b["w"][nu] = float(p[5]); b["g"][(ib, jb, nu)] = float(p[7])
    return blocks, order
def groups(vals, tol):
    """indices grouped into degenerate subspaces (sorted by value)."""
    items = sorted(vals.items(), key=lambda kv: kv[1]); gs = []
    for k, v in items:
        if gs and abs(v - gs[-1][1][-1]) <= tol: gs[-1][0].append(k); gs[-1][1].append(v)
        else: gs.append(([k], [v]))
    return [(g, float(np.mean(v))) for g, v in gs]
for spec in a.g_dirs.split(","):
    kname, ddir, edir = spec.split(":")
    if not (os.path.isdir(P(ddir)) and os.path.isdir(P(edir))): continue
    dfpt, order = read_dfpt_prt(P(f"{ddir}/ph.out"))
    # ph.x prints "q coord." in Cartesian 2pi/alat units (the input path must be given with q_in_cryst_coord = .true.)
    qc = np.array([np.array(q[:2]) * 2 * np.pi / alat for q in order]); s_d = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(qc, axis=0), axis=1))]); s_d /= s_d[-1]
    # EPW prtgkk files: <prefix>_g<ib><jb>_<nu>.dat, rows = q on the fine path (same kpt file)
    files = sorted(f for f in os.listdir(P(edir)) if f.startswith(a.g_prefix) and f.endswith(".dat")); pairs = sorted({tuple(int(c) for c in re.search(r"_g(\d)(\d)_", f).groups()) for f in files})
    s_e, nkf = path_coord_epw(P(a.epw_kpt))
    epw = {}
    for ib, jb in pairs:
        for nu in range(1, 7):
            f = P(f"{edir}/{a.g_prefix}_g{ib}{jb}_{nu}.dat")
            if os.path.exists(f): epw[(ib, jb, nu)] = np.loadtxt(f)        # cols: ib jb nu enk enkq omega |g_sym| |g| Re Im
    res = []
    for iq, q in enumerate(order):
        b = dfpt[q]; ie = int(np.argmin(np.abs(s_e - s_d[iq])))
        if abs(s_e[ie] - s_d[iq]) > 0.6 / nkf: continue
        # EPW energies/modes at this q from the first file
        f0 = epw[next(iter(epw))]; Ek_e = {ib: epw[(ib, jb, 1)][ie, 3] for ib, jb in pairs}; Ekq_e = {jb: epw[(ib, jb, 1)][ie, 4] for ib, jb in pairs}; w_e = {nu: epw[(pairs[0][0], pairs[0][1], nu)][ie, 5] for nu in range(1, 7)}
        # ph.x 'prt' prints energies with a q-block-dependent reference: align each block on the fixed-k spectrum (offset that
        # maps the most EPW k-energies onto DFPT k-energies within 5 meV; the k state does not depend on q)
        Ed_k = np.array(sorted(b["E"].values())); Ee_k = np.array(sorted(Ek_e.values())); best = (0, 0.0)
        for cand in (Ed_k[:, None] - Ee_k[None, :]).ravel():
            nmatch = sum(np.min(np.abs(Ed_k - (e + cand))) < 5e-3 for e in Ee_k)
            if nmatch > best[0]: best = (nmatch, cand)
        off = best[1]; b = {"E": {i: e - off for i, e in b["E"].items()}, "Eq": {i: e - off for i, e in b["Eq"].items()}, "w": b["w"], "g": b["g"]}
        gk_d, gkq_d, gw_d = groups(b["E"], a.tol_deg), groups(b["Eq"], a.tol_deg), groups(b["w"], a.tol_w)
        gk_e, gkq_e, gw_e = groups(Ek_e, a.tol_deg), groups(Ekq_e, a.tol_deg), groups(w_e, a.tol_w)
        for Ad, Ea in gk_d:
            Ae = [g for g, e in gk_e if abs(e - Ea) <= a.tol_match]
            if not Ae: continue
            for Bd, Eb in gkq_d:
                Be = [g for g, e in gkq_e if abs(e - Eb) <= a.tol_match]
                if not Be: continue
                for Md, wm in gw_d:
                    Me = [g for g, w in gw_e if abs(w - wm) <= a.tol_w]
                    if not Me: continue
                    G2d = sum(b["g"].get((i, j, nu), 0.0) ** 2 for i in Ad for j in Bd for nu in Md)
                    G2e = sum(epw[(i, j, nu)][ie, 7] ** 2 for i in Ae[0] for j in Be[0] for nu in Me[0] if (i, j, nu) in epw)
                    res.append((s_d[iq], Ea, Eb, wm, np.sqrt(G2d), np.sqrt(G2e)))
    res = np.array(res); out[f"g_{kname}"] = res
    big = res[:, 4] > 20; err = np.abs(res[:, 5] - res[:, 4])
    print(f"[|g| k={kname}] {len(res)} (q, subspace_k, subspace_k+q, mode-group) triples on {len(order)} q; DFPT G {res[:,4].min():.1f}..{res[:,4].max():.1f} meV; "
          f"max|EPW-DFPT| {err.max():.2f} meV, rms {np.sqrt((err**2).mean()):.2f} meV; rel (G>20 meV): max {np.max(err[big]/res[big,4]):.2%}, median {np.median(err[big]/res[big,4]):.2%}")
    # pi-only summary: subspaces within 1 eV of E_D at k and k+q
    EDg = out.get("bands_ED_epw", -4.2389); pi = (np.abs(res[:, 1] - EDg) < 1.5) & (np.abs(res[:, 2] - EDg) < 1.5) & big
    if pi.any(): print(f"   pi<->pi (|E-E_D|<1.5 eV, G>20 meV): {pi.sum()} entries, max rel {np.max(err[pi]/res[pi,4]):.2%}, median {np.median(err[pi]/res[pi,4]):.2%}")
    # fsthick window (P18 complement, 2026-09-21): only sums whose initial (k) AND final (k+q) subspaces lie within |E - E_D| <= fsthick
    # are relevant downstream (elecselfen/phonselfen use fsthick = 3.5 eV); q = M kept here (excluded only in the control figure)
    win = (np.abs(res[:, 1] - EDg) <= a.fsthick) & (np.abs(res[:, 2] - EDg) <= a.fsthick) & big; noM = win & ~np.isclose(res[:, 0], 0.634, atol=0.003)
    for lab, m in (("fsthick window", win), ("fsthick window, q=M excluded", noM)):
        if m.any(): print(f"   {lab} (|E_k-E_D|,|E_k+q-E_D| <= {a.fsthick} eV, G>20 meV): {m.sum()} sums, rel median {np.median(err[m]/res[m,4]):.2%}, max {np.max(err[m]/res[m,4]):.2%}, |dG| max {err[m].max():.2f} meV")
    out[f"g_{kname}_win_mask"] = win; out[f"g_{kname}_win_median"] = float(np.median(err[win] / res[win, 4])) if win.any() else np.nan; out[f"g_{kname}_win_max"] = float(np.max(err[win] / res[win, 4])) if win.any() else np.nan
    out[f"g_{kname}_win_noM_median"] = float(np.median(err[noM] / res[noM, 4])) if noM.any() else np.nan; out[f"g_{kname}_win_noM_max"] = float(np.max(err[noM] / res[noM, 4])) if noM.any() else np.nan; out["g_fsthick"] = a.fsthick
o = a.out or f"results/epw/validation_{a.tag}.npz"; os.makedirs(os.path.dirname(o), exist_ok=True); np.savez(o, units="energies eV; |g| meV; phonons cm^-1; decay Ry", tag=a.tag, **out); print(f"saved {o}")
