#!/usr/bin/env python
"""
epw_phself_post.py -- P6 post-processing of EPW phonselfen runs (read-only).
  --dir <run> [<run2>…] --tag <tag> -> results/epw/phself_<tag>.npz (several runs are merged and ordered along Γ–K–M–Γ) (omega, gamma HWHM & FWHM, lambda along the q list, both T) + key values at Gamma (E2g) and K (A1')
  --table <run dirs...>     -> convergence table of gamma(A1', K) HWHM (meV) per T, relative to --ref (substring of the run dir)
Convention (P7 calibration): EPW gamma___ is numerically the FWHM (2 x Dirac-cone Im Pi from EPW's own g); HWHM = gamma___/2. Energies meV; cm^-1 = meV x 8.06554.
"""
import argparse, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from electron_defect_interaction.electron_phonon import phself
ap = argparse.ArgumentParser(); ap.add_argument("--dir", nargs="*"); ap.add_argument("--tag"); ap.add_argument("--table", nargs="*"); ap.add_argument("--ref", default=None); a = ap.parse_args()
C = phself.MEV2CM

def key_values(R):
    sp = phself.special_points(R["q"]); rows = []
    for name, modes, label in (("G", (4, 5), "E2g(Gamma) LO/TO"), ("K", (2,), "A1'(K)")):
        for i in sp.get(name, []):
            for it, T in enumerate(R["T"]):
                g = R["gamma_fwhm"][it, i, list(modes)]; w = R["omega"][i, list(modes)]
                rows.append((label, i, T, w.mean(), g.mean(), g.max() - g.min()))
    return rows

if a.dir:
    R = phself.merge_runs([phself.read_phself(d) for d in a.dir]); tag = a.tag
    print(f"[phself {tag}] nkf {R['nkf']}x{R['nkf']} (full mesh), degaussw {R['degaussw']} eV, E_F {R['E_F']} eV, {len(R['q'])} q, T = {R['T']} K; gamma_epw = FWHM (P7 calibration), HWHM = gamma_epw/2")
    for label, i, T, w, g, spread in key_values(R):
        print(f"   {label:18s} iq={i+1:3d} T={T:6.1f} K: omega {w:8.3f} meV ({w*C:7.1f} cm^-1); gamma FWHM = EPW {g:8.4f} meV ({g*C:7.2f} cm^-1), HWHM {g/2:8.4f} meV; LO/TO spread {spread:.2e} meV")
    sp = phself.special_points(R["q"])
    for i in sp.get("G", []):
        print(f"   Gamma iq={i+1}: acoustic omega {R['omega'][i,:3]} meV, acoustic gamma(300K) {R['gamma_epw'][-1, i, :3]} meV (masked in the figure)")
    out = f"results/epw/phself_{tag}.npz"; os.makedirs("results/epw", exist_ok=True)
    np.savez(out, tag=tag, units="meV (omega, gamma); q crystal; s = path coordinate", convention="gamma_epw = raw EPW gamma___; P7 calibration (2026-09-11): gamma_fwhm = gamma_epw, gamma_hwhm = gamma_epw/2 (EPW value = 2 x Dirac-cone Im Pi from EPW's g; matches literature FWHM)",
             E_F=R["E_F"], E_D=-4.2389, nkf=R["nkf"], degaussw=R["degaussw"], fsthick=3.5, T=R["T"], q=R["q"], s=R["s"], omega=R["omega"], gamma_epw=R["gamma_epw"], gamma_hwhm=R["gamma_hwhm"], gamma_fwhm=R["gamma_fwhm"],
             lam=R["lam"] if R["lam"] is not None else np.array([]), key_values=np.array([(r[0], r[1], r[2], r[3], r[4]) for r in key_values(R)], dtype=object))
    print(f"   saved {out}")

if a.table is not None:
    rows = []
    for d in a.table:
        R = phself.read_phself(d); sp = phself.special_points(R["q"]); i = sp["K"][0]
        rows.append((os.path.basename(d.rstrip("/")), R["nkf"], R["degaussw"], [R["gamma_epw"][it, i, 2] for it in range(len(R["T"]))], R["T"]))
    ref = a.ref and next((r for r in rows if a.ref in r[0]), None)
    Ts = rows[0][4]; hdr = f"{'run':30s} {'nkf':>5s} {'degaussw':>9s}" + "".join(f"  gamma_epw(A1',K) {T:.0f}K (meV)" for T in Ts) + ("   rel. ref " + " / ".join(f"{T:.0f}K" for T in Ts) if ref else "")
    print(hdr)
    for r in rows:
        line = f"{r[0]:30s} {r[1]:5d} {r[2]:9.3f}" + "".join(f"  {g:22.4f}" for g in r[3])
        if ref: line += "   " + " / ".join(f"{g/gr-1:+8.2%}" if gr else "   n/a  " for g, gr in zip(r[3], ref[3]))
        print(line)
