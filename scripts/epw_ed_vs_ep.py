#!/usr/bin/env python
"""
epw_ed_vs_ep.py -- point-by-point ratio Gamma^ed(c) / Gamma^ep(T) on the energy axis eps - E_D (each curve relative to ITS OWN E_D),
Gamma^ed = c x Gamma_T (Lorentzian-averaged curve of resonance_<size>.npz, eta 0.02), Gamma^ep = Gamma_e of a selfen npz
(Lorentzian eta 0.02). Reports median / min / max of the ratio on |eps - E_D| <= window and the crossings (ratio = 1), plus
the window medians of both curves. Output: results/epw/ed_vs_ep_<tag>.npz. No new run.
Usage: python scripts/epw_ed_vs_ep.py --selfen results/epw/selfen_240_dg0.02_T300.npz --tag 24k24q [--resonance results/M/resonance_9x9.npz] [--window 3]
"""
import argparse, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--selfen", required=True); ap.add_argument("--tag", required=True)
ap.add_argument("--resonance", default="results/M/resonance_9x9.npz"); ap.add_argument("--window", type=float, default=3.0); a = ap.parse_args()
M = np.load(a.resonance, allow_pickle=True); S = np.load(a.selfen, allow_pickle=True); c = float(M["conc"])
xe = M["eg"] - float(M["E_D"]); Ged = c * M["Gamma_T"]; xp = S["eg"] - float(S["E_D"]); Gep = S["Gamma_e"]
m = np.abs(xp) <= a.window; x = xp[m]; ep = Gep[m]; ed = np.interp(x, xe, Ged); r = ed / ep
imin, imax = int(np.argmin(r)), int(np.argmax(r)); sgn = np.sign(r - 1); zc = np.where(np.diff(sgn) != 0)[0]
xc = np.array([x[i] + (x[i + 1] - x[i]) * (1 - r[i]) / (r[i + 1] - r[i]) for i in zc])
med_ed = float(np.median(ed)); med_ep = float(np.median(ep)); m12 = np.abs(x) <= 1.2
print(f"[ed/ep {a.tag}] c = {c*100:.0f} %, T = {float(S['T']):.0f} K, {int(S['n_mesh'])}^2, degaussw {float(S['degaussw'])} eV ; window +-{a.window} eV ({m.sum()} pts, 5 meV)")
print(f"   median Gamma^ed = {med_ed*1e3:.3f} meV (+-1.2 eV: {np.median(ed[m12])*1e3:.3f}), median Gamma^ep = {med_ep*1e3:.3f} meV (+-1.2 eV: {np.median(ep[m12])*1e3:.3f})")
print(f"   ratio Gamma^ed/Gamma^ep: median {np.median(r):.3f}, min {r[imin]:.3f} at {x[imin]:+.3f} eV, max {r[imax]:.3f} at {x[imax]:+.3f} eV ; crossings (ratio = 1) at {np.round(xc, 3).tolist()} eV")
out = f"results/epw/ed_vs_ep_{a.tag}.npz"; np.savez(out, x=x, Gamma_ed=ed, Gamma_ep=ep, ratio=r, c=c, T=float(S["T"]), median=np.median(r), min=r[imin], x_min=x[imin], max=r[imax], x_max=x[imax],
                                                  crossings=xc, median_ed=med_ed, median_ep=med_ep, selfen=a.selfen, resonance=a.resonance, units="eV; ratio dimensionless; energies rel. E_D of each chain")
print(f"saved {out}")
