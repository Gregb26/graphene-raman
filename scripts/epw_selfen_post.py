#!/usr/bin/env python
"""
epw_selfen_post.py -- P2/P4 post-processing of an EPW elecselfen run (read-only): parse epw.out, Gamma^ep_nk = 2 Im Sigma_nk,
irreducible-wedge weights of the nkf x nkf MP mesh, Lorentzian energy average (eta from config/production.json, 0.02 eV as in
chapter 4), medians of Gamma^ep(eps) over |eps - E_D| <= 3 and 1.2 eV. Energies eV, E_D explicit, no rescaling.
Usage: epw_selfen_post.py --dir <epw run dir> --tag <tag>            -> results/epw/selfen_<tag>.npz (one per temperature: _T300 ...)
       epw_selfen_post.py --table results/epw/selfen_*.npz            -> convergence table (medians, relative to a reference)
"""
import argparse, os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from electron_defect_interaction.config import load_production
from electron_defect_interaction.electron_phonon import selfen

ap = argparse.ArgumentParser(); ap.add_argument("--dir"); ap.add_argument("--tag"); ap.add_argument("--table", nargs="*"); ap.add_argument("--ref", default=None)
ap.add_argument("--E_D", type=float, default=-4.2389); a = ap.parse_args()
cfg = load_production(verbose=False); eta = cfg["eta_eV"]; E_D = a.E_D

if a.dir:
    runs = selfen.read_selfen(os.path.join(a.dir, "epw.out"))
    for r in runs:
        n = r["n_mesh"]; w = selfen.mp_weights(r["k"], n); ws = w[r["ik"]]
        inwin = np.abs(r["E"] - r["E_F"]) <= 3.5            # EPW: only eigenstates within fsthick are meaningful
        r = {**r, **{key: r[key][inwin] for key in ("ik", "ibnd", "E", "ReS", "ImS", "Gamma")}}; ws = ws[inwin]
        eg = np.arange(E_D - 3.5, E_D + 3.5 + 1e-9, 0.005)
        Ge = selfen.gamma_of_energy(r["E"], r["Gamma"], ws, eg, eta)
        m3, n3 = selfen.window_median(eg, Ge, E_D, 3.0); m12, n12 = selfen.window_median(eg, Ge, E_D, 1.2)
        s3 = np.abs(r["E"] - E_D) <= 3.0; s12 = np.abs(r["E"] - E_D) <= 1.2
        tag = f"{a.tag}_T{int(round(r['T']))}"
        print(f"[selfen {tag}] nkf {n}x{n}: {len(r['k'])} irreducible k (weights sum {w.sum():.6f}), {len(r['E'])} states in fsthick, T = {r['T']:.1f} K, "
              f"degaussw {r['degaussw']} eV, E_F(EPW input) {r['E_F']} eV, E_D {E_D} eV, eta {eta} eV")
        print(f"   Gamma^ep = 2 Im Sigma: states |eps-E_D|<=3 eV: n={s3.sum()}, median {np.median(r['Gamma'][s3])*1e3:.3f} meV, max {r['Gamma'][s3].max()*1e3:.2f} meV; "
              f"|eps-E_D|<=1.2 eV: n={s12.sum()}, median {np.median(r['Gamma'][s12])*1e3:.3f} meV")
        print(f"   Gamma^ep(eps) (Lorentzian eta={eta}): median over |eps-E_D|<=3 eV ({n3} pts) = {m3*1e3:.3f} meV; over <=1.2 eV ({n12} pts) = {m12*1e3:.3f} meV; "
              f"at E_D {np.interp(E_D, eg, Ge)*1e3:.3f} meV; min {Ge.min()*1e3:.3f} meV at {eg[np.argmin(Ge)]-E_D:+.3f} eV")
        out = f"results/epw/selfen_{tag}.npz"; os.makedirs("results/epw", exist_ok=True)
        np.savez(out, tag=tag, units="eV", convention="Gamma^ep = 2 Im Sigma (total width); Gamma^ed = -2 Im Sigma (ch. 4)", E_D=E_D, E_F_epw=r["E_F"], T=r["T"],
                 degaussw=r["degaussw"], n_mesh=n, eta=eta, fsthick=3.5, k=r["k"], w=w, ik=r["ik"], ibnd=r["ibnd"], E=r["E"], ReS=r["ReS"], ImS=r["ImS"], Gamma=r["Gamma"],
                 eg=eg, Gamma_e=Ge, median_3eV=m3, median_1p2eV=m12, n_states_3eV=int(s3.sum()), n_states_1p2eV=int(s12.sum()))
        print(f"   saved {out}")

if a.table is not None:
    rows = []
    for f in a.table:
        R = np.load(f); rows.append((os.path.basename(f), int(R["n_mesh"]), float(R["degaussw"]), float(R["T"]), float(R["median_3eV"]) * 1e3, float(R["median_1p2eV"]) * 1e3, int(R["n_states_3eV"])))
    ref = a.ref and next((r for r in rows if a.ref in r[0]), None)
    print(f"{'run':38s} {'nkf':>5s} {'degaussw':>9s} {'T(K)':>6s} {'med|3eV|(meV)':>14s} {'med|1.2eV|':>11s} {'n_st':>6s}" + ("   rel.ref(3eV)  rel.ref(1.2eV)" if ref else ""))
    for r in rows:
        line = f"{r[0]:38s} {r[1]:5d} {r[2]:9.3f} {r[3]:6.0f} {r[4]:14.3f} {r[5]:11.3f} {r[6]:6d}"
        if ref: line += f"   {r[4]/ref[4]-1:+12.2%}  {r[5]/ref[5]-1:+12.2%}"
        print(line)
