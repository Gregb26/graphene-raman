#!/usr/bin/env python
"""
potential_background.py -- C_N, in-plane background and vacuum offset of one supercell pair from its two pp.x filplots (R11, 2026-10-08).

The geometry (cell, atomic positions) is read from the filplot headers (qe_io.read_filplot): no .save, no scf.in needed. Works on any
pp.x 3D plot (Vks_* = plot_num 1, Vh_* = plot_num 11) as long as the pair shares its grid and cell. Measures (alignment.potential_background):
  C_N      Kumagai-Oba constant: spheres of --radius A around the atoms at true distance >= --frac r_max (config "alignment", R10 p1_size)
  plane    in-plane background: core-masked (--r-core) azimuthal profile of dV = V_d - V_p in the sheet plane (fig:Ved (c)), mean of the bins at
           r >= --frac |a1|/2, and its values at 2 and 3 A (sampling_table.csv)
  vacuum   plane-averaged dV for |z| > --z-vac A (tab:alignement (E), R10 c6)
Usage (from the repo root):
  .venv/bin/python scripts/m/potential_background.py <Vks_d> <Vks_p> [--json out.json] [--npz out.npz] [--gate 9x9]
--gate <size> compares C_N with config/production.json (alignment_C) and exits 3 beyond 1e-6 eV: the production Vks of that size must give
its C_N back (C_9 = -25.1437 meV). Energies in the json are in eV (printed in meV).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "src"))
from graphene_raman.config import load_production, alignment_C, HA2EV  # noqa: E402
from graphene_raman.io import qe_io  # noqa: E402
from graphene_raman.defects import alignment as al  # noqa: E402

BOHR = 0.529177210903
GATE_TOL = 1e-6     # eV


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    return x


def measure(path_d, path_p, radius=1.0, frac=al.PLATEAU_FRAC, r_core=0.5, dr=0.05, z_vac=5.0):
    """Read the two filplots (eV) and run alignment.potential_background. Returns (result, header_d, header_p)."""
    hd = qe_io.read_filplot(path_d)
    hp = qe_io.read_filplot(path_p)
    if hd["ngfft"] != hp["ngfft"]:
        raise SystemExit(f"grids differ: d {hd['ngfft']} vs p {hp['ngfft']}")
    if np.abs(hd["A_cols"] - hp["A_cols"]).max() > 1e-6:
        raise SystemExit(f"cells differ:\n{hd['A_cols']}\n{hp['A_cols']}")
    if hd["plot_num"] != hp["plot_num"]:
        raise SystemExit(f"plot_num differ: d {hd['plot_num']} vs p {hp['plot_num']}")
    V_d = hd.pop("V") * HA2EV
    V_p = hp.pop("V") * HA2EV
    res = al.potential_background(V_d, V_p, hd["x_red"], hp["x_red"], hd["A_cols"] * BOHR, radius, frac, r_core, dr, z_vac)
    return res, hd, hp


def summary(res, hd, hp, args):
    ko, rp, zp = res["kumagai_oba"], res["radial"], res["z"]
    return dict(files=dict(d=os.path.realpath(args.vks_d), p=os.path.realpath(args.vks_p)), plot_num=hd["plot_num"],
                n_atoms=dict(d=int(len(hd["x_red"])), p=int(len(hp["x_red"]))), ngfft=list(hd["ngfft"]),
                A_bohr_cols=hd["A_cols"], c_A=zp["c"], s_vac=res["s_vac"],
                params=dict(radius_A=args.radius, plateau_frac=args.frac, r_core_A=args.r_core, dr_A=args.dr, z_vac_A=args.z_vac),
                units="eV",
                C_N=dict(value=res["C_N"], rms=res["C_N_rms"], max_abs_dev=ko["max_abs_dev"], n=res["C_N_n"], r_max_A=res["r_max"], d_min_A=ko["d_min"]),
                plane=dict(value=res["plane"], rms=res["plane_rms"], n_bins=res["plane_n"], r_min_A=args.frac * rp["r_max"], r_max_A=rp["r_max"],
                           at=res["plane_at"], masked_fraction=rp["masked_fraction"], iz=rp["iz"]),
                vacuum=dict(value=res["vacuum"], std=res["vacuum_std"], min=zp["vacuum_min"], max=zp["vacuum_max"], n_planes=zp["n_vacuum"],
                            sheet_plane=zp["sheet_plane"]),
                mean3d=res["mean3d"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vks_d", help="filplot of the defective supercell")
    ap.add_argument("vks_p", help="filplot of the pristine supercell")
    ap.add_argument("--json", default=None, help="write the summary (eV) here")
    ap.add_argument("--npz", default=None, help="write the profiles and per-atom shifts here")
    ap.add_argument("--gate", default=None, help="size (e.g. 9x9): C_N must equal config alignment_C to 1e-6 eV, else exit 3")
    ap.add_argument("--radius", type=float, default=1.0, help="sphere radius, A (default 1.0)")
    ap.add_argument("--frac", type=float, default=al.PLATEAU_FRAC, help="plateau rule d >= frac * r_max (default 0.75)")
    ap.add_argument("--r-core", type=float, default=0.5, help="core mask radius of the in-plane profile, A (default 0.5)")
    ap.add_argument("--dr", type=float, default=0.05, help="radial bin width, A (default 0.05)")
    ap.add_argument("--z-vac", type=float, default=5.0, help="vacuum = planes with |z| > z_vac, A (default 5.0)")
    a = ap.parse_args()

    res, hd, hp = measure(a.vks_d, a.vks_p, a.radius, a.frac, a.r_core, a.dr, a.z_vac)
    s = summary(res, hd, hp, a)
    print(f"[background] {a.vks_d} - {a.vks_p} : plot_num {hd['plot_num']}, {s['n_atoms']['d']}/{s['n_atoms']['p']} atoms, grid {hd['ngfft']}, c = {s['c_A']:.3f} A")
    print(f"  C_N      = {res['C_N']*1e3:+.4f} meV  (rms {res['C_N_rms']*1e3:.2f}, {res['C_N_n']} atoms at d >= {s['C_N']['d_min_A']:.2f} A, spheres {a.radius} A)")
    print(f"  plane    = {res['plane']*1e3:+.2f} meV  (rms {res['plane_rms']*1e3:.2f}, {res['plane_n']} bins at r >= {s['plane']['r_min_A']:.2f} A ; "
          + ", ".join(f"r = {r:g} A: {v*1e3:+.1f}" for r, v in res["plane_at"].items()) + f" ; core mask {a.r_core} A)")
    print(f"  vacuum   = {res['vacuum']*1e3:+.3f} +- {res['vacuum_std']*1e3:.3f} meV  (|z| > {a.z_vac} A, {s['vacuum']['n_planes']} planes ; <dV>_3D {res['mean3d']*1e3:+.3f})")
    print(f"  C_N - vacuum = {(res['C_N'] - res['vacuum'])*1e3:+.2f} meV ; plane - vacuum = {(res['plane'] - res['vacuum'])*1e3:+.2f} meV")
    if a.gate:
        ref = alignment_C(load_production(verbose=False), a.gate)
        s["gate"] = dict(size=a.gate, config_eV=ref, diff_eV=res["C_N"] - ref, pass_=bool(abs(res["C_N"] - ref) < GATE_TOL))
        print(f"  gate {a.gate}: config {ref*1e3:+.4f} meV, diff {(res['C_N'] - ref)*1e3:.2e} meV -> {'PASS' if s['gate']['pass_'] else 'FAIL'}")
    if a.json:
        with open(a.json, "w") as fh:
            json.dump(jsonable(s), fh, indent=1)
        print(f"  wrote {a.json}")
    if a.npz:
        ko, rp, zp = res["kumagai_oba"], res["radial"], res["z"]
        np.savez(a.npz, rc_A=rp["rc"], rad_eV=rp["rad"], cnt=rp["cnt"], z_A=zp["z"], dV_xy_eV=zp["profile"],
                 dist_A=ko["shifts"]["dist"], shift_eV=ko["shifts"]["shift"], plateau_mask=ko["mask"], s_vac=res["s_vac"])
        print(f"  wrote {a.npz}")
    if a.gate and not s["gate"]["pass_"]:
        sys.exit(3)


if __name__ == "__main__":
    main()
