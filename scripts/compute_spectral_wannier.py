#!/usr/bin/env python
"""
compute_spectral_wannier.py
    Level 1 of the two-level convergence protocol (fixed supercell N): interpolate the coarse M to a
    dense output k-grid via Wannier functions and evaluate the per-defect on-shell rate with the EXACT
    local t-matrix. Sweeps (R_cut shell) x (output grid) x (eta) and reports the Gamma(eta, grid) MAP
    plus the joint-plateau check and the vacancy-resonance position.

    HARD gauge gate: the run REFUSES to start unless wannier_provenance.load_wannier_checked passes
    (matching sha256 for tb/u/u_dis from one wannier run). No warn-and-continue.

    Usage:
        python scripts/compute_spectral_wannier.py --size 7x7 --manifest <seed>_wann.json \
            --grids 60,120,240 --etas 0.05,0.02,0.01 --rcut 0,1,2 --nk-int 300
"""
import argparse
import numpy as np

from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import (
    Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order)
from electron_defect_interaction.defects.many_body import local_tmatrix as lt

HA2EV = 27.211386245988


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--size", required=True)
    p.add_argument("--manifest", required=True, help="wannier provenance manifest (hard gauge gate)")
    p.add_argument("--grids", default="60,120,240", help="output k-grid densities N (NxN)")
    p.add_argument("--etas", default="0.05,0.02,0.01", help="broadenings (eV)")
    p.add_argument("--rcut", default="0,1,2", help="R-shell cutoffs (max |R| in reduced coords)")
    p.add_argument("--nk-int", type=int, default=300, help="internal k-grid density for g0")
    p.add_argument("--out", default=None)
    p.add_argument("--dense", action="store_true",
                   help="use the zero-padded dense M (M_dense_<size>.npy) and the dense primitive .save / dense wannierization")
    return p.parse_args()


def main():
    args = parse_args()

    # --- HARD gauge gate: refuse to run unless tb/u/u_dis are the same, unmodified wannier run ---
    paths = wannier_provenance.load_wannier_checked(args.manifest)   # raises on any failure
    print(f"[gauge] provenance OK: {args.manifest}", flush=True)

    PF = {"5x5": (5, 25), "7x7": (4, 28), "8x8": (4, 32), "9x9": (3, 27)}
    if args.dense:
        D = PF[args.size][1]
        uc = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"
        mfile = f"results/M/M_dense_{args.size}.npy"
    else:
        uc = f"data/graphene/unit_cell/qe/defect_{args.size}.save"
        mfile = f"results/M/M_ed_{args.size}.npy"
    # NORMALIZATION CONTRACT (validated by test_local_tmatrix_real.py to 1e-13):
    #   * the LOCAL Wannier t-matrix needs the INTENSIVE real-space potential V_loc = <wR|V|w'R'>,
    #     i.e. Mwr built from the unit-cell-normalized M_raw (bloch_norm='unit_cell');
    #   * the DENSE Bloch T-matrix (single_defect.compute_T) needs M/N_cells ('supercell').
    #   Feeding M_norm here silently suppresses V_loc by 1/N_cells (Born limit, Gamma ~ 0).
    M = matrix_io.load_M_checked(mfile, require_bloch_norm=matrix_io.UNIT_CELL)
    k_coarse = qe_io.get_k_red(uc)

    U, k_U = read_w90_mat(paths["u"])
    U = U[_match_kpoint_order(k_U, k_coarse)]
    U_dis = None
    if paths.get("u_dis"):
        U_dis, k_Ud = read_w90_mat(paths["u_dis"])
        U_dis = U_dis[_match_kpoint_order(k_Ud, k_coarse)]
    Hwr, Rw, ndegen = read_w90_HR(paths["tb"])

    # interpolate coarse M -> Wannier real space once; locality guardrail (hard)
    Mwk = Mbk_to_Mwk(M, U, U_dis)
    MP = _infer_mp_grid(k_coarse)
    Mwr, R_mwr = Mwk_to_Mwr(Mwk, k_coarse, MP)
    R_mwr, R_d = lt.recenter_mwr(Mwr, R_mwr, MP)                     # defect -> origin, applied ONCE
    print(f"[recenter] defect site detected at R_d={R_d.tolist()} (supercell cell index); labels shifted so R0=0", flush=True)
    dist, wt = lt.mwr_locality(Mwr, R_mwr)                            # guardrail a (raises if off-center)
    print("[locality] ||Mwr(R,R0)|| (eV) vs |R-R0|:", [(float(d), round(float(w)*HA2EV, 4)) for d, w in zip(dist[:12], wt[:12])], "...", flush=True)

    grids = [int(x) for x in args.grids.split(",")]
    etas = [float(x) for x in args.etas.split(",")]
    rcuts = [float(x) for x in args.rcut.split(",")]
    k_int = lt.mp_grid(args.nk_int, args.nk_int, 1)

    # Dirac point from the Wannier bands (min pi/pi* gap on a fine grid); energies from _tb.dat are eV
    _, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=ndegen)
    gap = E_ref[:, 4] - E_ref[:, 3]
    iD = int(np.argmin(gap)); E_dirac = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
    win = (E_dirac - 3.0, E_dirac + 3.0)
    print(f"[dirac] E_Dirac = {E_dirac:.4f} eV (Wannier, min gap {gap[iD]*1e3:.1f} meV); window {win}", flush=True)

    print(f"\n{'Rcut':>5} {'grid':>6} {'eta(eV)':>9} {'median G*Ncells(meV)':>22} {'resonance E-ED(eV)':>19}")
    results = {}
    for rc in rcuts:
        Rloc = R_mwr[np.linalg.norm(R_mwr, axis=1) <= rc + 1e-9]
        V_loc, _ = lt.extract_V_loc(Mwr, R_mwr, Rloc)                # guardrail b
        for N in grids:
            k_out = lt.mp_grid(N, N, 1)
            _, E_out, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)
            for eta in etas:
                gamma = lt.scattering_rate_fast(Hwr, Rw, ndegen, V_loc, Rloc, k_out, eta,
                                                k_int=k_int, e_window=win, ne_per_eta=8)
                med = float(np.nanmedian(np.abs(gamma))) * 1e3
                # resonance: energy (rel. Dirac) of max on-shell rate within +-1.5 eV of Dirac
                E = E_out.T                                           # (nw, nk) eV
                m = np.isfinite(gamma) & (np.abs(E - E_dirac) <= 1.5)
                e_res = float(E[m][np.argmax(np.abs(gamma)[m])] - E_dirac) if m.any() else np.nan
                results[(rc, N, eta)] = (med, e_res)
                print(f"{rc:>5.0f} {N:>6} {eta:>9.3f} {med:>22.2f} {e_res:>19.3f}", flush=True)

    # joint plateau: region stable both as eta decreases AND grid refines (<=5%)
    print("\n[Level 1] joint plateau = a (Rcut, grid, eta) box where median G*Ncells varies <= 5% "
          "when eta halves AND grid doubles. Inspect the map above.")
    if args.out:
        np.savez(args.out, results=np.array([(rc, N, e, m, r) for (rc, N, e), (m, r) in results.items()]))
        print(f"saved {args.out}")


if __name__ == "__main__":
    main()
