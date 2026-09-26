#!/usr/bin/env python
"""
r6_states_identity.py -- R6 étape 1.2, 9x9 : porte A.2 sur états purs et sur les 30 états QE de R5 avec les M2 de results/M2/
(grossier 16 bandes : M_{L,NL,ed}_9x9 ; dense 20 bandes restreint aux 81 k : M_{L,NL}_dense_9x9, M_dense_9x9).
Identité (R5-A.2) : c†H_in c = (e_QE - Lu) c†c - <Psi_in|dV|Psi_out> + s c†c, s attendu = -décalage Lu (+24,65 meV).
Réutilise le pilote R5 (load_bases, qe_states_9x9, deltav_grid_and_projector, A2_states_loop) ; journal r6_log.txt ; sorties states/.
"""
import os
import sys
import json
import time

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
R5DIR = os.path.join(GQ, "defects", "R5_base_vs_M"); R4DIR = os.path.join(GQ, "defects", "R4_quasi_lie")
sys.path.insert(0, os.path.join(PROJ, "src")); sys.path.insert(0, R4DIR); sys.path.insert(0, R5DIR)
import r5_driver as r5  # noqa: E402  (chdir PROJ, sys.path src)
from electron_defect_interaction.io import matrix_io  # noqa: E402
from electron_defect_interaction.config import HA2EV  # noqa: E402
from electron_defect_interaction.defects import deltav_pw as dv  # noqa: E402  (module promu = r5_deltav_pw)
import r5_sc_projection as sp  # noqa: E402

M2 = os.path.join(PROJ, "results", "M2"); OUT = os.path.join(WORK, "states"); os.makedirs(OUT, exist_ok=True)
NT = r5.NT; N_CELLS = r5.N_CELLS


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"; print(line, flush=True)
    with open(os.path.join(WORK, "r6_log.txt"), "a") as f:
        f.write(line + "\n")


r5.log = log; r5.r4.log = log


def v2(path):
    meta = matrix_io.read_manifest(path)
    if meta is None or meta.get("bloch_norm") != matrix_io.UNIT_CELL or meta.get("units") != matrix_io.HARTREE or not str(meta.get("M_normalization", "")).startswith("v2"):
        raise ValueError(f"{path}: sidecar v2 unit_cell/hartree requis, trouvé {meta}")
    return path


def coarse_M2(part):
    return np.load(v2(os.path.join(M2, {"tot": "M_ed_9x9.npy", "L": "M_L_9x9.npy", "NL": "M_NL_9x9.npy"}[part]))) * HA2EV


def dense_M2(part, idx81):
    p = v2(os.path.join(M2, {"tot": "M_dense_9x9.npy", "L": "M_L_dense_9x9.npy", "NL": "M_NL_dense_9x9.npy"}[part]))
    Mm = np.load(p, mmap_mode="r"); nb = Mm.shape[0]
    M = np.array(Mm[np.ix_(np.arange(nb), idx81, np.arange(nb), idx81)]) * HA2EV; del Mm
    return M


def main():
    t_start = time.time(); out = dict(M2_dir=M2)
    B = r5.load_bases(with_dense=True); ng = B["ng"]; Om = B["Om"]
    Q = r5.qe_states_9x9(); ns = len(Q["e"]); hl = r5.highlighted(Q)
    log(f"[states] {ns} états QE (bandes {Q['bands1'][0]}-{Q['bands1'][-1]}), décalage Lu {Q['shift']*1e3:+.2f} meV, E_D(SC) {Q['E_D_SC']:.4f} eV")
    Ms = {16: {p: coarse_M2(p) for p in ("tot", "L", "NL")}, 20: {p: dense_M2(p, B["idx81"]) for p in ("tot", "L", "NL")}}
    out["M2_linearity"] = {str(nb): float(np.abs(Ms[nb]["tot"] - Ms[nb]["L"] - Ms[nb]["NL"]).max()) for nb in Ms}
    for nb in Ms:
        log(f"[states] M2 base {nb} : max|tot - L - NL| = {out['M2_linearity'][str(nb)]:.1e} eV ; max|L| {np.abs(Ms[nb]['L']).max():.4f}, max|NL| {np.abs(Ms[nb]['NL']).max():.4f} eV")
    # recouvrements c_nk
    c = {16: np.zeros((ns, 16, 81), complex), 20: np.zeros((ns, 20, 81), complex)}; t0 = time.time()
    for s in range(ns):
        A = sp.sc_state_grid(Q["Cw"][s], Q["mill"], ng, gamma_only=Q["gamma_only"])
        for nb in (16, 20):
            b = B["bases"][nb]; c[nb][s] = sp.bloch_overlaps(A, b["C"], b["nG"], b["flat"])
    log(f"[states] recouvrements des {ns} états sur les bases 16 et 20 en {time.time()-t0:.0f} s")
    dV, proj, flat_full, mill_full, geo = r5.deltav_grid_and_projector(B, Q["mill"]); out["deltaV"] = geo
    A2 = dict(states={}, states_Lx81={}); out["pure"] = {}; out["states"] = {}
    for nb in (16, 20):
        b = B["bases"][nb]; iK, iKp = b["iK"], b["iKp"]
        pairs = [(3, iK, 3, iK), (3, iK, 4, iK), (3, iK, 3, iKp), (0, 0, 0, 0), (2, 5, 7, 40), (10, 17, 1, 60)] + ([(17, iK, 19, iKp)] if nb == 20 else [])
        t0 = time.time(); res = dv.check_pure_bloch(pairs, b["C"], b["nG"], b["flat"], ng, Om, dV, Ms[nb]["L"], Ms[nb]["NL"], N_CELLS, proj, flat_full, workers=NT)
        dL = max(r["dL"] for r in res); dNL = max(r["dNL"] for r in res); ok = bool(max(dL, dNL) <= r5.TOL_PURE)
        ratios = [abs(r["VL_direct"]) / abs(r["ML_over_N"]) for r in res if abs(r["ML_over_N"]) > 1e-4]
        out["pure"][str(nb)] = dict(ok=ok, max_dL=float(dL), max_dNL=float(dNL), tol=r5.TOL_PURE, n_pairs=len(pairs), ratio_L=[float(x) for x in ratios], seconds=time.time() - t0,
                                    pairs=[dict(pair=list(r["pair"]), VL_direct=r["VL_direct"].real, ML_over_N=r["ML_over_N"].real, dL=r["dL"], VNL_direct=r["VNL_direct"].real, MNL_over_N=r["MNL_over_N"].real, dNL=r["dNL"]) for r in res])
        log(f"[A.2 M2] base {nb} états purs : max|direct - M2/81| = {max(dL, dNL):.2e} eV sur {len(pairs)} paires (L {dL:.1e}, NL {dNL:.1e}) -> {'OK' if ok else 'REFUSÉ'} ; rapport L direct/(M2_L/81) : {[round(x, 6) for x in ratios]}")
        if not ok:
            continue
        rows = []; Mf = {p: Ms[nb][p].reshape(nb * 81, nb * 81) for p in ("tot", "L", "NL")}
        r5.A2_states_loop("M2", nb, b, Mf, rows, c, Q, ns, ng, Om, dV, proj, flat_full, hl, A2, geo)
        out["states"][str(nb)] = rows
    with open(os.path.join(OUT, "states_results.json"), "w") as f:
        json.dump(r5.r4.jsonable(out), f, indent=1, ensure_ascii=False)
    # tableau markdown (états en évidence + agrégats)
    L = [f"# R6 1.2 — 9x9 : porte A.2 (états purs) et identité des 30 états QE avec M2 ({time.strftime('%Y-%m-%d %H:%M')})\n"]
    for nb in (16, 20):
        p = out["pure"][str(nb)]
        L.append(f"\nBase {nb} : états purs max|Δ| = {max(p['max_dL'], p['max_dNL']):.2e} eV (L {p['max_dL']:.1e}, NL {p['max_dNL']:.1e} ; seuil {p['tol']:.0e}) → {'OK' if p['ok'] else 'REFUSÉ'} ; rapport L direct/(M2_L/81) : {', '.join(f'{x:.6f}' for x in p['ratio_L'])}.\n")
        if str(nb) not in out["states"]:
            continue
        rows = out["states"][str(nb)]; sal = np.array([r["s_align_eV"] for r in rows])
        L.append(f"\n30 états, base {nb} : max|(i) c†(M2/81)c − (ii) direct| in-in L {max(r['dL_in_in'] for r in rows):.1e}, NL {max(r['dNL_in_in'] for r in rows):.1e} eV ; "
                 f"s_align médiane {np.median(sal)*1e3:+.2f} meV (min {sal.min()*1e3:+.2f}, max {sal.max()*1e3:+.2f}) ; décalage Lu {Q['shift']*1e3:+.2f} meV.\n")
        L.append("\n| bande | parité | ε − E_D QE | Σ\\|c\\|² | c†(M2/81)c tot (L ; NL) | ⟨in\\|ΔV^L\\|out⟩ | ⟨in\\|ΔV^NL\\|out⟩ | c†H_in c | (e_QE − Lu) c†c − ⟨in\\|ΔV\\|out⟩ | s (meV) |\n|---|---|---|---|---|---|---|---|---|---|")
        for lab in hl:
            r = rows[hl[lab]]
            L.append(f"| {r['band']} | {r['parity']} | {r['x_qe']:+.3f} | {r['norm_in']:.4f} | {r['cMc_tot']:+.4f} ({r['cMc_L']:+.4f} ; {r['cMc_NL']:+.4f}) | {r['VL_in_out']:+.4f} | {r['VNL_in_out']:+.4f} | {r['lhs_cHc']:+.4f} | {r['rhs']:+.4f} | {r['s_align_eV']*1e3:+.2f} |")
    with open(os.path.join(OUT, "states_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log(f"[states] terminé en {time.time()-t_start:.0f} s ; states/states_results.json, states/states_tables.md")


if __name__ == "__main__":
    main()
