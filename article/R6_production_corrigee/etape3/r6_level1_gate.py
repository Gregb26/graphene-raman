#!/usr/bin/env python
"""
r6_level1_gate.py -- R6 étape 3.3 : porte de niveau 1 sur les paramètres gelés (config/production.json : R_cut 3, grille 240, eta 0,02, N_k^int 300),
même critère que celui qui a figé la config le 2026-09-05 (NOTES_TGAMMA §3, C10 et C11 ; seuil 5 %) :
  C11 : à R_cut gelé, la médiane de |Gamma|·N_cells varie de <= 5 % quand la grille double (120 -> 240) ET quand eta est divisé par 2 (0,02 -> 0,01) ;
  C10 : à (grille, eta) gelées, <= 5 % quand R_cut passe de 3 à 4 (carte niveau 1 R_cut 0..4 et rcut_resigma) ;
  N_k^int (P13, §6c, même seuil 5 % sur la médiane de la fenêtre) : 300 -> 600 (et 150, 450 rapportés) ; Gamma_T(E_D) rapporté à part (v1 : +5,4 %).
Lit results_dir(config) : specwd_<S>_prod.npz (toutes tailles, verdict sur la taille de référence), resigma_9x9_rc*.npz, resigma_9x9_rc3_nk*.npz.
Écrit etape3/level1_gate.md ; code de sortie 3 si un paramètre gelé ne satisfait plus le critère (STOP avant 3.4).
"""
import os, sys, glob, json
import numpy as np
WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src")); os.chdir(PROJ)
from electron_defect_interaction.config import load_production, results_dir
cfg = load_production(verbose=False); RES = results_dir(cfg); RC, G, ETA, NKI, REF = cfg["R_cut"], cfg["grid"], cfg["eta_eV"], cfg["nk_int"], cfg["reference_size"]
THR = 0.05
def load_map(S):
    f = f"{RES}/specwd_{S}_prod.npz"
    if not os.path.exists(f): return None
    return {(int(x[0]), int(x[1]), round(float(x[2]), 4)): (float(x[3]), float(x[4])) for x in np.load(f)["results"]}
def rel(a, b): return abs(a - b) / abs(b)
L = ["# R6 3.3 — porte de niveau 1 sur les paramètres gelés (seuil 5 %, critère du 2026-09-05)\n",
     f"Config : R_cut {RC}, grille {G}, η {ETA} eV, N_k^int {NKI}, référence {REF}. Médianes de |Γ|·N_cells (meV) sur la fenêtre ±{cfg['e_window_eV']} eV.\n"]
fails = []; rows = []
for S in ("5x5", "6x6", "7x7", "8x8", "9x9", "12x12"):
    m = load_map(S)
    if m is None:
        rows.append(f"| {S} | carte absente | | | | | |")
        if S == REF: fails.append(f"{S}: carte de niveau 1 absente (specwd_{S}_prod.npz) : porte non évaluable")
        continue
    key = (RC, G, ETA); med, er = m[key]
    c11g = rel(med, m[(RC, G // 2, ETA)][0]) if (RC, G // 2, ETA) in m else float("nan")
    c11e = rel(med, m[(RC, G, round(ETA / 2, 4))][0]) if (RC, G, round(ETA / 2, 4)) in m else float("nan")
    c10 = rel(med, m[(RC + 1, G, ETA)][0]) if (RC + 1, G, ETA) in m else float("nan")
    ok = all(v <= THR for v in (c11g, c11e, c10) if np.isfinite(v))
    rows.append(f"| {S} | {med:.2f} | {er:+.3f} | {c11g*100:.2f} % | {c11e*100:.2f} % | {c10*100:.2f} % | {'OK' if ok else 'ÉCHEC'} |")
    if S == REF and not ok: fails.append(f"{S}: C11 grille {c11g*100:.2f} %, C11 eta {c11e*100:.2f} %, C10 R_cut {c10*100:.2f} %")
L += ["| taille | médiane (R_cut 3, 240², η 0,02) | E_res (eV) | C11 grille 120 → 240 | C11 η 0,02 → 0,01 | C10 R_cut 3 → 4 | verdict |", "|---|---|---|---|---|---|---|"] + rows
# N_k^int (référence) : fichiers resigma_<REF>_rc3_nk<N>.npz de rcut_resigma.py (Sigma_rc3 = <nk|t|nk>, Gamma = -2 Im Sigma, états de la fenêtre)
def med_gamma(z, rc):
    Sg = z[f"Sigma_rc{rc}"]; m = np.isfinite(Sg.real); G = -2 * Sg.imag; E = z["E_out"]; ED = float(z["E_D"])
    mm = m & (np.abs(E - ED) <= 1.5); e_res = float(E[mm][np.argmax(G[mm])] - ED) if mm.any() else float("nan")
    atK = m & (np.abs(E - ED) < 1e-6)
    return float(np.median(G[m]) * 1e3), e_res, (float(G[atK].mean() * 1e3) if atK.any() else float("nan")), int(m.sum())
nk = {}
for f in glob.glob(f"{RES}/resigma_{REF}_rc{RC}_nk*.npz"):
    n = int(os.path.basename(f).split("nk")[-1][:-4]); nk[n] = np.load(f)
if nk:
    ns = sorted(nk); refv = med_gamma(nk[max(ns)], RC); basev = med_gamma(nk[NKI], RC) if NKI in nk else None
    L.append(f"\nN_k^int ({REF}, R_cut {RC}, grille {G}, η {ETA}) : médiane de |Γ|·N_cells sur les états de la fenêtre ; référence N_k^int = {max(ns)}")
    L.append("| N_k^int | n états | médiane (meV) | écart vs max | écart vs 300 | E_res (eV) | Γ_T(E_D) (meV, états à K) | écart Γ(E_D) vs max |"); L.append("|---|---|---|---|---|---|---|---|")
    for n in ns:
        v = med_gamma(nk[n], RC)
        L.append(f"| {n} | {v[3]} | {v[0]:.2f} | {rel(v[0], refv[0])*100:+.2f} % | {(rel(v[0], basev[0])*100 if basev else float('nan')):+.2f} % | {v[1]:+.3f} | {v[2]:.2f} | {rel(v[2], refv[2])*100:+.2f} % |")
    if basev is not None and rel(basev[0], refv[0]) > THR: fails.append(f"N_k^int {NKI} vs {max(ns)} : médiane fenêtre {rel(basev[0], refv[0])*100:.2f} %")
    if basev is not None and np.isfinite(basev[2]) and rel(basev[2], refv[2]) > THR:
        L.append(f"\nNote : Γ_T(E_D) (états à K) {NKI} → {max(ns)} = {rel(basev[2], refv[2])*100:+.2f} % (> 5 %, comme en v1 : +5,4 % ; hors critère du gel, rapporté)")
else:
    L.append("\nN_k^int : aucun fichier resigma_*_nk*.npz (balayage non encore fait)"); fails.append("N_k^int : balayage absent, porte non évaluable")
# C10 par rcut_resigma (R_cut 0..4 à grille/eta gelées)
rz = {}
for f in glob.glob(f"{RES}/resigma_{REF}_rc[0-9]*.npz"):
    if "_nk" in os.path.basename(f): continue
    z = np.load(f)
    for k in z.files:
        if k.startswith("Sigma_rc"): rz[int(k[8:])] = med_gamma(z, int(k[8:]))
if rz:
    rcs = sorted(rz); L.append(f"\nR_cut par `rcut_resigma.py` ({REF}, grille {G}, η {ETA}, N_k^int {NKI}) :")
    L.append("| R_cut | médiane (meV) | écart vs R_cut max | E_res (eV) |"); L.append("|---|---|---|---|")
    for r in rcs: L.append(f"| {r} | {rz[r][0]:.2f} | {rel(rz[r][0], rz[max(rcs)][0])*100:+.2f} % | {rz[r][1]:+.3f} |")
    if RC in rz and (RC + 1) in rz and rel(rz[RC][0], rz[RC + 1][0]) > THR: fails.append(f"C10 (rcut_resigma) R_cut {RC} vs {RC+1} : {rel(rz[RC][0], rz[RC+1][0])*100:.2f} %")
L.append(f"\n**Verdict** : {'OK — paramètres gelés conformes au critère du 2026-09-05 avec M2 ; 3.4 peut démarrer' if not fails else 'ÉCHEC — STOP avant 3.4 : ' + ' ; '.join(fails)}")
os.makedirs(os.path.join(WORK, "etape3"), exist_ok=True)
open(os.path.join(WORK, "etape3", "level1_gate.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
sys.exit(3 if fails else 0)
