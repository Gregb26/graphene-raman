#!/usr/bin/env python
"""
r6_compare_v1_v2.py -- R6 étape 3.7 : table de correspondance ancien (results/M, v1 gelé) -> nouveau (results_dir = results/M2) pour les chiffres
cités dans NOTES_TGAMMA §2 et dans les tableaux/figures du mémoire (tab:rcut_M, tab:L_NL, tab:tests_M, fig_spectral, fig_plateau, fig_level2,
fig_locality, fig_rcut, fig_M_map, fig_M_scaling, fig_epw_vs_ed). Chiffres bruts ; « — » si le fichier v2 n'existe pas encore.
Écrit etape3/table_v1_v2.md.
"""
import os, sys, csv, json
import numpy as np
WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src")); os.chdir(PROJ)
from electron_defect_interaction.config import load_production, results_dir
cfg = load_production(verbose=False); V2 = results_dir(cfg); V1 = os.path.join(PROJ, cfg.get("results_dir_frozen", "results/M"))
RC, G, ETA, REF = cfg["R_cut"], cfg["grid"], cfg["eta_eV"], cfg["reference_size"]
rows = []
def add(q, a, b, src, where): rows.append((q, a, b, src, where))
def fmt(x, nd=3):
    if x is None: return "—"
    if isinstance(x, str): return x
    if isinstance(x, complex): return f"{x.real:+.{nd}f}{x.imag:+.{nd}f}i"
    try: return f"{float(x):.{nd}f}" if abs(float(x)) < 1e4 else f"{float(x):.4e}"
    except Exception: return str(x)
def csvrows(d, name):
    p = os.path.join(d, name)
    return list(csv.DictReader(open(p))) if os.path.exists(p) else None
def npz(d, name):
    p = os.path.join(d, name); return np.load(p, allow_pickle=True) if os.path.exists(p) else None
def both(fn, name, *a):
    r1 = fn(V1, name, *a); r2 = fn(V2, name, *a); return r1, r2

# ---- niveau 1 / niveau 2
def lvl1(d, S, rc, g, e):
    r = csvrows(d, "level1_summary.csv")
    if r is None: return None, None
    for x in r:
        if x["size"] == S and int(x["R_cut"]) == rc and int(x["grid"]) == g and abs(float(x["eta_eV"]) - e) < 1e-6:
            return float(x["median_Gamma_Ncells_meV"]), float(x["argmax_E_minus_ED_eV"])
    return None, None
for S in ("5x5", "6x6", "7x7", "8x8", "9x9", "12x12"):
    (m1, e1), (m2, e2) = lvl1(V1, S, RC, G, ETA), lvl1(V2, S, RC, G, ETA)
    add(f"médiane |Γ|·N_cells, {S}, R_cut {RC}, {G}², η {ETA} (meV)", fmt(m1, 2), fmt(m2, 2), "level1_summary.csv / level2_summary.csv", "fig_level2, fig_spectral (9×9), §2 NOTES_TGAMMA")
    add(f"E_res − E_D (argmax états ±1,5 eV), {S} (eV)", fmt(e1), fmt(e2), "level1_summary.csv", "fig_plateau, §2")
for (rc, g, e) in ((RC, G, 0.01), (RC, 120, ETA), (2, G, ETA), (4, G, ETA), (0, G, ETA)):
    (m1, e1), (m2, e2) = lvl1(V1, REF, rc, g, e), lvl1(V2, REF, rc, g, e)
    add(f"médiane |Γ|·N_cells, {REF}, R_cut {rc}, {g}², η {e} (meV)", fmt(m1, 2), fmt(m2, 2), "level1_summary.csv", "fig_plateau, fig_rcut, C10/C11")
# ---- niveau 2 familles
f1, f2 = both(csvrows, "level2_families.csv")
if f1:
    for x in f1:
        y = next((z for z in (f2 or []) if z["size"] == x["size"]), None)
        add(f"p_z–p_z sur le site de la lacune, {x['size']} dense (eV)", fmt(x["onsite_pz_vac_eV"]), fmt(y["onsite_pz_vac_eV"]) if y else "—", "level2_families.csv (mwr_locality.npz)", "fig_locality, tab familles")
        add(f"Re M^L / Re M^NL à K, paire π, {x['size']} (eV)", f"{x['ReML_K_eV']} / {x['ReMNL_K_eV']}", (f"{y['ReML_K_eV']} / {y['ReMNL_K_eV']}" if y else "—"), "level2_families.csv (M_analysis.npz)", "fig_M_map, §4.1.5")
# ---- L/NL Frobenius
l1, l2 = both(csvrows, "lnl_frobenius.csv")
if l1:
    for x in l1:
        y = next((z for z in (l2 or []) if z["size"] == x["size"]), None)
        add(f"⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, {x['size']} dense", fmt(x["ratio_full"]), fmt(y["ratio_full"]) if y else "—", "lnl_frobenius.csv", "tab:L_NL")
        add(f"⟨‖M^L‖_F⟩, {x['size']} (eV)", fmt(x["mean_fL_eV"], 4), fmt(y["mean_fL_eV"], 4) if y else "—", "lnl_frobenius.csv", "tab:L_NL")
# ---- localité
L1, L2 = both(npz, "mwr_locality.npz")
for S, tag in (("9x9", "9x9_dense"), ("12x12", "12x12_dense")):
    for key, lab in (("onsite_norm", "‖M_W(0,0)‖ (eV)"), ("onsite_pzvac", "p_z–p_z lacune (eV)")):
        a = float(L1[f"{tag}_{key}"]) if L1 is not None and f"{tag}_{key}" in L1.files else None
        b = float(L2[f"{tag}_{key}"]) if L2 is not None and f"{tag}_{key}" in L2.files else None
        add(f"localité M_W {S} dense : {lab}", fmt(a), fmt(b), "mwr_locality.npz", "fig_locality, §2")
    for Z, nm in ((L1, "v1"), (L2, "v2")):
        pass
    if L1 is not None and f"{tag}_w" in L1.files:
        w1 = L1[f"{tag}_w"]; w2 = L2[f"{tag}_w"] if (L2 is not None and f"{tag}_w" in L2.files) else None
        add(f"localité M_W {S} : ‖M_W(R,0)‖ première couronne (eV, 4 premières)", ", ".join(f"{v:.3f}" for v in w1[1:5]), (", ".join(f"{v:.3f}" for v in w2[1:5]) if w2 is not None else "—"), "mwr_locality.npz", "fig_locality")
# ---- R_cut (tab:rcut_M)
r1, r2 = both(csvrows, "m_rcut_convergence.csv")
if r1:
    for x in r1:
        y = next((z for z in (r2 or []) if z["size"] == x["size"] and z["R_cut"] == x["R_cut"]), None)
        add(f"tab:rcut_M {x['size']} R_cut {x['R_cut']} : max|ΔM|/max|M| (60²)", x["max_dM_over_maxM"], y["max_dM_over_maxM"] if y else "—", "m_rcut_convergence.csv", "tab:rcut_M")
# ---- N_k^int
n1, n2 = both(csvrows, "nkint_check_9x9.csv")
if n1:
    for x in n1:
        y = next((z for z in (n2 or []) if z["nk_int"] == x["nk_int"]), None)
        add(f"N_k^int {x['nk_int']} : médiane Γ (meV) / Re Σ (meV) / E_res / Γ(E_D)", f"{float(x['medG']):.2f} / {float(x['medRe']):.2f} / {float(x['E_res']):+.3f} / {float(x['G_ED']):.2f}",
            (f"{float(y['medG']):.2f} / {float(y['medRe']):.2f} / {float(y['E_res']):+.3f} / {float(y['G_ED']):.2f}" if y else "—"), "nkint_check_9x9.csv", "§6 NOTES_TGAMMA")
# ---- résonance 9x9
R1, R2 = both(npz, f"resonance_{REF}.npz")
def rstats(R):
    if R is None: return {}
    eg = R["eg"]; ED = float(R["E_D"]); m = np.abs(eg - ED) <= 3.0; GT = R["Gamma_T"]; GB = R["Gamma_Born"]
    d = dict(medGT=float(np.nanmedian(GT[m]) * 1e3), medGB=float(np.nanmedian(GB[m]) * 1e3), BoverT=float(np.nanmedian(GB[m] / GT[m])),
             BoverT_min=float(np.nanmin(GB[m] / GT[m])), BoverT_max=float(np.nanmax(GB[m] / GT[m])), MLmean=float(R["ML_diag_mean"]) * 1e3,
             ReT=float(R["ReTbar_at_ED"]), ImT=float(R["ImTbar_at_ED"]), zc=[float(v) for v in R["Tbar_zero_crossings"]])
    for k in ("peak_GT", "peak_GB", "peak_ratio", "peak_drho", "peak_rho_dis", "peak_absTbar", "peak_ImTbar", "min_absReTbar"):
        d[k] = float(R[k])
    for k in ("E_res_states", "median_GT_states_meV"):
        if k in R.files: d[k] = float(R[k])
    return d
s1, s2 = rstats(R1), rstats(R2)
for k, lab in (("medGT", "médiane de la courbe Γ_T sur ±3 eV (meV)"), ("medGB", "médiane Γ_Born (meV)"), ("BoverT", "Born/T médian"), ("BoverT_min", "Born/T min"), ("BoverT_max", "Born/T max"),
               ("MLmean", "moyenne diagonale de M^L (meV, C14 ancien)"), ("peak_GT", "pic Γ_T (eV)"), ("peak_GB", "pic Γ_Born"), ("peak_ratio", "pic Γ_T/ρ₀"), ("peak_drho", "pic δρ"), ("peak_rho_dis", "pic ρ_dis"),
               ("peak_absTbar", "pic |T̄|"), ("peak_ImTbar", "pic −Im T̄"), ("min_absReTbar", "min |Re T̄|"), ("ReT", "Re T̄(E_D) (eV)"), ("ImT", "Im T̄(E_D) (eV)"), ("E_res_states", "E_res états ±1,5 eV (resonance_metrics)"),
               ("median_GT_states_meV", "médiane |Γ_T| états (meV, resonance_metrics)")):
    add(f"résonance {REF} : {lab}", fmt(s1.get(k), 3 if k not in ("medGT", "medGB", "MLmean") else 2), fmt(s2.get(k), 3 if k not in ("medGT", "medGB", "MLmean") else 2), f"resonance_{REF}.npz", "fig_spectral, fig_epw_vs_ed (ch. 5), §2")
add(f"résonance {REF} : zéros de Re T̄", ", ".join(fmt(v) for v in s1.get("zc", [])) or "aucun", (", ".join(fmt(v) for v in s2.get("zc", [])) or "aucun") if s2 else "—", f"resonance_{REF}.npz", "§2")
# ---- critères 9x9
C1, C2 = both(npz, f"resonance_criteria_{REF}.npz")
def cstats(C, suf=""):
    if C is None or f"logdet_rel{suf}" not in C.files: return {}
    eg = C["eg"]; ED = float(C["E_D"]); ld = C[f"logdet_rel{suf}"]; ml = C[f"minlam{suf}"]; lm = C[f"lam_min{suf}"]
    jd = int(np.argmin(ld)); jl = int(np.argmin(ml))
    d = dict(det=float(np.exp(ld[jd])), det_at=float(eg[jd] - ED), lam=complex(lm[jl]), lam_at=float(eg[jl] - ED),
             sr=float(C[f"sumrule{suf}"]), srl=float(C[f"sumrule_lloyd{suf}"]), srw=float(C[f"sumrule_window{suf}"]))
    if suf == "" and "Gamma_c" in C.files:
        x = C["x_c"]; g = C["Gamma_c"] * 1e3; d.update(Gc_min=float(np.nanmin(g)), Gc_min_at=float(x[np.nanargmin(g)]), Gc_max=float(np.nanmax(g)), Gc_max_at=float(x[np.nanargmax(g)]),
                                                     Gc_ED=float(np.interp(0, x, g)), tau_m=658.2 / float(np.interp(-0.3, x, g)), tau_p=658.2 / float(np.interp(0.3, x, g)))
    return d
for suf, blab in (("", "matrice complète"), ("_pi", "bloc π"), ("_sigma", "bloc σ")):
    c1, c2 = cstats(C1, suf), cstats(C2, suf)
    add(f"critère {REF} ({blab}) : min |det|/max (ε − E_D)", (f"{c1['det']:.3e} ({c1['det_at']:+.3f})" if c1 else "—"), (f"{c2['det']:.3e} ({c2['det_at']:+.3f})" if c2 else "—"), f"resonance_criteria_{REF}.npz", "§2, C16")
    add(f"critère {REF} ({blab}) : λ_min (ε − E_D)", (f"{fmt(c1['lam'], 4)} ({c1['lam_at']:+.3f})" if c1 else "—"), (f"{fmt(c2['lam'], 4)} ({c2['lam_at']:+.3f})" if c2 else "—"), f"resonance_criteria_{REF}.npz", "§2, C16")
    add(f"Friedel {REF} ({blab}) : ∫δρ bande / Lloyd / ±3 eV (états)", (f"{c1['sr']:+.4f} / {c1['srl']:+.4f} / {c1['srw']:+.3f}" if c1 else "—"), (f"{c2['sr']:+.4f} / {c2['srl']:+.4f} / {c2['srw']:+.3f}" if c2 else "—"), f"resonance_criteria_{REF}.npz", "§2, C15")
c1, c2 = cstats(C1), cstats(C2)
if c1.get("Gc_min") is not None or c2.get("Gc_min") is not None:
    add(f"Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV)", (f"{c1['Gc_min']:.2f} ({c1['Gc_min_at']:+.2f}) / {c1['Gc_max']:.2f} ({c1['Gc_max_at']:+.2f}) / {c1['Gc_ED']:.2f}" if c1.get("Gc_min") is not None else "—"),
        (f"{c2['Gc_min']:.2f} ({c2['Gc_min_at']:+.2f}) / {c2['Gc_max']:.2f} ({c2['Gc_max_at']:+.2f}) / {c2['Gc_ED']:.2f}" if c2.get("Gc_min") is not None else "—"), f"resonance_criteria_{REF}.npz", "§2 (Kaasbjerg)")
    add(f"ħ/Γ à ∓0,3 eV, c = 0,1 % (fs)", (f"{c1['tau_m']:.0f} / {c1['tau_p']:.0f}" if c1.get("Gc_min") is not None else "—"), (f"{c2['tau_m']:.0f} / {c2['tau_p']:.0f}" if c2.get("Gc_min") is not None else "—"), f"resonance_criteria_{REF}.npz", "§2")
# ---- tests de la chaîne
t1, t2 = both(csvrows, "M_tests_summary.csv")
if t1:
    for x in t1:
        y = next((z for z in (t2 or []) if z["test"] == x["test"]), None)
        add(f"tab:tests_M : {x['test']}", f"{x['valeur']} ({x['verdict']})", (f"{y['valeur']} ({y['verdict']})" if y else "—"), "M_tests_summary.csv", "tab:tests_M")
    for z in (t2 or []):
        if not any(x["test"] == z["test"] for x in t1): add(f"tab:tests_M : {z['test']}", "—", f"{z['valeur']} ({z['verdict']})", "M_tests_summary.csv", "tab:tests_M (nouvelle ligne)")
# ---- C14 redéfini (v2 seulement)
Rs = npz(V2, f"resonance_{REF}_shiftL.npz")
if Rs is not None:
    for k in sorted(Rs.files):
        if k.startswith(("E_res_states", "peak_GT_shift", "median_GT_states")) and "shift" in k or k in ("E_res_states", "median_GT_states_meV"):
            add(f"C14 redéfini (V_loc + C·1) : {k}", "—", fmt(float(Rs[k]), 3), f"resonance_{REF}_shiftL.npz", "3.5")
# ---- chapitre 5 : Γ^ed/Γ^ep (epw_ed_vs_ep.py, chaîne mv0.02 ; v1 dans results/epw/, gelé ; v2 dans results_dir)
E1, E2 = npz(os.path.join(PROJ, "results", "epw"), "ed_vs_ep_24k24q_mv0.02.npz"), npz(V2, "ed_vs_ep_24k24q_mv0.02.npz")
if E1 is not None and E2 is not None:
    w12 = lambda E, k: float(np.median(E[k][np.abs(E["x"]) <= 1.2])) * 1e3
    add("Γ^ed/Γ^ep (c = 1 %, 300 K, ±3 eV) : médiane", fmt(float(E1["median"]), 3), fmt(float(E2["median"]), 3), "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed, NOTES_EPW (ch. 5)")
    add("Γ^ed/Γ^ep : min (à ε − E_D, eV)", f"{float(E1['min']):.3f} ({float(E1['x_min']):+.3f})", f"{float(E2['min']):.3f} ({float(E2['x_min']):+.3f})", "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
    add("Γ^ed/Γ^ep : max (à ε − E_D, eV)", f"{float(E1['max']):.3f} ({float(E1['x_max']):+.3f})", f"{float(E2['max']):.3f} ({float(E2['x_max']):+.3f})", "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
    add("Γ^ed/Γ^ep : croisements (ratio = 1, eV)", ", ".join(f"{v:+.3f}" for v in E1["crossings"]), ", ".join(f"{v:+.3f}" for v in E2["crossings"]), "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
    add("médiane Γ^ed, c = 1 % (meV) : ±3 eV / ±1,2 eV", f"{float(E1['median_ed'])*1e3:.3f} / {w12(E1, 'Gamma_ed'):.3f}", f"{float(E2['median_ed'])*1e3:.3f} / {w12(E2, 'Gamma_ed'):.3f}", "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed (ch. 5)")
    add("médiane Γ^ep, 300 K (meV) : ±3 eV / ±1,2 eV (indépendant de M)", f"{float(E1['median_ep'])*1e3:.3f} / {w12(E1, 'Gamma_ep'):.3f}", f"{float(E2['median_ep'])*1e3:.3f} / {w12(E2, 'Gamma_ep'):.3f}", "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed (ch. 5)")
L = [f"# R6 3.7 — table de correspondance ancien (v1, `{os.path.relpath(V1, PROJ)}`) → nouveau (v2, `{os.path.relpath(V2, PROJ)}`)\n",
     "| grandeur | v1 | v2 (M2) | fichier source | figure / tableau / note |", "|---|---|---|---|---|"]
esc = lambda t: str(t).replace("|", "\\|")
L += [f"| {esc(q)} | {esc(a)} | {esc(b)} | `{s}` | {esc(w)} |" for q, a, b, s, w in rows]
os.makedirs(os.path.join(WORK, "etape3"), exist_ok=True)
open(os.path.join(WORK, "etape3", "table_v1_v2.md"), "w").write("\n".join(L) + "\n"); print(f"{len(rows)} lignes -> etape3/table_v1_v2.md")
