#!/usr/bin/env python
"""
r7_driver.py -- pilote de la campagne R7 (super-cellules 15x15 et 18x18 ; famille 3m a cinq points).

Sous-commandes : d1 [--sizes 6,9,12,15,18] [--out d1] [--no-gate]   (protocole D1 de R4/R5, route « quadruplet » seulement)
                 tables [--out d1]                                     (tables + figures depuis le json)
D1 par taille N (parfaite p et lacune d, .save de production en lecture seule) :
  E_D = moyenne des 4 etats de la parfaite les plus proches de E_F (r5_alignment_ext.dirac_quadruplet) ; alignement Lu :
  alignment.far_atom_alignment (spheres 0,5 et 1,0 A autour de l'atome le plus loin de la lacune ; decalage a 1,0 A applique au
  spectre d) ; fenetre [-3, +1] eV ; r5_driver.window_states_gamma (qe_gamma_io : parite miroir, w2 disque 2 A, w1 disque 1 A) ;
  seuil 3 <w2> de la parfaite de meme N ; pi = etat impair localise de plus grand w2 ; doublet sigma = deux etats pairs de plus
  grand w2 ; comptages. Porte de regression : N = 6, 9, 12 doivent reproduire R5 c/c_results.json a 1e-6 (sinon exit 3).
  Ajustements eps(N) = eps_inf + a/N^p (p = 1, 2) sur les tailles presentes (r5_alignment_ext.size_fits) ; figures.
Sorties : <out>/d1_results.json, <out>/window_NxN.npz (extrait de la fenetre : valeurs propres, alignement, parite, w2, w1),
<out>/D1_tables.md, fig/ ; journal r7_log.txt. Modules R4 (src/, commis 75c8656) et pilotes R4/R5 reutilises sans modification
(leurs journaux sont rediriges ici). Aucun fichier de production modifie.
"""
import os
import sys
import re
import json
import time
import argparse

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                    # .../graphene/qe
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
R4DIR = os.path.join(GQ, "defects", "R4_quasi_lie"); R5DIR = os.path.join(GQ, "defects", "R5_base_vs_M")
sys.path.insert(0, os.path.join(PROJ, "src")); sys.path.insert(0, R4DIR); sys.path.insert(0, R5DIR); sys.path.insert(0, WORK)
os.chdir(PROJ)

from electron_defect_interaction.config import HA2EV
from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io import qe_gamma_io as qg
from electron_defect_interaction.defects import alignment as al
import r4_driver as r4
import r5_driver as r5
import r5_alignment_ext as ax


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "r7_log.txt"), "a") as f:
        f.write(line + "\n")


r4.log = log; r5.log = log                        # les aides R4/R5 ecrivent dans r7_log.txt, jamais dans leurs journaux

BOHR = r4.BOHR
R_W2, R_W1 = r4.R_W2, r4.R_W1
RADII = r4.RADII_A
WIN = r5.WIN
R5_C = os.path.join(R5DIR, "c", "c_results.json")
GATE_TOL = 1e-6
SC = os.path.join(GQ, "defects", "super_cell")


def ensure(sub):
    d = os.path.join(WORK, sub); os.makedirs(d, exist_ok=True); return d


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(r4.jsonable(d), f, indent=1, ensure_ascii=False)
    log(f"saved {os.path.relpath(path, WORK)}")


def savefig(fig, name):
    d = ensure("fig")
    fig.savefig(os.path.join(d, name + ".pdf")); fig.savefig(os.path.join(d, name + ".png"), dpi=200)
    log(f"[fig] {name}.pdf/.png")


def scf_info(path):
    """Mur (s), iterations, E_F (eV), RAM estimee par rang (MB), 2D cutoff, JOB DONE -- depuis scf.out."""
    txt = open(path, errors="ignore").read()
    it = re.search(r"convergence has been achieved in\s+(\d+) iterations", txt)
    ef = re.search(r"the Fermi energy is\s+([-\d.]+) ev", txt)
    ram = re.search(r"Estimated max dynamical RAM per process >\s+([\d.]+) (\w+)", txt)
    ramv = float(ram.group(1)) * (1024 if ram.group(2) == "GB" else 1) if ram else None
    return dict(wall_s=r5.scf_wall(path), iterations=int(it.group(1)) if it else None, E_F_out=float(ef.group(1)) if ef else None,
                ram_per_rank_MB=ramv, cutoff_2d=("running with the 2D cutoff" in txt), job_done=("JOB DONE" in txt),
                pwscf_version=(re.search(r"Program PWSCF v\.(\S+)", txt) or [None, None])[1])


def analyze_size(N, out_dir):
    P = r5.sc_paths(N); t0 = time.time()
    for k in ("d", "p"):
        assert os.path.isdir(P[k]), P[k]
    for k in ("vd", "vp", "scf_d", "scf_p"):
        assert os.path.isfile(P[k]), P[k]
    r = dict(N=N, c=1.0 / (2 * N * N), paths={k: P[k] for k in ("d", "p", "vd", "vp")}, scf_d=scf_info(P["scf_d"]), scf_p=scf_info(P["scf_p"]))
    A_b, Om = qe_io.get_A_volume(P["p"]); A_A = A_b * BOHR; x_p = qe_io.get_x_red(P["p"]); x_d = qe_io.get_x_red(P["d"])
    ng = tuple(int(v) for v in qe_io.get_ngfft(P["d"])); assert ng == tuple(int(v) for v in qe_io.get_ngfft(P["p"])), (N, ng)
    s_vac, i_vac_p, _ = al.vacancy_site(x_p, x_d, A_A); z0 = float(np.mean(x_d[:, 2]))
    dd = al.min_image_dist(x_d, s_vac, A_A); nn1 = np.where(dd < 1.8)[0]
    eP, efP, _ = qg.get_eigenvalues_spin(P["p"]); eP = eP[0] * HA2EV; efP *= HA2EV
    eD, efD, _ = qg.get_eigenvalues_spin(P["d"]); eD = eD[0] * HA2EV; efD *= HA2EV
    import xml.etree.ElementTree as ET
    nel_p = float(ET.parse(os.path.join(P["p"], "data-file-schema.xml")).getroot().find(".//output/band_structure/nelec").text)
    nel_d = float(ET.parse(os.path.join(P["d"], "data-file-schema.xml")).getroot().find(".//output/band_structure/nelec").text)
    r.update(nat_d=int(len(x_d)), nat_p=int(len(x_p)), nbnd_p=int(len(eP)), nbnd_d=int(len(eD)), nelec_p=nel_p, nelec_d=nel_d,
             E_F_p=float(efP), E_F_d=float(efD), s_vac=s_vac.tolist(), i_vac_p_1based=i_vac_p + 1, ngfft=list(ng), z0_red=z0,
             A_sc_A2=float(abs(np.linalg.det(A_A[:2, :2]))), first_neighbours_1based=[int(i) + 1 for i in nn1], d_first_A=[float(dd[i]) for i in nn1],
             top_band_minus_EF=dict(p=float(eP.max() - efP), d=float(eD.max() - efD)))
    # E_D : quadruplet de la parfaite (N = 3m)
    E_D, qidx, spread = ax.dirac_quadruplet(eP, efP)
    n_occ = int(round(nel_p / 2)); eS = np.sort(eP)
    r["E_D_quadruplet"] = dict(E_D=E_D, bands_1based=[i + 1 for i in qidx], quadruplet=sorted(float(eP[i]) for i in qidx), spread_eV=spread, E_D_minus_EF=E_D - efP,
                               homo=float(eS[n_occ - 1]), lumo=float(eS[n_occ]), n_occ=n_occ, n_states_within_1eV_above=int(((eP > E_D) & (eP <= E_D + 1.0)).sum()),
                               top_band_minus_ED=float(eP.max() - E_D))
    r["E_D"] = E_D; r["E_D_source"] = "quadruplet"
    if N % 3 != 0 or spread > 1e-3:
        log(f"[D1] N={N}: ATTENTION quadruplet non degenere (etalement {spread:.2e} eV) ou N != 3m")
    # alignement Lu
    V_P = r4.load_pot_eV(P["vp"]); V_D = r4.load_pot_eV(P["vd"]); assert V_P.shape == ng and V_D.shape == ng, (V_P.shape, V_D.shape, ng)
    ra = al.far_atom_alignment(V_D, V_P, x_d, x_p, A_A, RADII); shift = float(ra["shifts"][1.0])
    r["alignment"] = dict(i_far_1based=ra["i_far_d"] + 1, dist_far_A=ra["dist_far"], shift_05=float(ra["shifts"][0.5]), shift_10=shift,
                          means={str(k): [float(v[0]), float(v[1]), int(v[2]), int(v[3])] for k, v in ra["means"].items()},
                          mean3d_diff_meV=float((V_D.mean() - V_P.mean()) * 1e3),
                          dV_site_eV=float((V_D - V_P)[tuple(np.rint(s_vac * np.array(ng)).astype(int) % np.array(ng))]), radius_used_A=1.0)
    del V_P, V_D
    # etats de la fenetre
    mask2 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W2); mask1 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W1)
    r["disc_area_fraction"] = dict(r2=float(mask2.mean()), r1=float(mask1.mean()))
    WP = r5.window_states_gamma(P["p"], E_D, 0.0, ng, mask2, mask1, z0=z0); WD = r5.window_states_gamma(P["d"], E_D, shift, ng, mask2, mask1, z0=z0)
    assert WP["n"] > 0 and WD["n"] > 0, (N, WP["n"], WD["n"])
    thr = 3.0 * float(WP["w2"].mean()); r["threshold"] = dict(thr=thr, mean_w2_P=float(WP["w2"].mean()), max_w2_P=float(WP["w2"].max()), rule="w2 > 3 <w2>(parfaite N x N)")
    for tag, W in (("P", WP), ("D", WD)):
        r[f"window_{tag}"] = dict(n=int(W["n"]), bands=[int(W["b0"]) + 1, int(W["b1"])], n_even=W["n_even"], n_odd=W["n_odd"], parity_dev=W["parity_dev"],
                                  x=[float(v) for v in W["x"]], parity=[float(v) for v in W["parity"]], w2=[float(v) for v in W["w2"]], w1=[float(v) for v in W["w1"]], read_s=W["read_s"])
    loc = [dict(band=int(WD["bands1"][j]), x=float(WD["x"][j]), parity=float(WD["parity"][j]), w2=float(WD["w2"][j]), w1=float(WD["w1"][j]),
                relEF=float(WD["e"][WD["b0"] + j] - efD)) for j in range(WD["n"]) if WD["w2"][j] > thr]
    r["localized"] = loc; odd = [q for q in loc if q["parity"] < 0]; even = [q for q in loc if q["parity"] > 0]
    r["n_localized"] = dict(even=len(even), odd=len(odd))
    r["pi_state"] = max(odd, key=lambda q: q["w2"]) if odd else None
    r["sigma_doublet"] = sorted(sorted(even, key=lambda q: -q["w2"])[:2], key=lambda q: q["x"]) if even else None
    if r["sigma_doublet"] is not None and len(r["sigma_doublet"]) == 2:
        r["sigma_splitting_meV"] = float((r["sigma_doublet"][1]["x"] - r["sigma_doublet"][0]["x"]) * 1e3)
    # extrait npz de la fenetre (manifeste)
    np.savez(os.path.join(out_dir, f"window_{N}x{N}.npz"), N=N, E_D=E_D, E_D_source="quadruplet", E_F_p=efP, E_F_d=efD, shift_Lu_1A=shift, shift_Lu_05A=float(ra["shifts"][0.5]),
             thr_w2=thr, s_vac=s_vac, i_vac_p_1based=i_vac_p + 1, ngfft=np.array(ng), A_A=A_A, window_eV=np.array(WIN), radii_A=np.array([R_W2, R_W1]),
             e_p_all=eP, e_d_all=eD, P_bands1=WP["bands1"], P_e=eP[WP["b0"]:WP["b1"]], P_x=WP["x"], P_parity=WP["parity"], P_w2=WP["w2"], P_w1=WP["w1"],
             D_bands1=WD["bands1"], D_e=eD[WD["b0"]:WD["b1"]], D_e_aligned=eD[WD["b0"]:WD["b1"]] - shift, D_x=WD["x"], D_parity=WD["parity"], D_w2=WD["w2"], D_w1=WD["w1"])
    r["seconds"] = time.time() - t0
    log(f"[D1] N={N}: nat {r['nat_p']}/{r['nat_d']}, nbnd {r['nbnd_p']}/{r['nbnd_d']} ; E_D {E_D:.5f} (etalement {spread:.1e}, E_F-E_D {efP-E_D:+.4f}) ; Lu {shift*1e3:+.2f} meV "
        f"(0,5 A {ra['shifts'][0.5]*1e3:+.2f}) ; <V_d>-<V_p> {r['alignment']['mean3d_diff_meV']:+.2f} meV ; fenetre P {WP['n']} ({WP['n_even']}s/{WP['n_odd']}p), D {WD['n']} "
        f"({WD['n_even']}s/{WD['n_odd']}p) ; parite max dev {max(WP['parity_dev'], WD['parity_dev']):.1e} ; seuil {thr:.4f} ; localises s {len(even)} / p {len(odd)} ; "
        f"pi {r['pi_state'] and (round(r['pi_state']['x'],4), round(r['pi_state']['w2'],3), round(r['pi_state']['w1'],3), r['pi_state']['band'])} ; "
        f"sigma {r['sigma_doublet'] and [(round(q['x'],4), round(q['w2'],3)) for q in r['sigma_doublet']]} ({time.time()-t0:.0f} s)")
    return r


def regression_gate(S):
    """N = 6, 9, 12 contre R5 C (memes fichiers, meme code) : ecart <= GATE_TOL sur les grandeurs rapportees."""
    if not os.path.isfile(R5_C):
        log(f"[gate] {R5_C} absent : porte non evaluee"); return None, []
    ref = json.load(open(R5_C))["sizes"]; rows = []; ok = True
    for N, b in S.items():
        a = ref.get(str(N))
        if a is None:
            continue
        sm = lambda q: float(np.mean([s["x"] for s in q["sigma_doublet"]])) if q.get("sigma_doublet") else np.nan
        sw = lambda q: float(np.mean([s["w2"] for s in q["sigma_doublet"]])) if q.get("sigma_doublet") else np.nan
        checks = [("E_D", a["E_D"], b["E_D"]), ("Lu_1.0", a["alignment"]["shift_10"], b["alignment"]["shift_10"]), ("Lu_0.5", a["alignment"]["shift_05"], b["alignment"]["shift_05"]),
                  ("thr", a["threshold"]["thr"], b["threshold"]["thr"]), ("window_P_n", a["window_P"]["n"], b["window_P"]["n"]), ("window_D_n", a["window_D"]["n"], b["window_D"]["n"]),
                  ("n_loc_even", a["n_localized"]["even"], b["n_localized"]["even"]), ("n_loc_odd", a["n_localized"]["odd"], b["n_localized"]["odd"]),
                  ("pi_x", a["pi_state"]["x"], b["pi_state"]["x"]), ("pi_w2", a["pi_state"]["w2"], b["pi_state"]["w2"]), ("pi_w1", a["pi_state"]["w1"], b["pi_state"]["w1"]),
                  ("pi_band", a["pi_state"]["band"], b["pi_state"]["band"]), ("sigma_x_mean", sm(a), sm(b)), ("sigma_w2_mean", sw(a), sw(b))]
        for name, va, vb in checks:
            d = abs(float(va) - float(vb)); good = d <= GATE_TOL; ok &= good
            rows.append(dict(N=N, quantity=name, R5=float(va), R7=float(vb), abs_diff=d, ok=good))
    log(f"[gate] regression 6/9/12 contre R5 C : {'PASS' if ok else 'FAIL'} (max ecart {max((q['abs_diff'] for q in rows), default=0):.2e}, tol {GATE_TOL:.0e})")
    return ok, rows


def two_term_fit(N, e):
    """eps = eps_inf + a/N + b/N^2 (3 parametres ; supplement, hors prompt)."""
    N = np.asarray(N, float); e = np.asarray(e, float)
    X = np.stack([np.ones_like(N), 1 / N, 1 / N ** 2], axis=1); coef, *_ = np.linalg.lstsq(X, e, rcond=None); res = e - X @ coef
    return dict(eps_inf=float(coef[0]), a=float(coef[1]), b=float(coef[2]), rms=float(np.sqrt(np.mean(res ** 2))), max_resid=float(np.abs(res).max()), n=int(len(N)))


def do_fits(S):
    fits = {}
    Ns_all = sorted(S)
    for what in ("pi", "sigma"):
        pts = [(N, S[N]["pi_state"]["x"]) for N in Ns_all if S[N].get("pi_state")] if what == "pi" else \
              [(N, float(np.mean([q["x"] for q in S[N]["sigma_doublet"]]))) for N in Ns_all if S[N].get("sigma_doublet")]
        for lab, sel in (("all", pts), ("6-9-12", [p for p in pts if p[0] in (6, 9, 12)])):
            if len(sel) >= 3:
                Nn = [p[0] for p in sel]; xx = [p[1] for p in sel]
                fits[f"{what}_{lab}"] = dict(N=Nn, x=xx, fits=ax.size_fits(Nn, xx), two_term=two_term_fit(Nn, xx) if len(sel) >= 4 else None)
                f = fits[f"{what}_{lab}"]["fits"]
                log(f"[fit] {what} N {Nn} x {[round(v, 4) for v in xx]} : 1/N eps_inf {f[1]['eps_inf']:+.4f} (rms {f[1]['rms']:.4f}) ; 1/N^2 eps_inf {f[2]['eps_inf']:+.4f} (rms {f[2]['rms']:.4f})"
                    + (f" ; 1/N + 1/N^2 eps_inf {fits[f'{what}_{lab}']['two_term']['eps_inf']:+.4f} (rms {fits[f'{what}_{lab}']['two_term']['rms']:.4f})" if fits[f"{what}_{lab}"]["two_term"] else ""))
    return fits


def _f(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:+.{nd}f}"


def write_tables(out, out_dir):
    S = {int(k): v for k, v in out["sizes"].items()}; L = []
    L.append("#### D1 — E_D, alignement et états localisés par taille (famille 3m ; ε − E_D en eV ; alignement Lu 1,0 Å ; seuil = 3 × ⟨w₂⟩ de la parfaite N×N)\n")
    L.append("| N | nat p/d | nbnd p/d | E_F p / d (eV) | E_D quadruplet (étalement ; E_F − E_D meV) | Lu 1,0 / 0,5 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ (meV) | seuil w₂ | état π (ε − E_D ; w₂ ; w₁ ; bande) | doublet σ (ε − E_D ; w₂) | localisés σ / π | fenêtre P (σ/π) ; D (σ/π) | mur scf d / p (s ; it) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for N in sorted(S):
        r = S[N]; q = r["E_D_quadruplet"]; pi = r["pi_state"]; sg = r["sigma_doublet"]; a = r["alignment"]
        L.append(f"| {N} | {r['nat_p']}/{r['nat_d']} | {r['nbnd_p']}/{r['nbnd_d']} | {r['E_F_p']:.4f} / {r['E_F_d']:.4f} | {q['E_D']:.5f} ({q['spread_eV']:.1e} ; {(r['E_F_p']-q['E_D'])*1e3:+.1f}) | "
                 f"{a['shift_10']*1e3:+.2f} / {a['shift_05']*1e3:+.2f} | {a['mean3d_diff_meV']:+.2f} | {r['threshold']['thr']:.4f} | "
                 + (f"{pi['x']:+.4f} ; {pi['w2']:.3f} ; {pi['w1']:.3f} ; {pi['band']}" if pi else "aucun") + " | " + (", ".join(f"{s['x']:+.4f} ; {s['w2']:.3f}" for s in sg) if sg else "aucun")
                 + f" | {r['n_localized']['even']} / {r['n_localized']['odd']} | {r['window_P']['n']} ({r['window_P']['n_even']}/{r['window_P']['n_odd']}) ; {r['window_D']['n']} ({r['window_D']['n_even']}/{r['window_D']['n_odd']}) | "
                 f"{_f(r['scf_d']['wall_s'], 0)} ({r['scf_d']['iterations']}) / {_f(r['scf_p']['wall_s'], 0)} ({r['scf_p']['iterations']}) |")
    L.append("\n#### D1 — tous les états localisés (w₂ > seuil) par taille\n")
    for N in sorted(S):
        L.append(f"- N = {N} : " + ("; ".join(f"b{q['band']} {'σ' if q['parity'] > 0 else 'π'} {q['x']:+.4f} (w₂ {q['w2']:.3f}, w₁ {q['w1']:.3f})" for q in S[N]["localized"]) or "aucun"))
    L.append("\n#### D1 — ajustements ε(N) = ε_∞ + a/N^p (eV ; π = état impair localisé de plus grand w₂ ; σ = moyenne du doublet)\n")
    L.append("| état | points | N | ε − E_D | p = 1 : ε_∞ ; a ; rms ; max résidu | p = 2 : ε_∞ ; a ; rms ; max résidu | supplément 1/N + 1/N² : ε_∞ ; a ; b ; rms |")
    L.append("|---|---|---|---|---|---|---|")
    for key, f in out["fits"].items():
        what, lab = key.split("_", 1); f1 = f["fits"].get(1, f["fits"].get("1")); f2 = f["fits"].get(2, f["fits"].get("2")); tt = f.get("two_term")
        L.append(f"| {what} | {lab} | {f['N']} | {', '.join(f'{x:+.4f}' for x in f['x'])} | {f1['eps_inf']:+.4f} ; {f1['a']:+.3f} ; {f1['rms']:.4f} ; {f1['max_resid']:.4f} | "
                 f"{f2['eps_inf']:+.4f} ; {f2['a']:+.3f} ; {f2['rms']:.4f} ; {f2['max_resid']:.4f} | " + (f"{tt['eps_inf']:+.4f} ; {tt['a']:+.3f} ; {tt['b']:+.2f} ; {tt['rms']:.4f}" if tt else "—") + " |")
    if out.get("gate_rows"):
        L.append(f"\n#### D1 — porte de régression contre R5 C (tolérance {GATE_TOL:.0e}) : **{'PASS' if out['gate_ok'] else 'FAIL'}**\n")
        L.append("| N | grandeur | R5 | R7 | écart |"); L.append("|---|---|---|---|---|")
        for q in out["gate_rows"]:
            L.append(f"| {q['N']} | {q['quantity']} | {q['R5']:.6f} | {q['R7']:.6f} | {q['abs_diff']:.1e}{'' if q['ok'] else ' **FAIL**'} |")
    L.append("\n#### D1 — contrôles des runs (scf.out)\n")
    L.append("| N | cellule | PWSCF | 2D cutoff | JOB DONE | itérations | mur (s) | E_F scf.out (eV) | RAM QE par rang (MB) | dernière bande − E_F (eV) | états dans ]E_D, E_D + 1 eV] (parfaite) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for N in sorted(S):
        r = S[N]
        for c, k in (("parfaite", "scf_p"), ("lacune", "scf_d")):
            s = r[k]
            L.append(f"| {N} | {c} | {s['pwscf_version']} | {s['cutoff_2d']} | {s['job_done']} | {s['iterations']} | {_f(s['wall_s'], 1)} | {_f(s['E_F_out'], 4)} | {_f(s['ram_per_rank_MB'], 0)} | "
                     f"{r['top_band_minus_EF']['p' if k == 'scf_p' else 'd']:+.3f} | {r['E_D_quadruplet']['n_states_within_1eV_above'] if k == 'scf_p' else '—'} |")
    with open(os.path.join(out_dir, "D1_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log(f"[tables] {os.path.relpath(os.path.join(out_dir, 'D1_tables.md'), WORK)}")


def figs(out, tag=""):
    plt, pal = r5.fig_style(); S = {int(k): v for k, v in out["sizes"].items()}; F = out["fits"]
    Nf = np.linspace(5.5, 200, 400)
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.2))
    for ax_, what, col, title in ((axs[0], "pi", pal.ORANGE, "(a) état π quasi-lié (impair)"), (axs[1], "sigma", pal.NAVY, "(b) doublet σ (moyenne, pair)")):
        f = F.get(f"{what}_all")
        if f is None:
            continue
        Nn = np.array(f["N"], float); xx = np.array(f["x"])
        ax_.scatter(1 / Nn, xx, color=col, s=20, zorder=3)
        for N, x in zip(Nn, xx):
            ax_.annotate(f"{int(N)}", (1 / N, x), textcoords="offset points", xytext=(4, 3), fontsize=6, color=pal.INK)
        f1 = f["fits"].get(1, f["fits"].get("1")); f2 = f["fits"].get(2, f["fits"].get("2"))
        ax_.plot(1 / Nf, f1["eps_inf"] + f1["a"] / Nf, "-", color=col, lw=0.9, label=f"$\\varepsilon_\\infty + a/N$ : $\\varepsilon_\\infty$ = {f1['eps_inf']:+.3f} eV (rms {f1['rms']:.3f})")
        ax_.plot(1 / Nf, f2["eps_inf"] + f2["a"] / Nf ** 2, "--", color=col, lw=0.9, label=f"$\\varepsilon_\\infty + a/N^2$ : $\\varepsilon_\\infty$ = {f2['eps_inf']:+.3f} eV (rms {f2['rms']:.3f})")
        ax_.axhline(0, color=pal.MUTED, lw=0.5); ax_.set_xlim(0, 0.19); ax_.set_xlabel("$1/N$"); ax_.set_title(title, fontsize=8); ax_.legend(fontsize=5.5, loc="best")
    axs[0].set_ylabel(r"$\varepsilon - E_D$ (eV)")
    fig.tight_layout(); savefig(fig, "size_3m" + tag)
    fig, ax_ = plt.subplots(figsize=(4.2, 3.4))
    for N in sorted(S):
        for q in S[N]["localized"]:
            ax_.scatter(1 / N, q["x"], s=6 + 80 * q["w2"], color=pal.NAVY if q["parity"] > 0 else pal.ORANGE, alpha=0.75, lw=0)
    ax_.axhline(0, color=pal.MUTED, lw=0.5); ax_.set_xlabel("$1/N$"); ax_.set_ylabel(r"$\varepsilon - E_D$ (eV)"); ax_.set_xlim(0, 0.19); ax_.set_ylim(-3, 1)
    ax_.set_title("états localisés ($w_2$ > seuil) : bleu σ (pair), orange π (impair) ; aire ∝ $w_2$", fontsize=6.5)
    fig.tight_layout(); savefig(fig, "localized_3m" + tag)


def cmd_d1(a):
    out_dir = ensure(a.out); sizes = [int(s) for s in a.sizes.split(",")]; t_start = time.time()
    log(f"[D1] tailles {sizes} ; fenetre {WIN} ; rayons w2/w1 {R_W2}/{R_W1} A ; alignement Lu rayons {RADII} ; sorties {os.path.relpath(out_dir, WORK)}")
    S = {}
    for N in sizes:
        S[N] = analyze_size(N, out_dir)
    out = dict(sizes={str(N): S[N] for N in S}, window=WIN, radii=dict(w2=R_W2, w1=R_W1, lu=RADII), fits=do_fits(S), timing_s=time.time() - t_start)
    if not a.no_gate:
        ok, rows = regression_gate(S); out["gate_ok"] = ok; out["gate_rows"] = rows
    save_json(os.path.join(out_dir, "d1_results.json"), out)
    write_tables(out, out_dir); figs(out, "" if a.out == "d1" else "_" + a.out)
    log(f"[D1] termine en {time.time()-t_start:.0f} s")
    if not a.no_gate and out.get("gate_ok") is False:
        log("[D1] PORTE DE REGRESSION EN ECHEC : exit 3"); sys.exit(3)


def cmd_tables(a):
    out_dir = ensure(a.out); out = json.load(open(os.path.join(out_dir, "d1_results.json")))
    write_tables(out, out_dir); figs(out, "" if a.out == "d1" else "_" + a.out)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cmd", choices=["d1", "tables"]); p.add_argument("--sizes", default="6,9,12,15,18"); p.add_argument("--out", default="d1"); p.add_argument("--no-gate", action="store_true")
    a = p.parse_args()
    {"d1": cmd_d1, "tables": cmd_tables}[a.cmd](a)


if __name__ == "__main__":
    main()
