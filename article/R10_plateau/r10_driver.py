#!/usr/bin/env python
"""
r10_driver.py -- pilote unique de la campagne R10 (base cohérente pour les chiffres du chapitre 4). Lancé par submit_r10.sh (une tâche par job).

Sous-commandes (GO 1) :
  a0       porte A.0 : P1 (alignment.atom_sphere_shifts) sur 5, 6, 7, 8, 9, 12 ; moyenne du plateau (i) (1,0 Å, d >= 0,75 r_max) = C_i_eV de R9 au bit ;
           atomes vraiment les plus loin et leurs décalages = audit de l'image minimale (far_true_1based, far_true_meV à <= 1e-6 meV) ; P1 à l'atome de
           Lu = far_atom_alignment au bit et = Lu publié (R5 C). Échec : STOP sur A (code 3).
  a1       A.1 : P1 sur 10, 11, 15, 18, 21, 24, 27 (1,0 et 0,5 Å) : plateau (moyenne, rms, max|écart|, n), valeur « un atome » vraie (atomes à r_max en
           vraie image, moyenne en cas d'égalité), Lu publié en regard ; a/a1_results.json (13 tailles, C_N_eV = moyenne du plateau (i), D14),
           a/A1_tables.md, fig/offset_profiles_13.
  a2       A.2 : niveaux QE réalignés sans recalcul, x_plateau = x_Lu + Lu − C_N (π quasi-lié, paire σ ; R5 C 5…12, R7 D1 6…27) ; a/A2_tables.md,
           fig/levels_vs_invN.
  a3       A.3 : colonne QE de R9 C.1 refaite avec le plateau, modèles repliés de R9 (c_results.json) tels quels ; écart chaîne − QE ; a/A3_tables.md.
  atables  tables et figures de A depuis les json (aucun calcul).
  b        B : --parts b0,b1,b2. b0 porte (R9 B.0, 300²) ; b1 E_res contre la grille de sortie (120², 240², 480², 960² ; N_k^int 900 ; tel quel et
           plateau (i) ; chaîne « res » de R9 B ; fenêtre complète ±1,5 eV ; 240² d'abord = porte contre R9 B à 900²) ; b2 familles 5…12 (XML + json).
  btables  tables et figures de B depuis les json (aucun calcul).
Sous-commandes (GO 2, partie C ; chaîne de production : defect_mwr, C_N de la config, matrices dans matrices_dir, produits dans results_dir) :
  c0       portes C.0 : (a) C_N = 0, étiquettes actuelles, 9×9 (niveau 1 par scattering_rate_fast, tab:rcut_M R_cut 3) = R9 ; (b) plateau (i), 9×9 et 12×12,
           mêmes grandeurs = R9 (a3_results_plateau), V_loc = R9 au bit ; (c) Wigner-Seitz : Mwr_to_Mwk(ws) = sans ws sur la grille MP (1e-12), carte D.1
           tel quel sans ws = R9 D.1, avec ws = « vraie | brut » de l'audit (1e-9 eV Å²) ; (d) (i) avec C = Lu contre R9 D.1 exact (rapporté). Échec : code 3.
  c1post   fin de C.1 : m_rcut_resigma.csv (porte : npz de M2 -> csv de M2 à 1e-9) ; tab:tests_M : test d'or et porte A.2 recopiés de M2 (D12).
  c2       tab:rcut_M à trois colonnes (tel quel / plateau, étiquettes actuelles / plateau, Wigner-Seitz ; porte : dernière colonne = csv de C.1) ;
           carte de Kaasbjerg D.1 et bloc D.2, plateau (i) + Wigner-Seitz (tel quel en regard) ; fig/kaasbjerg_plateau_ws.
  c3       contrôles de l'audit : P-c2 (anneaux, porte : appartenance identique, valeurs à 1e-13 eV sous 0,433 a_sc), P-b2 (abscisses et poids = audit, M2),
           tab:tests_M par famille ; code 3 si une porte refuse.
  c4       planches avant (figures/) / après (fig/) : fig/avant_apres/ ; carte de M brute / alignée (fig/M_map_brut_aligne).
  c5       table de correspondance results/M2 -> results/M2_plateau (c/table_v2_plateau.md) ; README.md et MD5SUMS de results_dir.
  c6       anomalie ⟨ΔV⟩_3D d'A.1 : état des SCF (15…27, 5…12 en regard) et profil ⟨ΔV⟩(z) moyenné dans le plan (13 tailles) ; fig/dV_z_profiles.
Données de production en lecture seule (results/M2, caches Wannier, .save, Vks, config/, article/R4…R9) ; sorties dans ce répertoire (a/, b/, fig/,
cache/) ; journal r10_log.txt. Paramètres de production (config v2) sauf l'alignement (B) et les grilles de sortie de B.1.
"""
import os
import sys
import json
import time
import argparse
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                   # .../graphene/qe
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src"))
sys.path.insert(0, os.path.join(PROJ, "scripts"))                              # _palette.py
os.chdir(PROJ)

from electron_defect_interaction.config import load_production, dense_paths, results_dir, alignment_C, HA2EV  # noqa: E402
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance  # noqa: E402
from electron_defect_interaction.io import qe_gamma_io as qg  # noqa: E402
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mwr_to_Mwk, Mwr_to_Mwk_pairs, ws_images, ws_phase  # noqa: E402  (R10 (a))
from electron_defect_interaction.defects.many_body import local_tmatrix as lt  # noqa: E402
from electron_defect_interaction.defects import alignment as al  # noqa: E402  (P1, R9)

BOHR = 0.529177210903
NW = 5
K_RED = np.array([2 / 3, 1 / 3, 0.0]); KP_RED = np.array([1 / 3, 2 / 3, 0.0])
SCRATCH = "/home/gregb26/links/scratch/qe_tmp"                                 # .save de production (R4, R5, R7, R9)
ART = os.path.join(PROJ, "article")
R5_C = os.path.join(ART, "R5_base_vs_M", "c", "c_results.json")
R7_D1 = os.path.join(ART, "R7_tailles_3m", "d1_8pts", "d1_results.json")
R9_A1 = os.path.join(ART, "R9_controles", "a", "a1_results.json")
R9_A3P = os.path.join(ART, "R9_controles", "a", "a3_results_plateau.json")
R9_B = os.path.join(ART, "R9_controles", "b", "b_results.json")
R9_BP = os.path.join(ART, "R9_controles", "b", "b_results_plateau.json")
R9_C = os.path.join(ART, "R9_controles", "c", "c_results.json")
R9_AUDIT = os.path.join(ART, "R9_controles", "audit", "audit_results.json")
R9_WORK = os.path.join(os.path.dirname(WORK), "R9_controles")                 # npz de R9 (courbes de B), lecture seule
GATE6 = ["5x5", "6x6", "7x7", "8x8", "9x9", "12x12"]
NEW7 = ["10x10", "11x11", "15x15", "18x18", "21x21", "24x24", "27x27"]
ALL13 = sorted(GATE6 + NEW7, key=lambda s: int(s.split("x")[0]))
FAM3 = {6, 9, 12, 15, 18, 21, 24, 27}
PLATEAU_FRAC = 0.75                                                            # décision R9 : atomes à d >= 0,75 r_max
TIE = 1e-6                                                                     # égalités de distance (Å), comme l'audit
B_GRIDS = [240, 120, 480, 960]                                                 # 240² d'abord : porte contre R9 B
B_NK = 900
LOGF = os.path.join(WORK, "r10_log.txt")
M2_DIR = os.path.join(PROJ, "results", "M2")                                   # v2 tel quel, gelé (lecture seule) : référence de A, B et C
FIGURES_DIR = os.path.join(PROJ, "figures")                                    # lecture seule (planches de C.4)
R9_A3 = os.path.join(ART, "R9_controles", "a", "a3_results.json")
R9_D = os.path.join(ART, "R9_controles", "d", "d_results.json")
SC_DIR = os.path.join(GQ, "defects", "super_cell")                              # scf.in / scf.out / pp.in de production (C.6)


# ----------------------------------------------------------------------------------------------- utilitaires
def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOGF, "a") as f:
        f.write(line + "\n")


def ensure(sub):
    d = os.path.join(WORK, sub)
    os.makedirs(d, exist_ok=True)
    return d


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    if isinstance(x, complex):
        return [x.real, x.imag]
    return x


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(jsonable(d), f, indent=1, ensure_ascii=False)


def load_json(path, default=None):
    return json.load(open(path)) if os.path.exists(path) else ({} if default is None else default)


def git_head():
    try:
        return subprocess.run(["git", "-C", PROJ, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "?"


def md5(path):
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def size_n(S):
    return int(S.split("x")[0])


def fr(x, nd=2):
    """Nombre au format français (virgule décimale), signe explicite si demandé par nd < 0."""
    s = f"{x:+.{-nd}f}" if nd < 0 else f"{x:.{nd}f}"
    return s.replace(".", ",").replace("-", "−")                             # virgule décimale, vrai signe moins


def fig_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.join(PROJ, "figures", "memoire.mplstyle"))
    import _palette as pal
    return plt, pal


def savefig(fig, name):
    d = ensure("fig")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(d, f"{name}.{ext}"), dpi=200 if ext == "png" else None)
    log(f"[fig] {name}.pdf/.png")


# ----------------------------------------------------------------------------------------------- A : potentiels et P1
def sc_paths(S):
    return dict(d=f"{SCRATCH}/defect_{S}_d/defect_{S}_d.save", p=f"{SCRATCH}/defect_{S}_p/defect_{S}_p.save",
                vd=f"{GQ}/defects/super_cell/{S}/defective/Vks_{S}_d", vp=f"{GQ}/defects/super_cell/{S}/pristine/Vks_{S}_p")


def load_pot_eV(path):
    V, _ = qe_io.get_pot(path, subtract_mean=False, to_hartree=True)
    return V.transpose(2, 1, 0) * HA2EV


def plateau_stats(values, mask):
    v = np.asarray(values)[mask]; mu = float(v.mean())                        # mêmes opérations que r9_driver.plateau_stats
    return dict(mean=mu, rms=float(np.sqrt(((v - mu) ** 2).mean())), max_abs_dev=float(np.abs(v - mu).max()), n=int(mask.sum()))


def lu_published():
    """Lu (1,0 Å) publié par taille : R5 C (5…12), R7 D1 huit points (15…27). eV, avec la source et l'atome de Lu."""
    r5 = load_json(R5_C)["sizes"]; r7 = load_json(R7_D1)["sizes"]; out = {}
    for S in ALL13:
        n = str(size_n(S)); src, rec = ("R5 C", r5[n]) if n in r5 else ("R7 D1", r7[n])
        out[S] = dict(value_eV=float(rec["alignment"]["shift_10"]), source=src, i_far_1based=int(rec["alignment"]["i_far_1based"]),
                      dist_far_axis_A=float(rec["alignment"]["dist_far_A"]))
    return out


def p1_size(S, pub):
    """P1 autour de tous les atomes (1,0 et 0,5 Å), plateau (i), atomes vraiment les plus loin, contrôle de Lu ; profil npz dans a/."""
    t0 = time.time(); P = sc_paths(S); n = size_n(S)
    A_b, _ = qe_io.get_A_volume(P["d"]); A_A = A_b * BOHR
    x_p = qe_io.get_x_red(P["p"]); x_d = qe_io.get_x_red(P["d"])
    ng = tuple(int(v) for v in qe_io.get_ngfft(P["d"]))
    Vp = load_pot_eV(P["vp"]); Vd = load_pot_eV(P["vd"]); t_read = time.time() - t0
    assert Vp.shape == ng and Vd.shape == ng, (S, Vp.shape, Vd.shape, ng)
    r10 = al.atom_sphere_shifts(Vd, Vp, x_d, x_p, A_A, 1.0); r05 = al.atom_sphere_shifts(Vd, Vp, x_d, x_p, A_A, 0.5)
    lu = al.far_atom_alignment(Vd, Vp, x_d, x_p, A_A, (0.5, 1.0))
    mean3d = float(Vd.mean() - Vp.mean()); del Vd, Vp
    rmax = float(r10["dist"].max()); mfar = r10["dist"] >= PLATEAU_FRAC * rmax
    pl10 = plateau_stats(r10["shift"], mfar); pl05 = plateau_stats(r05["shift"], mfar)
    far = np.where(r10["dist"] >= rmax - TIE)[0]; ifar = int(r10["i_far_axis"])
    lp = pub[S]
    rec = dict(N=n, n_atoms=int(len(x_d)), ngfft=list(ng), r_max_A=rmax, plateau_frac=PLATEAU_FRAC, d_min_plateau_A=PLATEAU_FRAC * rmax,
               plateau_i_10=pl10, plateau_i_05=pl05, C_N_eV=pl10["mean"],
               far_true=dict(index1=[int(i) + 1 for i in far], dist_A=[float(r10["dist"][i]) for i in far],
                             shift10_meV=[float(r10["shift"][i] * 1e3) for i in far], shift05_meV=[float(r05["shift"][i] * 1e3) for i in far],
                             mean10_meV=float(np.mean(r10["shift"][far]) * 1e3), mean05_meV=float(np.mean(r05["shift"][far]) * 1e3)),
               Lu=dict(axis_atom_1based=ifar + 1, dist_axis_A=float(r10["dist_axis"][ifar]), dist_true_A=float(r10["dist"][ifar]),
                       P1_meV=float(r10["shift"][ifar] * 1e3), P1_05_meV=float(r05["shift"][ifar] * 1e3),
                       far_atom_alignment_meV=float(lu["shifts"][1.0] * 1e3), far_atom_alignment_05_meV=float(lu["shifts"][0.5] * 1e3),
                       same_as_function=bool(abs(r10["shift"][ifar] - lu["shifts"][1.0]) < 1e-12),
                       published_meV=lp["value_eV"] * 1e3, published_source=lp["source"], published_atom_1based=lp["i_far_1based"],
                       dev_published_meV=float((r10["shift"][ifar] - lp["value_eV"]) * 1e3)),
               mean3d_meV=mean3d * 1e3,
               files={k: os.path.realpath(v) for k, v in P.items()}, read_s=t_read, elapsed_s=time.time() - t0)
    np.savez(os.path.join(ensure("a"), f"profiles_{S}.npz"), dist=r10["dist"], dist_axis=r10["dist_axis"], shift10=r10["shift"], shift05=r05["shift"],
             npts10=r10["npts_d"], r_max=rmax, far=far, i_far_axis=ifar)
    log(f"[A] {S} : {len(x_d)} atomes, grille {ng} ; plateau (i) 1,0 Å {pl10['mean']*1e3:+.4f} meV (rms {pl10['rms']*1e3:.2f}, max|écart| "
        f"{pl10['max_abs_dev']*1e3:.2f}, {pl10['n']} atomes à d >= {PLATEAU_FRAC*rmax:.2f} Å) ; 0,5 Å {pl05['mean']*1e3:+.2f} ; un atome vrai "
        f"{rec['far_true']['index1']} -> {np.round(rec['far_true']['shift10_meV'], 2).tolist()} (moyenne {rec['far_true']['mean10_meV']:+.2f}) ; Lu : atome "
        f"{ifar+1} P1 {rec['Lu']['P1_meV']:+.3f} (fonction {rec['Lu']['far_atom_alignment_meV']:+.3f}, publié {lp['value_eV']*1e3:+.3f} {lp['source']}, "
        f"atome publié {lp['i_far_1based']}) ; <ΔV>_3D {mean3d*1e3:+.2f} ; lecture {t_read:.0f} s, total {time.time()-t0:.0f} s")
    return rec


def a1_path():
    return os.path.join(ensure("a"), "a1_results.json")


def cmd_a0(a):
    """Porte A.0 : six tailles de R9."""
    load_production(verbose=True); pub = lu_published()
    r9 = load_json(R9_A1)["sizes"]; aud = load_json(R9_AUDIT)["sizes"]
    out = load_json(a1_path()); out.setdefault("sizes", {})
    out.update(rule=dict(radius_A=1.0, plateau_frac=PLATEAU_FRAC, tie_A=TIE, C_N="moyenne du plateau (i), sans critère d'arrêt (D14)"), head=git_head())
    gate = dict(rows={}, ok=True)
    for S in GATE6:
        rec = p1_size(S, pub); out["sizes"][S] = rec; save_json(a1_path(), out)
        ci9 = float(r9[S]["C_i_eV"]); far9 = aud[S]["far"]["far_true_1based"]; fm9 = aud[S]["Lu"]["far_true_meV"]
        same_far = rec["far_true"]["index1"] == far9
        dfar = float(np.max(np.abs(np.array(rec["far_true"]["shift10_meV"]) - np.array(fm9)))) if same_far else float("inf")
        row = dict(C_i_R10_eV=rec["C_N_eV"], C_i_R9_eV=ci9, C_bit=bool(rec["C_N_eV"] == ci9), far_R10=rec["far_true"]["index1"], far_audit=far9,
                   far_same=bool(same_far), far_max_dev_meV=dfar, mean_far_R10_meV=rec["far_true"]["mean10_meV"], mean_far_audit_meV=float(np.mean(fm9)),
                   Lu_same_as_function=rec["Lu"]["same_as_function"], Lu_dev_published_meV=rec["Lu"]["dev_published_meV"])
        row["ok"] = bool(row["C_bit"] and same_far and dfar <= 1e-6 and row["Lu_same_as_function"] and abs(row["Lu_dev_published_meV"]) <= 1e-6)
        gate["rows"][S] = row; gate["ok"] &= row["ok"]
        log(f"[A.0] {S} : C_i {rec['C_N_eV']!r} (R9 {ci9!r}, au bit : {row['C_bit']}) ; un atome vrai {rec['far_true']['index1']} (audit {far9}), écart max "
            f"{dfar:.1e} meV, moyenne {row['mean_far_R10_meV']:+.2f} (audit {row['mean_far_audit_meV']:+.2f}) ; Lu = fonction : {row['Lu_same_as_function']}, "
            f"écart au publié {row['Lu_dev_published_meV']:.1e} meV -> {'OK' if row['ok'] else 'ÉCHEC'}")
    out["gate_A0"] = gate; save_json(a1_path(), out)
    log(f"[A.0] porte : {'OK' if gate['ok'] else 'ÉCHEC — STOP sur la partie A'}")
    if not gate["ok"]:
        sys.exit(3)


def cmd_a1(a):
    """A.1 : sept tailles nouvelles ; tableau 13 tailles ; figure."""
    load_production(verbose=True); pub = lu_published()
    out = load_json(a1_path())
    if not out.get("gate_A0", {}).get("ok"):
        log("[A.1] porte A.0 absente ou en échec : STOP"); sys.exit(3)
    for S in NEW7:
        out["sizes"][S] = p1_size(S, pub); save_json(a1_path(), out)
    out["C_N_eV"] = {S: out["sizes"][S]["C_N_eV"] for S in ALL13 if S in out["sizes"]}
    out["head"] = git_head(); save_json(a1_path(), out)
    a1_tables(out); a1_figure(out)


def a1_tables(out):
    L = ["# R10 — A.1 : C_N = moyenne du plateau (i) pour les 13 tailles (tables générées par r10_driver.py)", "",
         f"Plateau : atomes de la cellule avec lacune à distance vraie (image minimale) ≥ {fr(PLATEAU_FRAC)} r_max de la lacune ; sphères de P1 "
         "(`alignment.atom_sphere_shifts`) ; ΔV = ⟨V_d⟩ − ⟨V_p⟩ ; « un atome vrai » : atomes à r_max en vraie image (égalités à 1e-6 Å), moyenne en cas d'égalité ; "
         "Lu publié : atome de `far_atom` (réduction axe par axe), R5 C (5…12) et R7 D1 (15…27). meV sauf mention.", ""]
    g = out.get("gate_A0")
    if g:
        L += ["Porte A.0 (six tailles de R9) : " + " ; ".join(f"{S.replace('x', '×')} C_i au bit {'oui' if r['C_bit'] else 'NON'}, un atome {r['far_R10']} "
                                                             f"(écart max à l'audit {r['far_max_dev_meV']:.1e} meV)" for S, r in g["rows"].items())
              + f" → {'OK' if g['ok'] else 'ÉCHEC'}", ""]
    L += ["| N | r_max (Å) | plateau (i) 1,0 Å : moyenne ; rms ; max\\|écart\\| ; n | plateau 0,5 Å : moyenne | un atome vrai : indices ; valeurs 1,0 Å ; moyenne | "
          "Lu publié : atome ; valeur | C_N − Lu | ⟨ΔV⟩_3D |", "|---|---|---|---|---|---|---|---|"]
    for S in ALL13:
        if S not in out["sizes"]:
            continue
        r = out["sizes"][S]; p, q, f, u = r["plateau_i_10"], r["plateau_i_05"], r["far_true"], r["Lu"]
        L.append(f"| {S.replace('x', '×')} | {fr(r['r_max_A'])} | {fr(p['mean']*1e3, -2)} ; {fr(p['rms']*1e3)} ; {fr(p['max_abs_dev']*1e3)} ; {p['n']} | "
                 f"{fr(q['mean']*1e3, -2)} | {', '.join(map(str, f['index1']))} ; {' ; '.join(fr(v, -2) for v in f['shift10_meV'])} ; {fr(f['mean10_meV'], -2)} | "
                 f"{u['axis_atom_1based']} ; {fr(u['published_meV'], -2)} | {fr(p['mean']*1e3 - u['published_meV'], -2)} | {fr(r['mean3d_meV'], -2)} |")
    L += ["", "Contrôles par taille : P1 à l'atome de Lu = `far_atom_alignment` au bit ; écart au Lu publié :", ""]
    for S in ALL13:
        if S in out["sizes"]:
            u = out["sizes"][S]["Lu"]
            L.append(f"- {S.replace('x', '×')} : atome {u['axis_atom_1based']} (publié {u['published_atom_1based']}, {u['published_source']}) à "
                     f"{fr(u['dist_axis_A'])} Å axe par axe, {fr(u['dist_true_A'])} Å vrai ; P1 = fonction : {'oui' if u['same_as_function'] else 'NON'} ; "
                     f"écart au publié {u['dev_published_meV']:.1e} meV")
    with open(os.path.join(ensure("a"), "A1_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def a1_figure(out):
    plt, pal = fig_style()
    sizes = [S for S in ALL13 if S in out["sizes"]]
    fig, axes = plt.subplots(4, 4, figsize=(6.5, 7.6), squeeze=False); nout = []
    for ax in axes.ravel()[len(sizes):]:
        ax.set_visible(False)
    for i, (ax, S) in enumerate(zip(axes.ravel(), sizes)):
        r = out["sizes"][S]; Z = np.load(os.path.join(WORK, "a", f"profiles_{S}.npz"))
        d, s = Z["dist"], Z["shift10"] * 1e3; C = r["C_N_eV"] * 1e3; far = Z["far"]
        ax.axvspan(PLATEAU_FRAC * r["r_max_A"], r["r_max_A"] * 1.03, color=pal.LIGHT, alpha=0.35, lw=0)
        ax.plot(d, s, "o", ms=1.6, color=pal.NAVY, label=r"sphères de 1,0 Å")
        ax.axhline(C, color=pal.NAVY, lw=0.8, label=r"$C_N$ (moyenne du plateau)")
        ax.axhline(r["Lu"]["published_meV"], color=pal.REF, lw=0.8, ls="--", label="Lu publié (un atome, axe par axe)")
        ax.plot(d[far], s[far], "D", ms=3.2, mfc="none", mec=pal.ORANGE, mew=0.9, label="atomes vraiment les plus loin")
        keep = d >= 2.0                                                         # cadre sans les premiers voisins (1,42 Å)
        lo = min(s[keep].min(), C, r["Lu"]["published_meV"]); hi = max(s[keep].max(), C, r["Lu"]["published_meV"]); pad = 0.08 * (hi - lo)
        lo, hi = lo - pad, hi + pad; nout.append(int(((s < lo) | (s > hi)).sum()))
        ax.set_ylim(lo, hi); ax.set_title(f"({chr(97 + i)}) {S.replace('x', '×')}", fontsize=7, loc="left")
        ax.tick_params(labelsize=6)
        if i % 4 == 0:
            ax.set_ylabel(r"$\Delta V$ (meV)", fontsize=7)
        if i >= len(sizes) - 4:
            ax.set_xlabel("Distance à la lacune (Å)", fontsize=7)
    h, l_ = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l_, loc="lower center", ncol=2, fontsize=6, frameon=False, bbox_to_anchor=(0.5, 0.012))
    fig.text(0.5, 0.004, f"Hors cadre : {', '.join(map(str, sorted(set(nout))))} point(s) par panneau, les premiers voisins de la lacune (1,42 Å). Zone grisée : "
             + r"$d \geq 0{,}75\, r_\mathrm{max}$.",
             ha="center", va="bottom", fontsize=5.5, color=pal.MUTED)
    fig.tight_layout(rect=(0, 0.065, 1, 1)); savefig(fig, "offset_profiles_13"); plt.close(fig)


# ----------------------------------------------------------------------------------------------- A.2 et A.3 (json seulement)
def cn_table():
    out = load_json(a1_path())
    return {S: out["sizes"][S]["C_N_eV"] for S in ALL13 if S in out.get("sizes", {})}


def cmd_a2(a):
    C = cn_table(); r5 = load_json(R5_C)["sizes"]; r7 = load_json(R7_D1)["sizes"]
    rows = {}; checks = {}
    for S in ALL13:
        n = size_n(S); key = str(n)
        srcs = [("R5 C", r5[key])] if key in r5 else []
        if key in r7:
            srcs.append(("R7 D1", r7[key]))
        if S not in C or not srcs:
            log(f"[A.2] {S} : C_N ou niveaux absents, sauté"); continue
        per = {}
        for lab, rec in srcs:
            Lu = float(rec["alignment"]["shift_10"]); dx = Lu - C[S]
            pi = rec["pi_state"]; sg = rec["sigma_doublet"]
            per[lab] = dict(E_D=float(rec["E_D"]), E_D_source=rec["E_D_source"], Lu_eV=Lu, C_N_eV=C[S], shift_eV=dx,
                            pi=dict(band=int(pi["band"]), x_Lu=float(pi["x"]), x_plateau=float(pi["x"]) + dx, w2=float(pi["w2"])),
                            sigma=[dict(band=int(q["band"]), x_Lu=float(q["x"]), x_plateau=float(q["x"]) + dx, w2=float(q["w2"])) for q in sg])
        if len(per) == 2:                                                      # 6, 9, 12 : R5 C et R7 D1 (porte de R7 : identiques)
            p5, p7 = per["R5 C"], per["R7 D1"]
            checks[S] = max(abs(p5["pi"]["x_Lu"] - p7["pi"]["x_Lu"]), abs(p5["Lu_eV"] - p7["Lu_eV"]),
                            max(abs(u["x_Lu"] - v["x_Lu"]) for u, v in zip(p5["sigma"], p7["sigma"])))
        rows[S] = per
        p = per[srcs[0][0]]
        log(f"[A.2] {S} : Lu {p['Lu_eV']*1e3:+.2f}, C_N {C[S]*1e3:+.2f}, Lu − C_N {p['shift_eV']*1e3:+.2f} meV ; π {p['pi']['x_Lu']:+.4f} -> {p['pi']['x_plateau']:+.4f} ; "
            f"σ {[round(q['x_Lu'], 4) for q in p['sigma']]} -> {[round(q['x_plateau'], 4) for q in p['sigma']]}")
    res = dict(rows=rows, checks_R5_vs_R7=checks, relation="x_plateau = x_Lu + Lu − C_N (r5_driver.window_states_gamma l. 605)", head=git_head())
    save_json(os.path.join(ensure("a"), "a2_results.json"), res)
    a2_tables(res); a2_figure(res)


def a2_tables(res):
    L = ["# R10 — A.2 : niveaux QE réalignés sans recalcul (tables générées par r10_driver.py)", "",
         "x = ε − C − E_D (eV), C = Lu (publié, un atome) ou C_N (moyenne du plateau (i)) ; x_plateau = x_Lu + Lu − C_N ; w₂ inchangé. Sources : R5 C (5…12), "
         "R7 D1 huit points (15…27 ; 6, 9, 12 aussi, écart R5/R7 en fin de table).", "",
         "| N | famille | E_D (eV) ; définition | Lu (meV) | C_N (meV) | Lu − C_N (meV) | π : bande ; x Lu → x plateau ; w₂ | σ : bandes ; x Lu → x plateau ; w₂ | source |",
         "|---|---|---|---|---|---|---|---|---|"]
    for S, per in res["rows"].items():
        lab = "R5 C" if "R5 C" in per else "R7 D1"; p = per[lab]; n = size_n(S)
        sg = p["sigma"]
        L.append(f"| {S.replace('x', '×')} | {'3m' if n in FAM3 else 'non-3m'} | {fr(p['E_D'], 5)} ; {p['E_D_source']} | {fr(p['Lu_eV']*1e3, -2)} | {fr(p['C_N_eV']*1e3, -2)} | "
                 f"{fr(p['shift_eV']*1e3, -2)} | {p['pi']['band']} ; {fr(p['pi']['x_Lu'], -4)} → {fr(p['pi']['x_plateau'], -4)} ; {fr(p['pi']['w2'], 3)} | "
                 f"{', '.join(str(q['band']) for q in sg)} ; {' / '.join(fr(q['x_Lu'], -4) for q in sg)} → {' / '.join(fr(q['x_plateau'], -4) for q in sg)} ; "
                 f"{fr(sg[0]['w2'], 3)} | {' + '.join(per)} |")
    L += ["", "Écart R5 C / R7 D1 sur x_Lu et Lu (6, 9, 12) : " + " ; ".join(f"{S.replace('x', '×')} {v:.1e} eV" for S, v in res["checks_R5_vs_R7"].items())]
    with open(os.path.join(ensure("a"), "A2_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def a2_figure(res):
    plt, pal = fig_style()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.3), sharex=True)
    for ax, what, ttl in zip(axes, ("pi", "sigma"), (r"(a) état $\pi$ quasi-lié", r"(b) paire $\sigma$ (moyenne)")):
        for fam, mk in (("3m", "o"), ("non-3m", "s")):
            xs, yl, yp = [], [], []
            for S, per in res["rows"].items():
                n = size_n(S)
                if ("3m" if n in FAM3 else "non-3m") != fam:
                    continue
                p = per["R5 C" if "R5 C" in per else "R7 D1"]
                if what == "pi":
                    yl.append(p["pi"]["x_Lu"]); yp.append(p["pi"]["x_plateau"])
                else:
                    yl.append(np.mean([q["x_Lu"] for q in p["sigma"]])); yp.append(np.mean([q["x_plateau"] for q in p["sigma"]]))
                xs.append(1.0 / n)
            ax.plot(xs, yl, mk, ms=4.5, mfc="none", mec=pal.ORANGE, mew=0.9, ls="none", label=f"Lu publié, {fam}")
            ax.plot(xs, yp, mk, ms=3.5, color=pal.NAVY, ls="none", label=f"plateau $C_N$, {fam}")
        ax.set_title(ttl, loc="left", fontsize=9); ax.set_xlabel(r"$1/N$"); ax.axhline(0, color=pal.MUTED, lw=0.6)
    axes[0].set_ylabel(r"$\varepsilon - E_D$ (eV)")
    h, l_ = axes[0].get_legend_handles_labels(); fig.legend(h, l_, loc="lower center", ncol=4, fontsize=6, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1)); savefig(fig, "levels_vs_invN"); plt.close(fig)


def cmd_a3(a):
    C = cn_table(); c = load_json(R9_C); rows = {}
    for key, e in c["C1"].items():
        S = f"{key}x{key}"
        if S not in C:
            log(f"[A.3] {S} : C_N absent, sauté"); continue
        q = e["QE"]; Lu = float(q["Lu_1A"]); dx = Lu - C[S]
        qpi_lu = float(q["pi_quasi_bound"]["x"]); qsg_lu = [float(s["x"]) for s in q["even"]]
        row = dict(N=int(key), Lu_eV=Lu, C_N_eV=C[S], QE_E_D=float(q["E_D"]), QE_pi_Lu=qpi_lu, QE_pi_plateau=qpi_lu + dx,
                   QE_sigma_Lu=qsg_lu, QE_sigma_plateau=[v + dx for v in qsg_lu], C9_model_eV=float(c["C9_eV"]), C9_model_label=c["C9_label"])
        for var in ("brut", "aligne"):
            m = e[var]; mpi = float(m["odd"]["quasi_bound"]["x"]); msg = float(m["even_window"][0]["x"]) if m["even_window"] else float("nan")
            row[var] = dict(pi=mpi, sigma=msg, d_pi_vs_QE_Lu=mpi - qpi_lu, d_pi_vs_QE_plateau=mpi - (qpi_lu + dx),
                            d_sigma_vs_QE_Lu=msg - float(np.mean(qsg_lu)), d_sigma_vs_QE_plateau=msg - float(np.mean(qsg_lu)) - dx)
        rows[S] = row
        log(f"[A.3] {S} : QE π {qpi_lu:+.4f} (Lu) -> {qpi_lu+dx:+.4f} (plateau) ; modèle tel quel {row['brut']['pi']:+.4f} (écart {row['brut']['d_pi_vs_QE_plateau']*1e3:+.1f} meV), "
            f"aligné {row['aligne']['pi']:+.4f} ({row['aligne']['d_pi_vs_QE_plateau']*1e3:+.1f} meV)")
    res = dict(rows=rows, source=R9_C, note="modèles de R9 C.1 non recalculés (aligné : C_9 = Lu de R9) ; QE réaligné avec C_N de R10 A.1", head=git_head())
    save_json(os.path.join(ensure("a"), "a3_results.json"), res); a3_tables(res)


def a3_tables(res):
    L = ["# R10 — A.3 : colonne QE de R9 C.1 avec le plateau (tables générées par r10_driver.py)", "",
         "QE (R7) : x = ε − C − E_D(quadruplet de la parfaite QE) ; modèle (R9 C.1, non recalculé) : ε − E_D(quadruplet de la parfaite repliée), « aligné » = V_loc − "
         "C_9·P_boîte avec C_9 = Lu de R9. Écart chaîne − QE en meV. σ : QE = moyenne de la paire, modèle = état pair de la fenêtre.", "",
         "| N | Lu (meV) | C_N (meV) | QE π : Lu → plateau | modèle tel quel π ; écart à QE Lu / plateau | modèle aligné π ; écart à QE Lu / plateau | QE σ : Lu → plateau | "
         "modèle tel quel σ ; écart plateau | modèle aligné σ ; écart plateau |", "|---|---|---|---|---|---|---|---|---|"]
    for S, r in res["rows"].items():
        b, g = r["brut"], r["aligne"]
        L.append(f"| {S.replace('x', '×')} | {fr(r['Lu_eV']*1e3, -2)} | {fr(r['C_N_eV']*1e3, -2)} | {fr(r['QE_pi_Lu'], -4)} → {fr(r['QE_pi_plateau'], -4)} | "
                 f"{fr(b['pi'], -4)} ; {fr(b['d_pi_vs_QE_Lu']*1e3, -1)} / {fr(b['d_pi_vs_QE_plateau']*1e3, -1)} | {fr(g['pi'], -4)} ; {fr(g['d_pi_vs_QE_Lu']*1e3, -1)} / "
                 f"{fr(g['d_pi_vs_QE_plateau']*1e3, -1)} | {fr(np.mean(r['QE_sigma_Lu']), -4)} → {fr(np.mean(r['QE_sigma_plateau']), -4)} | {fr(b['sigma'], -4)} ; "
                 f"{fr(b['d_sigma_vs_QE_plateau']*1e3, -1)} | {fr(g['sigma'], -4)} ; {fr(g['d_sigma_vs_QE_plateau']*1e3, -1)} |")
    with open(os.path.join(ensure("a"), "A3_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def cmd_atables(a):
    out = load_json(a1_path())
    if out.get("sizes"):
        a1_tables(out); a1_figure(out)
    f2 = os.path.join(WORK, "a", "a2_results.json")
    if os.path.exists(f2):
        r = load_json(f2); a2_tables(r); a2_figure(r)
    f3 = os.path.join(WORK, "a", "a3_results.json")
    if os.path.exists(f3):
        a3_tables(load_json(f3))


# ----------------------------------------------------------------------------------------------- B : chaîne « res » de R9 B
def load_wannier(dp, k):
    """U, U_dis réordonnés sur les k du .save dense, H(R) ; porte de jauge (refus si le manifeste ne passe pas)."""
    paths = wannier_provenance.load_wannier_checked(dp["manifest"])
    U, kU = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_kpoint_order(kUd, k)]
    Hwr, Rw, nd = read_w90_HR(paths["tb"])
    return U, Ud, Hwr, Rw, nd


def setup_size(cfg, S="9x9", rc=3):
    """M_W du M2 dense (V† M V, double TF, k du XML), recentrage, V_loc (R_cut), boîte N×N, E_D (grille 90²), grille « res » : chaîne de R9 B."""
    t0 = time.time(); dp = dense_paths(cfg, S); n = size_n(S)
    M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k)
    U, Ud, Hwr, Rw, nd = load_wannier(dp, k)
    Mwr, R = Mwk_to_Mwr(Mbk_to_Mwk(M, U, Ud), k, MP); del M
    Rn, R_d = lt.recenter_mwr(Mwr, R, MP); lt.mwr_locality(Mwr, Rn)
    inb = np.all(np.mod(np.asarray(R, int)[:, :2], np.asarray(MP[:2], int)) < n, axis=1)       # = box_geometry de r9_driver.py
    Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]
    V, herm = lt.extract_V_loc(Mwr, Rn, Rloc)
    iloc = np.array([int(np.where((Rn == r).all(axis=1))[0][0]) for r in Rloc])
    _, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
    gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
    eta, ew, npe = float(cfg["eta_eV"]), float(cfg["e_window_eV"]), int(cfg["ne_per_eta"]); de = eta / npe
    egrid = np.arange(E_D - ew - eta, E_D + ew + eta + de, de)                                  # resonance_metrics.py l. 62
    A_uc, _ = qe_io.get_A_volume(dp["uc"]); Bc = 2 * np.pi * np.linalg.inv(A_uc * BOHR).T
    st = dict(S=S, n=n, MP=MP, R_d=R_d, Rloc=Rloc, V=V, herm=herm, in_box=inb[iloc], Hwr=Hwr, Rw=Rw, nd=nd, E_D=E_D, eta=eta, ew=ew, de=de, egrid=egrid, Bc=Bc)
    log(f"[setup] {S} : R_d {R_d.tolist()}, R_cut {rc} : {len(Rloc)} sites (dim {len(Rloc)*NW}, {int((~st['in_box']).sum())} hors boîte), herm V_loc {herm:.1e} ; "
        f"E_D {E_D:.6f} eV ; grille « res » {len(egrid)} énergies, pas {de*1e3:.2f} meV ; {time.time()-t0:.0f} s")
    return st


def hwk_chunked(Hwr, Rw, nd, k, chunk=400_000):
    out = np.empty((len(k), Hwr.shape[1], Hwr.shape[1]), dtype=complex)                        # = r9_driver.hwk_chunked
    for s in range(0, len(k), chunk):
        out[s:s + chunk] = lt.Hwr_to_Hwk(Hwr, Rw, k[s:s + chunk], ndegen=nd)[0]
    return out


def g0_res(st, nk):
    """g₀ sur l'amas, grille « res », N_k^int = nk ; cache npy vérifié par ses métadonnées (comme r9_driver.g0_cached)."""
    egrid = st["egrid"]; f = os.path.join(ensure("cache"), f"g0_{st['S']}_nk{nk}_res.npy"); fj = f.replace(".npy", ".json")
    meta = dict(S=st["S"], nk_int=nk, E_D=st["E_D"], e0=float(egrid[0]), de=st["de"], nE=len(egrid), eta=st["eta"], Rloc=st["Rloc"].tolist())
    if os.path.exists(f) and os.path.exists(fj):
        old = json.load(open(fj))
        if all((old.get(k_) == v) if not isinstance(v, float) else abs(old.get(k_, 1e9) - v) < 1e-12 for k_, v in meta.items()):
            log(f"[g0] {st['S']} {nk}² : cache"); return np.load(f)
    t0 = time.time(); k_int = lt.mp_grid(nk, nk, 1); Hk = hwk_chunked(st["Hwr"], st["Rw"], st["nd"], k_int)
    g0 = lt.local_green_batch(Hk, k_int, st["Rloc"], egrid, st["eta"]); del Hk
    np.save(f, g0); save_json(fj, dict(meta, elapsed_s=time.time() - t0, head=git_head()))
    log(f"[g0] {st['S']} {nk}² « res » : {len(egrid)} énergies en {(time.time()-t0)/60:.1f} min")
    return g0


def t_matrix(V, g0):
    I = np.eye(V.shape[0])
    return np.array([V @ np.linalg.solve(I - g0[j] @ V, I) for j in range(len(g0))])          # = resonance_metrics.py


def out_grid(st, N):
    """États de la grille de sortie N² : ε, U, sélection ±e_window, φ = <wR|nk> sur l'amas (disposition L*nw + w), indice de K."""
    t0 = time.time(); k_out = lt.mp_grid(N, N, 1)
    if len(k_out) <= 250_000:
        _, E_out, U_out = lt.Hwr_to_Hwk(st["Hwr"], st["Rw"], k_out, ndegen=st["nd"])            # une seule tranche (identique à R9 pour 240²)
    else:
        E_out = np.empty((len(k_out), NW)); U_out = np.empty((len(k_out), NW, NW), complex)
        for s in range(0, len(k_out), 230_400):
            _, E_out[s:s + 230_400], U_out[s:s + 230_400] = lt.Hwr_to_Hwk(st["Hwr"], st["Rw"], k_out[s:s + 230_400], ndegen=st["nd"])
    sel = np.abs(E_out - st["E_D"]) <= st["ew"]
    nL = len(st["Rloc"])
    phi = np.einsum("kL,kwn->knLw", lt._phase(k_out, st["Rloc"]), U_out, optimize=True).reshape(len(k_out), NW, nL * NW)
    iK = int(np.argmin(np.linalg.norm(np.mod(k_out - K_RED + 0.5, 1) - 0.5, axis=1)))
    dK = float(np.linalg.norm(np.mod(k_out[iK] - K_RED + 0.5, 1) - 0.5))
    log(f"[grille] {N}² : {len(k_out)} k, {int(sel.sum())} états à ±{st['ew']} eV ; K = k_out[{iK}] (écart {dK:.1e}) ; {time.time()-t0:.0f} s")
    return dict(N=N, k_out=k_out, E_out=E_out, U_out=U_out, sel=sel, phi=phi, iK=iK)


def lor(x, eta):
    return (eta / np.pi) / (x * x + eta * eta)


def observables(st, g, t):
    """Mêmes formules que resonance_metrics.py / r9_driver.resonance_observables ; courbes en tranches d'énergie (sommes par ligne inchangées)."""
    egrid, de, eta, E_D = st["egrid"], st["de"], st["eta"], st["E_D"]; E_out, sel, phi = g["E_out"], g["sel"], g["phi"]
    j_all = np.rint((E_out - egrid[0]) / de).astype(int)
    G = np.full(E_out.shape, np.nan)
    for j in np.unique(j_all[sel]):
        m = sel & (j_all == j); P = phi[m]
        G[m] = -2.0 * np.einsum("mi,mi->m", P.conj() @ t[j], P).imag
    Es = E_out[sel]; gs = G[sel]

    def curve(eg, chunk=128):
        c = np.empty(len(eg))
        for s in range(0, len(eg), chunk):
            w = lor(eg[s:s + chunk, None] - Es[None, :], eta); c[s:s + chunk] = (w * gs[None, :]).sum(1) / w.sum(1)
        return c
    eg5 = egrid[::2]; c5 = curve(eg5); c25 = curve(egrid)
    PK = phi[g["iK"], 3:5]
    tr = np.array([0.5 * np.trace(PK.conj() @ t[j] @ PK.T) for j in range(len(egrid))])
    mm = np.isfinite(G) & (np.abs(E_out - E_D) <= 1.5)
    ik_all, n_all = np.nonzero(mm); a_ = int(np.argmax(np.abs(G[mm])))
    atED = sel & (np.abs(E_out - E_D) < 1e-6)
    return dict(curve_eg5=c5, curve_eg25=c25, tr=tr,
                peak_GT_prod=float(eg5[int(np.argmax(c5))] - E_D), peak_GT_fine=float(egrid[int(np.argmax(c25))] - E_D), max_GT_fine=float(c25.max()),
                peak_ImTbar=float(egrid[int(np.argmax(-tr.imag))] - E_D), max_mImTbar=float((-tr.imag).max()),
                ReTbar_at_ED=float(np.interp(E_D, egrid, tr.real)), G_ED_meV=float(np.abs(G[atED]).mean() * 1e3) if atED.any() else float("nan"),
                n_ED=int(atED.sum()), median_GT_states_meV=float(np.nanmedian(np.abs(G[sel])) * 1e3), n_states=int(sel.sum()),
                E_res=float(E_out[ik_all[a_], n_all[a_]] - E_D), E_res_state=(int(ik_all[a_]), int(n_all[a_])), G_res_meV=float(np.abs(G[mm])[a_] * 1e3))


def crowns(st, g, band):
    """
    Couronnes de la grille N² autour de K et K′ (norme de Löschian du décalage (i, j)/N à la vallée la plus proche) : x = i² + j² − ij pour la cellule QE
    (b₁·b₂ = −½|b|², vérifié) — forme entière de r9_driver.crown_table ; énergies des états de la bande `band` par couronne, triées par énergie moyenne.
    """
    k, N, Bc = g["k_out"], g["N"], st["Bc"]; E = g["E_out"][:, band] - st["E_D"]
    b1, b2 = Bc[:2, 0], Bc[:2, 1]; cosang = float(b1 @ b2 / (np.linalg.norm(b1) * np.linalg.norm(b2)))
    assert abs(cosang + 0.5) < 1e-9, f"b1·b2/(|b1||b2|) = {cosang}"
    xs = np.full(len(k), np.iinfo(np.int64).max, dtype=np.int64)
    for Kv in (K_RED, KP_RED):
        ij = np.rint((np.mod(k - Kv + 0.5, 1.0) - 0.5)[:, :2] * N).astype(np.int64)
        xs = np.minimum(xs, ij[:, 0] ** 2 + ij[:, 1] ** 2 - ij[:, 0] * ij[:, 1])
    distinct = np.unique(xs[xs > 0]); rank = {int(x): i + 1 for i, x in enumerate(distinct)}
    order = np.argsort(xs, kind="stable"); xs_s = xs[order]; E_s = E[order]
    bounds = np.searchsorted(xs_s, distinct, side="left"), np.searchsorted(xs_s, distinct, side="right")
    rows = []
    for x, lo_, hi_ in zip(distinct, *bounds):
        e = E_s[lo_:hi_]
        rows.append(dict(x=int(x), rank=rank[int(x)], n_states=int(hi_ - lo_), E_mean=float(e.mean()), E_min=float(e.min()), E_max=float(e.max())))
    rows.sort(key=lambda r: r["E_mean"])
    return rows, xs


def crown_info(st, g, o):
    ik, n = o["E_res_state"]; rows, xs = crowns(st, g, n); x = int(xs[ik])
    i = next(j for j, r in enumerate(rows) if r["x"] == x); r = rows[i]
    nb = {}
    for lab, j in (("below", i - 1), ("above", i + 1)):
        if 0 <= j < len(rows):
            q = rows[j]; nb[lab] = dict(q, dE_res_minus_mean=o["E_res"] - q["E_mean"])
    return dict(x=x, rank=r["rank"], E_mean=r["E_mean"], E_min=r["E_min"], E_max=r["E_max"], n_states=r["n_states"], band=int(n),
                k_red=[float(v) for v in g["k_out"][ik]], neighbours=nb)


def c9_plateau():
    """C_9 = moyenne du plateau (i) : R10 A.1 si disponible (égal au bit à R9 par la porte A.0), sinon R9 A.1 (C_i_eV)."""
    r9 = float(load_json(R9_A1)["sizes"]["9x9"]["C_i_eV"]); r10 = load_json(a1_path()).get("sizes", {}).get("9x9", {}).get("C_N_eV")
    if r10 is not None:
        assert r10 == r9, (r10, r9)
        return float(r10), "R10 A.1 (= R9 C_i_eV au bit)"
    return r9, "R9 A.1 C_i_eV (A.1 de R10 absente)"


def rel(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return float(np.abs(a - b).max() / np.abs(b).max())


def cmd_b(a):
    parts = a.parts.split(","); cfg = load_production(verbose=True); d = ensure("b"); fres = os.path.join(d, "b_results.json")
    res = load_json(fres); res.setdefault("meta", {}).update(head=git_head(), config=dict(R_cut=cfg["R_cut"], eta=cfg["eta_eV"], ne_per_eta=cfg["ne_per_eta"],
                                                                                       e_window=cfg["e_window_eV"], nk_int_B1=B_NK, grids=B_GRIDS))
    st = setup_size(cfg) if ("b0" in parts or "b1" in parts) else None
    C9, c9src = c9_plateau(); res["meta"]["C9_eV"] = C9; res["meta"]["C9_source"] = c9src
    if st is not None:
        Va = st["V"] - C9 * np.diag(np.repeat(st["in_box"].astype(float), NW))                 # approximation (i) (= V_variant « aligne » de R9)
        VAR = {"brut": st["V"], "plateau": Va}
    if "b0" in parts:
        g = out_grid(st, int(cfg["grid"])); g0 = g0_res(st, 300); o = observables(st, g, t_matrix(st["V"], g0)); del g0
        z = np.load(os.path.join(M2_DIR, "resonance_9x9.npz"))
        dc = rel(o["curve_eg5"], z["Gamma_T"]); dt = rel(o["tr"], z["Tbar_tr"])
        ok = abs(o["peak_GT_prod"] - (-0.180)) < 5e-4 and abs(o["peak_ImTbar"] - (-0.177)) < 1e-3 and dc < 1e-6 and dt < 1e-6
        res["gate_B0"] = dict(peak_GT_prod=o["peak_GT_prod"], peak_ImTbar=o["peak_ImTbar"], curve_rel_dev_vs_R6=dc, Tbar_rel_dev_vs_R6=dt,
                              median_states_meV=o["median_GT_states_meV"], expected_median_states_meV=float(z["median_GT_states_meV"]), ok=bool(ok))
        save_json(fres, res)
        log(f"[B.0] 9x9 tel quel, 300², 240² : pic Γ_T {o['peak_GT_prod']:+.4f} (R6 −0,180), pic −Im T̄(K) {o['peak_ImTbar']:+.4f} (R6 −0,177) ; courbe vs "
            f"resonance_9x9.npz {dc:.1e}, T̄ {dt:.1e} ; médiane {o['median_GT_states_meV']:.4f} (R6 {float(z['median_GT_states_meV']):.4f}) -> {'OK' if ok else 'ÉCHEC'}")
        del g
        if not ok:
            log("[B.0] ÉCHEC de la porte : STOP sur la partie B"); sys.exit(3)
    if "b1" in parts:
        if not res.get("gate_B0", {}).get("ok"):
            log("[B.1] porte B.0 absente ou en échec : STOP"); sys.exit(3)
        g0 = g0_res(st, B_NK); T = {v: t_matrix(V, g0) for v, V in VAR.items()}; del g0
        cf = os.path.join(d, "b1_curves.npz"); curves = dict(np.load(cf)) if os.path.exists(cf) else {}
        r9 = {"brut": (load_json(R9_B)["9x9"]["brut"][str(B_NK)], np.load(os.path.join(R9_WORK, "b", "b_curves_9x9.npz"))["brut_nk900_GT"]),
              "plateau": (load_json(R9_BP)["9x9"]["aligne_plateau"][str(B_NK)], np.load(os.path.join(R9_WORK, "b", "b_curves_9x9_plateau.npz"))["aligne_plateau_nk900_GT"])}
        res.setdefault("B1", {})
        for N in B_GRIDS:
            t0 = time.time(); g = out_grid(st, N)
            for var in VAR:
                t1 = time.time(); o = observables(st, g, T[var]); cr = crown_info(st, g, o)
                row = {k_: v for k_, v in o.items() if k_ not in ("curve_eg5", "curve_eg25", "tr")}
                row.update(crown=cr, d_Eres_peak_prod_meV=abs(o["E_res"] - o["peak_GT_prod"]) * 1e3, d_Eres_peak_fine_meV=abs(o["E_res"] - o["peak_GT_fine"]) * 1e3,
                           C_eV=0.0 if var == "brut" else C9, elapsed_s=time.time() - t1)
                res["B1"].setdefault(var, {})[str(N)] = row
                curves[f"{var}_N{N}_GT25"] = o["curve_eg25"]; curves[f"{var}_N{N}_GT5"] = o["curve_eg5"]
                nb = cr["neighbours"]
                log(f"[B.1] {N}² {var} : E_res {o['E_res']:+.5f} (état {o['E_res_state']}, couronne x = {cr['x']}, rang {cr['rank']}, [{cr['E_min']:+.4f}, {cr['E_max']:+.4f}] ; "
                    f"voisines " + ", ".join(f"x={q['x']} {q['E_mean']:+.4f}" for q in nb.values()) + f") ; pic Γ_T {o['peak_GT_prod']:+.4f} (5 meV) / {o['peak_GT_fine']:+.4f} "
                    f"(2,5 meV) ; |E_res − pic| {row['d_Eres_peak_prod_meV']:.1f} / {row['d_Eres_peak_fine_meV']:.1f} meV ; médiane {o['median_GT_states_meV']:.2f} meV ; "
                    f"{o['n_states']} états ; {time.time()-t1:.0f} s")
                if N == 240:                                                    # porte : R9 B à N_k^int 900 (b_results*.json, b_curves*.npz)
                    ref, cref = r9[var]
                    gt = dict(E_res_state_same=list(o["E_res_state"]) == list(ref["E_res_state"]), E_res_dev=abs(o["E_res"] - ref["E_res"]),
                              peak_fine_dev=abs(o["peak_GT_fine"] - ref["peak_GT_fine"]), peak_prod_dev=abs(o["peak_GT_prod"] - ref["peak_GT_prod"]),
                              median_rel_dev=abs(o["median_GT_states_meV"] - ref["median_GT_states_meV"]) / ref["median_GT_states_meV"],
                              curve_rel_dev=rel(o["curve_eg25"], cref), crown_x=cr["x"], crown_x_R9=ref["E_res_crown"]["x"])
                    gt["ok"] = bool(gt["E_res_state_same"] and gt["E_res_dev"] < 1e-12 and gt["peak_fine_dev"] < 1e-9 and gt["peak_prod_dev"] < 1e-9
                                    and gt["median_rel_dev"] < 1e-9 and gt["curve_rel_dev"] < 1e-9 and gt["crown_x"] == gt["crown_x_R9"])
                    res.setdefault("gate_B1_240", {})[var] = gt; save_json(fres, res)
                    log(f"[B.1] porte 240² {var} contre R9 B (900²) : état {gt['E_res_state_same']}, E_res {gt['E_res_dev']:.1e}, pics {gt['peak_fine_dev']:.1e} / "
                        f"{gt['peak_prod_dev']:.1e}, médiane {gt['median_rel_dev']:.1e}, courbe {gt['curve_rel_dev']:.1e}, couronne {gt['crown_x']} (R9 {gt['crown_x_R9']}) "
                        f"-> {'OK' if gt['ok'] else 'ÉCHEC'}")
                    if not gt["ok"]:
                        log("[B.1] ÉCHEC de la porte 240² : STOP sur B.1"); sys.exit(3)
                save_json(fres, res)
            np.savez(cf, egrid=st["egrid"], E_D=st["E_D"], **{k_: v for k_, v in curves.items() if k_ not in ("egrid", "E_D")})
            log(f"[B.1] {N}² terminé en {time.time()-t0:.0f} s"); del g
    if "b2" in parts:
        res["B2"] = families(); save_json(fres, res)
    save_json(fres, res); b_tables(res)
    if "B1" in res:
        b_figure(res)


def xml_occupations(save):
    root = ET.parse(os.path.join(save, "data-file-schema.xml")).getroot(); bs = root.find(".//output/band_structure")
    ks = bs.findall("ks_energies"); assert len(ks) == 1, save
    occ = np.array(ks[0].find("occupations").text.split(), float)
    sm = root.find(".//output/band_structure/smearing")
    nel = float(bs.find("nelec").text)
    return occ, dict(smearing=(sm.text.strip() if sm is not None else None), degauss_Ha=(float(sm.get("degauss")) if sm is not None else None), nelec=nel,
                     k_weight=float(ks[0].find("k_point").get("weight")))


def families():
    """B.2 : 5…12, lecture des XML (E_F, gap à Γ de la parfaite, occupations) et des json (R5 C, R10 A.1, M_analysis, R9 clôture)."""
    cfg = load_production(verbose=False); r5 = load_json(R5_C)["sizes"]; C = cn_table()
    Zm = np.load(os.path.join(M2_DIR, "M_analysis.npz"), allow_pickle=True); a3p = load_json(R9_A3P)["level1"]
    rows = {}; gate_ok = True
    for n in range(5, 13):
        S = f"{n}x{n}"; P = sc_paths(S); rec = r5[str(n)]
        eP, efP, _ = qg.get_eigenvalues_spin(P["p"]); eP = eP[0] * HA2EV; efP *= HA2EV
        eD, efD, _ = qg.get_eigenvalues_spin(P["d"]); eD = eD[0] * HA2EV; efD *= HA2EV
        occ, meta = xml_occupations(P["d"]); _, metaP = xml_occupations(P["p"])
        n_occ = int(round(metaP["nelec"] / 2)); eS = np.sort(eP); homo, lumo = float(eS[n_occ - 1]), float(eS[n_occ])
        bpi = int(rec["pi_state"]["band"]); bsg = [int(q["band"]) for q in rec["sigma_doublet"]]
        chk = dict(E_F_p=abs(efP - rec["E_F_p"]), E_F_d=abs(efD - rec["E_F_d"]), gap=abs((lumo - homo) - rec["gamma_gap_p"]["gap"]),
                   pi_relEF=abs((eD[bpi - 1] - efD) - rec["pi_state"]["relEF"]), sigma_relEF=max(abs((eD[b - 1] - efD) - q["relEF"]) for b, q in zip(bsg, rec["sigma_doublet"])))
        ok = max(chk.values()) < 1e-9 and len(occ) == len(eD); gate_ok &= ok
        lv = a3p.get(S, {}).get("rc", {}).get("3", {}).get("variants", {}).get("aligne_plateau")
        rows[S] = dict(N=n, family="3m" if n % 3 == 0 else "non-3m", gap_gamma_p_eV=lumo - homo, homo_p=homo, lumo_p=lumo, E_F_p=efP, E_F_d=efD, dE_F_meV=(efD - efP) * 1e3,
                       pi=dict(band=bpi, e_minus_EF_d=float(eD[bpi - 1] - efD), occ=float(occ[bpi - 1])),
                       sigma=[dict(band=b, e_minus_EF_d=float(eD[b - 1] - efD), occ=float(occ[b - 1])) for b in bsg],
                       C_N_eV=C.get(S), max_M_dense_eV=float(Zm[f"scale_{S}_dense"]) if f"scale_{S}_dense" in Zm else None,
                       peak_GT_plateau_fine=(lv["peak_GT_fine"] if lv else None), peak_GT_plateau_prod=(lv["peak_GT_prod"] if lv else None),
                       smearing=meta, checks_vs_R5C=chk, ok=bool(ok))
        log(f"[B.2] {S} : gap Γ parfaite {lumo-homo:.4f} eV ; ΔE_F {(efD-efP)*1e3:+.1f} meV ; π ε − E_F {eD[bpi-1]-efD:+.4f} occ {occ[bpi-1]:.4f} ; σ "
            f"{[(round(float(eD[b-1]-efD), 4), round(float(occ[b-1]), 4)) for b in bsg]} ; C_N {C.get(S, float('nan'))*1e3 if C.get(S) is not None else float('nan'):+.2f} meV ; "
            f"contrôles R5 C {max(chk.values()):.1e} -> {'OK' if ok else 'ÉCART'}")
    return dict(rows=rows, gate_ok=bool(gate_ok), occupations_note="XML <ks_energies><occupations> de la cellule avec lacune, un k (Γ, poids 2), valeurs par état dans [0, 1]")


def b_tables(res):
    L = ["# R10 — B : argument E_res et familles (tables générées par r10_driver.py)", ""]
    g = res.get("gate_B0")
    if g:
        L += [f"Porte B.0 (9×9 tel quel, N_k^int 300, 240²) : pic de la courbe Γ_T {fr(g['peak_GT_prod'], -4)} eV (R6 −0,180), pic de −Im T̄(K) {fr(g['peak_ImTbar'], -4)} "
              f"(R6 −0,177), écart relatif à `resonance_9x9.npz` {g['curve_rel_dev_vs_R6']:.1e} (T̄ {g['Tbar_rel_dev_vs_R6']:.1e}) → {'OK' if g['ok'] else 'ÉCHEC'}", ""]
    for var, gt in res.get("gate_B1_240", {}).items():
        L += [f"Porte 240² ({var}) contre R9 B à N_k^int 900 : même état {gt['E_res_state_same']}, écarts E_res {gt['E_res_dev']:.1e}, pics {gt['peak_fine_dev']:.1e} / "
              f"{gt['peak_prod_dev']:.1e}, médiane {gt['median_rel_dev']:.1e}, courbe {gt['curve_rel_dev']:.1e}, couronne x = {gt['crown_x']} (R9 {gt['crown_x_R9']}) → "
              f"{'OK' if gt['ok'] else 'ÉCHEC'}", ""]
    if "B1" in res:
        m = res["meta"]
        L += [f"## B.1 — 9×9, R_cut 3, η 0,02 eV, N_k^int {B_NK}, chaîne « res » (grille d'énergie de g₀ au pas de 2,5 meV), fenêtre des états ±3 eV, E_res = argmax |Γ| des états "
              f"à ±1,5 eV ; plateau (i) : C_9 = {fr(m['C9_eV']*1e3, -4)} meV ({m['C9_source']}). Énergies en eV relatives à E_D (Wannier).", "",
              "| variante | grille | états | E_res ; état (k, bande) | couronne de E_res : x (rang) ; [min, max] | couronne voisine en dessous : x ; moyenne ; E_res − moyenne | "
              "couronne voisine au-dessus : x ; moyenne ; E_res − moyenne | pic Γ_T 5 meV ; 2,5 meV | \\|E_res − pic\\| 5 ; 2,5 meV (meV) | médiane états (meV) |",
              "|---|---|---|---|---|---|---|---|---|---|"]
        for var in ("brut", "plateau"):
            for N, r in sorted(res["B1"].get(var, {}).items(), key=lambda kv: int(kv[0])):
                c = r["crown"]; nb = c["neighbours"]
                def nbs(lab):
                    q = nb.get(lab)
                    return f"{q['x']} ; {fr(q['E_mean'], -4)} ; {fr(q['dE_res_minus_mean']*1e3, -1)} meV" if q else "—"
                L.append(f"| {'tel quel' if var == 'brut' else 'plateau (i)'} | {N}² | {r['n_states']} | {fr(r['E_res'], -5)} ; ({r['E_res_state'][0]}, {r['E_res_state'][1]}) | "
                         f"{c['x']} ({c['rank']}) ; [{fr(c['E_min'], -4)} ; {fr(c['E_max'], -4)}] | {nbs('below')} | {nbs('above')} | {fr(r['peak_GT_prod'], -4)} ; "
                         f"{fr(r['peak_GT_fine'], -4)} | {fr(r['d_Eres_peak_prod_meV'], 1)} ; {fr(r['d_Eres_peak_fine_meV'], 1)} | {fr(r['median_GT_states_meV'], 2)} |")
        L.append("")
    if "B2" in res:
        L += ["## B.2 — familles, 5…12 (XML et json ; ε − E_F de la cellule avec lacune, non aligné ; occupations par état, poids de Γ 2)", "",
              f"Contrôle contre R5 C (E_F, gap, ε − E_F des états π et σ) : {'OK' if res['B2']['gate_ok'] else 'ÉCART'}.", "",
              "| N | famille | gap de la parfaite à Γ (eV) | ΔE_F = E_F(d) − E_F(p) (meV) | π : ε − E_F (eV) ; occupation | σ : ε − E_F (eV) ; occupations | C_N plateau (meV) | "
              "max\\|M\\| dense, bandes 1–16 (eV) | pic de Γ_T aligné plateau, 2,5 ; 5 meV (eV, R9 clôture) |", "|---|---|---|---|---|---|---|---|---|"]
        for S, r in res["B2"]["rows"].items():
            L.append(f"| {S.replace('x', '×')} | {r['family']} | {fr(r['gap_gamma_p_eV'], 4)} | {fr(r['dE_F_meV'], -1)} | {fr(r['pi']['e_minus_EF_d'], -4)} ; {fr(r['pi']['occ'], 4)} | "
                     f"{' / '.join(fr(q['e_minus_EF_d'], -4) for q in r['sigma'])} ; {' / '.join(fr(q['occ'], 4) for q in r['sigma'])} | "
                     f"{fr(r['C_N_eV']*1e3, -2) if r['C_N_eV'] is not None else '—'} | {fr(r['max_M_dense_eV'], 3) if r['max_M_dense_eV'] is not None else '—'} | "
                     + (f"{fr(r['peak_GT_plateau_fine'], -4)} ; {fr(r['peak_GT_plateau_prod'], -4)}" if r["peak_GT_plateau_fine"] is not None else "—") + " |")
    with open(os.path.join(ensure("b"), "B_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def b_figure(res):
    plt, pal = fig_style()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.3))
    for var, col, lab in (("brut", pal.ORANGE, "tel quel"), ("plateau", pal.NAVY, "plateau (i)")):
        rr = sorted(res["B1"].get(var, {}).items(), key=lambda kv: int(kv[0]))
        if not rr:
            continue
        x = [1.0 / int(N) for N, _ in rr]
        a1.plot(x, [r["E_res"] for _, r in rr], "o", ms=4, color=col, label=rf"$E_\mathrm{{res}}$, {lab}")
        a1.plot(x, [r["peak_GT_fine"] for _, r in rr], "-s", ms=2.5, lw=1.0, mfc="none", color=col, label=rf"pic de $\Gamma_T$ (2,5 meV), {lab}")
        a2.plot(x, [r["d_Eres_peak_fine_meV"] for _, r in rr], "o-", ms=3.5, lw=1.0, color=col, label=lab)
    a1.set_xlabel(r"$1/N_\mathrm{out}$"); a1.set_ylabel(r"Énergie $\varepsilon - E_D$ (eV)"); a1.set_title(r"(a) $E_\mathrm{res}$ et pic de la courbe", loc="left", fontsize=9)
    a2.set_xlabel(r"$1/N_\mathrm{out}$"); a2.set_ylabel(r"$|E_\mathrm{res} - \mathrm{pic}|$ (meV)"); a2.set_title("(b) écart au pic (2,5 meV)", loc="left", fontsize=9)
    a1.legend(fontsize=6); a2.legend(fontsize=6)
    fig.tight_layout(); savefig(fig, "eres_vs_grid"); plt.close(fig)


def cmd_btables(a):
    res = load_json(os.path.join(WORK, "b", "b_results.json")); b_tables(res)
    if "B1" in res:
        b_figure(res)


# ----------------------------------------------------------------------------------------------- C : rejeu de production aligné (GO 2)
# Chaîne de production : lt.defect_mwr (R10 (a)), C_N = config["alignment"] (alignment_C), matrices M2 brutes (matrices_dir), produits de C1 dans
# results_dir (results/M2_plateau). Références : results/M2 (v2 tel quel, gelé), json de R9 (A.3, D, audit). Sorties du pilote dans c/ et fig/.
C_SIZES = ["9x9", "12x12"]                                                     # portes C.0 (a, b) et tab:rcut_M


def lu9():
    """Lu (1,0 Å) de la 9×9 : C_retenu_eV de R9 A.1 (étiquette « sans plateau : valeur Lu »), C de R9 D.1 « exact » (porte C.0 d)."""
    r = load_json(R9_A1)["sizes"]["9x9"]
    assert "Lu" in r["etiquette"], r["etiquette"]
    return float(r["C_retenu_eV"])


def c_inputs(cfg, S):
    """
    Entrées de la chaîne de production d'une taille dense : M2 (eV, portes de norme de matrix_io), k du XML (compute_spectral_wannier, rcut_resigma,
    resonance_*) et k ramenés sur la grille MP (m_rcut_convergence, analyze_M), U et U_dis réordonnés pour chacun (manifeste : porte de jauge),
    H(R) (read_w90_tb), cellule (bohr).
    """
    t0 = time.time(); dp = dense_paths(cfg, S)
    M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k); kmp = np.round(k * np.asarray(MP)) / np.asarray(MP)
    U, Ud, Hwr, Rw, nd = load_wannier(dp, k); Ump, Udmp, _, _, _ = load_wannier(dp, kmp)
    A_uc, _ = qe_io.get_A_volume(dp["uc"])
    log(f"[C entrées] {S} : M2 {M.shape} ({os.path.relpath(dp['mfile'], PROJ)}), MP {tuple(int(x) for x in MP)}, max|k_XML − k_MP| {np.abs(k - kmp).max():.1e} ; "
        f"{time.time()-t0:.0f} s")
    return dict(S=S, n=size_n(S), dp=dp, M=M, k=k, kmp=kmp, MP=MP, U=U, Ud=Ud, Ump=Ump, Udmp=Udmp, Hwr=Hwr, Rw=Rw, nd=nd, A_uc=A_uc)


def dirac_point(Hwr, Rw, nd):
    """E_D de production : milieu de π/π* au minimum du gap sur la grille 90² (compute_spectral_wannier.py, rcut_resigma.py)."""
    _, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
    gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap))
    return float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))


def level1_prod(cfg, inp, Mwr, Rn, rc=None, grid=None, eta=None):
    """Niveau 1 par la fonction de production lt.scattering_rate_fast, entrées de compute_spectral_wannier.py (R_cut, 240², η 0,02, N_k^int 300 par défaut)."""
    rc = cfg["R_cut"] if rc is None else rc; grid = cfg["grid"] if grid is None else grid; eta = cfg["eta_eV"] if eta is None else eta
    t0 = time.time(); Hwr, Rw, nd = inp["Hwr"], inp["Rw"], inp["nd"]
    Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]; V_loc, herm = lt.extract_V_loc(Mwr, Rn, Rloc)
    E_D = dirac_point(Hwr, Rw, nd); win = (E_D - float(cfg["e_window_eV"]), E_D + float(cfg["e_window_eV"]))
    k_int = lt.mp_grid(int(cfg["nk_int"]), int(cfg["nk_int"]), 1); k_out = lt.mp_grid(grid, grid, 1)
    _, E_out, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=nd)
    gamma = lt.scattering_rate_fast(Hwr, Rw, nd, V_loc, Rloc, k_out, eta, k_int=k_int, e_window=win, ne_per_eta=int(cfg["ne_per_eta"]))
    med = float(np.nanmedian(np.abs(gamma))) * 1e3
    E = E_out.T; m = np.isfinite(gamma) & (np.abs(E - E_D) <= 1.5)
    e_res = float(E[m][np.argmax(np.abs(gamma)[m])] - E_D)
    return dict(median_meV=med, E_res=e_res, E_D=E_D, n_sites=int(len(Rloc)), herm=float(herm), V=V_loc, Rloc=Rloc, elapsed_s=time.time() - t0)


def fine_pi(Mwr, nw, nR, Upi, Pm, Pp, mask):
    """π/π* de M sur la grille fine avec R, R' dans `mask` (m_rcut_convergence.py, recopiée telle quelle)."""
    nk = Pm.shape[0]
    Wm = Mwr.copy(); keep = np.asarray(mask, bool); Wm[:, ~keep, :, :] = 0; Wm[:, :, :, ~keep] = 0
    A = Wm.transpose(1, 0, 3, 2).reshape(nR * nw, nR * nw)
    X = A.reshape(nR * nw, nR, nw); Y = np.einsum("aRw,kR->akw", X, Pp, optimize=True).reshape(nR, nw, nk, nw)
    Y = np.einsum("Rwkv,kvn->Rwkn", Y, Upi, optimize=True)
    Z = np.einsum("Rwkn,KR->Kwkn", Y, Pm, optimize=True)
    return np.einsum("Kwkn,Kwm->mKnk", Z, Upi.conj(), optimize=True)


def rcut_rows(inp, Mwr, Rn, rcuts=range(7), nf=60, ws=None):
    """
    tab:rcut_M (m_rcut_convergence.py : grille fine nf², paire π/π* de H_W(k), R_cut en norme réduite sur les étiquettes recentrées) ; ws=None :
    phases e^{∓2πi k·Rn} (« étiquettes actuelles ») ; sinon ws_phase (Wigner-Seitz, comme m_rcut_convergence.py de (b)). Lignes au format du csv.
    """
    nw, nR = Mwr.shape[0], Mwr.shape[1]
    kf = lt.mp_grid(nf, nf, 1); nk = len(kf)
    _, Ef, Uf = lt.Hwr_to_Hwk(inp["Hwr"], inp["Rw"], kf, ndegen=inp["nd"]); Upi = np.ascontiguousarray(Uf[:, :, 3:5])
    if ws is None:
        Pm = np.exp(-2j * np.pi * (kf @ Rn.T)); Pp = np.exp(+2j * np.pi * (kf @ Rn.T))
    else:
        Pm = ws_phase(kf, ws, nR, -1); Pp = ws_phase(kf, ws, nR, +1)
    dist = np.linalg.norm(Rn, axis=1)
    ref = fine_pi(Mwr, nw, nR, Upi, Pm, Pp, np.ones(nR, bool)); mref = np.abs(ref).max()
    d_ref = np.array([[ref[n, i, n, i] for i in range(nk)] for n in range(2)]); mdref = np.abs(d_ref).max()
    sv_ref = np.linalg.svd(ref.transpose(1, 3, 0, 2), compute_uv=False)
    rows = []
    for rc in rcuts:
        mask = dist <= rc + 1e-9; Mt = fine_pi(Mwr, nw, nR, Upi, Pm, Pp, mask)
        e_max = np.abs(Mt - ref).max() / mref
        sv = np.linalg.svd(Mt.transpose(1, 3, 0, 2), compute_uv=False); e_sv = (np.abs(sv - sv_ref).max(-1) / np.maximum(sv_ref[:, :, 0], 1e-12)).max()
        d_t = np.array([[Mt[n, i, n, i] for i in range(nk)] for n in range(2)]); e_diag = np.abs(d_t - d_ref).max() / mdref
        e_diag_sv = np.abs(np.abs(d_t) - np.abs(d_ref)).max() / mdref
        rows.append(dict(size=inp["S"], R_cut=int(rc), n_sites=int(mask.sum()), fine_grid=int(nf), max_dM_over_maxM=float(e_max), sv_mismatch_pi_blocks=float(e_sv),
                         diag_max_dM_over_max=float(e_diag), diag_abs_mismatch=float(e_diag_sv), maxM_pi_eV=float(mref), maxMdiag_eV=float(mdref)))
        del Mt
    return rows


def rcut_csv_strings(r):
    """Une ligne de rcut_rows au format exact de m_rcut_convergence.csv (.4e, .4f)."""
    return [r["size"], str(r["R_cut"]), str(r["n_sites"]), str(r["fine_grid"]), f"{r['max_dM_over_maxM']:.4e}", f"{r['sv_mismatch_pi_blocks']:.4e}",
            f"{r['diag_max_dM_over_max']:.4e}", f"{r['diag_abs_mismatch']:.4e}", f"{r['maxM_pi_eV']:.4f}", f"{r['maxMdiag_eV']:.4f}"]


def pi_pair_U(Uk):
    """Bandes π, π* (ordre des énergies) : les deux états propres de plus grand poids p_z (r9_driver.pi_pair_U)."""
    w = np.abs(Uk[..., 3, :]) ** 2 + np.abs(Uk[..., 4, :]) ** 2
    top = np.sort(np.argsort(-w, axis=-1)[..., :2], axis=-1)
    return top[..., 0], top[..., 1]


def d1_setup(inp):
    """Carte D.1 de R9 (cmd_d) : k = K + δx̂ (δ = 0,01 |b₁|), carte k' 240², K et K' exacts exclus, disques de rayon 0,05 |b₁| autour de K et K'."""
    A_A = inp["A_uc"] * BOHR; A_cell = float(np.linalg.norm(np.cross(A_A[:, 0], A_A[:, 1])))
    Bc = 2 * np.pi * np.linalg.inv(A_A).T; b1 = float(np.linalg.norm(Bc[:, 0]))
    delta = 0.01 * b1; xhat = np.array([1.0, 0.0, 0.0]); k = K_RED + np.linalg.solve(Bc, delta * xhat)
    kdist = lambda kk, k0: np.linalg.norm(np.mod(kk - np.asarray(k0) + 0.5, 1.0) - 0.5, axis=1)
    kp = lt.mp_grid(240, 240, 1); exK = kdist(kp, K_RED) < 1e-9; exKp = kdist(kp, KP_RED) < 1e-9
    Ek, Uk = lt.Hwr_to_Hwk(inp["Hwr"], inp["Rw"], k[None], ndegen=inp["nd"])[1:]
    Ep, Up = lt.Hwr_to_Hwk(inp["Hwr"], inp["Rw"], kp, ndegen=inp["nd"])[1:]
    pk_v, pk_c = pi_pair_U(Uk[0][None]); pv, pc_ = pi_pair_U(Up)
    kc = (Bc[:2, :2] @ kp[:, :2].T).T; imgs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)]) @ Bc[:2, :2].T

    def disk(K0):
        K0c = Bc[:2, :2] @ K0[:2]; dd = np.min(np.linalg.norm(kc[:, None, :] - K0c[None, None, :] + imgs[None], axis=2), axis=1)
        return dd <= 0.05 * b1
    return dict(A_cell=A_cell, Bc=Bc, b1=b1, delta=delta, k=k, kp=kp, exK=exK, exKp=exKp, Ek=Ek, Uk=Uk, Up=Up, pk_v=pk_v, pk_c=pk_c, pv=pv, pc=pc_,
                kc=kc, dK=disk(K_RED) & ~exK, dKp=disk(KP_RED) & ~exKp)


def d1_map(D1, Mw, R, ws=None):
    """Ṽ = A_cell |M^{nn}_{k'k}| (valence π, conduction π*) sur la carte ; moyennes, min, max par disque (r9_driver.cmd_d)."""
    kp = D1["kp"]
    Mk = Mwr_to_Mwk_pairs(Mw, R, kp, D1["k"][None], ws=ws)[:, :, :, 0]
    Mb = np.einsum("kwb,wkW,WB->bkB", D1["Up"].conj(), Mk, D1["Uk"][0], optimize=True)
    Vv = D1["A_cell"] * np.abs(Mb[D1["pv"], np.arange(len(kp)), D1["pk_v"][0]]); Vc = D1["A_cell"] * np.abs(Mb[D1["pc"], np.arange(len(kp)), D1["pk_c"][0]])
    Vv[D1["exK"] | D1["exKp"]] = np.nan; Vc[D1["exK"] | D1["exKp"]] = np.nan
    stats = {}
    for lab, m in (("K", D1["dK"]), ("K'", D1["dKp"])):
        for band, Vm in (("valence", Vv), ("conduction", Vc)):
            stats[f"{band} | {lab}"] = dict(mean=float(np.nanmean(Vm[m])), min=float(np.nanmin(Vm[m])), max=float(np.nanmax(Vm[m])), n=int(np.isfinite(Vm[m]).sum()))
    return dict(valence=Vv, conduction=Vc, stats=stats)


def max_mean_diff(sa, sb):
    return float(max(abs(sa[k]["mean"] - sb[k]["mean"]) for k in sb))


def cmd_c0(a):
    """
    Portes C.0 : (a) C_N = 0, étiquettes actuelles (9×9) ; (b) plateau (i), étiquettes actuelles (9×9, 12×12), V_loc aligné = R9 au bit ;
    (c) Wigner-Seitz : M inchangé sur la grille MP, D.1 tel quel = « vraie » de l'audit ; (d) approximation (i) avec C = Lu contre R9 D.1 exact (rapporté).
    Échec de (a), (b) ou (c) : code 3 (STOP sur C).
    """
    cfg = load_production(verbose=True); d = ensure("c"); f = os.path.join(d, "c0_results.json")
    a3 = load_json(R9_A3); a3p = load_json(R9_A3P); out = dict(head=git_head(), gates={}); ok = {}
    TOL_MED = 0.01                                                             # meV (0.7 : « 3 132,60 meV (0,01 meV) »)
    for S in C_SIZES:
        inp = c_inputs(cfg, S); n = inp["n"]; C = alignment_C(cfg, S); M = inp["M"]; res = dict(C_N_eV=C)
        # chaîne sans alignement recopiée (Mbk_to_Mwk -> Mwk_to_Mwr -> recenter_mwr, k du XML : mwr_cached de R9) contre defect_mwr
        Mw0, R0 = Mwk_to_Mwr(Mbk_to_Mwk(M, inp["U"], inp["Ud"]), inp["k"], inp["MP"]); Rn0, Rd0 = lt.recenter_mwr(Mw0, R0, inp["MP"])
        dz = lt.defect_mwr(M, inp["U"], inp["Ud"], inp["k"], inp["MP"], n, C_N=0.0)
        dpl = lt.defect_mwr(M, inp["U"], inp["Ud"], inp["k"], inp["MP"], n, C_N=C)
        res["defect_mwr_C0_bitwise"] = bool(np.array_equal(dz["Mwr"], Mw0) and np.array_equal(dz["Rn"], Rn0) and np.array_equal(dz["R_d"], Rd0))
        res["R_d_plateau_eq_C0"] = bool(np.array_equal(dpl["R_d"], Rd0) and np.array_equal(dpl["Rn"], Rn0))
        del Mw0
        # (a) / (b) : niveau 1 (scattering_rate_fast) et V_loc
        if S == "9x9":
            L0 = level1_prod(cfg, inp, dz["Mwr"], dz["Rn"]); r9 = a3["level1"][S]["rc"]["3"]["variants"]["brut"]
            res["a_level1"] = dict(median_meV=L0["median_meV"], E_res=L0["E_res"], R9_median_meV=r9["median_meV"], R9_E_res=r9["E_res"],
                                   d_median_meV=L0["median_meV"] - r9["median_meV"], d_E_res=L0["E_res"] - r9["E_res"], elapsed_s=L0["elapsed_s"])
            ok[f"a niveau 1 {S}"] = abs(L0["median_meV"] - r9["median_meV"]) <= TOL_MED and abs(L0["E_res"] - r9["E_res"]) <= 1e-6
            log(f"[C.0 a] {S} C_N = 0 : médiane {L0['median_meV']:.4f} meV (R9 {r9['median_meV']:.4f}), E_res {L0['E_res']:+.6f} (R9 {r9['E_res']:+.6f}) ; {L0['elapsed_s']:.0f} s")
            V0 = L0["V"]; Rloc0 = L0["Rloc"]; del L0
        else:
            Rloc0 = dz["Rn"][np.linalg.norm(dz["Rn"], axis=1) <= cfg["R_cut"] + 1e-9]; V0, _ = lt.extract_V_loc(dz["Mwr"], dz["Rn"], Rloc0)
        Lp = level1_prod(cfg, inp, dpl["Mwr"], dpl["Rn"]); r9 = a3p["level1"][S]["rc"]["3"]["variants"]["aligne_plateau"]
        ib = np.array([dpl["in_box"][np.where((dpl["Rn"] == r).all(1))[0][0]] for r in Lp["Rloc"]])
        v_bit = bool(np.array_equal(Lp["Rloc"], Rloc0) and np.array_equal(Lp["V"], V0 - C * np.diag(np.repeat(ib.astype(float), NW))))
        res["b_level1"] = dict(median_meV=Lp["median_meV"], E_res=Lp["E_res"], R9_median_meV=r9["median_meV"], R9_E_res=r9["E_res"],
                               d_median_meV=Lp["median_meV"] - r9["median_meV"], d_E_res=Lp["E_res"] - r9["E_res"], elapsed_s=Lp["elapsed_s"],
                               Vloc_eq_R9_bitwise=v_bit, n_sites=Lp["n_sites"], n_out_of_box=int((~ib).sum()))
        ok[f"b niveau 1 {S}"] = abs(Lp["median_meV"] - r9["median_meV"]) <= TOL_MED and abs(Lp["E_res"] - r9["E_res"]) <= 1e-6
        ok[f"b V_loc au bit {S}"] = v_bit and res["defect_mwr_C0_bitwise"] and res["R_d_plateau_eq_C0"]
        log(f"[C.0 b] {S} plateau (C_N {C*1e3:+.4f} meV) : médiane {Lp['median_meV']:.4f} meV (R9 {r9['median_meV']:.4f}), E_res {Lp['E_res']:+.6f} "
            f"(R9 {r9['E_res']:+.6f}) ; V_loc = R9 au bit : {v_bit} ; defect_mwr(C_N = 0) = chaîne recopiée au bit : {res['defect_mwr_C0_bitwise']} ; {Lp['elapsed_s']:.0f} s")
        del Lp, V0, dz, dpl
        # tab:rcut_M, R_cut 3, étiquettes actuelles (k ramenés sur la grille MP, comme m_rcut_convergence.py)
        for lab, Cv, key, ref in (("a", 0.0, "brut", a3.get(f"rcut_{S}", {}).get("rows", {}).get("brut")),
                                  ("b", C, "aligne_plateau", a3p[f"rcut_{S}"]["rows"]["aligne_plateau"])):
            if lab == "a" and S != "9x9":
                continue
            dm = lt.defect_mwr(M, inp["Ump"], inp["Udmp"], inp["kmp"], inp["MP"], n, C_N=Cv)
            r3 = rcut_rows(inp, dm["Mwr"], dm["Rn"], rcuts=[3])[0]; r9v = ref[3]["max_dM_over_maxM"]
            res[f"{lab}_rcut3"] = dict(max_dM_over_maxM=r3["max_dM_over_maxM"], R9=r9v, rel=abs(r3["max_dM_over_maxM"] - r9v) / r9v)
            ok[f"{lab} tab:rcut_M {S}"] = abs(r3["max_dM_over_maxM"] - r9v) / r9v <= 1e-9
            log(f"[C.0 {lab}] {S} tab:rcut_M R_cut 3 ({key}) : {r3['max_dM_over_maxM']:.6e} (R9 {r9v:.6e}, écart relatif {res[f'{lab}_rcut3']['rel']:.1e})")
            del dm
        out[S] = res; save_json(f, dict(out, gates=ok))
        if S != "9x9":
            del inp, M
            continue
        # (c) Wigner-Seitz, 9×9 : M inchangé sur la grille MP (étiquettes brutes autour de R_d et recentrées autour de 0)
        dm = lt.defect_mwr(M, inp["Ump"], inp["Udmp"], inp["kmp"], inp["MP"], n, C_N=C)
        ws_raw = ws_images(dm["R"], dm["R_d"], inp["MP"], inp["A_uc"]); ws_rec = ws_images(dm["Rn"], np.zeros(3, int), inp["MP"], inp["A_uc"])
        M_plain = Mwr_to_Mwk(dm["Mwr"], dm["R"], inp["kmp"]); sc = np.abs(M_plain).max()
        c_raw = float(np.abs(Mwr_to_Mwk(dm["Mwr"], dm["R"], inp["kmp"], ws=ws_raw) - M_plain).max() / sc)
        M_rec = Mwr_to_Mwk(dm["Mwr"], dm["Rn"], inp["kmp"])
        c_rec = float(np.abs(Mwr_to_Mwk(dm["Mwr"], dm["Rn"], inp["kmp"], ws=ws_rec) - M_rec).max() / sc)
        del M_plain, M_rec, dm
        out["c_onGrid"] = dict(raw_labels_rel=c_raw, recentred_labels_rel=c_rec, n_tie_labels=int((ws_raw["n_tie"] > 1).sum()), n_images=int(len(ws_raw["w"])))
        ok["c M sur la grille MP"] = c_raw <= 1e-12 and c_rec <= 1e-12
        log(f"[C.0 c] 9x9 : Mwr_to_Mwk avec ws = sans ws sur la grille MP : étiquettes brutes {c_raw:.1e}, recentrées {c_rec:.1e} ; "
            f"{out['c_onGrid']['n_tie_labels']} étiquettes à égalité, {out['c_onGrid']['n_images']} images")
        # D.1 (k du XML, étiquettes brutes de Mwk_to_Mwr : pilote D.1 de R9) : tel quel sans ws = R9 D.1 ; tel quel avec ws (autour de R_d) = « vraie | brut » de l'audit
        D1 = d1_setup(inp); dd = load_json(R9_D); aud = load_json(R9_AUDIT)["interp_D1"]["carte"]["disk_means"]
        k_ok = float(np.abs(D1["k"] - np.array(dd["D1_meta"]["k_red"])).max())
        dz = lt.defect_mwr(M, inp["U"], inp["Ud"], inp["k"], inp["MP"], n, C_N=0.0)
        ws_d = ws_images(dz["R"], dz["R_d"], inp["MP"], inp["A_uc"])
        m_plain = d1_map(D1, dz["Mwr"], dz["R"]); m_ws = d1_map(D1, dz["Mwr"], dz["R"], ws=ws_d)
        g_plain = max_mean_diff(m_plain["stats"], dd["D1"]["brut"]); g_ws = max_mean_diff(m_ws["stats"], aud["vraie | brut"])
        out["c_D1"] = dict(k_vs_R9=k_ok, sans_ws=m_plain["stats"], avec_ws=m_ws["stats"], gate_sans_ws_vs_R9_D1=g_plain, gate_ws_vs_audit_vraie=g_ws)
        ok["c D.1 sans ws = R9 D.1"] = g_plain <= 1e-9 and k_ok <= 1e-12
        ok["c D.1 ws = audit « vraie »"] = g_ws <= 1e-9
        log(f"[C.0 c] D.1 tel quel : sans ws, valence K {m_plain['stats']['valence | K']['mean']:.6f} eV Å² (R9 {dd['D1']['brut']['valence | K']['mean']:.6f}), "
            f"écart max des moyennes {g_plain:.1e} ; avec ws {m_ws['stats']['valence | K']['mean']:.6f} (audit « vraie » {aud['vraie | brut']['valence | K']['mean']:.6f}), "
            f"écart max {g_ws:.1e} ; |k − k(R9)| {k_ok:.1e}")
        # (d) approximation (i) avec C = Lu, étiquettes brutes, contre R9 D.1 exact (F_W) : écart rapporté
        Clu = lu9(); dl = lt.defect_mwr(M, inp["U"], inp["Ud"], inp["k"], inp["MP"], n, C_N=Clu)
        m_lu = d1_map(D1, dl["Mwr"], dl["R"])
        out["d_D1"] = dict(C_Lu_eV=Clu, approx_i=m_lu["stats"], R9_exact=dd["D1"]["exact"],
                           diff={k: m_lu["stats"][k]["mean"] - dd["D1"]["exact"][k]["mean"] for k in dd["D1"]["exact"]})
        log(f"[C.0 d] D.1, (i) avec C = Lu ({Clu*1e3:+.4f} meV), étiquettes brutes : " + " ; ".join(
            f"{k} {m_lu['stats'][k]['mean']:.4f} (R9 exact {dd['D1']['exact'][k]['mean']:.4f}, écart {out['d_D1']['diff'][k]:+.4f})" for k in dd["D1"]["exact"]) + " eV Å²")
        del dz, dl, D1, m_plain, m_ws, m_lu, inp, M
    out["gates"] = ok; out["all_ok"] = bool(all(ok.values())); save_json(f, out); c0_tables(out)
    log("[C.0] portes : " + " ; ".join(f"{k} {'OK' if v else 'ÉCHEC'}" for k, v in ok.items()))
    if not out["all_ok"]:
        log("[C.0] ÉCHEC d'une porte : STOP sur C (code 3)"); sys.exit(3)


def c0_tables(out):
    g = out["gates"]; L = ["# R10 — C.0 : portes (tables générées par r10_driver.py c0)", ""]
    L += ["| porte | résultat |", "|---|---|"] + [f"| {k} | {'OK' if v else 'ÉCHEC'} |" for k, v in g.items()] + [""]
    L += ["| taille | variante | médiane Γ·N_cells (meV) ; R9 | E_res (eV) ; R9 | tab:rcut_M R_cut 3 ; R9 (écart relatif) |", "|---|---|---|---|---|"]
    for S in C_SIZES:
        r = out.get(S, {})
        if "a_level1" in r:
            x = r["a_level1"]; y = r.get("a_rcut3", {})
            L.append(f"| {S} | C_N = 0 | {fr(x['median_meV'], 4)} ; {fr(x['R9_median_meV'], 4)} | {fr(x['E_res'], -6)} ; {fr(x['R9_E_res'], -6)} | "
                     f"{y.get('max_dM_over_maxM', float('nan')):.6e} ; {y.get('R9', float('nan')):.6e} ({y.get('rel', float('nan')):.1e}) |")
        if "b_level1" in r:
            x = r["b_level1"]; y = r.get("b_rcut3", {})
            L.append(f"| {S} | plateau (C_N {fr(r['C_N_eV']*1e3, -4)} meV) | {fr(x['median_meV'], 4)} ; {fr(x['R9_median_meV'], 4)} | {fr(x['E_res'], -6)} ; {fr(x['R9_E_res'], -6)} | "
                     f"{y.get('max_dM_over_maxM', float('nan')):.6e} ; {y.get('R9', float('nan')):.6e} ({y.get('rel', float('nan')):.1e}) |")
    if "c_onGrid" in out:
        c = out["c_onGrid"]; L += ["", f"(c) Wigner-Seitz sur la grille MP 27² : écart relatif {c['raw_labels_rel']:.1e} (étiquettes brutes), {c['recentred_labels_rel']:.1e} (recentrées) ; "
                                   f"{c['n_tie_labels']} étiquettes à égalité, {c['n_images']} images."]
    if "c_D1" in out:
        c = out["c_D1"]; L += ["", "| D.1 (eV Å²) | valence K | conduction K | valence K′ | conduction K′ |", "|---|---|---|---|---|"]
        keys = ["valence | K", "conduction | K", "valence | K'", "conduction | K'"]
        L.append("| tel quel, sans ws | " + " | ".join(f"{c['sans_ws'][k]['mean']:.6f}" for k in keys) + " |")
        L.append("| tel quel, Wigner-Seitz | " + " | ".join(f"{c['avec_ws'][k]['mean']:.6f}" for k in keys) + " |")
        if "d_D1" in out:
            dd = out["d_D1"]
            L.append(f"| (i) C = Lu ({fr(dd['C_Lu_eV']*1e3, -4)} meV), sans ws | " + " | ".join(f"{dd['approx_i'][k]['mean']:.6f}" for k in keys) + " |")
            L.append("| R9 D.1 exact (F_W) | " + " | ".join(f"{dd['R9_exact'][k]['mean']:.6f}" for k in keys) + " |")
            L.append("| écart (i) − exact | " + " | ".join(f"{dd['diff'][k]:+.6f}" for k in keys) + " |")
        L.append(f"\nÉcarts max des moyennes : sans ws − R9 D.1 {c['gate_sans_ws_vs_R9_D1']:.1e} ; ws − audit « vraie | brut » {c['gate_ws_vs_audit_vraie']:.1e} eV Å².")
    with open(os.path.join(WORK, "c", "C0_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")

# ----------------------------------------------------------------------------------------------- C.1 (fin) : m_rcut_resigma.csv, lignes de tab:tests_M
MRR_COLS = ["size", "R_cut", "nL", "dim", "grid", "eta_eV", "window_eV", "E_D_eV", "med_ReSigma_meV", "mean_ReSigma_meV", "med_Gamma_meV",
            "mean_Gamma_meV", "med_absReSigma_over_Gamma", "med_abs_dReSigma_vs_rc4_meV", "max_abs_dReSigma_vs_rc4_meV", "med_abs_dGamma_vs_rc4_meV",
            "max_abs_dGamma_vs_rc4_meV", "rel_med_dReSigma", "rel_med_dGamma", "E_res_minus_E_D_eV", "n_states"]


def mrr_table(dirpath, S="9x9"):
    """m_rcut_resigma.csv à partir de resigma_<S>_rc0123.npz et resigma_<S>_rc4.npz (article/R6_production_corrigee/etape3/r6_m_rcut_resigma.py, recopiée)."""
    A = np.load(f"{dirpath}/resigma_{S}_rc0123.npz"); B = np.load(f"{dirpath}/resigma_{S}_rc4.npz")
    Sg = {rc: (A if rc < 4 else B)[f"Sigma_rc{rc}"] for rc in range(5)}; nL = {rc: int((A if rc < 4 else B)[f"nL_rc{rc}"]) for rc in range(5)}
    E = A["E_out"]; ED = float(A["E_D"]); m = ~np.isnan(Sg[3]); nw = Sg[3].shape[0]
    assert np.array_equal(m, ~np.isnan(Sg[4])) and np.allclose(E, B["E_out"], equal_nan=True) and float(B["E_D"]) == ED
    G4, R4 = -2 * Sg[4].imag[m], Sg[4].real[m]; rows = []
    for rc in range(5):
        G, R = -2 * Sg[rc].imag[m], Sg[rc].real[m]; mm = np.abs(E[m] - ED) <= 1.5
        rows.append(dict(size=S, R_cut=rc, nL=nL[rc], dim=nL[rc] * nw, grid=int(A["grid"]), eta_eV=float(A["eta"]), window_eV=float(A["e_window"]), E_D_eV=ED,
                         med_ReSigma_meV=np.median(R) * 1e3, mean_ReSigma_meV=R.mean() * 1e3, med_Gamma_meV=np.median(G) * 1e3, mean_Gamma_meV=G.mean() * 1e3,
                         med_absReSigma_over_Gamma=np.median(np.abs(R) / G),
                         med_abs_dReSigma_vs_rc4_meV=np.median(np.abs(R - R4)) * 1e3, max_abs_dReSigma_vs_rc4_meV=np.abs(R - R4).max() * 1e3,
                         med_abs_dGamma_vs_rc4_meV=np.median(np.abs(G - G4)) * 1e3, max_abs_dGamma_vs_rc4_meV=np.abs(G - G4).max() * 1e3,
                         rel_med_dReSigma=np.median(np.abs(R - R4)) / np.median(R4), rel_med_dGamma=np.median(np.abs(G - G4)) / np.median(G4),
                         E_res_minus_E_D_eV=float(E[m][mm][np.argmax(G[mm])] - ED), n_states=int(m.sum())))
    return rows


def cmd_c1post(a):
    """
    Fin de C.1 : (1) m_rcut_resigma.csv — porte : les définitions redonnent results/M2/m_rcut_resigma.csv depuis les npz de M2 à 1e-9 (sinon rien n'est
    écrit, code 3) ; (2) M_tests_summary.csv : ligne du test d'or de M2 recopiée (D12, valeur de R6) et lignes de la porte A.2 de M2 (mêmes matrices).
    """
    import csv
    cfg = load_production(verbose=True); RES = results_dir(cfg); out = dict(head=git_head())
    ref = list(csv.DictReader(open(os.path.join(M2_DIR, "m_rcut_resigma.csv")))); v1 = mrr_table(M2_DIR); worst = 0.0
    for r, x in zip(ref, v1):
        for c in MRR_COLS[1:]:
            u, w = float(r[c]), float(x[c]); worst = max(worst, abs(u - w) / max(abs(u), 1e-12))
    out["mrr_gate_rel"] = worst; log(f"[C.1 mrr] recalcul de results/M2/m_rcut_resigma.csv depuis les npz de M2 : écart relatif max {worst:.2e} (seuil 1e-9)")
    if not (len(ref) == len(v1) and worst <= 1e-9):
        save_json(os.path.join(ensure("c"), "c1post_results.json"), out); log("[C.1 mrr] refus : définitions non reproduites, rien n'est écrit"); sys.exit(3)
    v2 = mrr_table(RES)
    with open(os.path.join(RES, "m_rcut_resigma.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MRR_COLS); w.writeheader(); w.writerows(v2)
    for x in v2:
        log(f"[C.1 mrr] R_cut {x['R_cut']} : méd. Γ {x['med_Gamma_meV']:.2f} meV, méd. Re Σ {x['med_ReSigma_meV']:.2f} ; écart médian vs R_cut 4 : Γ {100*x['rel_med_dGamma']:.2f} %, "
            f"Re Σ {100*x['rel_med_dReSigma']:.2f} % ; E_res {x['E_res_minus_E_D_eV']:+.3f}")
    out["mrr"] = v2
    p = os.path.join(RES, "M_tests_summary.csv"); rows = list(csv.DictReader(open(p))); m2 = list(csv.DictReader(open(os.path.join(M2_DIR, "M_tests_summary.csv"))))
    gold = [r for r in m2 if r["test"].startswith("test d'or")]; a2 = [r for r in m2 if r["test"].startswith("porte A.2")]
    assert len(gold) == 1 and len(a2) == 2, (len(gold), len(a2))
    rows = [r for r in rows if not r["test"].startswith("porte A.2")]
    rows = [gold[0] if r["test"].startswith("test d'or") else r for r in rows] + a2
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["test", "quantité", "valeur", "seuil", "verdict", "source"]); w.writeheader(); w.writerows(rows)
    out["tests_rows"] = dict(n=len(rows), golden=gold[0], porte_A2=a2)
    log(f"[C.1 tests] {os.path.relpath(p, PROJ)} : {len(rows)} lignes ; test d'or de M2 recopié (D12) ; 2 lignes « porte A.2 » de M2")
    save_json(os.path.join(ensure("c"), "c1post_results.json"), out)


# ----------------------------------------------------------------------------------------------- C.2 : tab:rcut_M à trois colonnes, carte de Kaasbjerg
def cmd_c2(a):
    """
    tab:rcut_M (9×9, 12×12, R_cut 0…6) : tel quel / plateau (i), étiquettes actuelles / plateau (i), Wigner-Seitz ; porte : dernière colonne = m_rcut_convergence.csv
    de C.1 (valeurs au format du csv, chaînes identiques), code 3 sinon. Carte de Kaasbjerg (D.1 de R9) et bloc D.2 : plateau (i) + Wigner-Seitz, tel quel en regard.
    """
    import csv
    cfg = load_production(verbose=True); d = ensure("c"); f = os.path.join(d, "c2_results.json"); out = dict(head=git_head(), rcut={}, gates={})
    a3 = load_json(R9_A3); a3p = load_json(R9_A3P)
    pc = os.path.join(results_dir(cfg), "m_rcut_convergence.csv")
    if not os.path.exists(pc):
        log(f"[C.2] {pc} absent (C.1 locality non terminé) : STOP sur C.2"); sys.exit(3)
    c1 = {}
    for r in csv.DictReader(open(pc)):
        c1[(r["size"], r["R_cut"])] = [r[k] for k in ("size", "R_cut", "n_sites", "fine_grid", "max_dM_over_maxM", "sv_mismatch_pi_blocks", "diag_max_dM_over_max",
                                                     "diag_abs_mismatch", "maxM_pi_eV", "maxMdiag_eV")]          # dernière occurrence si le csv a été complété
    for S in C_SIZES:
        inp = c_inputs(cfg, S); n = inp["n"]; C = alignment_C(cfg, S); M = inp["M"]; t0 = time.time()
        dz = lt.defect_mwr(M, inp["Ump"], inp["Udmp"], inp["kmp"], inp["MP"], n, C_N=0.0)
        rows_tq = rcut_rows(inp, dz["Mwr"], dz["Rn"]); del dz
        dp_ = lt.defect_mwr(M, inp["Ump"], inp["Udmp"], inp["kmp"], inp["MP"], n, C_N=C)
        rows_pl = rcut_rows(inp, dp_["Mwr"], dp_["Rn"])
        ws = ws_images(dp_["Rn"], np.zeros(3, int), inp["MP"], inp["A_uc"])
        rows_ws = rcut_rows(inp, dp_["Mwr"], dp_["Rn"], ws=ws); del dp_
        mism = [r["R_cut"] for r in rows_ws if c1.get((S, str(r["R_cut"]))) != rcut_csv_strings(r)]
        r9b = a3[f"rcut_{S}"]["rows"]["brut"]; r9p = a3p[f"rcut_{S}"]["rows"]["aligne_plateau"]
        reg_b = max(abs(x["max_dM_over_maxM"] - y["max_dM_over_maxM"]) / y["max_dM_over_maxM"] for x, y in zip(rows_tq, r9b))
        reg_p = max(abs(x["max_dM_over_maxM"] - y["max_dM_over_maxM"]) / y["max_dM_over_maxM"] for x, y in zip(rows_pl, r9p))
        out["rcut"][S] = dict(C_N_eV=C, tel_quel=rows_tq, plateau=rows_pl, plateau_ws=rows_ws, n_tie_labels=int((ws["n_tie"] > 1).sum()),
                              csv_mismatch_R_cut=mism, R9_rel_tel_quel=reg_b, R9_rel_plateau=reg_p, elapsed_s=time.time() - t0)
        out["gates"][f"tab:rcut_M {S} : colonne Wigner-Seitz = m_rcut_convergence.csv de C.1"] = (len(mism) == 0)
        log(f"[C.2 rcut] {S} : " + " ; ".join(f"R_cut {a_['R_cut']} {a_['max_dM_over_maxM']:.3e} / {b_['max_dM_over_maxM']:.3e} / {c_['max_dM_over_maxM']:.3e}"
                                              for a_, b_, c_ in zip(rows_tq, rows_pl, rows_ws))
            + f" ; csv de C.1 : {'identique' if not mism else 'ÉCART aux R_cut ' + str(mism)} ; R9 (tel quel, plateau) {reg_b:.1e}, {reg_p:.1e} ; {time.time()-t0:.0f} s")
        save_json(f, out)
        if S != "9x9":
            del inp, M
            continue
        # carte de Kaasbjerg (D.1 de R9, k du XML) : tel quel et plateau (i), étiquettes recentrées + Wigner-Seitz (convention de production, P-c4)
        D1 = d1_setup(inp); maps = {}
        for lab, Cv in (("tel quel", 0.0), ("plateau", C)):
            dm = lt.defect_mwr(M, inp["U"], inp["Ud"], inp["k"], inp["MP"], n, C_N=Cv)
            wsd = ws_images(dm["Rn"], np.zeros(3, int), inp["MP"], inp["A_uc"])
            maps[lab] = d1_map(D1, dm["Mwr"], dm["Rn"], ws=wsd); del dm
            log(f"[C.2 D.1] {lab} + Wigner-Seitz : " + " ; ".join(f"{k} {v['mean']:.3f} [{v['min']:.3f}, {v['max']:.3f}] ({v['n']})" for k, v in maps[lab]["stats"].items()) + " eV Å²")
        # D.2 : bloc 2×2 de la paire π/π* à (K, K) sur la grille dense (M2 : tot, L, NL) ; plateau (i) par la forme fermée de D6 (−C_N N² V V†, sur tot et L)
        dp9 = inp["dp"]; kd = inp["k"]; V9 = np.einsum("kbw,kwv->kbv", inp["Ud"], inp["U"])
        _, eps = qe_io.get_k_eigenvalues(dp9["uc"], False); eps = np.asarray(eps)
        if eps.shape[0] != len(kd):
            eps = eps.T
        eps = eps * HA2EV; iK = int(np.argmin(np.linalg.norm(np.mod(kd - K_RED + 0.5, 1.0) - 0.5, axis=1)))
        w = np.abs(V9[iK][:, 3]) ** 2 + np.abs(V9[iK][:, 4]) ** 2; pk = sorted(np.argsort(-w)[:2], key=lambda b: eps[iK, b])
        dKK = -C * n * n * (V9[iK] @ V9[iK].conj().T)[np.ix_(pk, pk)]                                       # D_N(0) = N² ; eV
        blocks = {}
        for part, path in (("tot", dp9["mfile"]), ("L", dp9["mfile"].replace("M_dense_", "M_L_dense_")), ("NL", dp9["mfile"].replace("M_dense_", "M_NL_dense_"))):
            matrix_io.check_manifest(path, require_normalization=matrix_io.M_NORM_V2)
            Bk = np.array(np.load(path, mmap_mode="r")[np.ix_(pk, [iK], pk, [iK])]).reshape(2, 2) * HA2EV
            for var, Bm in (("tel quel", Bk),) + ((("plateau", Bk + dKK),) if part in ("tot", "L") else ()):
                ev = np.linalg.eigvalsh(0.5 * (Bm + Bm.conj().T))
                blocks[f"{part} | {var}"] = dict(eigenvalues_eVA2=[float(x) * D1["A_cell"] for x in ev], half_trace_eVA2=float(0.5 * np.trace(Bm).real) * D1["A_cell"],
                                                 row_norms_eVA2=[float(np.linalg.norm(Bm[i])) * D1["A_cell"] for i in range(2)], frobenius_eVA2=float(np.linalg.norm(Bm)) * D1["A_cell"])
        out["kaasbjerg"] = dict(D1={lab: m["stats"] for lab, m in maps.items()}, D2=blocks, pairK=[int(x) for x in pk],
                                meta=dict(k_red=D1["k"].tolist(), delta_Ainv=D1["delta"], A_cell=D1["A_cell"], disk_radius_Ainv=0.05 * D1["b1"], C_N_eV=C,
                                          excluded_K_Kp=int(D1["exK"].sum() + D1["exKp"].sum())),
                                R9_D1=load_json(R9_D)["D1"], R9_D2=load_json(R9_D).get("D2", {}).get("blocks", {}))
        np.savez(os.path.join(d, "c2_kaasbjerg_maps.npz"), kc=D1["kc"], Bc=D1["Bc"], **{f"{lab.replace(' ', '_')}_{b}": maps[lab][b] for lab in maps for b in ("valence", "conduction")})
        c2_figure(D1, maps)
        del D1, maps, inp, M
    out["all_ok"] = bool(all(out["gates"].values())); save_json(f, out); c2_tables(out)
    if not out["all_ok"]:
        log("[C.2] ÉCHEC de la porte (colonne Wigner-Seitz contre le csv de C.1) : code 3"); sys.exit(3)


def c2_figure(D1, maps):
    plt, pal = fig_style(); Bc = D1["Bc"]; kc = D1["kc"]
    imgs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)]) @ Bc[:2, :2].T
    kf = np.array([kc[i] + imgs[np.argmin(np.linalg.norm(kc[i] + imgs, axis=1))] for i in range(len(kc))])
    corners = np.array([Bc[:2, :2] @ np.array(c) for c in [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (-2/3, -1/3), (-1/3, -2/3), (1/3, -1/3)]])
    corners = corners[np.argsort(np.arctan2(corners[:, 1], corners[:, 0]))]; Kc = Bc @ K_RED
    labs = list(maps); fig, axes = plt.subplots(len(labs), 2, figsize=(6.5, 3.0 * len(labs)), squeeze=False)
    vmax = max(np.nanmax(maps[v][b]) for v in labs for b in ("valence", "conduction")); vmin = min(np.nanmin(maps[v][b]) for v in labs for b in ("valence", "conduction"))
    tl = {"tel quel": "M2 tel quel, Wigner-Seitz", "plateau": "plateau (i), Wigner-Seitz"}
    for i, lab in enumerate(labs):
        for j, (band, bl) in enumerate((("valence", r"valence ($\pi$)"), ("conduction", r"conduction ($\pi^*$)"))):
            ax = axes[i][j]; z = maps[lab][band]; m = np.isfinite(z)
            sc = ax.scatter(kf[m, 0], kf[m, 1], c=z[m], s=0.6, cmap=pal.CMAP_SEQ, vmin=vmin, vmax=vmax, rasterized=True)
            ax.plot(np.r_[corners[:, 0], corners[0, 0]], np.r_[corners[:, 1], corners[0, 1]], color=pal.REF, lw=0.6)
            ax.plot(Kc[0], Kc[1], "x", color=pal.ORANGE, ms=4)
            ax.set_aspect("equal"); ax.set_xlabel(r"$k'_x$ (Å$^{-1}$)", fontsize=7); ax.set_ylabel(r"$k'_y$ (Å$^{-1}$)", fontsize=7)
            ax.set_title(f"({'abcd'[2*i+j]}) {bl}, {tl[lab]}", fontsize=8)
            fig.colorbar(sc, ax=ax, label=r"$A_\mathrm{cell}|M_{k'k}|$ (eV Å$^2$)")
    fig.tight_layout(); savefig(fig, "kaasbjerg_plateau_ws"); plt.close(fig)


def c2_tables(out):
    L = ["# R10 — C.2 : tab:rcut_M à trois colonnes et grandeur de Kaasbjerg (tables générées par r10_driver.py c2)", ""]
    L += ["| porte | résultat |", "|---|---|"] + [f"| {k} | {'OK' if v else 'ÉCHEC'} |" for k, v in out["gates"].items()] + [""]
    for S, r in out["rcut"].items():
        L += [f"## tab:rcut_M, {S} (max|ΔM|/max|M| de la paire π/π*, grille fine 60² ; C_N {fr(r['C_N_eV']*1e3, -4)} meV ; {r['n_tie_labels']} étiquettes à égalité)", "",
              "| R_cut | sites | tel quel | plateau (i), étiquettes actuelles | plateau (i), Wigner-Seitz | SV (WS) | diagonale (WS) |", "|---|---|---|---|---|---|---|"]
        for x, y, z in zip(r["tel_quel"], r["plateau"], r["plateau_ws"]):
            L.append(f"| {x['R_cut']} | {x['n_sites']} | {x['max_dM_over_maxM']:.4e} | {y['max_dM_over_maxM']:.4e} | {z['max_dM_over_maxM']:.4e} | {z['sv_mismatch_pi_blocks']:.4e} | {z['diag_max_dM_over_max']:.4e} |")
        L += ["", f"max|M_π| (eV) : tel quel {r['tel_quel'][0]['maxM_pi_eV']:.4f}, plateau {r['plateau'][0]['maxM_pi_eV']:.4f}, plateau WS {r['plateau_ws'][0]['maxM_pi_eV']:.4f} ; "
              f"écart relatif max à R9 : tel quel {r['R9_rel_tel_quel']:.1e}, plateau (étiquettes actuelles) {r['R9_rel_plateau']:.1e} ; "
              f"colonne WS contre le csv de C.1 : {'identique' if not r['csv_mismatch_R_cut'] else 'écart ' + str(r['csv_mismatch_R_cut'])}.", ""]
    if "kaasbjerg" in out:
        kb = out["kaasbjerg"]; m = kb["meta"]; keys = ["valence | K", "conduction | K", "valence | K'", "conduction | K'"]
        L += [f"## D.1 : Ṽ = A_cell·|M_k'k|, k = K + δx̂ (δ = {m['delta_Ainv']:.4f} Å⁻¹), carte 240², disques de rayon {m['disk_radius_Ainv']:.4f} Å⁻¹ ({m['excluded_K_Kp']} points exclus) ; eV Å²", "",
              "| variante | valence K | conduction K | valence K′ | conduction K′ |", "|---|---|---|---|---|"]
        for lab, s in kb["D1"].items():
            L.append(f"| {lab} + Wigner-Seitz | " + " | ".join(f"{s[k]['mean']:.3f} [{s[k]['min']:.2f} ; {s[k]['max']:.2f}]" for k in keys) + " |")
        for lab, s in kb["R9_D1"].items():
            L.append(f"| R9 {lab} (étiquettes brutes) | " + " | ".join(f"{s[k]['mean']:.3f}" for k in keys) + " |")
        L += ["", "## D.2 : bloc 2 × 2 de la paire π/π* à (K, K) (eV Å²)", "", "| partie | variante | valeurs propres | ½ Tr | normes des lignes | ‖·‖_F |", "|---|---|---|---|---|---|"]
        for k, v in kb["D2"].items():
            p_, var = k.split(" | ")
            L.append(f"| {p_} | {var} | {v['eigenvalues_eVA2'][0]:.2f} ; {v['eigenvalues_eVA2'][1]:.2f} | {v['half_trace_eVA2']:.2f} | {v['row_norms_eVA2'][0]:.2f} ; {v['row_norms_eVA2'][1]:.2f} | {v['frobenius_eVA2']:.2f} |")
        L += ["", "Kaasbjerg (PRB 101, 045433, Fig. 3) : ~70 eV Å² près de K et K′ ; V₀ ≈ 27 eV ; super-cellule 11×11."]
    with open(os.path.join(WORK, "c", "C2_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")

# ----------------------------------------------------------------------------------------------- C.3 : contrôles des sorties touchées par l'audit
VED_SIZES = ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]
LOC_SIZES = ["5x5", "6x6", "7x7", "8x8", "9x9", "12x12"]


def ved_plane_R(S):
    """Distances au site dans le plan de la grille d'analyze_Ved.py : réduction axe par axe (avant (b)) et vraie image (P-c2, true_min_image_dist)."""
    scp = f"data/graphene/supercell/qe/defect_{S}_p.save"; scd = f"data/graphene/supercell/qe/defect_{S}_d.save"
    A, _ = qe_io.get_A_volume(scd); xp = np.mod(qe_io.get_x_red(scp), 1.0); xd = np.mod(qe_io.get_x_red(scd), 1.0)
    dmin = np.array([np.min(np.linalg.norm(np.mod(xd - p + 0.5, 1) - 0.5, axis=1)) for p in xp]); s_vac = xp[int(np.argmax(dmin))]
    nr = np.array(qe_io.get_ngfft(scd), int); i = np.arange(nr[0]); j = np.arange(nr[1])
    s1 = (i / nr[0] - s_vac[0] + 0.5) % 1 - 0.5; s2 = (j / nr[1] - s_vac[1] + 0.5) % 1 - 0.5; S1, S2 = np.meshgrid(s1, s2, indexing="ij")
    Ro = np.sqrt(((S1 * A[0, 0] + S2 * A[0, 1]) * BOHR) ** 2 + ((S1 * A[1, 0] + S2 * A[1, 1]) * BOHR) ** 2)
    I1, I2 = np.meshgrid(i / nr[0], j / nr[1], indexing="ij")
    Rn = al.true_min_image_dist(np.stack([I1.ravel(), I2.ravel(), np.full(I1.size, s_vac[2])], axis=1), s_vac, A * BOHR).reshape(I1.shape)
    return dict(A=A, nr=nr, Ro=Ro, Rn=Rn, a_sc=float(np.linalg.norm(A[:, 0]) * BOHR))


def cmd_c3(a):
    """
    C.3 : (1) P-c2, anneaux d'analyze_Ved (porte : appartenance identique et valeurs à 1e-13 eV sur les anneaux entièrement sous 0,433 a_sc ; au-delà : contre
    l'audit) ; (2) P-b2, abscisses de mwr_locality (porte : distances vraies et nombre de points déplacés = audit, poids hors site = M2, sur site − M2 = −C_N) ;
    (3) tab:tests_M par famille (porte : max|M| = M2, écarts de famille recalculés = csv) ; (4) portes D6 d'analyze_M (lecture). Code 3 si une porte échoue.
    """
    import csv
    cfg = load_production(verbose=True); RES = results_dir(cfg); d = ensure("c"); out = dict(head=git_head(), gates={})
    aud = load_json(R9_AUDIT)
    # (1) P-c2
    Vn = np.load(os.path.join(RES, "ved_analysis.npz")); Vo = np.load(os.path.join(M2_DIR, "ved_analysis.npz")); ring = {}
    for S in VED_SIZES:
        g = ved_plane_R(S); assert tuple(Vn[f"{S}_map"].shape) == tuple(g["nr"][:2]), (S, Vn[f"{S}_map"].shape, g["nr"])
        rmax = 0.5 * g["a_sc"]; edges = np.arange(0, rmax + 0.05, 0.05); lim = np.sqrt(3) / 4 * g["a_sc"]
        jmax = int(np.searchsorted(edges, lim, side="right")) - 1                                              # anneaux 0 … jmax−1 : bord supérieur <= lim
        top = edges[jmax]; sel = (g["Ro"] < top) | (g["Rn"] < top)
        bo = np.digitize(g["Ro"][sel], edges); bn = np.digitize(g["Rn"][sel], edges); n_move = int((bo != bn).sum())
        r = dict(lim_A=float(lim), n_full_rings=jmax, top_edge_A=float(top), n_points=int(sel.sum()), n_points_changing_ring=n_move,
                 max_abs_R_diff_A=float(np.abs(g["Rn"][sel] - g["Ro"][sel]).max()))
        for key, rck in (("rad", "rc"), ("rad_masked", "rc_masked")):
            rc = Vn[f"{S}_{rck}"]; assert np.array_equal(rc, Vo[f"{S}_{rck}"]); m = np.arange(len(rc)) < jmax
            a_, b_ = Vn[f"{S}_{key}"][m], Vo[f"{S}_{key}"][m]; fin = np.isfinite(a_) & np.isfinite(b_)
            dmax = float(np.abs(a_[fin] - b_[fin]).max()); nan_eq = bool(np.array_equal(np.isfinite(a_), np.isfinite(b_)))
            ia = int(np.argmax(np.where(fin, np.abs(a_ - b_), -1)))
            beyond = (~m) & np.isfinite(Vn[f"{S}_{key}"]) & np.isfinite(Vo[f"{S}_{key}"])
            dch = np.abs(Vn[f"{S}_{key}"][beyond] - Vo[f"{S}_{key}"][beyond])
            ab = aud["sizes"][S]["radial_analyze_Ved"]["sans masque" if key == "rad" else "masqué"]
            aud_bins = {round(b["r_A"], 3): b for b in ab["bins"]}
            d_aud = [abs(Vn[f"{S}_{key}"][k] * 1e3 - aud_bins[round(float(rc[k]), 3)]["true_meV"]) for k in range(len(rc)) if round(float(rc[k]), 3) in aud_bins]
            d_aud = np.array(d_aud, float)
            r[key] = dict(max_abs_diff_full_rings_eV=dmax, at_r_A=float(rc[m][ia]), nan_pattern_equal=nan_eq, n_rings_over_1e13=int((np.abs(a_[fin] - b_[fin]) > 1e-13).sum()),
                          n_rings_changed_beyond=int((dch > 0).sum()), max_change_beyond_meV=float(dch.max() * 1e3) if dch.size else 0.0,
                          r_at_max_change_beyond_A=float(rc[beyond][int(np.argmax(dch))]) if dch.size else None,
                          audit_n_bins_count_differs=ab["n_bins_count_differs"], audit_r_first_bin_differs_A=ab["r_first_bin_differs_A"],
                          audit_max_abs_diff_meV=ab["max_abs_diff_meV"], max_abs_vs_audit_true_meV=float(np.nanmax(d_aud)) if d_aud.size else None, n_bins_vs_audit=int(np.isfinite(d_aud).sum()))
            out["gates"][f"P-c2 {S} {key} : anneaux sous 0,433 a_sc à 1e-13 eV"] = bool(dmax <= 1e-13 and nan_eq)
        out["gates"][f"P-c2 {S} : appartenance aux anneaux identique"] = (n_move == 0)
        if S == "9x9":
            cnt = np.histogram(g["Rn"][g["Rn"] < rmax], edges)[0]; rcs = 0.5 * (edges[1:] + edges[:-1])
            r["counts_9x9"] = {f"{x:.3f}": int(cnt[int(np.argmin(np.abs(rcs - x)))]) for x in (9.625, 11.075)}
        ring[S] = r
        log(f"[C.3 P-c2] {S} : {jmax} anneaux sous 0,433 a_sc = {lim:.3f} Å, {r['n_points']} points, {n_move} changent d'anneau ; écart max profil {r['rad']['max_abs_diff_full_rings_eV']:.2e}, "
            f"masqué {r['rad_masked']['max_abs_diff_full_rings_eV']:.2e} eV ; au-delà : {r['rad']['n_rings_changed_beyond']} anneaux changés (audit {r['rad']['audit_n_bins_count_differs']}), "
            f"max {r['rad']['max_change_beyond_meV']:.2f} meV (audit {r['rad']['audit_max_abs_diff_meV']:.2f}) ; écart à l'audit « vraie » {r['rad']['max_abs_vs_audit_true_meV']:.1e} meV")
    out["P_c2"] = ring
    # (2) P-b2
    Ln = np.load(os.path.join(RES, "mwr_locality.npz")); Lo = np.load(os.path.join(M2_DIR, "mwr_locality.npz")); A_UNIT = np.array([[4.0354919061, -2.3298923383], [4.0354919061, 2.3298923383]]) / 4.6597846766
    Acols = np.eye(3); Acols[:2, :2] = A_UNIT.T; loc = {}
    for S in LOC_SIZES:
        for grid in ("dense", "coarse"):
            tag = f"{S}_{grid}"
            if f"{tag}_dist" not in Ln.files:
                continue
            D = cfg["dense"][S]["D"] if grid == "dense" else size_n(S); ax_ = np.arange(D) - D // 2
            Rl = np.array([(i, j, 0) for i in ax_ for j in ax_], int); wsl = ws_images(Rl, np.zeros(3, int), (D, D, 1), Acols)
            dax = np.linalg.norm(Rl[:, :2] @ A_UNIT, axis=1); n_moved = int((np.abs(wsl["dist"] - dax) > 1e-9).sum())
            ref = aud["R_grids"][S][grid]; C = alignment_C(cfg, S)
            wn, wo = np.sort(Ln[f"{tag}_w"]), np.sort(Lo[f"{tag}_w"]); d_on = float(Ln[f"{tag}_onsite_pzvac"]) - float(Lo[f"{tag}_onsite_pzvac"])
            x = dict(D=D, max_abscissa_a=float(Ln[f"{tag}_dist"].max()), audit_max_true_dist_a=ref["max_true_dist_a"], n_moved=n_moved, audit_n_dist_differs=ref["n_dist_differs"],
                     offsite_weights_max_abs_diff=float(np.abs(wn[:-1] - wo[:-1]).max()), onsite_pzvac_diff_eV=d_on, minus_C_N_eV=-C,
                     onsite_norm_new=float(Ln[f"{tag}_onsite_norm"]), onsite_norm_M2=float(Lo[f"{tag}_onsite_norm"]))
            loc[tag] = x
            out["gates"][f"P-b2 {tag}"] = bool(abs(x["max_abscissa_a"] - ref["max_true_dist_a"]) <= 1e-9 and n_moved == ref["n_dist_differs"]
                                             and x["offsite_weights_max_abs_diff"] == 0.0 and abs(d_on + C) <= 1e-12)
            log(f"[C.3 P-b2] {tag} : abscisse max {x['max_abscissa_a']:.4f} a (audit {ref['max_true_dist_a']:.4f}), points déplacés {n_moved} (audit {ref['n_dist_differs']}), "
                f"poids hors site {x['offsite_weights_max_abs_diff']:.1e}, p_z–p_z sur site − M2 {d_on*1e3:+.4f} meV (−C_N {-C*1e3:+.4f})")
    out["P_b2"] = loc
    # (3) tab:tests_M par famille
    rows = list(csv.DictReader(open(os.path.join(RES, "M_tests_summary.csv")))); An = np.load(os.path.join(RES, "M_analysis.npz"), allow_pickle=True)
    Ao = np.load(os.path.join(M2_DIR, "M_analysis.npz"), allow_pickle=True); fam = {}
    sc_eq = all(float(An[f"scale_{S}_dense"]) == float(Ao[f"scale_{S}_dense"]) and float(An[f"scale_{S}_coarse"]) == float(Ao[f"scale_{S}_coarse"]) for S in LOC_SIZES)
    for famk, flab in (("3m", "3m"), ("non3m", "non-3m")):
        Sf = cfg["families"][famk]; v = np.array([float(An[f"scale_{S}_dense"]) for S in Sf]); sp = float((v.max() - v.min()) / v.mean())
        row = [r for r in rows if r["test"].startswith(f"convention intensive (cellule unitaire), famille {flab} ")]
        fam[flab] = dict(sizes=Sf, spread=sp, csv=row[0] if row else None)
        out["gates"][f"tab:tests_M famille {flab}"] = bool(sc_eq and len(row) == 1 and row[0]["valeur"] == f"{sp:.1e}")
    fam["scale_eq_M2"] = sc_eq; fam["six_sizes_line_absent"] = not any(r["test"] == "convention intensive (cellule unitaire)" for r in rows)
    out["gates"]["tab:tests_M : ancienne ligne six tailles retirée"] = fam["six_sizes_line_absent"]
    out["tests_M"] = fam
    out["D6_gates_analyze_M"] = {k: float(An[k]) for k in An.files if k.startswith("align_")}
    log(f"[C.3 tests_M] 3m {fam['3m']['spread']:.3e}, non-3m {fam['non-3m']['spread']:.3e} ; max|M| = M2 : {sc_eq} ; portes D6 d'analyze_M : {out['D6_gates_analyze_M']}")
    out["all_ok"] = bool(all(out["gates"].values())); save_json(os.path.join(d, "c3_results.json"), out); c3_tables(out)
    log("[C.3] portes : " + " ; ".join(f"{k} {'OK' if v else 'ÉCHEC'}" for k, v in out["gates"].items()))
    if not out["all_ok"]:
        log("[C.3] au moins une porte refusée : parties concernées en STOP (code 3 en fin de C.3)"); sys.exit(3)


def c3_tables(out):
    L = ["# R10 — C.3 : contrôles des sorties touchées par l'audit (tables générées par r10_driver.py c3)", "",
         "| porte | résultat |", "|---|---|"] + [f"| {k} | {'OK' if v else 'REFUSÉE'} |" for k, v in out["gates"].items()] + [""]
    L += ["## P-c2 : anneaux d'analyze_Ved (pas 0,05 Å)", "",
          "| N | 0,433 a_sc (Å) | anneaux entiers dessous | points ; changent d'anneau | écart max profil (eV) à r | écart max masqué (eV) à r | anneaux changés au-delà (audit) | max au-delà (meV) (audit) | écart à l'audit « vraie » (meV) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for S, r in out["P_c2"].items():
        a_, b_ = r["rad"], r["rad_masked"]
        L.append(f"| {S} | {fr(r['lim_A'], 3)} | {r['n_full_rings']} | {r['n_points']} ; {r['n_points_changing_ring']} | {a_['max_abs_diff_full_rings_eV']:.2e} à {fr(a_['at_r_A'], 3)} | "
                 f"{b_['max_abs_diff_full_rings_eV']:.2e} à {fr(b_['at_r_A'], 3)} | {a_['n_rings_changed_beyond']} ({a_['audit_n_bins_count_differs']}) | {fr(a_['max_change_beyond_meV'], 2)} ({fr(a_['audit_max_abs_diff_meV'], 2)}) | "
                 f"{a_['max_abs_vs_audit_true_meV']:.1e} |")
    if "counts_9x9" in out["P_c2"].get("9x9", {}):
        L += ["", "9×9, effectifs en vraie image : " + " ; ".join(f"anneau {k} Å : {v} points" for k, v in out["P_c2"]["9x9"]["counts_9x9"].items()) + " (audit : 534, 518)."]
    L += ["", "## P-b2 : abscisses de mwr_locality (unités de a)", "", "| taille | D | abscisse max (audit) | points déplacés (audit) | poids hors site − M2 | p_z–p_z sur site − M2 (meV) ; −C_N |",
          "|---|---|---|---|---|---|"]
    for t, x in out["P_b2"].items():
        L.append(f"| {t} | {x['D']} | {x['max_abscissa_a']:.4f} ({x['audit_max_true_dist_a']:.4f}) | {x['n_moved']} ({x['audit_n_dist_differs']}) | {x['offsite_weights_max_abs_diff']:.1e} | "
                 f"{fr(x['onsite_pzvac_diff_eV']*1e3, -4)} ; {fr(x['minus_C_N_eV']*1e3, -4)} |")
    f_ = out["tests_M"]
    L += ["", "## tab:tests_M par famille", ""] + [f"- famille {k} ({', '.join(v['sizes'])}) : (max − min)/moyenne = {v['spread']:.3e} ; csv « {v['csv']['valeur'] if v['csv'] else '—'} » ({v['csv']['verdict'] if v['csv'] else '—'})"
                                                  for k, v in f_.items() if isinstance(v, dict)]
    L += [f"- max|M| (denses et grossiers) = M2 : {f_['scale_eq_M2']} ; ligne « six tailles » absente : {f_['six_sizes_line_absent']}",
          "", "Portes D6 d'`analyze_M.py` (C.1) : " + " ; ".join(f"{k} {v:.1e}" for k, v in out["D6_gates_analyze_M"].items())]
    with open(os.path.join(WORK, "c", "C3_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------------- C.4 : planches avant/après
def cmd_c4(a):
    """Planches avant (figures/, lecture seule) / après (R10_plateau/fig/) pour chaque figure présente des deux côtés ; carte M brute / alignée (D6)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    cfg = load_production(verbose=False); dd = ensure(os.path.join("fig", "avant_apres")); out = dict(head=git_head(), planches={})
    names = sorted(n[:-4] for n in os.listdir(os.path.join(WORK, "fig")) if n.endswith(".png") and os.path.exists(os.path.join(FIGURES_DIR, n)))
    for name in names:
        A_ = mpimg.imread(os.path.join(FIGURES_DIR, name + ".png")); B_ = mpimg.imread(os.path.join(WORK, "fig", name + ".png"))
        same = A_.shape == B_.shape and float(np.abs(A_.astype(float) - B_.astype(float)).max()) == 0.0
        fig, axs = plt.subplots(1, 2, figsize=(13, 13 * max(A_.shape[0] / A_.shape[1], B_.shape[0] / B_.shape[1]) / 2 + 0.6))
        for ax, im, t in ((axs[0], A_, "avant : figures/ (v2 tel quel, results/M2)"), (axs[1], B_, "après : R10 (plateau (i), results/M2_plateau)")):
            ax.imshow(im); ax.set_axis_off(); ax.set_title(t, fontsize=10)
        fig.suptitle(f"{name}" + (" — identique au pixel" if same else ""), fontsize=11); fig.tight_layout()
        fig.savefig(os.path.join(dd, f"{name}.png"), dpi=110); plt.close(fig)
        out["planches"][name] = dict(identique_au_pixel=bool(same), shape_avant=list(A_.shape), shape_apres=list(B_.shape))
    log(f"[C.4] {len(names)} planches dans fig/avant_apres/ ; identiques au pixel : {[n for n, v in out['planches'].items() if v['identique_au_pixel']]}")
    # carte de M à k = K (fig_M_map) : brut / aligné (D6), même échelle par ligne
    plt_, pal = fig_style(); Z = np.load(os.path.join(results_dir(cfg), "M_analysis.npz"), allow_pickle=True)
    fig, axs = plt_.subplots(3, 2, figsize=(6.5, 8.4))
    for i, (key, lab) in enumerate((("map_Vpi", r"$\pi$ : $A_\mathrm{cell}|M|$"), ("map_Vpistar", r"$\pi^*$ : $A_\mathrm{cell}|M|$"), ("map_Lpar", r"$\pi$ : $M^L$ le long de $M$"))):
        lo = min(Z[key].min(), Z[key + "_aligned"].min()); hi = max(Z[key].max(), Z[key + "_aligned"].max())
        for j, (suf, t) in enumerate((("", "brut"), ("_aligned", "aligné (D6)"))):
            ax = axs[i][j]; sc = ax.scatter(Z["map_kx"], Z["map_ky"], c=Z[key + suf], cmap=pal.CMAP_SEQ, vmin=lo, vmax=hi, s=10, marker="h", linewidths=0)
            ax.plot(*Z["map_K"], "o", mfc="none", mec=pal.ORANGE, mew=1.0, ms=6); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f"({'abcdef'[2*i+j]}) {lab}, {t}", fontsize=8); fig.colorbar(sc, ax=ax, label=r"eV Å$^2$")
    fig.tight_layout(); savefig(fig, "M_map_brut_aligne"); plt_.close(fig)
    save_json(os.path.join(ensure("c"), "c4_results.json"), out)


# ----------------------------------------------------------------------------------------------- C.5 : table de correspondance, README et MD5 de results_dir
def cmd_c5(a):
    """
    Table de correspondance v2 tel quel (results/M2) → v2 plateau (results/M2_plateau), format de table_v1_v2.md (logique de r6_compare_v1_v2.py, recopiée ;
    clés absentes : « — »), + lignes R10 (C_N, D6, C14, tab:rcut_M à trois colonnes, Kaasbjerg) et ch. 5 ; README.md et MD5SUMS de results_dir.
    """
    import csv
    cfg = load_production(verbose=False); V1 = M2_DIR; V2 = results_dir(cfg); RC, G, REF = cfg["R_cut"], cfg["grid"], cfg["reference_size"]
    rows = []
    def add(q, x, y, src, where): rows.append((q, x, y, src, where))
    def fmt(x, nd=3):
        if x is None: return "—"
        if isinstance(x, str): return x
        try: return f"{float(x):.{nd}f}" if abs(float(x)) < 1e4 else f"{float(x):.4e}"
        except Exception: return str(x)
    def csvrows(dd, name):
        p = os.path.join(dd, name); return list(csv.DictReader(open(p))) if os.path.exists(p) else None
    def npz(dd, name):
        p = os.path.join(dd, name); return np.load(p, allow_pickle=True) if os.path.exists(p) else None
    def fnpz(Z, k):
        return float(Z[k]) if (Z is not None and k in Z.files) else None
    # C_N
    for S in ALL13:
        add(f"C_N plateau (i), {S} (meV)", "0 (tel quel)", fmt(alignment_C(cfg, S) * 1e3, 4), "config/production.json (alignment)", "D14, R10 A.1")
    # niveau 1
    def lvl1(dd, S, rc, g, e):
        r = csvrows(dd, "level1_summary.csv")
        for x in (r or []):
            if x["size"] == S and int(x["R_cut"]) == rc and int(x["grid"]) == g and abs(float(x["eta_eV"]) - e) < 1e-6:
                return float(x["median_Gamma_Ncells_meV"]), float(x["argmax_E_minus_ED_eV"])
        return None, None
    for S in LOC_SIZES:
        (m1, e1), (m2, e2) = lvl1(V1, S, RC, G, cfg["eta_eV"]), lvl1(V2, S, RC, G, cfg["eta_eV"])
        add(f"médiane |Γ|·N_cells, {S}, R_cut {RC}, {G}², η {cfg['eta_eV']} (meV)", fmt(m1, 2), fmt(m2, 2), "level1_summary.csv", "fig_level2, fig_spectral, §2 NOTES_TGAMMA")
        add(f"E_res − E_D (argmax états ±1,5 eV), {S} (eV)", fmt(e1), fmt(e2), "level1_summary.csv", "fig_plateau, §2")
    for (rc, g, e) in ((RC, G, 0.01), (RC, 120, cfg["eta_eV"]), (2, G, cfg["eta_eV"]), (4, G, cfg["eta_eV"]), (0, G, cfg["eta_eV"])):
        (m1, _), (m2, _) = lvl1(V1, REF, rc, g, e), lvl1(V2, REF, rc, g, e)
        add(f"médiane |Γ|·N_cells, {REF}, R_cut {rc}, {g}², η {e} (meV)", fmt(m1, 2), fmt(m2, 2), "level1_summary.csv", "fig_plateau, fig_rcut, C10/C11")
    # familles
    f1, f2 = csvrows(V1, "level2_families.csv"), csvrows(V2, "level2_families.csv")
    for x in (f1 or []):
        y = next((z for z in (f2 or []) if z["size"] == x["size"]), None)
        add(f"p_z–p_z sur le site de la lacune, {x['size']} dense (eV)", fmt(x["onsite_pz_vac_eV"]), fmt(y["onsite_pz_vac_eV"]) if y else "—", "level2_families.csv (mwr_locality.npz)", "fig_locality, tab familles")
        add(f"Re M^L / Re M^NL à K, paire π, {x['size']} (eV)", f"{x['ReML_K_eV']} / {x['ReMNL_K_eV']}", (f"{y['ReML_K_eV']} / {y['ReMNL_K_eV']}" if y else "—"), "level2_families.csv (M_analysis.npz)", "fig_M_map, §4.1.5")
        add(f"Re M^L à K aligné (D6), {x['size']} (eV)", "—", (y.get("ReML_K_aligned_eV", "—") if y else "—"), "level2_families.csv (M_analysis.npz)", "tab familles (D6)")
    # tab:L_NL
    l1, l2 = csvrows(V1, "lnl_frobenius.csv"), csvrows(V2, "lnl_frobenius.csv")
    for x in (l1 or []):
        y = next((z for z in (l2 or []) if z["size"] == x["size"]), None)
        add(f"⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩, {x['size']} dense : brut ; aligné (D6)", fmt(x["ratio_full"]), (f"{fmt(y['ratio_full'])} ; {fmt(y.get('ratio_full_aligned'))}" if y else "—"), "lnl_frobenius.csv", "tab:L_NL")
        add(f"⟨‖M^L‖_F⟩, {x['size']} (eV) : brut ; aligné", fmt(x["mean_fL_eV"], 4), (f"{fmt(y['mean_fL_eV'], 4)} ; {fmt(y.get('mean_fL_aligned_eV'), 4)}" if y else "—"), "lnl_frobenius.csv", "tab:L_NL")
    # carte de M (D6)
    A1, A2 = npz(V1, "M_analysis.npz"), npz(V2, "M_analysis.npz")
    if A2 is not None and "map_Vpi_aligned" in A2.files:
        iK = int(np.argmin(np.linalg.norm(np.stack([A2["map_kx"], A2["map_ky"]], 1) - A2["map_K"][None], axis=1)))
        for key, lab in (("map_Vpi", "Ṽ_π à K"), ("map_Vpistar", "Ṽ_π* à K"), ("map_Lpar", "M^L_∥ à K")):
            add(f"carte fig_M_map, {lab} (eV Å²) : brut ; aligné", fmt(float(A1[key][iK]), 2) if A1 is not None else "—", f"{fmt(float(A2[key][iK]), 2)} ; {fmt(float(A2[key + '_aligned'][iK]), 2)}", "M_analysis.npz", "fig_M_map (D6)")
        for k in sorted(x for x in A2.files if x.startswith("align_")):
            add(f"porte D6 {k}", "—", f"{float(A2[k]):.1e}", "M_analysis.npz", "D6")
    # localité
    L1, L2 = npz(V1, "mwr_locality.npz"), npz(V2, "mwr_locality.npz")
    for S in ("9x9", "12x12"):
        for key, lab in (("onsite_norm", "‖M_W(0,0)‖ (eV)"), ("onsite_pzvac", "p_z–p_z lacune (eV)")):
            add(f"localité M_W {S} dense : {lab}", fmt(fnpz(L1, f"{S}_dense_{key}")), fmt(fnpz(L2, f"{S}_dense_{key}")), "mwr_locality.npz", "fig_locality, §2")
        add(f"localité {S} dense : abscisse max (a)", fmt(float(L1[f"{S}_dense_dist"].max())) if L1 is not None else "—", fmt(float(L2[f"{S}_dense_dist"].max())) if L2 is not None else "—", "mwr_locality.npz", "fig_locality (P-b2)")
    # tab:rcut_M
    r1, r2 = csvrows(V1, "m_rcut_convergence.csv"), csvrows(V2, "m_rcut_convergence.csv"); c2 = load_json(os.path.join(WORK, "c", "c2_results.json")).get("rcut", {})
    for x in (r1 or []):
        y = next((z for z in (r2 or []) if z["size"] == x["size"] and z["R_cut"] == x["R_cut"]), None)
        extra = ""
        if x["size"] in c2:
            p_ = c2[x["size"]]["plateau"][int(x["R_cut"])]["max_dM_over_maxM"]; extra = f" (plateau, étiquettes actuelles : {p_:.4e})"
        add(f"tab:rcut_M {x['size']} R_cut {x['R_cut']} : max|ΔM|/max|M| (60²)", x["max_dM_over_maxM"], (y["max_dM_over_maxM"] if y else "—") + extra, "m_rcut_convergence.csv ; c/c2_results.json", "tab:rcut_M")
    # m_rcut_resigma, N_k^int
    s1, s2 = csvrows(V1, "m_rcut_resigma.csv"), csvrows(V2, "m_rcut_resigma.csv")
    for x in (s1 or []):
        y = next((z for z in (s2 or []) if z["R_cut"] == x["R_cut"]), None)
        add(f"Re Σ / Γ, 9x9, R_cut {x['R_cut']} : méd. Γ ; méd. Re Σ (meV) ; E_res", f"{float(x['med_Gamma_meV']):.2f} ; {float(x['med_ReSigma_meV']):.2f} ; {float(x['E_res_minus_E_D_eV']):+.3f}",
            (f"{float(y['med_Gamma_meV']):.2f} ; {float(y['med_ReSigma_meV']):.2f} ; {float(y['E_res_minus_E_D_eV']):+.3f}" if y else "—"), "m_rcut_resigma.csv", "C10, C11")
    n1, n2 = csvrows(V1, "nkint_check_9x9.csv"), csvrows(V2, "nkint_check_9x9.csv")
    for x in (n1 or []):
        y = next((z for z in (n2 or []) if z["nk_int"] == x["nk_int"]), None)
        add(f"N_k^int {x['nk_int']} : médiane Γ / Re Σ (meV) / E_res / Γ(E_D)", f"{float(x['medG']):.2f} / {float(x['medRe']):.2f} / {float(x['E_res']):+.3f} / {float(x['G_ED']):.2f}",
            (f"{float(y['medG']):.2f} / {float(y['medRe']):.2f} / {float(y['E_res']):+.3f} / {float(y['G_ED']):.2f}" if y else "—"), "nkint_check_9x9.csv", "§6 NOTES_TGAMMA")
    # résonances 9×9, 6×6, 12×12
    for S in (REF, "6x6", "12x12"):
        R1, R2 = npz(V1, f"resonance_{S}.npz"), npz(V2, f"resonance_{S}.npz")
        def rst(R):
            if R is None: return {}
            eg = R["eg"]; ED = float(R["E_D"]); m = np.abs(eg - ED) <= 3.0; GT = R["Gamma_T"]; GB = R["Gamma_Born"]
            dd = dict(medGT=float(np.nanmedian(GT[m]) * 1e3), medGB=float(np.nanmedian(GB[m]) * 1e3), BoverT=float(np.nanmedian(GB[m] / GT[m])), ReT=float(R["ReTbar_at_ED"]), ImT=float(R["ImTbar_at_ED"]))
            for k in ("peak_GT", "peak_GB", "peak_ratio", "peak_drho", "peak_rho_dis", "peak_absTbar", "peak_ImTbar", "min_absReTbar", "E_res_states", "median_GT_states_meV"):
                if k in R.files: dd[k] = float(R[k])
            dd["zc"] = ", ".join(f"{v:+.3f}" for v in R["Tbar_zero_crossings"]) or "aucun"
            return dd
        a_, b_ = rst(R1), rst(R2)
        for k, lab in (("medGT", "médiane courbe Γ_T ±3 eV (meV)"), ("medGB", "médiane Γ_Born (meV)"), ("BoverT", "Born/T médian"), ("peak_GT", "pic Γ_T (eV)"), ("peak_ratio", "pic Γ_T/ρ₀"),
                       ("peak_drho", "pic δρ"), ("peak_absTbar", "pic |T̄|"), ("peak_ImTbar", "pic −Im T̄"), ("min_absReTbar", "min |Re T̄|"), ("ReT", "Re T̄(E_D) (eV)"), ("ImT", "Im T̄(E_D) (eV)"),
                       ("E_res_states", "E_res états ±1,5 eV"), ("median_GT_states_meV", "médiane |Γ_T| états (meV)"), ("zc", "zéros de Re T̄")):
            add(f"résonance {S} : {lab}", fmt(a_.get(k), 2 if k in ("medGT", "medGB", "median_GT_states_meV") else 3), fmt(b_.get(k), 2 if k in ("medGT", "medGB", "median_GT_states_meV") else 3),
                f"resonance_{S}.npz", "fig_spectral, fig_epw_vs_ed (ch. 5), §2" if S == REF else "R6 3.4 (D9)")
        C1_, C2_ = npz(V1, f"resonance_criteria_{S}.npz"), npz(V2, f"resonance_criteria_{S}.npz")
        for suf, bl in (("", "complet"), ("_pi", "bloc π"), ("_sigma", "bloc σ")):
            def cs(Cz):
                if Cz is None or f"logdet_rel{suf}" not in Cz.files: return "—"
                eg = Cz["eg"]; ED = float(Cz["E_D"]); ld = Cz[f"logdet_rel{suf}"]; ml = Cz[f"minlam{suf}"]; jd = int(np.argmin(ld)); jl = int(np.argmin(ml))
                return f"{np.exp(ld[jd]):.3e} ({eg[jd]-ED:+.3f}) ; |λ| {ml[jl]:.4f} ({eg[jl]-ED:+.3f}) ; Friedel {float(Cz[f'sumrule{suf}']):+.4f}"
            add(f"critère {S} ({bl}) : min|det|/max ; min|λ| ; ∫δρ", cs(C1_), cs(C2_), f"resonance_criteria_{S}.npz", "§2, C15, C16")
    # C14
    Rs1, Rs2 = npz(V1, f"resonance_{REF}_shiftL.npz"), npz(V2, f"resonance_{REF}_shiftL.npz")
    for Rz, lab in ((Rs1, "M2 (±25 meV, V_loc tel quel)"), (Rs2, "plateau (±9,05 meV, V_loc aligné)")):
        if Rz is None: continue
        for k in sorted(Rz.files):
            if k.startswith(("E_res_states", "peak_GT_shift", "median_GT_states")):
                add(f"C14 {lab} : {k}", fmt(float(Rz[k]), 3) if Rz is Rs1 else "—", fmt(float(Rz[k]), 3) if Rz is Rs2 else "—", f"resonance_{REF}_shiftL.npz", "C14 (D10)")
    # tab:tests_M
    t1, t2 = csvrows(V1, "M_tests_summary.csv"), csvrows(V2, "M_tests_summary.csv")
    for x in (t1 or []):
        y = next((z for z in (t2 or []) if z["test"] == x["test"]), None)
        add(f"tab:tests_M : {x['test']}", f"{x['valeur']} ({x['verdict']})", (f"{y['valeur']} ({y['verdict']})" if y else "— (ligne retirée)"), "M_tests_summary.csv", "tab:tests_M")
    for z in (t2 or []):
        if not any(x["test"] == z["test"] for x in (t1 or [])):
            add(f"tab:tests_M : {z['test']}", "—", f"{z['valeur']} ({z['verdict']})", "M_tests_summary.csv", "tab:tests_M (nouvelle ligne)")
    # Kaasbjerg (C.2)
    kb = load_json(os.path.join(WORK, "c", "c2_results.json")).get("kaasbjerg")
    if kb:
        for k in ("valence | K", "conduction | K", "valence | K'", "conduction | K'"):
            add(f"D.1 Kaasbjerg 9×9, {k} (eV Å²) : R9 brut → tel quel WS ; plateau WS", fmt(kb["R9_D1"]["brut"][k]["mean"], 3), f"{fmt(kb['D1']['tel quel'][k]['mean'], 3)} ; {fmt(kb['D1']['plateau'][k]['mean'], 3)}", "c/c2_results.json", "carte de Kaasbjerg")
    # chapitre 5
    E1, E2 = npz(V1, "ed_vs_ep_24k24q_mv0.02.npz"), npz(V2, "ed_vs_ep_24k24q_mv0.02.npz")
    if E1 is not None and E2 is not None:
        w12 = lambda E, k: float(np.median(E[k][np.abs(E["x"]) <= 1.2])) * 1e3
        add("Γ^ed/Γ^ep (c = 1 %, 300 K, ±3 eV) : médiane", fmt(float(E1["median"])), fmt(float(E2["median"])), "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed, NOTES_EPW (ch. 5)")
        add("Γ^ed/Γ^ep : min (à ε − E_D)", f"{float(E1['min']):.3f} ({float(E1['x_min']):+.3f})", f"{float(E2['min']):.3f} ({float(E2['x_min']):+.3f})", "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
        add("Γ^ed/Γ^ep : max (à ε − E_D)", f"{float(E1['max']):.3f} ({float(E1['x_max']):+.3f})", f"{float(E2['max']):.3f} ({float(E2['x_max']):+.3f})", "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
        add("Γ^ed/Γ^ep : croisements (ratio = 1)", ", ".join(f"{v:+.3f}" for v in E1["crossings"]), ", ".join(f"{v:+.3f}" for v in E2["crossings"]), "ed_vs_ep_24k24q_mv0.02.npz", "NOTES_EPW (ch. 5)")
        add("médiane Γ^ed, c = 1 % (meV) : ±3 eV / ±1,2 eV", f"{float(E1['median_ed'])*1e3:.3f} / {w12(E1, 'Gamma_ed'):.3f}", f"{float(E2['median_ed'])*1e3:.3f} / {w12(E2, 'Gamma_ed'):.3f}", "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed (ch. 5)")
        add("médiane Γ^ep, 300 K (meV) : ±3 eV / ±1,2 eV (indépendant de M)", f"{float(E1['median_ep'])*1e3:.3f} / {w12(E1, 'Gamma_ep'):.3f}", f"{float(E2['median_ep'])*1e3:.3f} / {w12(E2, 'Gamma_ep'):.3f}", "ed_vs_ep_24k24q_mv0.02.npz", "fig_epw_vs_ed (ch. 5)")
    esc = lambda t: str(t).replace("|", "\\|")
    L = [f"# R10 C.5 — table de correspondance v2 tel quel (`{os.path.relpath(V1, PROJ)}`) → v2 plateau (i) (`{os.path.relpath(V2, PROJ)}`)\n",
         "| grandeur | v2 tel quel (M2) | v2 plateau (M2_plateau) | fichier source | figure / tableau / note |", "|---|---|---|---|---|"]
    L += [f"| {esc(q)} | {esc(x)} | {esc(y)} | `{s}` | {esc(w)} |" for q, x, y, s, w in rows]
    open(os.path.join(ensure("c"), "table_v2_plateau.md"), "w").write("\n".join(L) + "\n"); log(f"[C.5] {len(rows)} lignes -> c/table_v2_plateau.md")
    # README et MD5SUMS de results_dir (produits seulement ; journaux exclus)
    files = sorted(fn for fn in os.listdir(V2) if os.path.isfile(os.path.join(V2, fn)) and not fn.startswith(("README", "MD5SUMS")))
    stamp = time.strftime("%Y-%m-%d")
    with open(os.path.join(V2, f"MD5SUMS_{stamp}.txt"), "w") as fh:
        for fn in files:
            fh.write(f"{md5(os.path.join(V2, fn))}  {fn}\n")
    txt = [f"# {os.path.relpath(V2, PROJ)} — produits de production du chapitre 4, base plateau (i) (R10, {stamp})", "",
           "- Base : ΔV aligné par la moyenne du plateau (i) (config/production.json, bloc `alignment` ; C_N des 13 tailles, source "
           "`article/R10_plateau/a/a1_results.json`) ; M_W(R, R) − C_N·𝕀₅ sur les mailles de la boîte N×N (approximation (i), `local_tmatrix.defect_mwr`) ; "
           "interpolation hors grille par les images de Wigner-Seitz (`ws_images`, `ws_phase`) ; corrections de l'audit de l'image minimale (P-b2, P-c2, P-c4).",
           f"- Matrices : `{cfg['matrices_dir']}` (M2 brutes, lecture seule) ; ce répertoire ne contient que des produits (npz, csv) et les journaux (`logs/`, non versionnés).",
           f"- Rejeu : campagne R10 (`graphene/qe/defects/R10_plateau/`, copie `article/R10_plateau/`), lanceur `submit_r10.sh` (tâches c0, c1, c1f, c2, c3) ; HEAD {git_head()}.",
           "- Référence v2 tel quel : `results/M2` (gelé) ; table de correspondance : `article/R10_plateau/c/table_v2_plateau.md`.",
           "- `ks_reconstruction.npz` copié de `results/M2` (indépendant de M ; md5 identique), `ved_analysis.npz` refait (P-c2).",
           f"- md5 : `MD5SUMS_{stamp}.txt` ({len(files)} fichiers)."]
    open(os.path.join(V2, "README.md"), "w").write("\n".join(txt) + "\n"); log(f"[C.5] README.md et MD5SUMS_{stamp}.txt ({len(files)} fichiers) dans {os.path.relpath(V2, PROJ)}")


# ----------------------------------------------------------------------------------------------- C.6 : anomalie ⟨ΔV⟩_3D d'A.1
def scf_state(path):
    """État d'un scf.out de pw.x : convergence, nombre d'itérations, dernière « estimated scf accuracy » (Ry), JOB DONE ; paramètres lus dans scf.in."""
    txt = open(path, errors="replace").read(); import re
    m = re.search(r"convergence has been achieved in\s+(\d+)\s+iterations", txt); acc = re.findall(r"estimated scf accuracy\s+<\s+([0-9.Ee+-]+)\s+Ry", txt)
    it = len(re.findall(r"^\s*iteration #\s*\d+", txt, flags=re.M)); en = re.findall(r"^!\s+total energy\s+=\s+([0-9.Ee+-]+)\s+Ry", txt, flags=re.M)
    fe = re.findall(r"the Fermi energy is\s+([0-9.Ee+-]+)\s+ev", txt)
    inp = open(path.replace("scf.out", "scf.in"), errors="replace").read()
    par = {k: (re.search(rf"^\s*{k}\s*=\s*([^\s!,]+)", inp, flags=re.M | re.I).group(1) if re.search(rf"^\s*{k}\s*=\s*([^\s!,]+)", inp, flags=re.M | re.I) else None)
           for k in ("conv_thr", "mixing_beta", "electron_maxstep", "ecutwfc", "degauss", "occupations", "nspin", "nat", "assume_isolated", "tot_charge")}
    kp = re.search(r"K_POINTS\s*[({]?\s*(\w+)", inp)
    return dict(converged=bool(m), n_iter_reported=int(m.group(1)) if m else None, n_iterations_listed=it, last_scf_accuracy_Ry=float(acc[-1]) if acc else None,
                job_done="JOB DONE" in txt, not_converged="convergence NOT achieved" in txt, total_energy_Ry=float(en[-1]) if en else None,
                fermi_eV=float(fe[-1]) if fe else None, params=par, k_points=kp.group(1) if kp else None)


def cmd_c6(a):
    """
    C.6 : ⟨ΔV⟩_3D d'A.1 (N²⟨ΔV⟩ = 715–750 meV pour 5…12, écart pour 15…27) : état des SCF (d et p) de 15…27 (5…12 en regard) et profil ⟨ΔV⟩(z) moyenné dans
    le plan pour les 13 tailles (feuillet contre vide). Lecture seule ; chiffres bruts.
    """
    d = ensure("c"); out = dict(head=git_head(), sizes={}); a1 = load_json(a1_path()).get("sizes", {})
    for S in ALL13:
        t0 = time.time(); n = size_n(S); r = dict()
        for t, lab in (("defective", "d"), ("pristine", "p")):
            r[f"scf_{lab}"] = scf_state(os.path.join(SC_DIR, S, t, "scf.out"))
            pin = open(os.path.join(SC_DIR, S, t, "pp.in"), errors="replace").read(); import re
            r[f"pp_{lab}"] = {k: (re.search(rf"^\s*{k}\s*=\s*([^\s!,]+)", pin, flags=re.M | re.I).group(1) if re.search(rf"^\s*{k}\s*=\s*([^\s!,]+)", pin, flags=re.M | re.I) else None)
                              for k in ("plot_num", "filplot", "output_format", "iflag")}
        p = sc_paths(S); Vd = load_pot_eV(p["vd"]); Vp = load_pot_eV(p["vp"]); dV = Vd - Vp; del Vd, Vp
        A, _ = qe_io.get_A_volume(p["d"]); c = float(A[2, 2] * BOHR); nz = dV.shape[2]
        prof = dV.mean(axis=(0, 1)); del dV                                                                     # eV, par plan z
        x_d = qe_io.get_x_red(p["d"]); z_sheet = float(np.mod(np.median(x_d[:, 2]), 1.0))
        z = ((np.arange(nz) / nz - z_sheet + 0.5) % 1.0 - 0.5) * c                                             # Å, relatif au feuillet
        mean3d = float(prof.mean()); contrib = {}
        for zc in (1.0, 2.0, 3.0, 5.0):
            m = np.abs(z) <= zc; contrib[f"|z|<={zc:g}A"] = dict(sum_over_nz_meV=float(prof[m].sum() / nz * 1e3), n_planes=int(m.sum()))
        vac = np.abs(z) > 5.0; edge = int(np.argmax(np.abs(z)))
        r.update(N=n, c_A=c, nz=nz, z_sheet_red=z_sheet, mean3d_meV=mean3d * 1e3, N2_mean3d_meV=n * n * mean3d * 1e3, A1_mean3d_meV=a1.get(S, {}).get("mean3d_meV"),
                 contributions=contrib, vacuum_mean_meV=float(prof[vac].mean() * 1e3), vacuum_std_meV=float(prof[vac].std() * 1e3),
                 vacuum_min_meV=float(prof[vac].min() * 1e3), vacuum_max_meV=float(prof[vac].max() * 1e3), at_cell_edge_meV=float(prof[edge] * 1e3),
                 sheet_plane_meV=float(prof[int(np.argmin(np.abs(z)))] * 1e3), elapsed_s=time.time() - t0)
        np.savez(os.path.join(d, f"c6_profile_{S}.npz"), z_A=z, dV_xy_eV=prof)
        out["sizes"][S] = r; save_json(os.path.join(d, "c6_results.json"), out)
        log(f"[C.6] {S} : SCF d {'convergé' if r['scf_d']['converged'] else 'NON convergé'} en {r['scf_d']['n_iter_reported']} it. (précision {r['scf_d']['last_scf_accuracy_Ry']:.1e} Ry), "
            f"p {'convergé' if r['scf_p']['converged'] else 'NON convergé'} en {r['scf_p']['n_iter_reported']} it. ({r['scf_p']['last_scf_accuracy_Ry']:.1e} Ry) ; ⟨ΔV⟩_3D {r['mean3d_meV']:+.3f} meV "
            f"(A.1 {r['A1_mean3d_meV']}), N²⟨ΔV⟩ {r['N2_mean3d_meV']:+.1f} meV ; vide (|z| > 5 Å) {r['vacuum_mean_meV']:+.3f} ± {r['vacuum_std_meV']:.3f} meV ; feuillet {r['sheet_plane_meV']:+.1f} meV ; "
            f"{r['elapsed_s']:.0f} s")
    c6_tables(out); c6_figure(out)


def c6_tables(out):
    L = ["# R10 — C.6 : anomalie ⟨ΔV⟩_3D d'A.1 (tables générées par r10_driver.py c6)", "",
         "## SCF (scf.out, scf.in)", "", "| N | d : convergé ; itérations ; dernière précision (Ry) | p : convergé ; itérations ; dernière précision (Ry) | conv_thr ; mixing_beta ; K_POINTS (d) | E_F d ; p (eV) |",
         "|---|---|---|---|---|"]
    for S, r in out["sizes"].items():
        sd, sp = r["scf_d"], r["scf_p"]
        L.append(f"| {S} | {'oui' if sd['converged'] else 'NON'} ; {sd['n_iter_reported']} ; {sd['last_scf_accuracy_Ry']:.1e} | {'oui' if sp['converged'] else 'NON'} ; {sp['n_iter_reported']} ; "
                 f"{sp['last_scf_accuracy_Ry']:.1e} | {sd['params']['conv_thr']} ; {sd['params']['mixing_beta']} ; {sd['k_points']} | {sd['fermi_eV']} ; {sp['fermi_eV']} |")
    L += ["", "## ⟨ΔV⟩(z) moyenné dans le plan (ΔV = V_d − V_p, eV ; z relatif au feuillet)", "",
          "| N | c (Å) | nz | ⟨ΔV⟩_3D (meV) ; A.1 | N²⟨ΔV⟩ (meV) | part |z| ≤ 1 / 2 / 3 / 5 Å (meV) | vide |z| > 5 Å : moyenne ± écart-type ; [min, max] (meV) | bord de cellule (meV) | plan du feuillet (meV) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for S, r in out["sizes"].items():
        c = r["contributions"]
        L.append(f"| {S} | {fr(r['c_A'], 3)} | {r['nz']} | {fr(r['mean3d_meV'], -3)} ; {fr(r['A1_mean3d_meV'], -2) if r['A1_mean3d_meV'] is not None else '—'} | {fr(r['N2_mean3d_meV'], -1)} | "
                 + " / ".join(fr(c[k]["sum_over_nz_meV"], -3) for k in c) + f" | {fr(r['vacuum_mean_meV'], -3)} ± {fr(r['vacuum_std_meV'], 3)} ; [{fr(r['vacuum_min_meV'], -3)} ; {fr(r['vacuum_max_meV'], -3)}] | "
                 f"{fr(r['at_cell_edge_meV'], -3)} | {fr(r['sheet_plane_meV'], -1)} |")
    L += ["", "« part |z| ≤ z_c » = Σ_{plans |z| ≤ z_c} ⟨ΔV⟩_xy / nz : contribution de la tranche à ⟨ΔV⟩_3D (la somme sur tous les plans redonne ⟨ΔV⟩_3D)."]
    with open(os.path.join(WORK, "c", "C6_tables.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def c6_figure(out):
    plt, pal = fig_style(); fig, axs = plt.subplots(2, 1, figsize=(6.5, 6.4))
    cols = [pal.NAVY, pal.ORANGE, pal.GREEN, pal.PINK, pal.GOLD]
    for S in out["sizes"]:
        Z = np.load(os.path.join(WORK, "c", f"c6_profile_{S}.npz")); o = np.argsort(Z["z_A"]); n = size_n(S)
        big = n >= 15; col = cols[[15, 18, 21, 24, 27].index(n)] if big else pal.LIGHT
        for ax in axs:
            ax.plot(Z["z_A"][o], Z["dV_xy_eV"][o] * 1e3, color=col, lw=1.0 if big else 0.6, label=(S if big else None))
    axs[0].set_yscale("symlog", linthresh=1.0); axs[0].set_ylabel(r"$\langle\Delta V\rangle_{xy}$ (meV)"); axs[0].set_title("(a) profil complet (5 à 12 en gris)", loc="left", fontsize=9)
    axs[1].set_ylim(-5, 5); axs[1].set_ylabel(r"$\langle\Delta V\rangle_{xy}$ (meV)"); axs[1].set_xlabel(r"$z - z_\mathrm{feuillet}$ (Å)")
    axs[1].set_title(r"(b) zoom sur $\pm 5$ meV", loc="left", fontsize=9); axs[0].legend(ncol=5, fontsize=7, loc="upper right")
    fig.tight_layout(); savefig(fig, "dV_z_profiles"); plt.close(fig)


# ----------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="R10 : pilote unique (voir la docstring)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("a0", "a1", "a2", "a3", "atables", "btables", "c0", "c1post", "c2", "c3", "c4", "c5", "c6"):
        sub.add_parser(c)
    p = sub.add_parser("b"); p.add_argument("--parts", default="b0,b1,b2")
    a = ap.parse_args()
    log(f"=== {a.cmd} {vars(a)} ; HEAD {git_head()} ; pilote md5 {md5(os.path.abspath(__file__))[:12]} ; job {os.environ.get('SLURM_JOB_ID', '-')}")
    {"a0": cmd_a0, "a1": cmd_a1, "a2": cmd_a2, "a3": cmd_a3, "atables": cmd_atables, "b": cmd_b, "btables": cmd_btables,
     "c0": cmd_c0, "c1post": cmd_c1post, "c2": cmd_c2, "c3": cmd_c3, "c4": cmd_c4, "c5": cmd_c5, "c6": cmd_c6}[a.cmd](a)


if __name__ == "__main__":
    main()
