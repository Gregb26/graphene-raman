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

from electron_defect_interaction.config import load_production, dense_paths, results_dir, HA2EV  # noqa: E402
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance  # noqa: E402
from electron_defect_interaction.io import qe_gamma_io as qg  # noqa: E402
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order  # noqa: E402
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
        z = np.load(os.path.join(results_dir(cfg), "resonance_9x9.npz"))
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
    Zm = np.load(os.path.join(results_dir(cfg), "M_analysis.npz"), allow_pickle=True); a3p = load_json(R9_A3P)["level1"]
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


# ----------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="R10 : pilote unique (voir la docstring)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for c in ("a0", "a1", "a2", "a3", "atables", "btables"):
        sub.add_parser(c)
    p = sub.add_parser("b"); p.add_argument("--parts", default="b0,b1,b2")
    a = ap.parse_args()
    log(f"=== {a.cmd} {vars(a)} ; HEAD {git_head()} ; pilote md5 {md5(os.path.abspath(__file__))[:12]} ; job {os.environ.get('SLURM_JOB_ID', '-')}")
    {"a0": cmd_a0, "a1": cmd_a1, "a2": cmd_a2, "a3": cmd_a3, "atables": cmd_atables, "b": cmd_b, "btables": cmd_btables}[a.cmd](a)


if __name__ == "__main__":
    main()
