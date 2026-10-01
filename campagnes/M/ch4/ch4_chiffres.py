#!/usr/bin/env python
"""Chiffres du chapitre 4 — source unique pour la réécriture (prompt de Greg du 2026-09-30).

Sous-commandes (toutes locales, quelques secondes, aucun calcul QE ni t/Γ) :

  check    porte de la partie 2 (les 13 C_N de config/production.json redonnés par la règle ≥ 0,75 r_max à ≤ 1e-9 eV
           depuis campagnes/R/R10_plateau/a/profiles_<S>.npz) et md5 des sources d'en-tête ;
  table    partie 1 : campagnes/M/ch4/table_v1_final.md (v1 → non aligné → final, statut mécanique, chiffres nouveaux a–l,
           lignes retirées) ;
  regions  partie 2 : campagnes/M/ch4/alignement_regions.{md,npz,pdf,png} (régions d'échantillonnage (A) ≥ 0,75 r_max,
           (B) Kumagai–Oba 2D ≥ N·a/2, (C) Kumagai–Oba 3D min(N·a/2, c/2), (D) site unique) ;
  notes    partie 3 : vérifie que les valeurs finales citées dans NOTES_TGAMMA.md sont celles des sources ; écrit le diff
           non commité de NOTES_TGAMMA.md dans NOTES_TGAMMA_partie5.diff (NOTES_TGAMMA.diff = diff de la partie 3) ;
  tex      parties 4 et 5.5 : campagnes/M/ch4/defauts_nombres.md (tous les nombres du chapitre 4 du mémoire, appariés à la table
           principale et aux compléments ; pourcentages × 100, notation scientifique entière, entiers courts en tabular) ;
           le chapitre est lu dans le dépôt du mémoire (memoire/chapitres/défauts.tex ; dépôt voisin
           msc-graphene-raman-defects, ou $MSC_THESIS, ou --tex CHEMIN) : ce dépôt n'en garde plus de copie ;
  (la sous-commande table ajoute la section « Compléments » de la partie 5 : 5.1 tab:rcut_M, 5.2 niveau 1,
   5.3 tab:échantillonnage, 5.4 scalaires lus sur les courbes des npz, avec porte sur les valeurs publiées v1 / v2)
  all      table, regions, tex, notes.

Nomenclature (donnée une fois dans table_v1_final.md) : « non aligné » = tel quel (results/M2) ; « alignement à site
unique » = Lu ; « alignement de Kumagai–Oba » = plateau (i) ; « final » = v2 + Kumagai–Oba (results/M2_plateau) ;
ΔV_PA^(N) = C_N = alignment.C_N_eV de config/production.json.

Règles : tout chiffre est copié d'un fichier source (csv, npz, json, table de rapport) ; les seuls calculs sont ceux de la
partie 2 et les colonnes arithmétiques demandées (différences, produits) ; aucun arrondi (format de la source) ; deux
sources différentes pour la même grandeur → les deux ; erreur d'archive → rapportée, jamais corrigée. Lecture seule sur
article/, results/, config/, src/, scripts/, figures/. Racine du dépôt déduite de __file__ (jamais codée en dur).
"""
import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import unicodedata
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]   # campagnes/M/ch4/ -> racine
HERE = Path(__file__).resolve().parent
ART = ROOT / "campagnes" / "R"
R6 = ART / "R6_production_corrigee"
R10 = ART / "R10_plateau"
R9 = ART / "R9_controles"
R8 = ART / "R8_kaasbjerg"
R5 = ART / "R5_base_vs_M"
R4 = ART / "R4_quasi_lie"
M2 = ROOT / "results" / "M2"
M2P = ROOT / "results" / "M2_plateau"
CFG = ROOT / "config" / "production.json"
TAB_R6 = R6 / "etape3" / "table_v1_v2.md"
TAB_R10 = R10 / "c" / "table_v2_plateau.md"
MD5SUMS = M2P / "MD5SUMS_2026-09-30.txt"
THESIS = Path(os.environ.get("MSC_THESIS") or ROOT.parent / "msc-graphene-raman-defects")  # dépôt du mémoire (lecture seule)
TEX = THESIS / "memoire" / "chapitres" / "défauts.tex"
NOTES = ROOT / "notes" / "NOTES_TGAMMA.md"  # à la racine du dépôt jusqu'au 2026-09-30
SIZES13 = ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12", "15x15", "18x18", "21x21", "24x24", "27x27"]
BOHR_A = 0.529177210903

NUM_RE = re.compile(r"[-+−]?\d+(?:[.,]\d+)?(?:e[-+]?\d+)?", re.I)


# ----------------------------------------------------------------------------------------------------------------------
# utilitaires
# ----------------------------------------------------------------------------------------------------------------------
def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def head(repo=ROOT):
    try:
        return subprocess.check_output(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:  # noqa: BLE001
        return "?"


def rel(p):
    return os.path.relpath(str(p), str(ROOT))


def jload(p):
    return json.load(open(p, encoding="utf-8"))


def npz(p):
    return np.load(p, allow_pickle=True)


def fnum(x):
    """Format d'un flottant lu dans un npz/json : repr (plus court aller-retour), sans arrondi."""
    if isinstance(x, (np.floating, float)):
        return repr(float(x))
    if isinstance(x, (np.integer, int)):
        return str(int(x))
    return str(x)


def parse_numbers(cell):
    """Tous les nombres d'une cellule (ordre d'apparition), '—' → []."""
    if cell is None:
        return []
    c = cell.replace("−", "-").replace("\u202f", "").replace("\u00a0", "")
    if c.strip() in ("—", "-", ""):
        return []
    out = []
    for m in NUM_RE.finditer(c):
        s = m.group(0).replace(",", ".")
        try:
            out.append(float(s))
        except ValueError:
            pass
    return out


def decimals(cell):
    d = 0
    for m in re.finditer(r"\d+[.,](\d+)", cell or ""):
        d = max(d, len(m.group(1)))
    return d


def md_table(header, rows):
    out = ["| " + " | ".join(re.sub(r"(?<!\\)\|", "\\|", str(h)) for h in header) + " |", "|" + "---|" * len(header)]
    for r in rows:
        out.append("| " + " | ".join(re.sub(r"(?<!\\)\|", "\\|", str(c)) for c in r) + " |")
    return "\n".join(out)


def esc(s):
    return str(s)  # l'échappement des « | » est fait une seule fois par md_table


# ----------------------------------------------------------------------------------------------------------------------
# check (porte de la partie 2)
# ----------------------------------------------------------------------------------------------------------------------
def region_A(dist, r_max):
    return dist >= 0.75 * r_max


def cmd_check(args=None, quiet=False):
    cfg = jload(CFG)
    C = cfg["alignment"]["C_N_eV"]
    a1 = jload(R10 / "a" / "a1_results.json")
    worst = 0.0
    lines = []
    for S in SIZES13:
        d = npz(R10 / "a" / f"profiles_{S}.npz")
        m = region_A(d["dist"], float(d["r_max"]))
        cn = float(d["shift10"][m].mean())
        dev = cn - C[S]
        worst = max(worst, abs(dev))
        lines.append((S, int(m.sum()), a1["sizes"][S]["plateau_i_10"]["n"], cn, C[S], dev))
    ok = worst <= 1e-9
    if not quiet:
        print(f"[check] porte 0.4 : max|C_N(règle) − C_N(config)| = {worst:.3e} eV → {'PASS' if ok else 'FAIL'}")
        for S, n, nj, cn, cc, dev in lines:
            print(f"   {S:6s} n={n:4d} (json {nj:4d})  C={cn:+.16f}  config={cc:+.16f}  écart={dev:+.1e}")
        print(f"[check] md5 a1_results.json {md5(R10 / 'a' / 'a1_results.json')} (config source_md5 {cfg['alignment']['source_md5']})")
        print(f"[check] md5 table_v1_v2.md {md5(TAB_R6)} ; table_v2_plateau.md {md5(TAB_R10)} ; MD5SUMS {md5(MD5SUMS)}")
    if not ok:
        raise SystemExit("porte 0.4 refusée")
    return lines


# ----------------------------------------------------------------------------------------------------------------------
# partie 1 : table v1 → final
# ----------------------------------------------------------------------------------------------------------------------
def parse_table(p):
    rows = []
    for i, l in enumerate(open(p, encoding="utf-8"), 1):
        if not l.startswith("|") or l.startswith("|---") or l.startswith("| grandeur"):
            continue
        cells = [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", l.strip())[1:-1]]
        if len(cells) != 5:
            raise SystemExit(f"{p} l.{i} : {len(cells)} cellules")
        rows.append({"line": i, "g": cells[0], "c1": cells[1], "c2": cells[2], "src": cells[3], "lieu": cells[4]})
    return rows


def norm(s):
    s = unicodedata.normalize("NFKC", s).lower().replace("×", "x").replace("’", "'")
    return re.sub(r"\s+", " ", s).strip()


# réécritures mécaniques d'étiquettes (R6 → R10), listées dans le rapport
REWRITE_R6 = [
    (" : brut ; aligné (D6)", ""), (" : brut ; aligné", ""),
    ("médiane γ (mev) / re σ (mev)", "médiane γ / re σ (mev)"),
    ("médiane de la courbe γ_t sur ±3 ev", "médiane courbe γ_t ±3 ev"),
    (" (resonance_metrics)", ""), ("(mev, resonance_metrics)", "(mev)"),
    ("(à ε − e_d, ev)", "(à ε − e_d)"), ("(ratio = 1, ev)", "(ratio = 1)"),
    ("c14 redéfini (v_loc + c·1) :", "c14 m2 (±25 mev, v_loc tel quel) :"),
    ("(matrice complète)", "(complet)"),
]
# critères : 3 lignes R6 (det, λ_min, Friedel) ↔ 1 ligne R10 combinée ; la table finale garde les 3 lignes R6
CRIT_R6 = re.compile(r"^(critère|friedel) 9x9 \((complet|bloc π|bloc σ)\) : (.*)$")


def key_r6(g):
    k = norm(g)
    for a, b in REWRITE_R6:
        k = k.replace(a, b)
    return k


def key_r10(g):
    k = norm(g)
    return k.replace(" : brut ; aligné (d6)", "").replace(" : brut ; aligné", "")


def split_r10_criteria(cell):
    """'1.261e-04 (-0.785) ; |λ| 0.0019 (-0.787) ; Friedel -1.0005' → (det, lam, friedel)."""
    parts = [p.strip() for p in cell.split(";")]
    if len(parts) != 3:
        return None
    return parts[0], parts[1], parts[2].replace("Friedel", "").strip()


# lignes v1 seulement dont la valeur finale est lue dans un npz final (D1)
D1_NPZ = {
    48: (M2P / "mwr_locality.npz", "9x9_dense_w", lambda d: ", ".join(f"{v:.3f}" for v in d["9x9_dense_w"][1:5])),
    51: (M2P / "mwr_locality.npz", "12x12_dense_w", lambda d: ", ".join(f"{v:.3f}" for v in d["12x12_dense_w"][1:5])),
    77: (M2P / "resonance_9x9.npz", "peak_GB", lambda d: f"{float(d['peak_GB']):.3f}"),
    80: (M2P / "resonance_9x9.npz", "peak_rho_dis", lambda d: f"{float(d['peak_rho_dis']):.3f}"),
}
# lignes « retiré (décision) » (numéro de ligne R6 → référence) ; les E_res sont repérés par leur étiquette
RETIRE_R6 = {75: "C14 à 67 meV (R10 C.1 : C14 = ±9,05 meV)", 89: "zéro de det à −2,530 eV = doublet σ (R4 D3) ; critère par bloc (R10 C.1)",
             90: "idem (R4 D3 ; R10 C.1)", 112: "test « convention intensive » six tailles (R10 C.3 : lignes par famille)"}
for _l in range(117, 125):
    RETIRE_R6[_l] = "C14 à ±25 meV (R10 C.1 : C14 = ±9,05 meV)"
RETIRE_R10 = {180: "test « convention intensive » six tailles (R10 C.3)"}
for _l in range(152, 160):
    RETIRE_R10[_l] = "C14 à ±25 meV (R10 C.1)"


def is_eres(g):
    k = norm(g)
    return k.startswith("e_res") or ": e_res états" in k or ": e_res_states" in k


def statut(v1, final, v2=None):
    n1, nf = parse_numbers(v1), parse_numbers(final)
    if not n1 and not nf:
        return "sans valeur", "—"
    if not n1:
        return "nouveau", "—"
    if not nf:
        return "sans équivalent final", "—"
    if len(n1) != len(nf):
        return "remplacé", "—"
    same = all(abs(a - b) <= 1e-6 * max(abs(a), abs(b), 1e-300) for a, b in zip(n1, nf))
    d = max(decimals(v1), decimals(final))
    diff = " / ".join(f"{b - a:+.{d}f}" for a, b in zip(n1, nf))
    return ("inchangé" if same else "remplacé"), diff


def build_main_rows():
    A = parse_table(TAB_R6)
    B = parse_table(TAB_R10)
    kb = {}
    for r in B:
        kb.setdefault(key_r10(r["g"]), r)
    used = set()
    rows, notes = [], []
    for a in A:
        k = key_r6(a["g"])
        b = kb.get(k)
        fin, src_f, lieu, v2_r10, r10_line = None, None, a["lieu"], None, None
        m = CRIT_R6.match(k)
        if b is None and m:
            kind, block, what = m.groups()
            k10 = f"critère 9x9 ({block}) : min|det|/max ; min|λ| ; ∫δρ"
            b10 = kb.get(k10)
            if b10 is not None:
                sp = split_r10_criteria(b10["c2"])
                sp2 = split_r10_criteria(b10["c1"])
                idx = 0 if kind == "critère" and "det" in what else (1 if kind == "critère" else 2)
                fin = sp[idx]
                v2_r10 = sp2[idx]
                src_f = f"{b10['src']} ; `table_v2_plateau.md` l.{b10['line']} (cellule combinée, élément {idx + 1})"
                used.add(b10["line"])
                r10_line = b10["line"]
                if idx == 1:
                    notes.append(f"R6 l.{a['line']} : v1/v2 = λ_min complexe, final = |λ| (forme de la table R10)")
                if idx == 2:
                    d = npz(M2P / f"resonance_criteria_9x9.npz")
                    suf = {"complet": "", "bloc π": "_pi", "bloc σ": "_sigma"}[block]
                    fin = f"{float(d['sumrule' + suf]):.4f} / {float(d['sumrule_lloyd' + suf]):.4f} / {float(d['sumrule_window' + suf]):+.3f}"
                    src_f = f"`results/M2_plateau/resonance_criteria_9x9.npz` : sumrule{suf}, sumrule_lloyd{suf}, sumrule_window{suf} (D1 ; table R10 l.{b10['line']} : {sp[2]})"
        elif b is not None:
            fin, src_f, v2_r10, r10_line = b["c2"], f"{b['src']} ; `table_v2_plateau.md` l.{b['line']}", b["c1"], b["line"]
            used.add(b["line"])
            if norm(b["src"]) != norm(a["src"]):
                notes.append(f"R6 l.{a['line']} / R10 l.{b['line']} : fichier source différent ({a['src']} vs {b['src']})")
        elif a["line"] in D1_NPZ:
            p, key, f = D1_NPZ[a["line"]]
            fin = f(npz(p))
            src_f = f"`{rel(p)}` : {key} (D1, hors table R10)"
        if v2_r10 is not None and norm(v2_r10) != norm(a["c2"]):
            notes.append(f"R6 l.{a['line']} v2 « {a['c2']} » ≠ R10 v2 « {v2_r10} »")
        st, diff = statut(a["c1"], fin)
        if a["line"] in RETIRE_R6:
            st = f"retiré (décision) — {RETIRE_R6[a['line']]}"
        elif is_eres(a["g"]):
            st = f"retiré (décision) — E_res (R10 B.1) ; mécanique : {st}"
        v2 = a["c2"] if v2_r10 is None or norm(v2_r10) == norm(a["c2"]) else f"{a['c2']} [R10 : {v2_r10}]"
        rows.append({"g": a["g"], "v1": a["c1"], "v2": v2, "final": fin if fin is not None else "—", "diff": diff, "statut": st,
                     "src": src_f or f"— (v1 : `table_v1_v2.md` l.{a['line']}, {a['src']})", "lieu": lieu, "r6": a["line"],
                     "r10": r10_line})
    n_matched = sum(1 for r in rows if r["r10"])
    for b in B:
        if b["line"] in used:
            continue
        st, diff = statut("—", b["c2"])
        if b["line"] in RETIRE_R10:
            st = f"retiré (décision) — {RETIRE_R10[b['line']]}"
        elif is_eres(b["g"]):
            st = f"retiré (décision) — E_res (R10 B.1) ; mécanique : {st}"
        rows.append({"g": b["g"], "v1": "—", "v2": b["c1"], "final": b["c2"], "diff": diff, "statut": st,
                     "src": f"{b['src']} ; `table_v2_plateau.md` l.{b['line']}", "lieu": b["lieu"], "r6": None, "r10": b["line"]})
    counts = {"R6": len(A), "R10": len(B), "appariées": n_matched, "v1 seulement": sum(1 for r in rows if r["r6"] and not r["r10"]),
              "final seulement": sum(1 for r in rows if r["r6"] is None), "total": len(rows)}
    return rows, notes, counts


# ---- chiffres nouveaux a–l -------------------------------------------------------------------------------------------
def report_table_lines(path, first, last):
    """Lignes d'une table markdown d'un rapport (table de rapport = source autorisée), avec numéro de ligne."""
    lines = open(path, encoding="utf-8").read().split("\n")
    return [(i, lines[i - 1]) for i in range(first, last + 1) if lines[i - 1].startswith("|") and not lines[i - 1].startswith("|---")]


def cells(line):
    return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", line.strip())[1:-1]]


def loc(x, w2=None):
    return f"{x:.4f}" + (f" ({w2:.3f})" if w2 is not None else "")


def sec_a():
    out = ["### a. Porte A.2 (ΔV appliqué directement aux états de Bloch purs contre M/N_cells)", "",
           "Alignement : aucun (M bruts, non alignés) ; E_D sans objet (éléments de matrice). Sources : R5 A.2 = table de `campagnes/R/R5_base_vs_M/R5_rapport.md` "
           "(l.374–377 ; détail par paire dans `a/a_results.json` : A2.pure[\"16\"/\"20\"]) ; R6 1.3 = `campagnes/R/R6_production_corrigee/etape1/gate/gate_table.md` "
           "(json `gate_<S>_<niveau>.json`) ; R6 1.2 = `etape1/assemble_summary.jsonl`.", "", "R5 A.2 (9×9, convention de production v1 = M^L en norme super-cellule) :", ""]
    tl = report_table_lines(R5 / "R5_rapport.md", 374, 377)
    hdr = ["ligne"] + cells(tl[0][1])
    out.append(md_table(hdr, [[f"R5_rapport.md l.{i}"] + cells(l) for i, l in tl[1:]]))
    out += ["", "R6 1.3, porte A.2 sur chaque M2 (seuil 1e-6 eV ; `gate_table.md`) :", ""]
    gl = [(i, l) for i, l in enumerate(open(R6 / "etape1" / "gate" / "gate_table.md", encoding="utf-8"), 1) if l.startswith("|") and not l.startswith("|---")]
    hdr = ["ligne"] + cells(gl[0][1])
    out.append(md_table(hdr, [[f"gate_table.md l.{i}"] + cells(l) for i, l in gl[1:]]))
    out += ["", "R6 1.2, réassemblage M2 = N_cells·M^L + M^NL (`assemble_summary.jsonl`, une ligne par fichier ; max|·| en Ha) :", ""]
    rows = []
    for i, l in enumerate(open(R6 / "etape1" / "assemble_summary.jsonl", encoding="utf-8"), 1):
        if "{" not in l:
            continue
        j = json.loads(l[l.index("{"):])
        name = (j.get("out_m") or j.get("out_l") or j.get("out") or "?").split("/")[-1]
        ma, ck = j.get("max_abs", {}), j.get("checks", {})
        rows.append([f"l.{i}", name, j.get("n_cells", "—"), fnum(ma.get("L", "—")), fnum(ma.get("NL", "—")), fnum(ma.get("M", "—")), fnum(ck.get("herm_rel", "—")), j.get("ok", "—")])
    out.append(md_table(["ligne", "fichier", "N_cells", "max|N_cells·M^L|", "max|M^NL|", "max|M2|", "herm. rel.", "ok"], rows))
    return "\n".join(out)


def sec_b():
    d4 = jload(R6 / "etape2" / "d4" / "d4_results.json")
    a2 = jload(R10 / "a" / "a2_results.json")["rows"]["9x9"]["R5 C"]
    a3 = jload(R10 / "a" / "a3_results.json")
    ED = d4["E_D"]
    align_qe = f"site unique, Lu = {fnum(a2['Lu_eV'])} eV (R5 C)"
    order = [("QE_D1", "QE (D1, R6 2.1 ; E_D = E_D.SC_P, quadruplet de la parfaite)", ED["SC_P"], align_qe), ("a1_tot", "(a1) M2 16 bandes, total", ED["uc16_at_K"], "aucun (non aligné)"),
             ("a1_L", "(a1) M2^L seul", ED["uc16_at_K"], "aucun"), ("a1_NL", "(a1) M^NL seul", ED["uc16_at_K"], "aucun"),
             ("a2_16", "(a2) M2 dense ⊂ 81 k, 16 b.", ED["uc16_at_K"], "aucun"), ("a3_20", "(a3) idem 20 b.", ED["uc20dense_at_K"], "aucun"),
             ("b_5wf", "(b) 5 WF, V†εV ⊕ M2_W/81", ED["wannier_at_K"], "aucun"), ("c_all", "(c-all) H(R) + M2_W replié", ED["wannier_at_K"], "aucun"),
             ("c_3", "(c-3) idem, R_cut 3", ED["wannier_at_K"], "aucun")]
    rows = []
    for key, label, ed, al in order:
        v = d4["variants"][key]
        ev = "; ".join(loc(x, w) for x, w in v["blocks"]["even"]["localized"]) or "aucun"
        od = "; ".join(loc(x, w) for x, w in v["blocks"]["odd"]["localized"]) or "aucun"
        rs = d4["rigid_shift_fit"].get(key)
        rows.append([label, fnum(v.get("E_D", ed)), al, (f"{rs['shift_eV'] * 1e3:+.1f} ({rs['n_states']} états)" if rs else "—"), ev, od, f"`d4_results.json` : variants.{key}"])
    rows.append(["QE, alignement de Kumagai–Oba (R10 A.2)", fnum(a2["E_D"]), f"C_9 = {fnum(a2['C_N_eV'])} eV", "—",
                 "; ".join(loc(s["x_plateau"], s["w2"]) for s in a2["sigma"]), loc(a2["pi"]["x_plateau"], a2["pi"]["w2"]), "`R10_plateau/a/a2_results.json` : rows.9x9.\"R5 C\".*.x_plateau"])
    r9 = a3["rows"]["9x9"]
    rows.append([f"(c-3) « aligné » de R10 A.3 = modèle de R9 C.1 avec C_9 à site unique ({fnum(r9['C9_model_eV'])} eV), **pas Kumagai–Oba** ; aucune valeur K–O",
                 fnum(jload(R9 / "c" / "c_results.json")["E_D_wannier_K"]), f"site unique ({a3['note']})", "—", loc(r9["aligne"]["sigma"]), loc(r9["aligne"]["pi"]),
                 "`R10_plateau/a/a3_results.json` : rows.9x9.aligne ; `R9_controles/c/c_results.json`"])
    rows.append(["(c-3) « brut » de R10 A.3 (= R6 (c-3), non aligné)", "idem", "aucun", "—", loc(r9["brut"]["sigma"]), loc(r9["brut"]["pi"]), "`a3_results.json` : rows.9x9.brut"])
    out = ["### b. Escalier D4 9×9 (états localisés, ε − E_D en eV ; (w₂) ; R6 2.1)", "",
           "E_D par marche = clé `E_D` de chaque variante (`d4_results.json` : E_D.{SC_P, uc16_at_K, uc20dense_at_K, wannier_at_K}) ; « décalage rigide » = `rigid_shift_fit.shift_eV` "
           "(269 états < E_D − 4 eV, QE aligné site unique contre modèle sur son E_D). Portes de R6 2.1 : gate1 " + fnum(d4["gates"]["gate1_fold_HR"]["max_dev_eV"]) +
           " eV, gate2 (k = m/27) " + fnum(d4["gates"]["gate2_fold_Mwr"]["max_dev_eV_snapped_k"]) + " eV.", "",
           md_table(["marche", "E_D (eV)", "alignement", "décalage rigide (meV)", "σ (pairs) localisés", "π (impairs) localisés", "source"], rows)]
    return "\n".join(out)


def sec_c():
    b6 = jload(R6 / "etape2" / "b3" / "b3_results_M2.json")
    b5 = jload(R5 / "b" / "b3_results.json")
    rows = []
    for n, (x, w) in b6["pi_vs_n"]:
        p6 = b6["per_n"][str(n)]
        sig6 = "; ".join(loc(a, b) for a, b in p6["blocks"]["even"]["localized"]) or "aucun"
        p5 = b5["per_n"][str(n)]
        x5, w5 = [t for t in b5["pi_vs_n"] if t[0] == n][0][1:]
        sig5 = "; ".join(loc(a, b) for a, b in p5["blocks"]["even"]["localized"]) or "aucun"
        x81 = [t for t in b5["pi_vs_n_Lx81"] if t[0] == n][0][1]
        rows.append([n, p6["dim"], loc(x, w), sig6, loc(x5, w5), sig5, f"{x81:.4f}"])
    out = ["### c. Convergence en bandes, 9×9, (a1) à n bandes (R6 2.2 ; R5 B.3)", "",
           f"Non aligné (M2 nb128 réassemblé, `M_ed_9x9_nb128`) ; E_D = {fnum(b6['E_D_128'])} eV (`b3_results_M2.json` : E_D_128) ; premier n avec une paire σ à < 0,3 eV de +0,101 : "
           f"{b6['first_n_sigma_near_0p101']} (v2) / {b5['first_n_sigma_near_0p101']} (v1, convention de production) / {b5['first_n_sigma_near_0p101_Lx81']} (v1 diagnostique M^L × 81). QE : π −0,737 (0,246), σ +0,101 ×2 (0,713) (R6 2.2).", "",
           md_table(["n", "dim", "π v2 (ε − E_D ; w₂)", "σ v2", "π v1 (production)", "σ v1", "π v1 M^L × 81"], rows),
           "", "Sources : `R6_production_corrigee/etape2/b3/b3_results_M2.json` (pi_vs_n, per_n[n].blocks.even.localized) ; `R5_base_vs_M/b/b3_results.json` (pi_vs_n, per_n, pi_vs_n_Lx81)."]
    return "\n".join(out)


def families_from_A2():
    fam = {}
    for l in open(R10 / "a" / "A2_tables.md", encoding="utf-8"):
        if l.startswith("| ") and "×" in l.split("|")[1]:
            c = cells(l)
            fam[c[0].replace("×", "x")] = c[1]
    return fam


def sec_d():
    a2 = jload(R10 / "a" / "a2_results.json")["rows"]
    fam = families_from_A2()
    rows = []
    for S in SIZES13:
        src = "R7 D1" if "R7 D1" in a2[S] else "R5 C"
        r = a2[S][src]
        rows.append([S, fam.get(S, "?"), fnum(r["E_D"]), f"{r['Lu_eV'] * 1e3:.2f}", f"{r['C_N_eV'] * 1e3:.4f}", r["pi"]["band"], f"{r['pi']['x_Lu']:.4f} → {r['pi']['x_plateau']:.4f}", f"{r['pi']['w2']:.3f}",
                     ", ".join(str(s["band"]) for s in r["sigma"]), " / ".join(f"{s['x_Lu']:.4f} → {s['x_plateau']:.4f}" for s in r["sigma"]), " / ".join(f"{s['w2']:.3f}" for s in r["sigma"]), src])
    return "\n".join(["### d. Niveaux QE π et σ par taille, site unique → Kumagai–Oba (R10 A.2)", "",
                      "x = ε − C − E_D (eV) ; x_plateau = x_Lu + Lu − C_N ; w₂ inchangé ; E_D par taille (quadruplet pour 3m, maille pour N ≠ 3m) ; famille lue dans `a/A2_tables.md`.", "",
                      md_table(["N", "famille", "E_D (eV)", "Lu (meV)", "C_N (meV)", "π bande", "π : x site unique → x K–O", "w₂", "σ bandes", "σ : x site unique → x K–O", "w₂", "source"], rows),
                      "", "Source : `campagnes/R/R10_plateau/a/a2_results.json` : rows[S][source].{E_D, Lu_eV, C_N_eV, pi, sigma}."])


def sec_e():
    a3 = jload(R10 / "a" / "a3_results.json")
    rows = []
    for S, r in a3["rows"].items():
        rows.append([S, f"{r['QE_pi_Lu']:.4f} → {r['QE_pi_plateau']:.4f}", f"{r['brut']['pi']:.4f}", f"{r['brut']['d_pi_vs_QE_Lu'] * 1e3:+.1f} / {r['brut']['d_pi_vs_QE_plateau'] * 1e3:+.1f}",
                     f"{r['aligne']['pi']:.4f}", f"{r['aligne']['d_pi_vs_QE_Lu'] * 1e3:+.1f} / {r['aligne']['d_pi_vs_QE_plateau'] * 1e3:+.1f}", f"{r['brut']['sigma']:.4f} / {r['aligne']['sigma']:.4f}",
                     " / ".join(f"{v:.4f}" for v in r["QE_sigma_plateau"])])
    C9 = a3["rows"]["9x9"]["C9_model_eV"]
    return "\n".join(["### e. Chaîne repliée − QE, N = 6…27 (R10 A.3)", "",
                      f"Modèles de R9 C.1 non recalculés. **La colonne « aligné » utilise C_9 à site unique = {fnum(C9)} eV (−24,65 meV), et non C_9 de Kumagai–Oba (−25,14 meV)** "
                      f"(`a3_results.json` : C9_model_label = « {a3['rows']['9x9']['C9_model_label']} » ; note : « {a3['note']} »). QE : E_D par taille (A.2), colonne « plateau » = Kumagai–Oba. Écarts en meV, chaîne − QE.", "",
                      md_table(["N", "QE π : site unique → K–O", "modèle non aligné π", "écart vs QE site unique / K–O", "modèle aligné (site unique) π", "écart", "σ modèle non aligné / aligné", "σ QE K–O"], rows),
                      "", "Source : `campagnes/R/R10_plateau/a/a3_results.json` : rows[S].{QE_pi_Lu, QE_pi_plateau, brut, aligne, QE_sigma_plateau}."])


def sec_f():
    b = jload(R10 / "b" / "b_results.json")
    rows = []
    for var, lab in (("brut", "non aligné"), ("plateau", "final (Kumagai–Oba)")):
        for g in ("120", "240", "480", "960"):
            r = b["B1"][var][g]
            rows.append([lab, f"{g}²", r["n_states"], fnum(r["E_res"]), f"{r['crown']['x']} ({r['crown']['rank']})", fnum(r["peak_GT_prod"]), fnum(r["peak_GT_fine"]), f"{r['median_GT_states_meV']:.2f}", fnum(r["C_eV"])])
    return "\n".join(["### f. E_res et pic de la courbe Γ_T contre la grille de sortie (R10 B.1 ; 9×9, R_cut 3, η 0,02, N_k^int 900, E_D Wannier)", "",
                      "E_res = argmax |Γ| des états à ±1,5 eV (grandeur retirée du ch. 4, R10 B.1 ; ici pour mémoire) ; pics de Γ_T au pas 5 meV (prod) et 2,5 meV (fine).", "",
                      md_table(["variante", "grille", "états", "E_res − E_D (eV)", "couronne x (rang)", "pic Γ_T 5 meV", "pic Γ_T 2,5 meV", "médiane états (meV)", "C (eV)"], rows),
                      "", "Source : `campagnes/R/R10_plateau/b/b_results.json` : B1[brut|plateau][grille].{E_res, crown, peak_GT_prod, peak_GT_fine, median_GT_states_meV, C_eV}."])


def sec_g():
    b = jload(R10 / "b" / "b_results.json")["B2"]["rows"]
    rows = []
    for S in ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]:
        r = b[S]
        rows.append([S, r["family"], fnum(r["gap_gamma_p_eV"]), f"{r['dE_F_meV']:.1f}", f"{r['pi']['e_minus_EF_d']:.4f} ; {r['pi']['occ']:.4f}",
                     " / ".join(f"{s['e_minus_EF_d']:.4f} ; {s['occ']:.4f}" for s in r["sigma"]), f"{r['C_N_eV'] * 1e3:.2f}", (fnum(r["max_M_dense_eV"]) if r.get("max_M_dense_eV") is not None else "—")])
    out = ["### g. Familles N mod 3 (R10 B.2, C.1, C.3 ; R9 R.2)", "",
           "B.2 : ε − E_F de la cellule avec lacune (non aligné) ; occupations XML par état (Γ, poids 2) ; gap de la parfaite à Γ.", "",
           md_table(["N", "famille", "gap parfaite Γ (eV)", "ΔE_F (meV)", "π : ε − E_F ; occ.", "σ : ε − E_F ; occ. (chacun)", "C_N (meV)", "max|M| dense (eV)"], rows),
           "", "Source : `campagnes/R/R10_plateau/b/b_results.json` : B2.rows[S].", "",
           "Médianes de niveau 1 par taille et par famille (final, 240², η 0,02, R_cut 3, N_k^int 300 ; `results/M2_plateau/level2_families.csv`, lignes dans l'ordre du fichier) :", ""]
    rows = []
    with open(M2P / "level2_families.csv", encoding="utf-8") as f:
        for i, r in enumerate(csv.DictReader(f), 2):
            rows.append([f"l.{i}", r["size"], r["family"], r["vacancy_sublattice"], r["median_Gamma_Ncells_meV"], r["onsite_pz_vac_eV"], r["ReML_K_eV"], r["ReMNL_K_eV"], r["ReML_K_aligned_eV"]])
    out.append(md_table(["ligne", "N", "famille", "sous-réseau", "médiane Γ·N_cells (meV)", "p_z–p_z lacune (eV)", "Re M^L à K", "Re M^NL à K", "Re M^L à K aligné (D6)"], rows))
    c3 = jload(R10 / "c" / "c3_results.json")["tests_M"]
    out += ["", "**Test C.3 de R10** (son vrai nom : « convention intensive (cellule unitaire) », (max − min)/moyenne de max|M| des M denses complets, bandes 1–16, par famille ; `c/c3_results.json` : tests_M) :", "",
            md_table(["famille", "tailles", "écart relatif (valeur pleine)", "csv (`M_tests_summary.csv`)", "seuil", "verdict"],
                     [[k, ", ".join(v["sizes"]), fnum(v["spread"]), v["csv"]["valeur"], v["csv"]["seuil"], v["csv"]["verdict"]] for k, v in c3.items() if isinstance(v, dict)]),
            "", "**Test R9 R.2** (son vrai nom : `offdiag_intensive`, max|M| hors k = k′, trois variantes tel quel / site unique / Kumagai–Oba ; table de `campagnes/R/R9_controles/R9_rapport.md` l.733–737, valeurs de la clôture de R9) :", ""]
    rows = [[f"R9_rapport.md l.{i}"] + cells(l) for i, l in report_table_lines(R9 / "R9_rapport.md", 733, 737)]
    out.append(md_table(["ligne"] + rows[0][1:], rows[1:]))
    out.append("")
    out.append("(R9 l.739 : « (max − min)/moyenne : grossiers 7,98e-2 (trois variantes) ; denses six tailles 7,71e-2 / 7,61e-2 / 7,69e-2 ».)")
    return "\n".join(out)


def sec_h():
    c2 = jload(R10 / "c" / "c2_results.json")["kaasbjerg"]
    rows = []
    for lab, d in (("non aligné + Wigner-Seitz (final tel quel)", c2["D1"]["tel quel"]), ("final (Kumagai–Oba) + Wigner-Seitz", c2["D1"]["plateau"]),
                   ("R9 brut (étiquettes brutes)", c2["R9_D1"]["brut"]), ("R9 exact (F_W, C = site unique)", c2["R9_D1"]["exact"])):
        rows.append([lab] + [f"{d[k]['mean']:.3f} [{d[k]['min']:.2f} ; {d[k]['max']:.2f}]" for k in ("valence | K'", "conduction | K'", "valence | K", "conduction | K")])
    out = ["### h. Kaasbjerg (R10 C.2 : Fig. 3 ; R8 : Fig. 13–14)", "",
           "D.1 (9×9, k = K + δx̂, disques de 0,1471 Å⁻¹ sur la carte 240², eV Å², moyenne [min ; max]) ; **K′ en tête** (insensible à l'alignement). E_D Wannier ; variantes = alignement de M_W.", "",
           md_table(["variante", "valence | K′", "conduction | K′", "valence | K", "conduction | K"], rows), "",
           "D.2 (bloc 2×2 π/π* à (K, K), eV Å²) :", ""]
    rows = [[k, ", ".join(f"{v:.2f}" for v in d["eigenvalues_eVA2"]), f"{d['half_trace_eVA2']:.2f}", ", ".join(f"{v:.2f}" for v in d["row_norms_eVA2"])] for k, d in c2["D2"].items()]
    out.append(md_table(["bloc", "valeurs propres", "½ Tr", "normes de ligne"], rows))
    out.append("\nSource : `campagnes/R/R10_plateau/c/c2_results.json` : kaasbjerg.{D1, R9_D1, D2}.")
    dos = jload(R8 / "out" / "dos" / "dos_results.json")["metrics"]
    spec = jload(R8 / "out" / "spec" / "spec_results.json")["runs"]
    a7 = jload(R8 / "out" / "7a" / "7a_results.json")
    out += ["", "Fig. 13 (DOS, états/eV/maille/spin, ρ − ρ₀ sur [−1, 0] eV, grille 600², η_G 50 meV, N_k^int 900, η_t 20 meV ; R8 étape 2) et Fig. 14 (A_k, η_G 25 meV, R8 étape 3) ; alignement = variante de V_loc, E_D Wannier :", ""]
    rows = []
    for var, lab in (("tel_quel", "non aligné"), ("aligne", "final (Kumagai–Oba)")):
        for c in ("0.001", "0.01"):
            m = dos[f"{var}_600_c{c}"]
            s = spec[f"{var}_c{c}"]
            rows.append([lab, f"{float(c) * 100:g} %", fnum(m["drho_max_pos_eV"]), fnum(m["drho_max"]), fnum(m["drho_fwhm_eV"]), fnum(m["rho_at_ED"]),
                         "; ".join(f"{x:.4f} ({h:.1f})" for x, h in s["K_two_highest"]), (fnum(s["K_gap_eV"]) if s["K_gap_eV"] is not None else "— (un seul maximum)")])
    for c in ("0.001", "0.01"):
        m = dos[f"eta_unique_600_c{c}"]
        s = spec[f"eta_unique_c{c}"]
        rows.append(["η_t = η_G (R8 `eta_unique` ; **signalée, non utilisée**)", f"{float(c) * 100:g} %", fnum(m["drho_max_pos_eV"]), fnum(m["drho_max"]), fnum(m["drho_fwhm_eV"]), fnum(m["rho_at_ED"]),
                     "; ".join(f"{x:.4f} ({h:.1f})" for x, h in s["K_two_highest"]), (fnum(s["K_gap_eV"]) if s["K_gap_eV"] is not None else "—")])
    for c in ("0.1", "1.0"):
        m = a7["fig13_metrics"][c]
        dots = [x["eps_eV"] for x in a7["fig14"]["1.0"]["dots"] if x["col"] == 0] if c == "1.0" else []
        rows.append(["Kaasbjerg, valeurs lues (R8 7a)", f"{c} %", fnum(m["drho_max_pos_eV"]), fnum(m["drho_max"]), fnum(m["drho_fwhm_eV"]), fnum(m["rho_at_0"]),
                     ("; ".join(f"{e:.4f}" for e in dots) + " (points blancs, colonne K)" if dots else "colonne K saturée, non lisible"), ("0,100 (R8_rapport.md, 7a)" if c == "1.0" else "—")])
    out.append(md_table(["variante", "c_i", "position du max de ρ − ρ₀ (eV)", "hauteur", "largeur à mi-hauteur (eV)", "ρ(E_D)", "maxima de A_K (ε ; hauteur)", "gap à K (eV)"], rows))
    out.append("\nSources : `campagnes/R/R8_kaasbjerg/out/dos/dos_results.json` : metrics[<variante>_600_c<c>] ; `out/spec/spec_results.json` : runs[<variante>_c<c>].{K_two_highest, K_gap_eV} ; "
               "`out/7a/7a_results.json` : fig13_metrics[c], fig14[\"1.0\"].dots (col 0) ; gap 0,100 : `campagnes/R/R8_kaasbjerg/R8_rapport.md`, section 7a.")
    return "\n".join(out)


def sec_i():
    d = npz(M2P / "resonance_9x9_shiftL.npz")
    d2 = npz(M2 / "resonance_9x9_shiftL.npz")
    keys = ["median_GT_states_meV", "median_GT_states_shiftp9_05_meV", "median_GT_states_shiftm9_05_meV", "E_res_states", "E_res_states_shiftp9_05", "E_res_states_shiftm9_05", "peak_GT", "peak_GT_shiftp9_05", "peak_GT_shiftm9_05", "shifts_meV", "C_N_eV"]
    rows = [[k, (fnum(d[k]) if d[k].shape == () else ", ".join(fnum(v) for v in d[k]))] for k in keys if k in d.files]
    keys2 = ["median_GT_states_meV", "median_GT_states_shiftp25_meV", "median_GT_states_shiftm25_meV", "E_res_states", "E_res_states_shiftp25", "E_res_states_shiftm25", "peak_GT_shiftp25", "peak_GT_shiftm25", "shifts_meV", "ML_diag_mean"]
    rows2 = [[k, (fnum(d2[k]) if d2[k].shape == () else ", ".join(fnum(v) for v in d2[k]))] for k in keys2 if k in d2.files]
    return "\n".join(["### i. C14 à ±9,05 meV (R10 C.1 ; 9×9, V_loc aligné Kumagai–Oba, ±rms du plateau)", "",
                      "E_D Wannier ; `results/M2_plateau/resonance_9x9_shiftL.npz` :", "", md_table(["clé", "valeur"], rows), "",
                      "Pour mémoire (retiré) : C14 à ±25 meV autour du V_loc non aligné, `results/M2/resonance_9x9_shiftL.npz` :", "", md_table(["clé", "valeur"], rows2)])


def sec_j():
    def read(p):
        with open(p, encoding="utf-8") as f:
            return {int(r["R_cut"]): r for r in csv.DictReader(f)}
    a, b = read(M2 / "m_rcut_resigma.csv"), read(M2P / "m_rcut_resigma.csv")
    C9 = jload(CFG)["alignment"]["C_N_eV"]["9x9"] * 1e3
    rows, prev = [], None
    for rc in sorted(b):
        nL = int(b[rc]["nL"])
        ra, rb = float(a[rc]["med_ReSigma_meV"]), float(b[rc]["med_ReSigma_meV"])
        da = db = dn = cn = "—"
        if prev is not None:
            da, db = f"{ra - prev[0]:+.4f}", f"{rb - prev[1]:+.4f}"
            dn = nL - prev[2]
            cn = f"{C9 * dn:+.4f}"
        rows.append([rc, nL, a[rc]["med_ReSigma_meV"], b[rc]["med_ReSigma_meV"], da, db, dn, cn, a[rc]["med_Gamma_meV"], b[rc]["med_Gamma_meV"]])
        prev = (ra, rb, nL)
    return "\n".join(["### j. Re Σ médian contre R_cut (9×9, 240², η 0,02, N_k^int 300, fenêtre ±3 eV ; `m_rcut_resigma.csv`)", "",
                      f"Colonnes arithmétiques : Δ = valeur(R_cut) − valeur(R_cut − 1) ; C_9 × Δn = {C9:.4f} meV × (n_L(R_cut) − n_L(R_cut − 1)), n_L = 1, 5, 13, 29, 49.", "",
                      md_table(["R_cut", "n_L", "Re Σ méd. non aligné (meV)", "Re Σ méd. final (meV)", "Δ non aligné", "Δ final", "Δn_L", "C_9 × Δn_L (meV)", "Γ méd. non aligné (meV)", "Γ méd. final (meV)"], rows),
                      "", "Sources : `results/M2/m_rcut_resigma.csv` et `results/M2_plateau/m_rcut_resigma.csv` : med_ReSigma_meV, med_Gamma_meV, nL ; C_9 : `config/production.json` alignment.C_N_eV.9x9."])


def sec_k():
    c6 = jload(R10 / "c" / "c6_results.json")["sizes"]
    rows = []
    for S in SIZES13:
        r = c6[S]
        part = r["contributions"]["|z|<=3A"]["sum_over_nz_meV"]
        rows.append([S, fnum(r["mean3d_meV"]), fnum(r["N2_mean3d_meV"]), fnum(part), f"{r['N'] ** 2 * part:.4f}", f"{r['vacuum_mean_meV']:.3f} ± {r['vacuum_std_meV']:.3f}", fnum(r["at_cell_edge_meV"]), r["c_A"], r["nz"]])
    return "\n".join(["### k. ⟨ΔV⟩_3D et profil en z (R10 C.6)", "",
                      "Colonne arithmétique : N² × part |z| ≤ 3 Å. Aucun alignement (potentiels bruts).", "",
                      md_table(["N", "⟨ΔV⟩_3D (meV)", "N²⟨ΔV⟩ (meV)", "part |z| ≤ 3 Å (meV)", "N² × part", "vide |z| > 5 Å : moyenne ± écart-type", "bord de cellule (meV)", "c (Å)", "nz"], rows),
                      "", "Source : `campagnes/R/R10_plateau/c/c6_results.json` : sizes[S].{mean3d_meV, N2_mean3d_meV, contributions[\"|z|<=3A\"].sum_over_nz_meV, vacuum_mean_meV, vacuum_std_meV, at_cell_edge_meV, c_A, nz}."])


def sec_l():
    rows = []
    for lab, p in (("non aligné", M2 / "ed_vs_ep_24k24q_mv0.02.npz"), ("final", M2P / "ed_vs_ep_24k24q_mv0.02.npz")):
        d = npz(p)
        rows.append([lab, fnum(d["median"]), f"{fnum(d['min'])} (à {fnum(d['x_min'])})", f"{fnum(d['max'])} (à {fnum(d['x_max'])})", ", ".join(fnum(v) for v in d["crossings"]), fnum(d["median_ed"]), fnum(d["median_ep"]), fnum(d["c"]), fnum(d["T"])])
    t = [l for l in open(TAB_R10, encoding="utf-8") if l.startswith("| médiane Γ^ed")][0]
    return "\n".join(["### l. Chapitre 5 : Γ^ed/Γ^ep et Γ^ed médians (R10 C.1 ; c = 1 %, 300 K, ±3 eV, eV pour Γ)", "",
                      md_table(["variante", "médiane Γ^ed/Γ^ep", "min (à ε − E_D)", "max (à ε − E_D)", "croisements", "médiane Γ^ed (eV)", "médiane Γ^ep (eV)", "c", "T"], rows),
                      "", "Médianes sur ±1,2 eV (pas de clé scalaire dans le npz) : table R10 l.195 : « " + cells(t)[0] + " | " + cells(t)[1] + " | " + cells(t)[2] + " ».",
                      "", "Sources : `results/M2/ed_vs_ep_24k24q_mv0.02.npz`, `results/M2_plateau/ed_vs_ep_24k24q_mv0.02.npz` : median, min, x_min, max, x_max, crossings, median_ed, median_ep."])


def sec_retires():
    out = ["## Grandeurs retirées (décision), pour la traçabilité", ""]
    # E_res : déjà dans la table principale (statut) ; ici les 6 tailles final/non aligné
    def l1(p):
        with open(p, encoding="utf-8") as f:
            return {(r["size"], int(r["R_cut"]), int(r["grid"]), float(r["eta_eV"])): r for r in csv.DictReader(f)}
    a, b = l1(M2 / "level1_summary.csv"), l1(M2P / "level1_summary.csv")
    rows = [[S, a[(S, 3, 240, 0.02)]["argmax_E_minus_ED_eV"], b[(S, 3, 240, 0.02)]["argmax_E_minus_ED_eV"]] for S in ["5x5", "6x6", "7x7", "8x8", "9x9", "12x12"]]
    out += ["1. **E_res à toutes les tailles** (R10 B.1) — `level1_summary.csv` (R_cut 3, 240², η 0,02), clé argmax_E_minus_ED_eV :", "", md_table(["N", "non aligné", "final"], rows), ""]
    c = jload(R9 / "c" / "c_results.json")
    rows = [[k, fnum(v["eps_inf"]), fnum(v.get("a", "—")), fnum(v["rms"])] for k, v in c["fits"].items()] + [[k + " (9…27)", fnum(v["eps_inf"]), fnum(v.get("a", "—")), fnum(v["rms"])] for k, v in c["fits_9_27"].items()]
    out += ["2. **ε_∞ des extrapolations 1/N et 1/N²** (R9 C) — `campagnes/R/R9_controles/c/c_results.json` : fits, fits_9_27 :", "", md_table(["série | grandeur | loi", "ε_∞ (eV)", "a", "rms"], rows), ""]
    a1 = jload(R10 / "a" / "a1_results.json")["sizes"]
    rows = [[S, a1[S]["Lu"]["published_atom_1based"], fnum(a1[S]["Lu"]["published_meV"]), a1[S]["Lu"]["published_source"], fnum(a1[S]["Lu"]["dist_true_A"]), fnum(a1[S]["far_true"]["mean10_meV"]), fnum(a1[S]["C_N_eV"] * 1e3)] for S in SIZES13]
    out += ["3. **Valeurs d'alignement à site unique** (R10 A) — `a1_results.json` : sizes[S].Lu.{published_atom_1based, published_meV, dist_true_A}, far_true.mean10_meV, C_N_eV :", "",
            md_table(["N", "atome (1-based)", "site unique (meV)", "source publiée", "distance vraie (Å)", "site vraiment le plus loin (meV)", "C_N K–O (meV)"], rows), ""]
    out += ["4. **C14 à 67 meV et à ±25 meV** : table principale, R6 l.75 et l.117–124 (statut « retiré ») ; remplacé par C14 à ±9,05 meV (section i).",
            "5. **Zéro de det à −2,530 eV présenté comme état π** : table principale, R6 l.89–90 (v1) ; R4 D3 (`campagnes/R/R4_quasi_lie/R4_rapport.md` l.669–675) : vecteur propre 100 % bloc σ, doublet E des trois sp² de l'atome retiré ; remplacé par le pôle σ du critère par bloc (final −0,785/−0,787 eV).",
            "6. **Test « convention intensive » six tailles** : table principale, R6 l.112 / R10 l.180 (v1 2,7e-2, v2 7,7e-2 « À VOIR ») ; remplacé par les lignes par famille (section g)."]
    return "\n".join(out)


def cmd_table(args=None):
    rows, notes, counts = build_main_rows()
    hdr = ["grandeur", "v1 (mémoire actuel)", "non aligné (v2)", "final", "final − v1", "statut", "source finale (fichier : clé ou ligne)", "lieu"]
    body = [[esc(r["g"]), esc(r["v1"]), esc(r["v2"]), esc(r["final"]), r["diff"], r["statut"], r["src"], esc(r["lieu"])] for r in rows]
    from collections import Counter
    cs = Counter(r["statut"].split(" —")[0].split(" ;")[0] for r in rows)
    txt = [f"# Chiffres du chapitre 4 — table v1 → final", "",
           f"Date : {date.today().isoformat()} ; HEAD `{head()}` ; md5 `campagnes/R/R6_production_corrigee/etape3/table_v1_v2.md` {md5(TAB_R6)} ; "
           f"`campagnes/R/R10_plateau/c/table_v2_plateau.md` {md5(TAB_R10)} ; `results/M2_plateau/MD5SUMS_2026-09-30.txt` {md5(MD5SUMS)}. Généré par `campagnes/M/ch4/ch4_chiffres.py table` (ne pas éditer à la main).", "",
           "## Nomenclature", "",
           "- **non aligné** (= « tel quel » des archives) : ΔV brut, v2 (M2 = N_cells·M^L + M^NL), `results/M2` ;",
           "- **alignement à site unique** (= « Lu » des archives) : potentiel au site le plus éloigné de la lacune (Lu et al. 2019) ;",
           "- **alignement de Kumagai–Oba** (= « plateau (i) » des archives) : moyenne des potentiels de site (sphères de 1 Å) sur les atomes à ≥ 0,75 r_max de la lacune (Kumagai et Oba, PRB 89, 195205, 2014 ; région élargie, voir `alignement_regions.md`) ;",
           "- symbole du mémoire ΔV_PA^(N) = C_N = `alignment.C_N_eV` de `config/production.json` ;",
           "- **final** = v2 + alignement de Kumagai–Oba, `results/M2_plateau` (R10). Les archives (rapports R4–R10, tables sources) gardent leurs étiquettes.", "",
           "## Règles", "",
           "Appariement : « grandeur » normalisée (NFKC, minuscules, × → x, espaces réduits) après les réécritures d'étiquettes suivantes (R6 → R10) : " +
           " ; ".join(f"« {a} » → « {b} »" for a, b in REWRITE_R6) + " ; les 9 lignes « critère/Friedel 9x9 » de R6 (det, λ_min, Friedel × complet/π/σ) s'apparient aux 3 lignes combinées de R10 (cellule découpée sur « ; ») ; "
           "les lignes v1 seulement 48, 51, 77, 80 prennent leur valeur finale dans un npz de `results/M2_plateau` (D1). "
           "Statut mécanique : tous les nombres d'une cellule sont extraits dans l'ordre ; « inchangé » si même compte et chaque |final − v1| ≤ 1e-6 relatif ; « remplacé » sinon ; « nouveau » sans v1 ; « sans équivalent final » sans final ; "
           "« retiré (décision) » pour la liste du prompt seulement (E_res, C14 67/±25 meV, zéro de det −2,530 eV, convention intensive six tailles), avec la référence, le statut mécanique étant conservé après « ; ». "
           "« final − v1 » : élément par élément, au nombre de décimales le plus grand des deux cellules ; « — » si les comptes diffèrent. Valeurs copiées telles quelles (format des tables sources ; npz : repr).", "",
           f"Comptes : R6 {counts['R6']} lignes, R10 {counts['R10']} lignes ; appariées {counts['appariées']} ; v1 seulement {counts['v1 seulement']} ; final seulement {counts['final seulement']} ; total {counts['total']}. "
           "Statuts : " + ", ".join(f"{k} {v}" for k, v in sorted(cs.items())) + ".", "",
           "## Table principale (ordre : lignes de R6, puis lignes propres à R10)", "", md_table(hdr, body), "",
           "## Écarts entre les deux tables sources (rapportés, non corrigés)", ""] + [f"- {n}" for n in notes] + ["", "## Chiffres nouveaux", "",
           "Chaque section nomme la source (fichier : clé ou ligne), l'alignement et l'E_D utilisés.", ""]
    txt += [sec_a(), "", sec_b(), "", sec_c(), "", sec_d(), "", sec_e(), "", sec_f(), "", sec_g(), "", sec_h(), "", sec_i(), "", sec_j(), "", sec_k(), "", sec_l(), "", sec_retires(), ""]
    ctxt, _, ok4, _ = complements()
    txt += [ctxt]
    print(f"[table] compléments 5.1–5.4 ajoutés ; porte 5.4 : {'PASS' if ok4 else 'REFUSÉE (STOP sur 5.4)'}")
    (HERE / "table_v1_final.md").write_text("\n".join(txt), encoding="utf-8")
    print(f"[table] {HERE / 'table_v1_final.md'} : {counts} ; statuts {dict(cs)} ; {len(notes)} écarts entre tables")
    return rows, counts


# ----------------------------------------------------------------------------------------------------------------------
# partie 2 : régions d'alignement
# ----------------------------------------------------------------------------------------------------------------------
def lattice_from_scf(S, a1):
    """a et c de la super-cellule lus dans scf.in (CELL_PARAMETERS bohr), répertoire des Vks de a1_results.json."""
    p = Path(a1["sizes"][S]["files"]["vd"]).parent / "scf.in"
    lines = open(p, encoding="utf-8", errors="replace").read().split("\n")
    i = [k for k, l in enumerate(lines) if l.strip().startswith("CELL_PARAMETERS")][0]
    unit = lines[i].split()[1].lower() if len(lines[i].split()) > 1 else "bohr"
    vec = [np.array([float(x) for x in lines[i + j].split()[:3]]) for j in (1, 2, 3)]
    fac = BOHR_A if unit.startswith("bohr") else 1.0
    N = int(S.split("x")[0])
    return float(np.linalg.norm(vec[0]) * fac / N), float(vec[2][2] * fac), str(p)


def stats(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return 0, np.nan, np.nan, np.nan
    m = x.mean()
    rms = np.sqrt(((x - m) ** 2).mean())
    return int(x.size), float(m), float(rms), float(rms / np.sqrt(x.size))


def cmd_regions(args=None):
    a1 = jload(R10 / "a" / "a1_results.json")
    fam = families_from_A2()
    wout = open(ROOT / "results" / "wannier" / "27x27" / "wannier.wout", encoding="utf-8", errors="replace").read().split("\n")
    wa = [l for l in wout if l.strip().startswith("a_1")][0].split()
    wc = [l for l in wout if l.strip().startswith("a_3")][0].split()
    a_w, c_w = float(np.hypot(float(wa[1]), float(wa[2]))), float(wc[3])
    rows, store = [], {}
    for S in SIZES13:
        d = npz(R10 / "a" / f"profiles_{S}.npz")
        dist, sh, r_max = d["dist"], d["shift10"], float(d["r_max"])
        N = int(S.split("x")[0])
        a, c, scf = lattice_from_scf(S, a1)
        rB = N * a / 2.0
        rC = min(rB, c / 2.0)
        mA, mB, mC = dist >= 0.75 * r_max, dist >= rB - 1e-6, dist >= rC - 1e-6
        n_circle = int((np.abs(dist - rB) < 1e-6).sum())  # sites sur le cercle inscrit : inclus dans (B) (tolérance 1e-6 Å)
        sA, sB, sC = stats(sh[mA]), stats(sh[mB]), stats(sh[mC])
        lu = a1["sizes"][S]["Lu"]
        ft = a1["sizes"][S]["far_true"]
        dBA = sB[1] - sA[1]
        rows.append([S, fam.get(S, "?"), f"{a:.6f}", f"{c:.6f}", f"{r_max:.4f}", f"{rB:.4f}", f"{rC:.4f}",
                     f"{sA[0]} ; {sA[1] * 1e3:.4f} ; {sA[2] * 1e3:.4f} ; {sA[3] * 1e3:.4f}", f"{sB[0]} ; {sB[1] * 1e3:.4f} ; {sB[2] * 1e3:.4f} ; {sB[3] * 1e3:.4f}", n_circle,
                     f"{sC[0]} ; {sC[1] * 1e3:.4f} ; {sC[2] * 1e3:.4f} ; {sC[3] * 1e3:.4f}", f"{lu['published_meV']:.4f} (atome {lu['published_atom_1based']}, {lu['dist_true_A']:.3f} Å)",
                     f"{ft['mean10_meV']:.4f} (atomes {', '.join(map(str, ft['index1']))}, {r_max:.3f} Å)", f"{dBA * 1e3:+.4f}", f"{abs(dBA) / sA[2]:.4f}"])
        store[S] = dict(N=N, a=a, c=c, r_max=r_max, r_B=rB, r_C=rC, nA=sA[0], meanA=sA[1], rmsA=sA[2], nB=sB[0], meanB=sB[1], rmsB=sB[2], n_circle=n_circle, nC=sC[0], meanC=sC[1], rmsC=sC[2],
                        Lu=lu["published_meV"] * 1e-3, far_true=ft["mean10_meV"] * 1e-3, family=fam.get(S, "?"), scf=scf)
    np.savez(HERE / "alignement_regions.npz", **{f"{S}_{k}": v for S, dd in store.items() for k, v in dd.items()}, sizes=np.array(SIZES13), a_wout=a_w, c_wout=c_w)
    hdr = ["N", "famille", "a (Å)", "c (Å)", "r_max (Å)", "N·a/2 (Å)", "min(N·a/2, c/2) (Å)", "(A) ≥ 0,75 r_max : n ; moyenne ; rms ; rms/√n (meV)", "(B) K–O 2D ≥ N·a/2 : n ; moyenne ; rms ; rms/√n",
           "sites à < 1e-6 Å du cercle", "(C) K–O 3D : n ; moyenne ; rms ; rms/√n", "(D) site unique publié (meV)", "(D) site vraiment le plus loin (meV)", "(B) − (A) (meV)", "|(B) − (A)| / rms(A)"]
    txt = [f"# Régions d'échantillonnage de l'alignement, 13 tailles", "",
           f"Date : {date.today().isoformat()} ; HEAD `{head()}`. Généré par `ch4_chiffres.py regions` (ne pas éditer). Profils par atome de R10 A (`campagnes/R/R10_plateau/a/profiles_<S>.npz` : `dist` = distance vraie à la lacune en image minimale, "
           "`shift10` = potentiel de site, sphère de 1 Å, `r_max`). Statistiques : moyenne, rms = écart-type de population (ddof 0, comme `a1_results.json`), rms/√n.", "",
           "Régions : (A) nous, ≥ 0,75 r_max (porte : redonne les 13 C_N de `config/production.json` au bit) ; (B) Kumagai–Oba 2D : sites de la cellule de Wigner-Seitz plane centrée sur la lacune hors du cercle inscrit, ≥ N·a/2 (= (√3/2) r_max ; les sites à moins de 1e-6 Å du cercle sont inclus et comptés) ; "
           "(C) Kumagai–Oba 3D tel quel : hors de la sphère inscrite, rayon min(N·a/2, c/2) ; (D) site unique : valeur publiée (`a1_results.json` : Lu.published_meV) et site le plus éloigné en image minimale vraie (far_true.mean10_meV).", "",
           f"Réseau : a et c lus dans le `scf.in` de chaque taille (CELL_PARAMETERS bohr, a = |a₁|/N, × {BOHR_A} Å/bohr ; chemin = répertoire des `Vks` de `a1_results.json`) ; contrôle `results/wannier/27x27/wannier.wout` : a = {a_w:.6f} Å, c = {c_w:.6f} Å.", "",
           md_table(hdr, rows), "", "Figure : `alignement_regions.{pdf,png}` (moyenne ± rms contre N pour (A), (B), (C), (D) ; familles 3m / non-3m distinguées par le marqueur ; non installée dans `figures/`). Données : `alignement_regions.npz`.",
           "", "Sources des scf.in : " + " ; ".join(f"{S} `{store[S]['scf']}`" for S in SIZES13[:1]) + " (même motif pour les autres tailles)."]
    (HERE / "alignement_regions.md").write_text("\n".join(txt), encoding="utf-8")
    make_regions_figure(store)
    print(f"[regions] {HERE / 'alignement_regions.md'} ; npz ; figure")
    return store


def make_regions_figure(store):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        from graphene_raman.plotting.palette import CYCLE  # noqa: F401
        pal = CYCLE
    except Exception:  # noqa: BLE001
        pal = ["#000080", "#eb6834", "#2ca02c", "#e8c547", "#d94f90", "#66b3ff"]
    from graphene_raman.plotting.palette import use_style; use_style()
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    Ns = np.array([store[S]["N"] for S in SIZES13])
    fam3 = np.array([store[S]["family"] == "3m" for S in SIZES13])
    series = [("A", "(A) $\\geq 0{,}75\\,r_\\mathrm{max}$ (nous)", "meanA", "rmsA", pal[0]), ("B", "(B) Kumagai--Oba 2D, $\\geq N a/2$", "meanB", "rmsB", pal[1]),
              ("C", "(C) Kumagai--Oba 3D, $\\min(Na/2, c/2)$", "meanC", "rmsC", pal[2])]
    for i, (k, lab, mk, rk, col) in enumerate(series):
        m = np.array([store[S][mk] for S in SIZES13]) * 1e3
        r = np.array([store[S][rk] for S in SIZES13]) * 1e3
        off = (i - 1) * 0.12
        for sel, mark in ((fam3, "o"), (~fam3, "s")):
            ax.errorbar(Ns[sel] + off, m[sel], yerr=r[sel], fmt=mark, color=col, ms=4, capsize=2, lw=1, label=lab if mark == "o" else None)
    lu = np.array([store[S]["Lu"] for S in SIZES13]) * 1e3
    for sel, mark in ((fam3, "o"), (~fam3, "s")):
        ax.plot(Ns[sel] + 0.24, lu[sel], mark, mfc="none", color=pal[3], ms=4, label="(D) site unique" if mark == "o" else None)
    ax.plot([], [], "o", color="gray", label="famille $N = 3m$")
    ax.plot([], [], "s", color="gray", label="famille $N \\neq 3m$")
    ax.set_xlabel("Taille de la super-cellule $N$")
    ax.set_ylabel("$\\Delta V_\\mathrm{PA}^{(N)}$ (meV), moyenne $\\pm$ rms")
    ax.set_xticks(Ns)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(HERE / "alignement_regions.pdf")
    fig.savefig(HERE / "alignement_regions.png", dpi=200)


# ----------------------------------------------------------------------------------------------------------------------
# partie 3 : diff de NOTES_TGAMMA.md et vérification des valeurs finales citées
# ----------------------------------------------------------------------------------------------------------------------
def cmd_notes(args=None):
    # NOTES_TGAMMA.diff = diff de la partie 3 (commité, jamais réécrit) ; les modifications non commitées vont dans NOTES_TGAMMA_partie5.diff
    diff = subprocess.run(["git", "-C", str(ROOT), "diff", "--no-color", "--", str(NOTES.relative_to(ROOT))], capture_output=True, text=True).stdout
    dfile = HERE / "NOTES_TGAMMA_partie5.diff"
    if diff.strip():
        dfile.write_text(diff, encoding="utf-8")
    txt = open(NOTES, encoding="utf-8").read()
    checks = []
    def l1(p):
        with open(p, encoding="utf-8") as f:
            return {(r["size"], int(r["R_cut"]), int(r["grid"]), float(r["eta_eV"])): r for r in csv.DictReader(f)}
    L = l1(M2P / "level1_summary.csv")
    checks.append(("médiane 9×9 R_cut 3 (level1_summary.csv)", f"{float(L[('9x9', 3, 240, 0.02)]['median_Gamma_Ncells_meV']):.2f}"))
    with open(M2P / "m_rcut_resigma.csv", encoding="utf-8") as f:
        rr = {int(r["R_cut"]): r for r in csv.DictReader(f)}
    checks.append(("Re Σ médian R_cut 3 (m_rcut_resigma.csv)", f"{float(rr[3]['med_ReSigma_meV']):.2f}"))
    d = npz(M2P / "resonance_9x9.npz")
    for k, fmt in (("ReTbar_at_ED", "{:.3f}"), ("ImTbar_at_ED", "{:.3f}"), ("peak_GT", "{:.3f}"), ("peak_drho", "{:.3f}"), ("min_absReTbar", "{:.3f}")):
        checks.append((f"resonance_9x9.npz : {k}", fmt.format(float(d[k]))))
    c = npz(M2P / "resonance_criteria_9x9.npz")
    checks.append(("resonance_criteria_9x9.npz : sumrule", f"{float(c['sumrule']):.4f}"))
    checks.append(("resonance_criteria_9x9.npz : sigma_flag_det_rel", f"{float(c['sigma_flag_det_rel']):.3e}".replace("e-0", "e-")))
    cfg = jload(CFG)
    checks.append(("config alignment.C_N_eV.9x9 (meV)", f"{cfg['alignment']['C_N_eV']['9x9'] * 1e3:.4f}"))
    with open(M2P / "nkint_check_9x9.csv", encoding="utf-8") as f:
        nk = {int(r["nk_int"]): r for r in csv.DictReader(f)}
    checks.append(("nkint_check_9x9.csv : G_ED nk 300", f"{float(nk[300]['G_ED']):.2f}"))
    s = npz(M2P / "resonance_9x9_shiftL.npz")
    checks.append(("shiftL : median_GT_states_shiftp9_05_meV", f"{float(s['median_GT_states_shiftp9_05_meV']):.2f}"))
    _, _, ok4, sc = complements()
    if ok4:
        checks.append(("5.4 Born/T min (final)", sc["final"]["BT_min"]))
        checks.append(("5.4 Born/T max (final)", sc["final"]["BT_max"]))
        checks.append(("5.4 ħ/Γ (final)", sc["final"]["tau"]))
    bad = []
    for lab, val in checks:
        v = val.replace(".", ",")
        if v.startswith("-"):
            v = "−" + v[1:]  # vrai signe moins en tête (les exposants gardent le trait d'union)
        v2 = v
        if "," in v and len(v.split(",")[0]) > 3:  # espace fine des milliers dans NOTES (« 3 189,01 »)
            ip, fp = v.split(",")
            v2 = f"{ip[:-3]} {ip[-3:]},{fp}"
        ok = (v in txt) or (v2 in txt) or (val in txt)
        print(f"[notes] {'OK ' if ok else 'ABSENT'} {lab} = {val}")
        if not ok:
            bad.append(lab)
    print(f"[notes] diff non commité : {dfile if diff.strip() else 'aucun'} ({len(diff.splitlines())} lignes) ; {len(checks) - len(bad)}/{len(checks)} valeurs finales retrouvées dans NOTES_TGAMMA.md")
    return bad


# ----------------------------------------------------------------------------------------------------------------------
# partie 5 : compléments (5.1 tab:rcut_M, 5.2 niveau 1, 5.3 tab:échantillonnage, 5.4 scalaires lus sur les courbes)
# ----------------------------------------------------------------------------------------------------------------------
RES3 = (("v1", ROOT / "results" / "M"), ("non aligné", M2), ("final", M2P))


def pct(s):
    """Fraction du csv → pourcentage exact (décalage de la virgule, aucun arrondi) : '4.0851e-01' → '40.851'."""
    from decimal import Decimal
    return format(Decimal(s).scaleb(2), "f")


def compl_5_1():
    data = {}
    dups = {}
    for lab, d in RES3:
        with open(d / "m_rcut_convergence.csv", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        from collections import Counter
        c = Counter((r["size"], r["R_cut"]) for r in rows)
        dups[lab] = (len(rows), [k for k, v in c.items() if v > 1])
        data[lab] = {(r["size"], int(r["R_cut"])): r for r in rows}
    rows, idx = [], []
    for S in ("9x9", "12x12"):
        for rc in range(7):
            r = {lab: data[lab].get((S, rc)) for lab, _ in RES3}
            sv = [pct(r[lab]["sv_mismatch_pi_blocks"]) if r[lab] else "—" for lab, _ in RES3]
            dg = [pct(r[lab]["diag_max_dM_over_max"]) if r[lab] else "—" for lab, _ in RES3]
            rows.append([S, rc, r["final"]["n_sites"]] + sv + dg)
            idx.append({"g": f"5.1 tab:rcut_M {S} R_cut {rc} π–π* (%)", "v1": sv[0], "v2": sv[1], "final": sv[2]})
            idx.append({"g": f"5.1 tab:rcut_M {S} R_cut {rc} diag (%)", "v1": dg[0], "v2": dg[1], "final": dg[2]})
    eq = [f"{S} R_cut {rc}" for (S, rc), r in sorted(data["final"].items()) if r["max_dM_over_maxM"] == r["diag_max_dM_over_max"]]
    out = ["### 5.1 tab:rcut_M du mémoire : π–π* (`sv_mismatch_pi_blocks`) et diag (`diag_max_dM_over_max`), en %", "",
           "Sources : `results/M/m_rcut_convergence.csv` (v1), `results/M2/m_rcut_convergence.csv` (non aligné), `results/M2_plateau/m_rcut_convergence.csv` (final) ; grille fine 60² ; "
           "% = fraction du csv × 100 (décalage de la virgule, chiffres du csv conservés).", "",
           md_table(["taille", "R_cut", "cellules", "π–π* v1", "π–π* non aligné", "π–π* final", "diag v1", "diag non aligné", "diag final"], rows), "",
           "Lignes en double : " + " ; ".join(f"{lab} : {n} lignes, couples (taille, R_cut) répétés : {d if d else 'aucun'}" for lab, (n, d) in dups.items()) +
           ". **Aucune ligne en double dans les trois csv** (le csv v1 est seulement écrit en deux blocs : 9×9 et 12×12 R_cut 0–3, puis 9×9 et 12×12 R_cut 4–6). "
           "Colonnes égales dans le csv final (`max_dM_over_maxM` = `diag_max_dM_over_max`) : " + (", ".join(eq) or "aucune") +
           " ; `diag_max_dM_over_max` = `diag_abs_mismatch` sur toutes les lignes des trois csv."]
    return "\n".join(out), idx


def compl_5_2():
    def l1(p):
        with open(p, encoding="utf-8") as f:
            return {(r["size"], int(r["R_cut"])): r["median_Gamma_Ncells_meV"] for r in csv.DictReader(f) if int(r["grid"]) == 240 and float(r["eta_eV"]) == 0.02}
    data = {lab: l1(d / "level1_summary.csv") for lab, d in RES3}
    rows, idx = [], []
    for S in ("6x6", "9x9", "12x12", "5x5", "7x7", "8x8"):
        N = int(S.split("x")[0])
        for rc in range(5):
            v = [data[lab].get((S, rc), "—") for lab, _ in RES3]
            cite = "oui" if (rc <= 3 or S == "9x9") else "non (« --- » dans le mémoire)"
            rows.append([S, N % 3, rc] + v + [cite])
            idx.append({"g": f"5.2 niveau 1 {S} R_cut {rc} médiane Γ·N_cells (meV)", "v1": v[0], "v2": v[1], "final": v[2]})
    extra = "—"
    p = ROOT / "results" / "M" / "m_rcut_resigma.csv"
    if p.exists():
        with open(p, encoding="utf-8") as f:
            rr = {int(r["R_cut"]): r["med_Gamma_meV"] for r in csv.DictReader(f)}
        extra = rr.get(4, "—")
        idx.append({"g": "5.2 niveau 1 9x9 R_cut 4 (v1, m_rcut_resigma.csv)", "v1": extra, "v2": "—", "final": "—"})
    out = ["### 5.2 Tableau niveau 1 du mémoire (tab:convergence_gamma de `défauts.tex`, dépôt du mémoire) : médiane Γ·N_cells (meV), 240², η 0,02, N_k^int 300", "",
           "Sources : `level1_summary.csv` de `results/M` (v1), `results/M2` (non aligné), `results/M2_plateau` (final), colonne `median_Gamma_Ncells_meV` ; ordre des tailles = celui du mémoire.", "",
           md_table(["taille", "N mod 3", "R_cut", "v1", "non aligné", "final", "cité dans le mémoire"], rows), "",
           f"`results/M/level1_summary.csv` n'a aucune ligne R_cut 4 (carte v1 à R_cut 0–3). Le point 9×9, R_cut 4 du mémoire (2468) est dans `results/M/m_rcut_resigma.csv` : med_Gamma_meV = {extra}."]
    return "\n".join(out), idx


def compl_5_3():
    from decimal import Decimal
    C = jload(CFG)["alignment"]["C_N_eV"]
    rows, idx = [], []
    with open(M2P / "sampling_table.csv", encoding="utf-8") as f:
        tab = {int(r["N"]): r for r in csv.DictReader(f)}
    for N in (6, 9, 12, 5, 7, 8, 10, 11):
        r = tab[N]
        cn = C[f"{N}x{N}"] * 1e3
        raw = r["Ved_radial_1.42A_meV"]
        al = float(Decimal(raw)) - cn
        rows.append([N, r["N_mod_3"], r["dE_F_meV"], raw, f"{cn:.4f}", f"{al:+.4f}"])
        idx.append({"g": f"5.3 tab:échantillonnage N = {N} ΔE_F (meV)", "v1": r["dE_F_meV"], "v2": r["dE_F_meV"], "final": r["dE_F_meV"]})
        idx.append({"g": f"5.3 tab:échantillonnage N = {N} V̄^L(a_CC) (meV)", "v1": raw, "v2": raw, "final": f"{al:+.4f}"})
    same = md5(M2P / "sampling_table.csv") == md5(M2 / "sampling_table.csv") == md5(ROOT / "results" / "M" / "sampling_table.csv")
    out = ["### 5.3 tab:échantillonnage : ΔE_F et V̄^L_ed(a_CC), 8 tailles (meV)", "",
           f"Sources : `results/M2_plateau/sampling_table.csv` (colonnes `dE_F_meV`, `Ved_radial_1.42A_meV` ; fichier identique dans `results/M` et `results/M2` : {'md5 égaux' if same else 'md5 DIFFÉRENTS'}) ; "
           "C_N : `config/production.json` alignment.C_N_eV (× 10³, 4 décimales). Colonne arithmétique : aligné = non aligné − C_N. ΔE_F ne dépend pas de l'alignement.", "",
           md_table(["N", "N mod 3", "ΔE_F", "V̄^L(a_CC) non aligné", "C_N", "V̄^L(a_CC) aligné = non aligné − C_N"], rows)]
    return "\n".join(out), idx


def curve_scalars(res_path, crit_path):
    """Mêmes définitions que `r6_compare_v1_v2.py` (rstats, cstats) pour les lignes 73, 74, 98, 99 de table_v1_v2.md."""
    R, Cr = npz(res_path), npz(crit_path)
    eg, ED = R["eg"], float(R["E_D"])
    m = np.abs(eg - ED) <= 3.0
    ratio = R["Gamma_Born"][m] / R["Gamma_T"][m]
    x, g = Cr["x_c"], Cr["Gamma_c"] * 1e3
    return {"BT_min": f"{float(np.nanmin(ratio)):.3f}", "BT_max": f"{float(np.nanmax(ratio)):.3f}",
            "Gc": f"{float(np.nanmin(g)):.2f} ({float(x[np.nanargmin(g)]):+.2f}) / {float(np.nanmax(g)):.2f} ({float(x[np.nanargmax(g)]):+.2f}) / {float(np.interp(0, x, g)):.2f}",
            "tau": f"{658.2 / float(np.interp(-0.3, x, g)):.0f} / {658.2 / float(np.interp(0.3, x, g)):.0f}"}


def compl_5_4():
    pub = {r["line"]: r for r in parse_table(TAB_R6)}
    keys = (("BT_min", 73, "Born/T min sur ±3 eV"), ("BT_max", 74, "Born/T max sur ±3 eV"),
            ("Gc", 98, "Γ_T à c = 0,1 % sur ±1 eV : min (à) / max (à) / E_D (meV)"), ("tau", 99, "ħ/Γ à ∓0,3 eV, c = 0,1 % (fs)"))
    sc = {lab: curve_scalars(d / "resonance_9x9.npz", d / "resonance_criteria_9x9.npz") for lab, d in RES3}
    gate, ok = [], True
    for k, line, glab in keys:
        for lab, col in (("v1", "c1"), ("non aligné", "c2")):
            same = sc[lab][k] == pub[line][col]
            ok &= same
            gate.append([f"R6 l.{line}", glab, lab, pub[line][col], sc[lab][k], "OK" if same else "DIFFÉRENT"])
    out = ["### 5.4 Scalaires lus sur les courbes des npz (définitions de `table_v1_v2.md` l.73, 74, 98, 99)", "",
           "Lecture : Born/T = Γ_Born/Γ_T de `resonance_9x9.npz` (`eg`, `Gamma_Born`, `Gamma_T`) sur |ε − E_D| ≤ 3 eV, min et max ; Γ_T à c = 0,1 % = `Gamma_c` × 10³ de `resonance_criteria_9x9.npz` "
           "(`x_c` sur ±1 eV) : min (position), max (position), interpolation linéaire à 0 ; ħ/Γ = 658,2 meV·fs / Γ interpolé à −0,3 et +0,3 eV — code de `campagnes/R/R6_production_corrigee/etape3/r6_compare_v1_v2.py` (rstats, cstats).", "",
           f"**Porte** (la même lecture sur les npz v1 et non alignés redonne les valeurs publiées à la dernière décimale) : **{'PASS' if ok else 'REFUSÉE'}**.", "",
           md_table(["ligne publiée", "grandeur", "npz", "publié", "relu", "verdict"], gate)]
    idx = []
    if not ok:
        out += ["", "**STOP sur 5.4** : la porte est refusée, aucune valeur finale n'est donnée."]
        return "\n".join(out), idx, ok, sc
    rows = []
    for k, line, glab in keys:
        st, diff = statut(sc["v1"][k], sc["final"][k])
        src = "`results/M2_plateau/resonance_9x9.npz` : eg, Gamma_Born, Gamma_T" if k.startswith("BT") else "`results/M2_plateau/resonance_criteria_9x9.npz` : x_c, Gamma_c"
        rows.append([glab, sc["v1"][k], sc["non aligné"][k], sc["final"][k], diff, st, src, f"R6 l.{line} (sans équivalent final dans la table principale)"])
        idx.append({"g": f"5.4 {glab}", "v1": sc["v1"][k], "v2": sc["non aligné"][k], "final": sc["final"][k]})
    out += ["", "Valeurs finales (même lecture sur les npz de `results/M2_plateau`) :", "",
            md_table(["grandeur", "v1", "non aligné (v2)", "final", "final − v1", "statut", "source finale", "lieu"], rows)]
    return "\n".join(out), idx, ok, sc


def complements():
    t1, i1 = compl_5_1()
    t2, i2 = compl_5_2()
    t3, i3 = compl_5_3()
    t4, i4, ok4, sc = compl_5_4()
    txt = ["## Compléments (partie 5, 2026-09-30)", "",
           "Mêmes règles : valeurs copiées ou lues sur les fichiers de `results/` ; colonnes arithmétiques seulement là où c'est dit.", "", t1, "", t2, "", t3, "", t4, ""]
    return "\n".join(txt), i1 + i2 + i3 + i4, ok4, sc


# ----------------------------------------------------------------------------------------------------------------------
# partie 4 (appariement refait en 5.5) : nombres de défauts.tex
# ----------------------------------------------------------------------------------------------------------------------
TEX_STRIP = [r"\\(?:label|ref|eqref|pageref|cite[tp]?|autoref|cref|Cref|nameref|hyperref)\*?(?:\[[^\]]*\])?\{[^}]*\}", r"\\includegraphics(?:\[[^\]]*\])?\{[^}]*\}",
             r"\\(?:begin|end)\{[^}]*\}(?:\{[^}]*\})?", r"\\(?:h|v)space\*?\{[^}]*\}", r"\\(?:input|include)\{[^}]*\}", r"\\multicolumn\{\d+\}\{[^}]*\}"]
PCT = "％"  # marqueur interne de « \% »
TEX_SCI = r"(?:\d+(?:[.,]\d+)?\s*(?:\\times|\\cdot|×)\s*)?10\^\{?\s*[-+−]?\d+\s*\}?"
TEX_PLAIN = r"\d+(?:[.,]\d+)?(?:e[-+]?\d+)?"
TEX_NUM = re.compile(r"(?<![\w\^_{\\.,])[-+−]?(?:" + TEX_SCI + "|" + TEX_PLAIN + ")")


def tex_numbers(path):
    """(ligne, nombre, contexte, est_pourcentage, dans_tabular) pour chaque nombre du fichier."""
    out, in_tab = [], False
    for i, raw in enumerate(open(path, encoding="utf-8"), 1):
        if "\\begin{tabular" in raw:
            in_tab = True
        l = "" if raw.lstrip().startswith("%") else re.split(r"(?<!\\)%", raw)[0]
        l = l.replace("\\%", PCT)
        l = l.replace("---", " — ").replace("--", " – ")  # tirets LaTeX : « 2--32 » n'est pas un nombre négatif
        for pat in TEX_STRIP:
            l = re.sub(pat, " ", l)
        for m in TEX_NUM.finditer(l):
            s = m.group(0).strip()
            if re.fullmatch(r"(19|20)\d\d", s):  # année
                continue
            is_pct = re.match(r"^\s*\$?\s*(?:~|\\,|\;|\s)*" + PCT, l[m.end():]) is not None
            pre = l[max(0, m.start() - 60):m.start()].replace(PCT, "\\%").split()[-8:]
            post = l[m.end():m.end() + 60].replace(PCT, "\\%").split()[:8]
            out.append((i, s, " ".join(pre) + " ⟦" + s + "⟧ " + " ".join(post), is_pct, in_tab))
        if "\\end{tabular" in raw:
            in_tab = False
    return out


def tex_value(s):
    """Valeur, tolérance (½ unité de la dernière décimale, à l'échelle de l'exposant), nombre de décimales, notation scientifique ?"""
    t = s.replace("−", "-").replace(",", ".").replace(" ", "").replace("{", "").replace("}", "")
    m = re.match(r"^([-+]?)(?:(\d+(?:\.\d+)?)(?:\\times|\\cdot|×))?10\^([-+]?\d+)$", t)
    if m:
        man = m.group(2) or "1"
        e = int(m.group(3))
        d = len(man.split(".")[1]) if "." in man else 0
        v = float(man) * 10.0 ** e * (-1 if m.group(1) == "-" else 1)
        return v, 0.5 * 10.0 ** (e - d), d, True
    m = re.match(r"^([-+]?\d+(?:\.\d+)?)(?:e([-+]?\d+))?$", t)
    if not m:
        return None, 0.0, 0, False
    e = int(m.group(2)) if m.group(2) else 0
    d = len(m.group(1).split(".")[1]) if "." in m.group(1) else 0
    return float(m.group(1)) * 10.0 ** e, 0.5 * 10.0 ** (e - d), d, bool(m.group(2))


def cmd_tex(args=None):
    tex = Path(getattr(args, "tex", None) or TEX).expanduser().resolve()
    if not tex.is_file():
        sys.exit(f"[tex] chapitre introuvable : {tex} (dépôt du mémoire attendu à côté de celui-ci ; sinon $MSC_THESIS ou --tex CHEMIN)")
    rows, _, _ = build_main_rows()
    _, crow, _, _ = complements()
    index = []  # (valeur, étiquette de ligne, colonne)
    for n, r in enumerate(rows, 1):
        for col in ("v1", "final", "v2"):
            for v in parse_numbers(r[col]):
                index.append((v, f"table l.{n}", col))
    for r in crow:
        for col in ("v1", "final", "v2"):
            for v in parse_numbers(r[col]):
                index.append((v, "compl. " + r["g"], col))
    nums = tex_numbers(tex)
    cited = set()
    out_rows, n_ok, n_short, unmatched = [], 0, 0, []
    for line, s, ctx, is_pct, in_tab in nums:
        v, tol, d, sci = tex_value(s)
        if v is not None and d == 0 and not sci and abs(v) < 100 and not in_tab:
            n_short += 1
            out_rows.append([line, s + (" %" if is_pct else ""), ctx, "non comparé (entier court < 100, hors tabular)"])
            continue
        hit = []
        if v is not None:
            rtol = 1e-6 * abs(v)
            for val, lab, col in index:
                if abs(val - v) <= tol + 1e-15 or abs(val - v) <= rtol:
                    hit.append((lab, col))
                elif is_pct and (abs(val * 100 - v) <= tol + 1e-15):
                    hit.append((lab, col + " (×100)"))
                elif v != 0 and (abs(val - v * 1e3) <= 1e-6 * abs(v * 1e3) or abs(val - v * 1e-3) <= 1e-6 * abs(v * 1e-3)):
                    hit.append((lab, col + " (×10³)"))
        hit = sorted(set(hit), key=lambda h: (not h[0].startswith("compl."), h))
        if hit:
            n_ok += 1
            cited.update(h[0] for h in hit)
        else:
            unmatched.append((line, s + (" %" if is_pct else ""), ctx))
        out_rows.append([line, s + (" %" if is_pct else ""), ctx, (" ; ".join(f"{lab} ({col})" for lab, col in hit[:6]) + (f" … (+{len(hit) - 6})" if len(hit) > 6 else "")) if hit else "non apparié"])
    never = [n for n in range(1, len(rows) + 1) if f"table l.{n}" not in cited]
    n_un = len(unmatched)
    txt = [f"# Nombres cités dans `défauts.tex` du mémoire (lecture seule ; {len(nums)} nombres)", "",
           f"Date : {date.today().isoformat()} ; HEAD `{head()}` ; chapitre lu dans le dépôt du mémoire, `{'/'.join(tex.parts[-4:])}` (HEAD `{head(tex.parent)}`), md5 {md5(tex)}. Généré par `ch4_chiffres.py tex` (appariement refait en partie 5.5). Exclusions : commentaires, \\label, \\ref, \\eqref, \\cite, "
           "\\includegraphics, \\multicolumn, années (19xx/20xx), indices et exposants (précédés de ^ ou _). "
           "Appariement : le nombre du texte (virgule = point décimal) est comparé à tous les nombres des colonnes v1, final, non aligné de la table principale de `table_v1_final.md` (« table l.N » = rang dans cette table) "
           "et de sa section « Compléments » (« compl. 5.x … ») ; tolérance ½ unité de la dernière décimale du texte. "
           "(i) Un nombre suivi de « \\% » est aussi comparé aux fractions de la table × 100, marqué « (×100) ». "
           "(ii) La notation scientifique est lue en entier (mantisse × 10^exposant, « 10^{-12} » = 1e-12) et la tolérance suit l'exposant. "
           "(iii) Les entiers sans décimale inférieurs à 100 sont comparés quand ils sont dans un environnement `tabular`, et non comparés ailleurs. "
           "Les correspondances à un facteur 10³ (meV ↔ eV) sont marquées « (×10³) ». Un nombre peut s'apparier à plusieurs lignes par coïncidence (N, N mod 3, R_cut, 0.5, 3.10…) : correspondances mécaniques, pas des identifications ; "
           "les lignes des compléments sont listées en premier.", "",
           f"Comptes : appariés {n_ok} ; non appariés {n_un} ; entiers courts hors tabular non comparés {n_short} ; lignes de la table principale jamais citées {len(never)} / {len(rows)}.", "",
           md_table(["ligne tex", "nombre", "contexte (± 8 mots)", "ligne(s) appariée(s)"], out_rows), "",
           f"## Non appariés restants ({n_un})", "", md_table(["ligne tex", "nombre", "contexte (± 8 mots)"], [list(u) for u in unmatched]), "",
           "## Lignes de la table principale jamais citées", "", ", ".join(f"l.{n} ({rows[n - 1]['g'][:60]})" for n in never)]
    (HERE / "defauts_nombres.md").write_text("\n".join(txt), encoding="utf-8")
    print(f"[tex] {len(nums)} nombres ; appariés {n_ok} ; non appariés {n_un} ; entiers courts hors tabular non comparés {n_short} ; lignes jamais citées {len(never)}/{len(rows)}")



def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["check", "table", "regions", "notes", "tex", "all"])
    ap.add_argument("--tex", help="chemin de défauts.tex (défaut : dépôt du mémoire voisin, ou $MSC_THESIS)")
    a = ap.parse_args()
    cmd_check(a, quiet=(a.cmd != "check"))
    if a.cmd in ("table", "all"):
        cmd_table(a)
    if a.cmd in ("regions", "all"):
        cmd_regions(a)
    if a.cmd in ("tex", "all"):
        cmd_tex(a)
    if a.cmd in ("notes", "all"):
        cmd_notes(a)


if __name__ == "__main__":
    main()
