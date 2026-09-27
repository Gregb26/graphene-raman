#!/usr/bin/env python
"""R6 3.6 : ajoute à <results_dir>/M_tests_summary.csv (écrit par analyze_M.py) les lignes de la porte A.2 (etape1/gate/gate_*.json, J3b) ;
idempotent (remplace les lignes existantes 'porte A.2')."""
import os, sys, csv, glob, json
WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src")); os.chdir(PROJ)
from electron_defect_interaction.config import load_production, results_dir
RES = results_dir(load_production(verbose=False)); p = os.path.join(RES, "M_tests_summary.csv")
rows = list(csv.DictReader(open(p))); rows = [r for r in rows if not r["test"].startswith("porte A.2")]
for r in rows:   # sbatch --export coupe GOLDEN_RESULT à la virgule : la source du test d'or est rétablie ici (job r6golden 21852238, valeur lue 1.80e-13)
    if r["test"].startswith("test d'or") and "non défini" in r["source"]:
        r["source"] = "r6golden_21852238 (submit_golden_dense.sh 5x5, M2)"
# les json de la porte sont dans <R6>/gate/ (J3b) ; etape1/gate/ n'existe que dans la copie article/ (1re version : chemin faux, 0 ligne ajoutée)
G = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(WORK, "gate", "gate_*.json"))) if "v1_" not in os.path.basename(f)]
if not G: sys.exit(f"[porte A.2] aucun gate_*.json dans {os.path.join(WORK, 'gate')} : rien ajouté")
for x in G: x["tag"] = x["size"] + (f" ({x['nb']} b)" if x["level"] == "coarse" and x["nb"] != 16 else "")
for lev, lab in (("coarse", "grossiers"), ("dense", "denses")):
    g = [x for x in G if x["level"] == lev]
    if not g: continue
    rows.append(dict(test=f"porte A.2 : ΔV appliqué aux états de Bloch purs contre M2/N_cells, {len(g)} fichiers {lab} (N = {', '.join(sorted({x['tag'] for x in g}, key=lambda s: (int(s.split('x')[0]), s)))})",
                     quantité="max|⟨nk|ΔV^L|n′k′⟩ − M2^L/N_cells|, idem NL (eV)", valeur=f"L {max(x['max_dL'] for x in g):.1e} ; NL {max(x['max_dNL'] for x in g):.1e}", seuil="1e-6 eV",
                     verdict="OK" if all(x["ok"] for x in g) else "REFUSÉ", source="gate_M_normalization.py (R6 J3b 21820491)"))
with open(p, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["test", "quantité", "valeur", "seuil", "verdict", "source"]); w.writeheader(); w.writerows(rows)
print(f"{p}: {len(rows)} lignes, dont {sum(r['test'].startswith('porte A.2') for r in rows)} « porte A.2 » ({len(G)} json)")
