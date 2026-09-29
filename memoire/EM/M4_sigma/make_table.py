"""
Tableau d'EM3 (chiffres du §2.5) à partir des sorties de m4_prod.py, sans recopie à la main :
  em_table.csv : format long, une ligne par (laser, variante) : σ_xx, σ_yy (σ/σ₀), ⟨|e_x,y·ħv_cv|²⟩ sur l'anneau
                 (unités (ħv_F)²/2), nœuds δ_x, δ_y (degrés), rapport min–max de |ħv_cv| au mode complet ;
  em_table.tex : `tabular` (booktabs) prêt à envelopper dans un `table` du mémoire, point décimal.
Relancer : .venv/bin/python memoire/EM/M4_sigma/make_table.py
"""

import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
scal = json.loads((HERE / "em_scalars.json").read_text())
stats = json.loads((HERE / "em_ring_stats.json").read_text())
sig = scal["sigma_over_sigma0_at_lasers"]   # variante -> laser -> {xx, yy, xy}
LASERS = {"1.96": 633, "2.33": 532, "2.54": 488}   # eV : nm
VARIANTS = ("full", "centres_only", "no_berry")

# CSV, format long
with open(HERE / "em_table.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["hw_eV", "lambda_nm", "variant", "sigma_xx", "sigma_yy", "avg_x", "avg_y", "node_x_deg", "node_y_deg",
                "ratio_min", "ratio_max"])
    for e, nm in LASERS.items():
        for v in VARIANTS:
            r = stats[e][v]
            w.writerow([e, nm, v, f"{sig[v][e]['xx']:.5f}", f"{sig[v][e]['yy']:.5f}", f"{r['avg'][0]:.5f}", f"{r['avg'][1]:.5f}",
                        f"{r['node'][0]:.1f}", f"{r['node'][1]:.1f}", f"{r['ratio'][0]:.4f}", f"{r['ratio'][1]:.4f}"])

# LaTeX : une ligne par laser
f4 = lambda x: f"{x:.4f}"
deg = lambda x: f"${x:+.1f}^\\circ$" if abs(x) > 1e-9 else "$0$"
rows = []
for e, nm in LASERS.items():
    s, r = {v: sig[v][e] for v in VARIANTS}, stats[e]
    rows.append(" & ".join([str(nm), e, f4(s["full"]["xx"]), f4(s["centres_only"]["xx"]), f4(s["no_berry"]["xx"]), f4(s["no_berry"]["yy"]),
                            f4(r["full"]["avg"][0]), f4(r["centres_only"]["avg"][0]), deg(r["full"]["node"][0]), deg(r["no_berry"]["node"][0])])
                + r" \\")
tex = "\n".join([
    r"% généré par memoire/EM/M4_sigma/make_table.py (ne pas éditer à la main) ; nécessite booktabs",
    r"{\setlength{\tabcolsep}{4pt}%",
    r"\begin{tabular}{rr cccc cc cc}",
    r"\toprule",
    r" & & \multicolumn{4}{c}{$\sigma/\sigma_0$} & \multicolumn{2}{c}{$\langle|\hbar v^x_{cv}|^2\rangle_\theta$ ($(\hbar v_F)^2/2$)} & \multicolumn{2}{c}{nœud $\delta_x$} \\",
    r"\cmidrule(lr){3-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
    r"$\lambda$ (nm) & $\hbar\omega$ (eV) & complet & centres & sans Berry, $xx$ & sans Berry, $yy$ & complet & centres & complet & sans Berry \\",
    r"\midrule",
    *rows,
    r"\bottomrule",
    r"\end{tabular}}",
    (f"% ħv_F = {scal['hv_F_eV_A']['inter']:.3f} eV·Å (v_F = {scal['v_F_m_s'] / 1e5:.2f}×10⁵ m/s) ; ħω_froz = {scal['hw_froz_eV']:.2f} eV ; "
     f"N = {scal['N']}, η = {scal['eta_eV']} eV ; ⟨|ħv^x_cv|²⟩ en unités de (ħv_F)²/2 ; commit {scal['provenance']['commit']}"),
])
(HERE / "em_table.tex").write_text(tex + "\n")
print(tex)
