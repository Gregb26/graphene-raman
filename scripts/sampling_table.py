"""Sampling table: N, N mod 3, dE_F = E_F(d)-E_F(p) (meV), KS dense reconstruction (max/mean meV), max|V_ed^L| boundary (meV, raw)."""
import csv, os, numpy as np, xml.etree.ElementTree as ET
HA = 27.211386245988; rows = []
K = np.load("results/M/ks_reconstruction.npz") if os.path.exists("results/M/ks_reconstruction.npz") else {}
V = np.load("results/M/ved_analysis.npz") if os.path.exists("results/M/ved_analysis.npz") else {}
for N in (5, 6, 7, 8, 9, 10, 11, 12):
    S = f"{N}x{N}"; ef = {}
    for k in "pd":
        f = f"data/graphene/supercell/qe/defect_{S}_{k}.save/data-file-schema.xml"
        if os.path.exists(f): ef[k] = float(ET.parse(f).getroot().find("output/band_structure/fermi_energy").text)
    dEf = (ef["d"] - ef["p"]) * HA * 1e3 if len(ef) == 2 else np.nan
    ksmax = float(K[f"{S}_dense_max_meV"]) if f"{S}_dense_max_meV" in K else np.nan; ksmean = float(K[f"{S}_dense_mean_meV"]) if f"{S}_dense_mean_meV" in K else np.nan
    b = float(V[f"{S}_b_all"]) * 1e3 if f"{S}_b_all" in V else np.nan
    rows.append([N, N % 3, f"{dEf:+.1f}", f"{ksmax:.3f}", f"{ksmean:.3f}", f"{b:.1f}"])
with open("results/M/sampling_table.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["N", "N_mod_3", "dE_F_meV", "KS_dense_max_meV", "KS_dense_mean_meV", "Ved_boundary_max_meV"]); w.writerows(rows)
print(f"{'N':>3} {'N%3':>4} {'dE_F(meV)':>10} {'KS max':>8} {'KS mean':>8} {'|V|bord':>8}")
for r in rows: print(f"{r[0]:>3} {r[1]:>4} {r[2]:>10} {r[3]:>8} {r[4]:>8} {r[5]:>8}")
