"""Level 2 by family (N = 3m vs others): median Gamma*N_cells at the frozen parameters, dense on-site pz-pz of the vacancy
sublattice, Re M^L / Re M^NL at K. Writes results/M/level2_families.csv."""
import csv, os, numpy as np
from electron_defect_interaction.config import load_production
cfg = load_production(verbose=False); RC, G, E = cfg["R_cut"], cfg["grid"], cfg["eta_eV"]
L = np.load("results/M/mwr_locality.npz"); A = np.load("results/M/M_analysis.npz", allow_pickle=True)
rows = []
for S in ["6x6", "9x9", "12x12", "5x5", "7x7", "8x8"]:
    f = f"results/M/specwd_{S}_prod.npz"
    if not os.path.exists(f): continue
    r = np.load(f)["results"]; m = {(int(x[0]), int(x[1]), round(float(x[2]), 4)): float(x[3]) for x in r}
    g = m.get((RC, G, E), np.nan); N = int(S.split("x")[0])
    on = float(L[f"{S}_dense_onsite_pzvac"]) if f"{S}_dense_onsite_pzvac" in L else np.nan; sub = str(L[f"{S}_dense_vac_sublattice"]) if f"{S}_dense_vac_sublattice" in L else "?"
    reL = float(A[f"lnl_{S}_ReL"]) if f"lnl_{S}_ReL" in A else np.nan; reN = float(A[f"lnl_{S}_ReNL"]) if f"lnl_{S}_ReNL" in A else np.nan; dk = float(A[f"lnl_{S}_dK"]) if f"lnl_{S}_dK" in A else np.nan
    rows.append([S, N, "3m" if N % 3 == 0 else "non-3m", sub, f"{g:.2f}", f"{on:.3f}", f"{reL:+.4f}", f"{reN:+.4f}", f"{dk:.3f}"])
with open("results/M/level2_families.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["size", "N", "family", "vacancy_sublattice", "median_Gamma_Ncells_meV", "onsite_pz_vac_eV", "ReML_K_eV", "ReMNL_K_eV", "dk_to_K_Ainv"]); w.writerows(rows)
print(f"{'size':>6} {'fam':>7} {'sub':>3} {'G*Nc(meV)':>10} {'on-site':>8} {'ReM^L(K)':>9} {'ReM^NL(K)':>10} {'|dk|':>6}")
for r in rows: print(f"{r[0]:>6} {r[2]:>7} {r[3]:>3} {r[4]:>10} {r[5]:>8} {r[6]:>9} {r[7]:>10} {r[8]:>6}")
for fam in ("3m", "non-3m"):
    g = np.array([float(r[4]) for r in rows if r[2] == fam]); print(f"family {fam}: Gamma*N_cells mean {g.mean():.1f} meV, spread (max-min)/mean {100*(g.max()-g.min())/g.mean():.2f} %")
