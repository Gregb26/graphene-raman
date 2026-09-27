#!/usr/bin/env python
"""R6 : <results_dir>/m_rcut_resigma.csv à partir de resigma_9x9_rc0123.npz et resigma_9x9_rc4.npz (rcut_resigma.py), mêmes colonnes et
mêmes définitions que results/M/m_rcut_resigma.csv (v1, commit a9855af, produit hors de scripts/). Contrôle : appliqué à results/M, le
script doit redonner le csv v1 (écart relatif max imprimé) ; il n'écrit results_dir qu'après ce contrôle.
Usage : python r6_m_rcut_resigma.py [--size 9x9] [--tol 1e-9]"""
import os, sys, csv, argparse
import numpy as np
WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src")); os.chdir(PROJ)
from electron_defect_interaction.config import load_production, results_dir
ap = argparse.ArgumentParser(); ap.add_argument("--size", default="9x9"); ap.add_argument("--tol", type=float, default=1e-9); a = ap.parse_args()
RES = results_dir(load_production(verbose=False)); FROZEN = "results/M"
COLS = ["size", "R_cut", "nL", "dim", "grid", "eta_eV", "window_eV", "E_D_eV", "med_ReSigma_meV", "mean_ReSigma_meV", "med_Gamma_meV",
        "mean_Gamma_meV", "med_absReSigma_over_Gamma", "med_abs_dReSigma_vs_rc4_meV", "max_abs_dReSigma_vs_rc4_meV", "med_abs_dGamma_vs_rc4_meV",
        "max_abs_dGamma_vs_rc4_meV", "rel_med_dReSigma", "rel_med_dGamma", "E_res_minus_E_D_eV", "n_states"]

def table(d):
    A = np.load(f"{d}/resigma_{a.size}_rc0123.npz"); B = np.load(f"{d}/resigma_{a.size}_rc4.npz")
    S = {rc: (A if rc < 4 else B)[f"Sigma_rc{rc}"] for rc in range(5)}; nL = {rc: int((A if rc < 4 else B)[f"nL_rc{rc}"]) for rc in range(5)}
    E = A["E_out"]; ED = float(A["E_D"]); m = ~np.isnan(S[3]); nw = S[3].shape[0]
    assert np.array_equal(m, ~np.isnan(S[4])) and np.allclose(E, B["E_out"], equal_nan=True) and float(B["E_D"]) == ED
    G4, R4 = -2 * S[4].imag[m], S[4].real[m]; rows = []
    for rc in range(5):
        G, R = -2 * S[rc].imag[m], S[rc].real[m]; mm = np.abs(E[m] - ED) <= 1.5
        rows.append(dict(size=a.size, R_cut=rc, nL=nL[rc], dim=nL[rc] * nw, grid=int(A["grid"]), eta_eV=float(A["eta"]), window_eV=float(A["e_window"]), E_D_eV=ED,
                         med_ReSigma_meV=np.median(R) * 1e3, mean_ReSigma_meV=R.mean() * 1e3, med_Gamma_meV=np.median(G) * 1e3, mean_Gamma_meV=G.mean() * 1e3,
                         med_absReSigma_over_Gamma=np.median(np.abs(R) / G),
                         med_abs_dReSigma_vs_rc4_meV=np.median(np.abs(R - R4)) * 1e3, max_abs_dReSigma_vs_rc4_meV=np.abs(R - R4).max() * 1e3,
                         med_abs_dGamma_vs_rc4_meV=np.median(np.abs(G - G4)) * 1e3, max_abs_dGamma_vs_rc4_meV=np.abs(G - G4).max() * 1e3,
                         rel_med_dReSigma=np.median(np.abs(R - R4)) / np.median(R4), rel_med_dGamma=np.median(np.abs(G - G4)) / np.median(G4),
                         E_res_minus_E_D_eV=float(E[m][mm][np.argmax(G[mm])] - ED), n_states=int(m.sum())))
    return rows

ref = list(csv.DictReader(open(f"{FROZEN}/m_rcut_resigma.csv"))); v1 = table(FROZEN); worst = 0.0
for r, x in zip(ref, v1):
    for c in COLS[1:]:
        u, w = float(r[c]), float(x[c]); worst = max(worst, abs(u - w) / max(abs(u), 1e-12))
print(f"[contrôle] recalcul de {FROZEN}/m_rcut_resigma.csv depuis les npz v1 : écart relatif max {worst:.2e} (seuil {a.tol:g})")
if worst > a.tol: sys.exit("[refus] définitions non reproduites : rien n'est écrit")
v2 = table(RES); out = f"{RES}/m_rcut_resigma.csv"
with open(out, "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=COLS); wr.writeheader(); wr.writerows(v2)
for x in v2:
    print(f"  R_cut {x['R_cut']} : med Γ {x['med_Gamma_meV']:.2f} meV, med Re Σ {x['med_ReSigma_meV']:.2f}, écart médian par état vs R_cut 4 : Γ {100*x['rel_med_dGamma']:.2f} %, "
          f"Re Σ {100*x['rel_med_dReSigma']:.2f} % ; |Re Σ|/Γ méd. {x['med_absReSigma_over_Gamma']:.3f} ; E_res {x['E_res_minus_E_D_eV']:+.3f}")
print(f"écrit {out}")
