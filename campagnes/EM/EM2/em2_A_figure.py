"""
EM2, control figure of A: (a) in-plane |hbar v_cv| on the 2.33 eV ring, QE (bands.x lp, 48 k) against the three
Wannier variants (dense ring, 720 angles); (b) relative deviation of the full velocity from QE, k by k, on the three
laser rings. Style and palette of the thesis (graphene_raman.plotting: memoire.mplstyle, palette.py); French text.
Reads em2_A.npz (em2_A_compare.py). Writes em2_A_anneaux.{pdf,png}.
"""

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
from graphene_raman.plotting.palette import use_style; use_style()
from graphene_raman.plotting.palette import NAVY, ORANGE, GREEN, GOLD, REF, SKY
from graphene_raman.electron_photon import centres_only, compute_velocity, make_grid_tb, make_wannier_tb, ring

E_D = -4.238895
A = np.load(HERE / "em2_A.npz")
tb = make_wannier_tb(REPO / "results" / "wannier" / "27x27" / "wannier_tb.dat")
K = make_grid_tb(tb, 3).K

k, th, _, _ = ring(tb, K, 2.33, E_D, ntheta=720)
dense = {}
for name, (model, mode) in {"full": (tb, "berry"), "centres_only": (centres_only(tb), "berry"),
                            "no_berry": (tb, "no_berry")}.items():
    _, eps, _, hv = compute_velocity(model, k, mode)
    n = np.sum(eps < E_D, axis=1); ii = np.arange(len(k))
    dense[name] = np.sqrt(np.abs(hv[ii, 0, n, n - 1])**2 + np.abs(hv[ii, 1, n, n - 1])**2)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.9))
deg = np.degrees(th)
a1.plot(deg, dense["full"], color=NAVY, label="Wannier complet")
a1.plot(deg, dense["centres_only"], color=GREEN, label="centres seuls")
a1.plot(deg, dense["no_berry"], color=ORANGE, ls="--", label="sans Berry")
m = (A["sets"] == "ring") & (A["hw"] == 2.33)
a1.plot(np.degrees(A["theta"][m]), A["hv_qe"][m], ls="none", marker="o", ms=3.5, mfc="none", color=REF,
        label=r"DFT (\texttt{bands.x}, $\hat p$)")
a1.set_xlabel(r"Angle $\theta$ autour de $K$ (degrés)")
a1.set_ylabel(r"$|\hbar v_{cv}|$ dans le plan (eV\,\AA)")
a1.set_xticks([0, 60, 120, 180, 240, 300, 360]); a1.set_xlim(0, 360); a1.set_ylim(2, 11.2)
a1.legend(fontsize=7, loc="upper center", ncol=2)
a1.set_title(r"(a) anneau à 2.33 eV", loc="left")

for w, c in zip([1.96, 2.33, 2.54], [SKY, NAVY, GOLD]):
    m = (A["sets"] == "ring") & (A["hw"] == w)
    rel = 100 * (A["hv_full"][m] / A["hv_qe"][m] - 1)
    a2.plot(np.degrees(A["theta"][m]), rel, marker="o", ms=3, color=c, label=f"{w:.2f} eV")
a2.axhline(0, color=REF, lw=0.6)
a2.set_xlabel(r"Angle $\theta$ autour de $K$ (degrés)")
a2.set_ylabel(r"Écart Wannier complet $-$ DFT (\%)")
a2.set_xticks([0, 60, 120, 180, 240, 300, 360]); a2.set_xlim(0, 360)
a2.legend(fontsize=7)
a2.set_title(r"(b) écart relatif sur $|\hbar v_{cv}|$", loc="left")

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"em2_A_anneaux.{ext}")
print("wrote em2_A_anneaux.{pdf,png}")
