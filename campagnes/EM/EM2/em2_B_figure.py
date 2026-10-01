"""
EM2, control figure of B: (a) sigma_xx/sigma_0 of M4 (full, N = 1200, eta = 0.04 eV) and of postw90 3.1.0 kubo
(berry_kmesh 1201, default settings and transl_inv = true), g_s = 2 applied to postw90; (b) |difference| with M4 on a
log scale, split into its causes: prefactor 1/Delta (our code, em2_B_ours.py), postw90 default, postw90 transl_inv,
and what is left of postw90 (transl_inv, no ws_distance) against our sigma with the same grid and prefactor.
Style and palette of the thesis, French text. Writes em2_B_sigma.{pdf,png}.
"""

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
sys.path.insert(0, str(HERE))
from graphene_raman.plotting.palette import use_style; use_style()
from graphene_raman.plotting.palette import NAVY, ORANGE, GREEN, GOLD, PINK, REF
from em2_B_compare import read_postw90   # noqa: E402  (reads postw90/<dir>/wannier-kubo_S_*.dat, g_s = 2 applied)

m4 = np.load(REPO / "campagnes" / "EM" / "M4_sigma" / "em_sigma_full_N1200_eta0.04.npz")
hw, s_m4 = m4["hw"], m4["sigma"][:, 0, 0]
ours = np.load(HERE / "em2_B_ours.npz")
s_def = read_postw90("k1201")[1][:, 0, 0]
s_ti = read_postw90("k1201_ti")[1][:, 0, 0]
s_ti_nows = read_postw90("k1201_ti_nows")[1][:, 0, 0]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 2.9))
a1.plot(hw, s_m4, color=NAVY, label=r"M4 (ce travail)")
a1.plot(hw[::15], s_def[::15], ls="none", marker="o", ms=3, mfc="none", color=REF, label=r"\texttt{postw90}")
a1.plot(hw[7::15], s_ti[7::15], ls="none", marker="s", ms=2.5, mfc="none", color=GOLD,
        label=r"\texttt{postw90}, \texttt{transl\_inv}")
a1.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)")
a1.set_ylabel(r"$\sigma_{xx}/\sigma_0$")
a1.set_xlim(0, 6)
a1.legend(fontsize=7, loc="upper left")
a1.set_title(r"(a) conductivité optique", loc="left")

curves = [
    (np.abs(ours["shift_delta"][:, 0, 0] - s_m4), ORANGE, "-", r"préfacteur $1/\Delta\varepsilon$ seul"),
    (np.abs(s_def - s_m4), REF, "-", r"\texttt{postw90} $-$ M4"),
    (np.abs(s_ti - s_m4), GOLD, "-", r"\texttt{postw90}, \texttt{transl\_inv} $-$ M4"),
    (np.abs(s_ti_nows - ours["gam1201"][:, 0, 0]), GREEN, "-", r"reste (même grille, même préfacteur)"),
]
for y, c, ls, lab in curves:
    a2.semilogy(hw, np.maximum(y, 1e-10), color=c, ls=ls, lw=1.0, label=lab)
a2.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)")
a2.set_ylabel(r"$|\Delta\sigma_{xx}|/\sigma_0$")
a2.set_xlim(0, 6); a2.set_ylim(1e-10, 1e-1)
a2.legend(fontsize=6.5, loc="lower right")
a2.set_title(r"(b) écarts à M4", loc="left")

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(HERE / f"em2_B_sigma.{ext}")
print("wrote em2_B_sigma.{pdf,png}")
