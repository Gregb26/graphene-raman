#!/usr/bin/env python
"""
Figure du §2.5 (couplage électron-photon, série EM), style figures/memoire.mplstyle, français. Données :
campagnes/EM/EM3/ (em3_ring_2p33.npz, em2_postw90_sigma_ti.npz : make_em3_data.py), campagnes/EM/EM2/em2_A.npz (DFT directe)
et la production M4 (campagnes/EM/M4_sigma/ : em_sigma_*_N1200_eta0.04.npz, em_scalars.json).
  fig_em_coupling : (a) |ħv_cv| dans le plan sur l'anneau de 2.33 eV (532 nm) en fonction de θ, trois variantes de la
                    vitesse (complète, centres seuls, sans Berry) et les 48 points de la DFT directe (EM2) ;
                    (b) σ_xx(ω)/σ₀ des trois variantes, postw90 avec transl_inv (EM2, décision 7 d'EM.md), ħω_froz,
                    médaillon sur les trois lasers. σ_yy sans Berry = σ_yy centres seuls à 1e-14 (τ_B − τ_A selon x).
La carte autour de K (F19, M4_sigma/em_map_K.npz) reste calculée mais n'est plus montrée (décision 10).
Usage : python scripts/fig/make_figures_em.py [--outdir figures/electron_photon]
"""
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]            # scripts/fig/ -> racine
plt.style.use(ROOT / "figures" / "memoire.mplstyle")
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # scripts/ : _palette, _bands, _paths
from _palette import NAVY, ORANGE, GREEN, GOLD, REF, INK, MUTED

ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default=str(ROOT / "figures" / "electron_photon")); a = ap.parse_args()
EM = ROOT / "campagnes" / "EM"
LASERS = {1.96: "633", 2.33: "532", 2.54: "488"}        # eV : nm
VARIANTS = ("full", "centres_only", "no_berry")
# une couleur et un trait par variante de la vitesse, les mêmes dans les deux panneaux ; références en marqueurs gris ouverts
# (DFT : cercles, postw90 : carrés), posés sur le mode complet qu'ils valident
STYLE = {"full": dict(color=NAVY, lw=1.5, ls="-", label="complet"),
         "centres_only": dict(color=GREEN, lw=1.2, ls="-.", label="centres seuls"),
         "no_berry": dict(color=ORANGE, lw=1.3, ls=":", label="sans connexion de Berry")}
C_LASER = GOLD


def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())


def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    w, h = fig.get_size_inches(); plt.close(fig); print(f"écrit {a.outdir}/{name}.pdf/.png ({w:.2f} × {h:.2f} po)")


scal = json.loads((EM / "M4_sigma" / "em_scalars.json").read_text())
hw_froz = scal["hw_froz_eV"]
RING = np.load(EM / "EM3" / "em3_ring_2p33.npz")
A = np.load(EM / "EM2" / "em2_A.npz"); dft = (A["sets"] == "ring") & np.isclose(A["hw"], float(RING["hw"]))
SIG = {v: np.load(EM / "M4_sigma" / f"em_sigma_{v}_N1200_eta0.04.npz") for v in VARIANTS}
P = np.load(EM / "EM3" / "em2_postw90_sigma_ti.npz")
hw = SIG["full"]["hw"]; s = {v: z["sigma"][:, 0, 0] for v, z in SIG.items()}

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.3))

# ---------------- (a) |ħv_cv| sur l'anneau de 2.33 eV
th = np.degrees(RING["theta"])
for v in VARIANTS:
    a1.plot(th, RING[f"hv_{v}"], **STYLE[v])
a1.plot(np.degrees(A["theta"][dft]), A["hv_qe"][dft], ls="none", marker="o", ms=3.5, mfc="none", mew=0.8, color=REF,
        label="DFT", zorder=3)
hv_F = float(RING["hv_F"])
a1.axhline(hv_F, color=MUTED, lw=0.6, zorder=0); a1.text(270, hv_F + 0.08, r"$\hbar v_F$", fontsize=8, color=INK, ha="center", va="bottom")
a1.text(0.03, 0.97, rf"$\hbar\omega$ = {float(RING['hw']):.2f} eV ({LASERS[2.33]} nm)", transform=a1.transAxes, fontsize=8,
        color=INK, ha="left", va="top")
a1.set_xlim(0, 360); a1.set_xticks(range(0, 361, 60)); a1.set_ylim(2, 10)
a1.set_xlabel(r"Angle $\theta$ sur l'anneau ($^\circ$)"); a1.set_ylabel(r"$|\hbar\mathbf{v}_{cv}|$ (eV$\,$\AA)")
panel(a1, "a")

# ---------------- (b) conductivité optique
for v in VARIANTS:
    a2.plot(hw, s[v], **{**STYLE[v], "label": None})
REF_P = dict(ls="none", marker="s", mfc="none", mew=0.7, color=REF, zorder=3)
a2.plot(P["hw"][5::10], P["sigma"][5::10, 0, 0], ms=2.8, label="postw90", **REF_P)   # tous les 0.1 eV, pic de van Hove (4.05) compris
a2.axhline(1, color=MUTED, lw=0.6, zorder=0)
a2.axvline(hw_froz, color=MUTED, lw=0.8, ls="-.", zorder=0)
a2.text(hw_froz + 0.07, 3.0, "limite de la fenêtre gelée", rotation=90, fontsize=7, color=INK, ha="left", va="bottom")
for e in LASERS: a2.axvline(e, ymax=0.04, color=C_LASER, lw=1.0)
a2.set_xlim(0, 6); a2.set_ylim(0, 8.5)
a2.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)"); a2.set_ylabel(r"$\sigma_{xx}(\omega)/\sigma_0$")

# médaillon : les trois lasers
ins = a2.inset_axes([0.14, 0.52, 0.42, 0.40])
for v in VARIANTS:
    ins.plot(hw, s[v], **{**STYLE[v], "lw": STYLE[v]["lw"] - 0.2, "label": None})
ins.plot(P["hw"][::5], P["sigma"][::5, 0, 0], ms=2.8, **REF_P)
for e, nm in LASERS.items():
    ins.axvline(e, color=C_LASER, lw=0.8, zorder=0); last = e == max(LASERS)   # « nm » une seule fois, sur la dernière, alignée à gauche pour ne pas toucher la précédente
    ins.text(e - (0.02 if last else 0), 1.81, nm + (" nm" if last else ""), fontsize=6, ha="left" if last else "center", va="bottom", color=INK, clip_on=False)
ins.set_xlim(1.8, 2.7); ins.set_ylim(1.1, 1.8); ins.tick_params(labelsize=7)
panel(a2, "b")

h1, l1 = a1.get_legend_handles_labels(); h2, l2 = a2.get_legend_handles_labels()
fig.legend(h1 + h2, l1 + l2, loc="lower center", ncol=len(l1 + l2), fontsize=8, frameon=False, handlelength=2.4, columnspacing=1.4)
fig.tight_layout(rect=(0, 0.07, 1, 1)); save(fig, "fig_em_coupling")
