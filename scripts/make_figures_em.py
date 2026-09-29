#!/usr/bin/env python
"""
Figure du §2.5 (couplage électron-photon, série EM), style figures/memoire.mplstyle, français. Données de la production M4
(memoire/EM/M4_sigma/ : em_map_K.npz, em_sigma_*_N1200_eta0.04.npz, em_scalars.json) et, s'il existe, σ(ω) de postw90
(memoire/EM/EM2/em2_postw90_sigma.npz : hw, sigma (nw, 3, 3) en σ/σ₀), ajouté au panneau (b).
  fig_em_coupling : (a) |ħv^x_cv|²/(ħv_F)² autour de K (vitesse complète) et anneaux résonants aux trois lasers ;
                    (b) σ_xx(ω)/σ₀ des trois variantes de la vitesse (complète, centres seuls, sans Berry), ħω_froz.
                    σ_yy sans Berry = σ_yy centres seuls à 1e-14 (τ_B − τ_A selon x) : une seule courbe.
Usage : python scripts/make_figures_em.py [--outdir figures]
"""
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
plt.style.use(ROOT / "figures" / "memoire.mplstyle")
from _palette import NAVY, ORANGE, GREEN, GOLD, REF, INK, MUTED, CMAP_SEQ

ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default=str(ROOT / "figures")); a = ap.parse_args()
M4 = ROOT / "memoire" / "EM" / "M4_sigma"
EM2 = ROOT / "memoire" / "EM" / "EM2" / "em2_postw90_sigma.npz"
LASERS = {1.96: "633", 2.33: "532", 2.54: "488"}        # eV : nm
# couleurs : une par variante de la vitesse ; les lasers (anneaux en (a), repères en (b)) forment une catégorie à part
C_FULL, C_CENTRES, C_NOBERRY, C_LASER, C_POSTW90 = NAVY, GREEN, ORANGE, GOLD, REF


def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())


def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    w, h = fig.get_size_inches(); plt.close(fig); print(f"écrit {a.outdir}/{name}.pdf/.png ({w:.2f} × {h:.2f} po)")


scal = json.loads((M4 / "em_scalars.json").read_text())
hv_F, hw_froz = scal["hv_F_eV_A"]["inter"], scal["hw_froz_eV"]
MAP = np.load(M4 / "em_map_K.npz")
SIG = {v: np.load(M4 / f"em_sigma_{v}_N1200_eta0.04.npz") for v in ("full", "centres_only", "no_berry")}
hw = SIG["full"]["hw"]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.45), gridspec_kw={"width_ratios": [1, 1.25]})

# ---------------- (a) carte autour de K
q = MAP["q"]; Z = MAP["P_full"][..., 0] / hv_F**2
im = a1.pcolormesh(q, q, Z, cmap=CMAP_SEQ, shading="auto", rasterized=True, vmin=0)
cs = a1.contour(q, q, MAP["deps"], levels=sorted(LASERS), colors=C_LASER, linewidths=1.0)
a1.text(0.5, 0.04, r"anneaux : $\hbar\omega$ = 1.96, 2.33, 2.54 eV", transform=a1.transAxes, ha="center", va="bottom",
        fontsize=7, color=INK, bbox=dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85))
a1.axhline(0, color=MUTED, lw=0.6, ls=":")
a1.plot(0, 0, marker="+", color=INK, ms=6, mew=0.9); a1.text(0.012, 0.012, "K", fontsize=8, color=INK)
a1.set_aspect("equal"); a1.set_xlim(q[0], q[-1]); a1.set_ylim(q[0], q[-1])
a1.set_xlabel(r"$q_x$ (\AA$^{-1}$)"); a1.set_ylabel(r"$q_y$ (\AA$^{-1}$)")
cb = fig.colorbar(im, ax=a1, fraction=0.046, pad=0.03); cb.set_label(r"$|\hbar v^x_{cv}|^2/(\hbar v_F)^2$", fontsize=9)
cb.ax.tick_params(labelsize=8)
panel(a1, "a")

# ---------------- (b) conductivité optique
s = {v: z["sigma"] for v, z in SIG.items()}
a2.plot(hw, s["full"][:, 0, 0], color=C_FULL, lw=1.5, label="complet")
a2.plot(hw, s["centres_only"][:, 0, 0], color=C_CENTRES, lw=1.2, ls="-.", label=r"centres seuls ($=$ sans Berry, $yy$)")
a2.plot(hw, s["no_berry"][:, 0, 0], color=C_NOBERRY, lw=1.2, ls="--", label=r"sans Berry, $xx$")
if EM2.exists():
    P = np.load(EM2); a2.plot(P["hw"][::10], P["sigma"][::10, 0, 0], ls="none", marker="o", ms=3, mfc="none", color=C_POSTW90, label="postw90")
for e in LASERS: a2.axvline(e, color=C_LASER, lw=0.8, zorder=0)
a2.axvline(hw_froz, color=MUTED, lw=0.8, ls="-.", zorder=0)
a2.text(hw_froz + 0.06, 7.3, r"$\hbar\omega_\mathrm{froz}$", fontsize=8, color=INK, va="top")
a2.axhline(1, color=MUTED, lw=0.6, ls=":", zorder=0)
a2.set_xlim(0, 6); a2.set_ylim(0, 7.5)
a2.set_xlabel(r"Énergie du photon $\hbar\omega$ (eV)"); a2.set_ylabel(r"$\sigma(\omega)/\sigma_0$")

# médaillon : les trois lasers
ins = a2.inset_axes([0.14, 0.52, 0.42, 0.40])
for v, kw in (("full", dict(color=C_FULL, lw=1.3)), ("centres_only", dict(color=C_CENTRES, lw=1.0, ls="-.")),
              ("no_berry", dict(color=C_NOBERRY, lw=1.0, ls="--"))):
    ins.plot(hw, s[v][:, 0, 0], **kw)
for e, nm in LASERS.items():
    ins.axvline(e, color=C_LASER, lw=0.8, zorder=0); last = e == max(LASERS)   # « nm » une seule fois, sur la dernière, alignée à gauche pour ne pas toucher la précédente
    ins.text(e - (0.02 if last else 0), 1.81, nm + (" nm" if last else ""), fontsize=6, ha="left" if last else "center", va="bottom", color=INK, clip_on=False)
ins.set_xlim(1.8, 2.7); ins.set_ylim(1.1, 1.8); ins.tick_params(labelsize=7)
panel(a2, "b")

h, l = a2.get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=len(l), fontsize=8, frameon=False, handlelength=2.2)
fig.tight_layout(rect=(0, 0.07, 1, 1)); save(fig, "fig_em_coupling")
