"""Figure du chapitre 2 (structure électronique) : parties locale et non locale du pseudo-potentiel ONCV du carbone
(fig:KB du mémoire, ex `KB_projectors_C.pdf`, ex `figures/electron/fig_kb_pseudo_C.pdf` dessinée à la main en juin 2026).
(a) V_PPL(r) contre le potentiel de Coulomb −Z_val/r (Z_val = 4) ; (b) partie radiale des projecteurs de Kleinman–Bylander
β_{iℓ}(r), tracée telle que le UPF la tabule (PP_BETA = r·β(r), nulle à l'origine), comme dans la figure d'origine du mémoire ;
les projecteurs s'annulent au-delà de r_c ≈ 1,2–1,3 bohr.
Données : le C.upf de la maille (`--upf`, défaut : scripts/validation/_paths.upf()). Style et palette : graphene_raman.plotting.
Usage : python scripts/fig/make_figures_electron.py [--upf C.upf] [--outdir figures/electron]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_raman.io.pseudo_io import read_upf
from graphene_raman.plotting.palette import use_style, save, NAVY, ORANGE, SKY, GOLD, REF, MUTED

ROOT = Path(__file__).resolve().parents[2]            # scripts/fig/ -> racine
sys.path.insert(0, str(ROOT / "scripts" / "validation"))   # _paths.py (chemins des .save locaux, EDI_DATA)
import _paths  # noqa: E402

Z_VAL = 4.0                                            # électrons de valence du carbone (2s2 2p2)
R_MAX = 6.0                                            # bohr, fenêtre tracée en (a)
R_MAX_B = 2.0                                          # bohr, fenêtre tracée en (b) (les projecteurs sont nuls au-delà de r_c ~ 1,3)
L_NAME = {0: "s", 1: "p", 2: "d"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--upf", default=str(ROOT / _paths.upf()))
    ap.add_argument("--outdir", default=str(ROOT / "figures" / "electron"))
    a = ap.parse_args()
    use_style()
    ekb_li, fr_li, r, lmax, imax, V_L = read_upf(a.upf)          # Hartree, bohr ; fr = r * beta
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.5, 2.9))
    # (a) partie locale contre Coulomb
    sel = (r > 0) & (r <= R_MAX)
    ax1.plot(r[sel], V_L[sel], color=NAVY, label=r"$V_\mathrm{PPL}(r)$")
    ax1.plot(r[sel], -Z_VAL / r[sel], color=REF, ls="--", label=r"$-Z_\mathrm{val}/r$")
    ax1.set_ylim(-8, -0.5)
    ax1.set_xlim(0, R_MAX); ax1.set_xlabel(r"Rayon $r$ (bohr)"); ax1.set_ylabel(r"Potentiel local (Ha)")
    ax1.legend(loc="lower right"); ax1.set_title("(a)", loc="left")
    # (b) projecteurs radiaux beta_{i l}(r)
    colors = {(0, 0): NAVY, (0, 1): SKY, (1, 0): ORANGE, (1, 1): GOLD}
    selb = r <= R_MAX_B
    for l in range(lmax + 1):
        for i in range(imax):
            if not np.any(fr_li[l, i]):
                continue
            ax2.plot(r[selb], fr_li[l, i][selb], color=colors.get((l, i), MUTED), label=rf"$\beta_{{{i + 1}{L_NAME.get(l, l)}}}(r)$")
    ax2.set_xlim(0, R_MAX_B); ax2.set_xlabel(r"Rayon $r$ (bohr)"); ax2.set_ylabel(r"$r\,\beta_{i\ell}(r)$ (bohr$^{-1/2}$)")
    ax2.legend(loc="center right", fontsize=8); ax2.set_title("(b)", loc="left")
    fig.tight_layout()
    save(fig, "fig_kb_pseudo_C", a.outdir)
    print(f"écrit {a.outdir}/fig_kb_pseudo_C.pdf/.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
