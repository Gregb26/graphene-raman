"""Figure du chapitre 2 (phonons du graphène parfait), style graphene_raman.plotting, texte en français :

  fig_phfreq_phdos   (a) dispersion matdyn ω_νq sur Γ–K–M–Γ, (b) DOS de phonons matdyn sur grille dense ; cm⁻¹.
                     Même contenu que `fig_epw_phonons` (section 7 de make_figures_epw.py), extraite ici parce qu'elle ne relève que des
                     phonons (ex `make_fig_phonons.py`, rapatrié le 2026-09-30).

Données : <phonon_dir>/validation_<tag>.npz (ph_s, ph_F_matdyn) et <phonon_dir>/phdos_<tag>.npz, <phonon_dir> = config.phonon_dir()
(results/phonon ; copies des npz de results/epw de la chaîne de production 24k-24q_mv0.02).
Sortie : figures/electron_phonon/fig_phfreq_phdos.{pdf,png}, reproductible au bit (palette.save).
Usage : python scripts/fig/make_figures_phonons.py [--tag 24k24q_mv0.02] [--outdir figures/electron_phonon]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_raman.config import load_production, phonon_dir
from graphene_raman.plotting.palette import use_style, save, NAVY, MUTED

ROOT = Path(__file__).resolve().parents[2]            # scripts/fig/ -> racine
PHONON = phonon_dir(load_production(verbose=False))   # results/phonon (config)

# chemin Γ–K–M–Γ : longueurs |ΓK| = 2/3, |KM| = 1/3, |MΓ| = 1/√3 (unités 2π/a)
L = np.array([2 / 3, 1 / 3, 1 / np.sqrt(3)]); TICKS = np.concatenate([[0], np.cumsum(L)]) / L.sum(); TLAB = [r"$\Gamma$", "K", "M", r"$\Gamma$"]


def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())


def path_axis(ax):
    ax.set_xticks(TICKS); ax.set_xticklabels(TLAB); ax.set_xlim(0, 1)
    for t in TICKS[1:-1]:
        ax.axvline(t, color=MUTED, lw=0.6)


def fig_phfreq_phdos(a):
    V = np.load(os.path.join(PHONON, f"validation_{a.tag}.npz"), allow_pickle=True)
    P = np.load(os.path.join(PHONON, f"phdos_{a.tag}.npz"), allow_pickle=True)
    Fm = V["ph_F_matdyn"]; sq = V["ph_s"]
    fig, (d1, d2) = plt.subplots(1, 2, figsize=(6.5, 3.6), sharey=True, gridspec_kw={"width_ratios": [3, 1]})
    for m in range(Fm.shape[1]):
        d1.plot(sq, Fm[:, m], color=NAVY, lw=1.1)
    nk = P["nk_scf"]; nq = P["nq_ph"]; nd = P["nq_dos"]
    d1.set_ylabel(r"$\omega_{\nu\mathbf{q}}$ (cm$^{-1}$)"); d1.set_ylim(0, 1705); path_axis(d1)
    d1.set_title(rf"{nk[0]}$\times${nk[1]} $\mathbf{{k}}$, {nq[0]}$\times${nq[1]} $\mathbf{{q}}$, $\sigma_\mathrm{{MV}}$ = {float(P['sigma_mv_Ry']):g} Ry", loc="left", fontsize=9); panel(d1, "a")
    d2.plot(P["dos"], P["freq_cm"], color=NAVY, lw=1.0); d2.fill_betweenx(P["freq_cm"], 0, P["dos"], color=NAVY, alpha=0.12, lw=0)
    d2.set_xlim(0, None); d2.set_xticks([]); d2.set_xlabel(r"DOS (états/cm$^{-1}$)", fontsize=9, labelpad=6)
    d2.text(0.95, 0.66, rf"{nd[0]}$\times${nd[1]} $\mathbf{{q}}$" + "\n" + rf"$\Delta E$ = {float(P['deltaE_cm']):g} cm$^{{-1}}$" + "\n" + rf"$\sigma$ = {float(P['degauss_dos_cm']):g} cm$^{{-1}}$",
            transform=d2.transAxes, ha="right", va="center", fontsize=7); panel(d2, "b")
    d2.tick_params(axis="y", labelleft=False)
    fig.tight_layout(w_pad=0.4)
    save(fig, "fig_phfreq_phdos", a.outdir); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="24k24q_mv0.02", help="validation_<tag>.npz et phdos_<tag>.npz (chaîne de production)")
    ap.add_argument("--outdir", default=str(ROOT / "figures" / "electron_phonon"))
    a = ap.parse_args()
    use_style()
    fig_phfreq_phdos(a); print(f"écrit {a.outdir}/fig_phfreq_phdos.pdf/.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
