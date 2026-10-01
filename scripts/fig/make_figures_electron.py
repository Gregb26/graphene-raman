"""Figures du chapitre 2 (structure électronique du graphène parfait), style graphene_raman.plotting, texte en français :

  fig_kb_pseudo_C            parties locale et non locale du pseudo-potentiel ONCV du carbone (fig:KB ; redessinée le 2026-09-30) :
                             (a) V_PPL(r) contre −Z_val/r (Z_val = 4) ; (b) projecteurs de Kleinman–Bylander tels que le UPF les tabule
                             (PP_BETA = r·β(r), nuls à l'origine, nuls au-delà de r_c ≈ 1,2–1,3 bohr). Données : <electron_dir>/C.upf (copie du pseudo du `.save` de la maille ; `--upf`).
  fig_electron_convergence   tests de convergence SCF (ex `qe_pp/plot_convergence.py`, fonction plot_all, rapatriée le 2026-09-30) :
                             (a) E_cut à grille modérée, (b) grille k, (c) E_cut à grille optimale, (d) élargissement par type de smearing ;
                             ΔE = (E − E_réf)/N_at en meV/atome, seuil 0,1 meV/atome. Données : <electron_dir>/{ecut,kpoint,ecut_cross,smearing}.dat.
  fig_ebands_edos            (ex `qe_pp/plot_ebands_edos.py`) : (a) structure de bandes sur Γ–K–M–Γ, (b) DOS électronique, énergies
                             relatives à E_F. Données : <electron_dir>/{ebands,edos}.dat (bands.x, dos.x).

<electron_dir> = config.electron_dir() (results/electron). Sortie : figures/electron/<nom>.{pdf,png}, reproductible au bit (palette.save).
Usage : python scripts/fig/make_figures_electron.py [--only kb_pseudo_C,electron_convergence,ebands_edos] [--upf C.upf] [--outdir figures/electron]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from graphene_raman.config import load_production, electron_dir
from graphene_raman.io.pseudo_io import read_upf
from graphene_raman.plotting.palette import use_style, save, NAVY, ORANGE, SKY, GOLD, GREEN, PINK, REF, MUTED, INK

ROOT = Path(__file__).resolve().parents[2]            # scripts/fig/ -> racine

Z_VAL = 4.0                                            # électrons de valence du carbone (2s2 2p2)
R_MAX = 6.0                                            # bohr, fenêtre tracée en (a)
R_MAX_B = 2.0                                          # bohr, fenêtre tracée en (b) (les projecteurs sont nuls au-delà de r_c ~ 1,3)
L_NAME = {0: "s", 1: "p", 2: "d"}
ELECTRON = electron_dir(load_production(verbose=False))   # results/electron (config)

# ---- convergence SCF (ex plot_convergence.py de qe_pp, contenu identique)
RY_TO_MEV = 13605.693_122_994
NAT = 2                                                # atomes par maille
THRESH = 0.1                                           # seuil de convergence, meV/atome
SMEARING_LABELS = {"gaussian": "gaussienne", "mp": "Methfessel-Paxton", "mv": "Marzari-Vanderbilt", "fd": "Fermi-Dirac"}
SMEARING_COLORS = {"gaussian": GREEN, "mp": PINK, "mv": GOLD, "fd": SKY}
SMEARING_MARKERS = {"gaussian": "o", "mp": "s", "mv": "^", "fd": "D"}


def fig_kb_pseudo_C(a):
    """(a) V_PPL(r) contre −Z_val/r ; (b) r·β_{iℓ}(r) des projecteurs KB, tels que tabulés dans le UPF."""
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




def load_ecut(fname):
    """ecut.dat ou ecut_cross.dat : ΔE par atome (meV) référencé à la dernière ligne (la plus convergée), qui est retirée."""
    data = np.loadtxt(os.path.join(ELECTRON, fname), skiprows=1)
    ecut, etot = data[:, 0], data[:, 1]
    dE = (etot - etot[-1]) / NAT * RY_TO_MEV
    return ecut[:-1], dE[:-1]


def load_kpoint(fname="kpoint.dat"):
    data = np.loadtxt(os.path.join(ELECTRON, fname), skiprows=1)
    K, etot = data[:, 0], data[:, 1]
    dE = (etot - etot[-1]) / NAT * RY_TO_MEV
    return K[:-1], dE[:-1]


def load_smearing(fname="smearing.dat"):
    """{type : (degauss, ΔE)} ; ΔE référencé au plus petit degauss de chaque type (retiré)."""
    raw = {}
    with open(os.path.join(ELECTRON, fname)) as f:
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            raw.setdefault(parts[0], []).append((float(parts[1]), float(parts[2])))
    out = {}
    for sm, vals in raw.items():
        vals = sorted(vals); dg = np.array([v[0] for v in vals]); e = np.array([v[1] for v in vals])
        out[sm] = (dg[1:], (e[1:] - e[0]) / NAT * RY_TO_MEV)
    return out


def _finalize(ax, xlabel, ylabel=r"$\Delta E$ (meV/atome)"):
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel); ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{x:g}"))
    ax.legend(framealpha=0.9, fontsize=6); ax.grid(True, which="both", ls=":", alpha=0.4)


def fig_electron_convergence(a, k_label=r"$12\times 12$", ecut_label="100 Ry", k_opt_label=r"$27\times27$"):
    """Quatre panneaux de convergence SCF (contenu de `plot_convergence.plot_all`)."""
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.2)); ax1, ax2, ax3, ax4 = axes.flatten()
    thr = dict(color=REF, lw=0.8, ls="--", alpha=0.8, label=f"Seuil {THRESH} meV/atome")
    ecut, dE = load_ecut("ecut.dat")
    ax1.plot(ecut, dE, "o-", color=NAVY, lw=1.2, ms=4, label=f"Grille {k_label}"); ax1.axhline(THRESH, **thr)
    ax1.set_title(r"(a) Convergence en coupure $E_\mathrm{cut}$", fontsize=8, loc="left"); ax1.set_xlim(32); _finalize(ax1, r"$E_\mathrm{cut}$ (Ry)")
    K, dE = load_kpoint()
    ax2.plot(K, np.abs(dE), "s-", color=ORANGE, lw=1.2, ms=4, label=rf"$E_\mathrm{{cut}}$ = {ecut_label}"); ax2.axhline(THRESH, **thr)
    ax2.set_title("(b) Convergence en densité de la grille", fontsize=8, loc="left"); ax2.set_xticks([9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39]); ax2.set_xlim(6, 40)
    _finalize(ax2, r"Grille $N \times N$"); ax2.set_ylim(bottom=0.001, top=10)
    ecut1, dE1 = load_ecut("ecut.dat"); ecut2, dE2 = load_ecut("ecut_cross.dat")
    ax3.plot(ecut2, dE2, "o-", color=NAVY, lw=1.2, ms=4, label=f"Grille optimale {k_opt_label}")
    ax3.plot(ecut1, dE1, "x--", color=ORANGE, lw=1.0, ms=4, alpha=0.8, label=f"Grille modérée {k_label}"); ax3.axhline(THRESH, **thr)
    ax3.set_title(r"(c) Convergence croisée en coupure $E_\mathrm{cut}$", fontsize=8, loc="left"); ax3.set_xlim(32); _finalize(ax3, r"$E_\mathrm{cut}$ (Ry)")
    for sm, (dg, dE) in load_smearing().items():
        ax4.plot(dg, np.abs(dE), SMEARING_MARKERS.get(sm, "o") + "-", color=SMEARING_COLORS.get(sm, INK), lw=1.2, ms=4, label=SMEARING_LABELS.get(sm, sm))
    ax4.axhline(THRESH, **thr)
    ax4.set_title(r"(d) Convergence en occupation et en élargissement $\sigma$", fontsize=8, loc="left"); ax4.set_xlim(left=0)
    _finalize(ax4, r"$\sigma$ (Ry)"); ax4.set_ylim(bottom=0.002); ax4.legend(loc="lower right", fontsize=6)
    fig.tight_layout()
    save(fig, "fig_electron_convergence", a.outdir); plt.close(fig)


def load_ebands(fname="ebands.dat"):
    data = np.loadtxt(os.path.join(ELECTRON, fname)); k = np.unique(data[:, 0])
    ebands = np.reshape(data[:, 1], (-1, len(k)))
    return k, ebands


def load_edos(fname="edos.dat"):
    energy, edos, idos = np.loadtxt(os.path.join(ELECTRON, fname), unpack=True)
    return energy, edos, idos


def fig_ebands_edos(a, high_sym_id=(0, 0.6667, 1.0000, 1.5774), high_sym_labels=(r"$\Gamma$", "K", "M", r"$\Gamma$"),
                    ebands_fermi=-1.8584, edos_fermi=-4.2352):
    """(a) bandes bands.x sur Γ–K–M–Γ, (b) DOS dos.x ; chacune référencée à son propre E_F (contenu de `plot_ebands_edos`)."""
    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(6.5, 4.0), gridspec_kw={"width_ratios": [3, 1]})
    k, ebands = load_ebands(); ebands = ebands - ebands_fermi
    for n in range(len(ebands)):
        ax[0].plot(k, ebands[n], lw=0.9, alpha=0.9, color=NAVY)
    for x in high_sym_id:
        ax[0].axvline(x, lw=0.6, color=MUTED)
    ax[0].set_xlim(min(k), max(k)); ax[0].set_xticks(list(high_sym_id)); ax[0].set_xticklabels(list(high_sym_labels))
    ax[0].set_ylabel(r"Énergie $\varepsilon - E_F$ (eV)"); ax[0].set_xlabel(r"Chemin dans l'espace $\mathbf{k}$"); ax[0].set_ylim(-20, 21)
    ax[0].axhline(0, linestyle=(0, (5, 5)), lw=0.75, color=ORANGE); ax[0].text(1.27, -1.4, r"Énergie de Fermi $E_F$", fontsize=7, color=ORANGE)
    ax[0].text(0.02, 0.98, "(a)", transform=ax[0].transAxes, fontsize=9, va="top", ha="left")
    energy, edos, _ = load_edos(); energy = energy - edos_fermi
    ax[1].plot(edos, energy, lw=0.75, color=NAVY); ax[1].set_xticks([]); ax[1].set_xlabel(r"DOS (états/eV)", fontsize=8); ax[1].set_xlim(0, None)
    ax[1].fill_betweenx(energy, 0, edos, where=(energy < 0), facecolor=SKY, alpha=0.35, lw=0)
    ax[1].axhline(0, linestyle=(0, (5, 5)), lw=0.75, color=ORANGE)
    ax[1].text(0.02, 0.98, "(b)", transform=ax[1].transAxes, fontsize=9, va="top", ha="left")
    fig.tight_layout()
    save(fig, "fig_ebands_edos", a.outdir); plt.close(fig)


FIGURES = {"kb_pseudo_C": fig_kb_pseudo_C, "electron_convergence": fig_electron_convergence, "ebands_edos": fig_ebands_edos}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default="", help="noms séparés par des virgules (sans le préfixe fig_)")
    ap.add_argument("--upf", default=os.path.join(ELECTRON, "C.upf"))
    ap.add_argument("--outdir", default=str(ROOT / "figures" / "electron"))
    a = ap.parse_args()
    use_style()
    for n in [n for n in a.only.split(",") if n] or list(FIGURES):
        FIGURES[n](a); print(f"écrit {a.outdir}/fig_{n}.pdf/.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
