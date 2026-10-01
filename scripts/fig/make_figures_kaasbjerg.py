"""Les cinq figures de la campagne R8 (DOS et A_k moyennées sur le désordre, Kaasbjerg PRB 101, 045433) retenues pour le chapitre 4,
redessinées au style graphene_raman.plotting depuis les sorties de `campagnes/R/R8_kaasbjerg/r8_driver.py` (aucun calcul).
Tracé extrait le 2026-09-30 de `r8_driver.cmd_fig` (quatre premières) et de la fin de `r8_driver.cmd_sigeff` (`sigma_K`, qui lit
`sigma_K_9x9.npz` au lieu de recalculer Σ), au contenu identique :

  fig_kaasbjerg_dos_c          (a) ρ(ε) tel quel, ρ₀ 1 200² et c_i = 0,1 et 1 % ; (b), (c) ρ − ρ₀ des trois variantes (600² plein, 300² pointillé)
  fig_kaasbjerg_spectral_GKM   A_k(ε) sur Γ–K–M, 2 concentrations × 3 variantes, branches où A ≥ 1 % du maximum
  fig_kaasbjerg_sensibilites   max de ρ − ρ₀ contre N_k^int, R_cut, taille, η_G ; ligne : valeur de Kaasbjerg (Fig. 13, c_i = 1 %)
  fig_kaasbjerg_superposition  nos courbes (ρ₀ 1 200², η 15 meV) sur celles extraites de Kaasbjerg, ×1 (par spin) et ×2
  fig_kaasbjerg_sigma_K        Σ_nn et Σ^eff_nK à k = K, c_i = 1 %, tel quel et aligné (étape 4)

Données : <campaigns_results_dir>/R8_kaasbjerg/ = results/campagnes/R8_kaasbjerg/ (config ; copies de `out/` de la campagne) :
dos/dos_9x9.npz, spec/spectral_GKM_9x9.npz, sens/sens_results.json, 7a/fig13_VA_interp_1meV.csv, 7a/7a_results.json, sigeff/sigma_K_9x9.npz, sigeff/sigeff_results.json.
Sortie : figures/electron_defect/fig_kaasbjerg_<nom>.{pdf,png} (reproductibles au bit, palette.save).
Usage : python scripts/fig/make_figures_kaasbjerg.py [--only nom,nom] [--outdir figures/electron_defect]
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

from graphene_raman.plotting.palette import use_style, save, NAVY, ORANGE, GREEN, SKY, GOLD, REF, MUTED, CMAP_SEQ
from graphene_raman.config import load_production, campaigns_dir

ROOT = Path(__file__).resolve().parents[2]            # scripts/fig/ -> racine
R8 = campaigns_dir(load_production(verbose=False), "R8_kaasbjerg")

# constantes de la campagne (r8_driver.py, 2026-09-28) : variantes, concentrations, grilles
VARIANTS = ("tel_quel", "aligne")
C_LIST = (0.001, 0.01)                   # c_i, défauts par maille
GRID_SENS, NK_MAIN = 600, 900
P0 = 1200                                # ρ₀ des figures : 1 200², η 15 meV (Greg, 2026-09-29)
col = {"tel_quel": NAVY, "aligne": ORANGE, "eta_unique": GREEN}
lab = {"tel_quel": "non aligné", "aligne": "Kumagai–Oba", "eta_unique": r"$\eta_t = \eta_G$"}


def jload(*p):
    return json.load(open(os.path.join(R8, *p), encoding="utf-8"))


def dos_inputs():
    D = np.load(os.path.join(R8, "dos", "dos_9x9.npz"))
    Kx = np.loadtxt(os.path.join(R8, "7a", "fig13_VA_interp_1meV.csv"), delimiter=",", skiprows=1)
    eK, rK = Kx[:, 0], {"pristine": Kx[:, 1], 0.001: Kx[:, 2], 0.01: Kx[:, 3]}
    K7 = jload("7a", "7a_results.json")["fig13_metrics"]
    egr = D["eg_rel"]
    rho0_fig = D[f"pristine_eta15_{P0}"]
    prot = lambda key: rho0_fig + (D[key] - D[key.rsplit("_c", 1)[0] + "_rho0"])          # rho0(1200^2, 15 meV) + delta rho
    return D, eK, rK, K7, egr, rho0_fig, prot


def fig_dos_c(a):
    D, eK, rK, K7, egr, rho0_fig, prot = dos_inputs()
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 2.5), sharex=True)
    ax[0].plot(egr, rho0_fig, color=REF, lw=0.8, label=r"$\rho_0$ : 1\,200$^2$, $\eta$ 15 meV")
    for cv, cc in ((0.001, SKY), (0.01, NAVY)):
        ax[0].plot(egr, prot(f"tel_quel_{GRID_SENS}_c{cv}"), color=cc, lw=0.9, label=rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%")
    for i, cv in enumerate(C_LIST):
        for name in ("tel_quel", "aligne", "eta_unique"):
            ax[i + 1].plot(egr, D[f"{name}_{GRID_SENS}_c{cv}"] - D[f"{name}_{GRID_SENS}_rho0"], color=col[name], lw=0.9, label=lab[name])
            ax[i + 1].plot(egr, D[f"{name}_300_c{cv}"] - D[f"{name}_300_rho0"], color=col[name], lw=0.5, ls=":")
        ax[i + 1].set_title(rf"({'bc'[i]}) $\rho - \rho_0$, $c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%", fontsize=9)
    ax[0].set_title(r"(a) $\rho(\varepsilon)$, non aligné", fontsize=9)
    ax[0].set_ylabel(r"DOS (états/eV/maille/spin)"); ax[1].legend(fontsize=6); ax[0].legend(fontsize=6)
    for x in ax:
        x.set_xlabel(r"Énergie $\varepsilon - E_D$ (eV)"); x.set_xlim(-1.2, 1.2)
    fig.tight_layout(); save(fig, "fig_kaasbjerg_dos_c", a.outdir); plt.close(fig)


def fig_spectral_GKM(a):
    Sp = np.load(os.path.join(R8, "spec", "spectral_GKM_9x9.npz"))
    names = ("tel_quel", "aligne", "eta_unique")
    fig, ax = plt.subplots(2, 3, figsize=(6.5, 4.4), sharex=True, sharey=True)
    s, idx, ege = Sp["s"], Sp["idx"], Sp["eg_rel"]
    for i, cv in enumerate(C_LIST):
        for j, name in enumerate(names):
            A = Sp[f"A_{name}_c{cv}"]; x = ax[i, j]
            x.imshow(A, origin="lower", aspect="auto", extent=(s[0], s[-1], ege[0], ege[-1]), cmap=CMAP_SEQ, vmin=0, vmax=np.percentile(A, 99.5))
            for br in ("low", "up"):                                 # plotted only where A >= 1 % of the map maximum (json keeps all)
                e_ = Sp[f"{br}_{name}_c{cv}"]; ok_ = np.isfinite(e_)
                h_ = np.zeros_like(e_); h_[ok_] = A[np.abs(ege[:, None] - e_[None, ok_]).argmin(0), np.where(ok_)[0]]
                vis = ok_ & (h_ >= 0.01 * A.max()); x.plot(s[vis], e_[vis], ".", ms=0.8, color=ORANGE)
            x.plot(s, Sp["E_pristine"], color=REF, lw=0.4, ls="--")
            x.set_xticks(s[idx]); x.set_xticklabels([r"$\Gamma$", r"$K$", r"$M$"]); x.set_ylim(-1.2, 1.2)
            x.set_title(f"({'abcdef'[3 * i + j]}) {lab[name]}, " + rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%", fontsize=8)
        ax[i, 0].set_ylabel(r"$\varepsilon - E_D$ (eV)")
    fig.tight_layout(); save(fig, "fig_kaasbjerg_spectral_GKM", a.outdir); plt.close(fig)


def fig_sensibilites(a):
    Se = jload("sens", "sens_results.json"); K7 = jload("7a", "7a_results.json")["fig13_metrics"]
    fig, ax = plt.subplots(1, 4, figsize=(6.5, 2.3), sharey=True)
    groups = (("N_k_int", "N_k_int", r"$N_k^\mathrm{int}$"), ("R_cut", "R_cut", r"$R_\mathrm{cut}$"), ("taille", "size", "taille"),
              ("eta_G", "eta_G", r"$\eta_G$ (eV)"))
    for x, (g, fld, xl) in zip(ax, groups):
        for var in VARIANTS:
            rows = [r for r in Se["rows"] if (r["group"] == g or (g == "R_cut" and r["group"] == "N_k_int" and r["N_k_int"] == NK_MAIN)
                                               or (g in ("taille", "eta_G") and r["group"] == "N_k_int" and r["N_k_int"] == NK_MAIN))
                    and r["variant"] == var]
            xs = [r[fld] if fld != "size" else (9 if r["size"] == "9x9" else 12) for r in rows]
            o = np.argsort(xs); x.plot(np.array(xs)[o], np.array([r["drho_max_pos_eV"] for r in rows])[o], "o-", ms=3, lw=0.8, color=col[var],
                                       label=lab[var])
        x.axhline(K7["1.0"]["drho_max_pos_eV"], color=MUTED, lw=0.6, ls="--")
        x.set_xlabel(xl)
    ax[0].set_ylabel(r"max de $\rho - \rho_0$ (eV)"); ax[0].legend(fontsize=6)
    fig.tight_layout(); save(fig, "fig_kaasbjerg_sensibilites", a.outdir); plt.close(fig)


def fig_superposition(a):
    D, eK, rK, K7, egr, rho0_fig, prot = dos_inputs()
    fig, ax = plt.subplots(2, 2, figsize=(6.5, 4.2), sharex=True)
    for i, f_ in enumerate((1, 2)):
        for j, cv in enumerate(C_LIST):
            x = ax[i, j]
            x.plot(eK, rK["pristine"], color=MUTED, lw=0.6, ls="--", label="Kaasbjerg, parfait")
            x.plot(eK, rK[cv], color=REF, lw=1.2, label="Kaasbjerg")
            for name in ("tel_quel", "aligne", "eta_unique"):
                x.plot(egr, f_ * prot(f"{name}_{GRID_SENS}_c{cv}"), color=col[name], lw=0.8, label=lab[name])
            x.set_xlim(-1.2, 1.2)
            x.set_title(f"({'abcd'[2 * i + j]}) " + rf"$c_i$ = {cv * 100:g}".replace(".", ",") + r"\,\%, " + ("par spin" if f_ == 1 else r"$\times 2$ (spin compté)"),
                        fontsize=8)
        ax[i, 0].set_ylabel(r"DOS (eV$^{-1}$)")
    for x in ax[1]:
        x.set_xlabel(r"Énergie $\varepsilon - E_D$ (eV)")
    ax[0, 0].plot([], [], " ", label=r"(nous : $\rho_0$ 1\,200$^2$, $\eta$ 15 meV)")      # legend line only
    ax[0, 0].legend(fontsize=6)
    fig.tight_layout(); save(fig, "fig_kaasbjerg_superposition", a.outdir); plt.close(fig)


def fig_sigma_K(a):
    Z = np.load(os.path.join(R8, "sigeff", "sigma_K_9x9.npz")); runs = jload("sigeff", "sigeff_results.json")["runs"]
    eg = Z["eg_rel"]
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 2.8), sharey=True)
    for iv, var in enumerate(VARIANTS):
        x = axs[iv]; Se = Z[f"Seff_{var}_K"]; Sb = Z[f"Sigma_{var}_K"]
        x.plot(eg, eg - runs[f"{var} K"]["eps_k"][0], color=REF, lw=0.7, ls="--", label=r"$\varepsilon - \varepsilon_{nK}$")
        x.plot(eg, Sb[:, 0, 0].real, color=SKY, lw=0.7, label=r"Re $\Sigma_{nn}$")
        x.plot(eg, Sb[:, 0, 0].imag, color=GOLD, lw=0.7, label=r"Im $\Sigma_{nn}$")
        x.plot(eg, Se[:, 0].real, color=NAVY, lw=1.0, label=r"Re $\Sigma^\mathrm{eff}_{nK}$")
        x.plot(eg, Se[:, 0].imag, color=ORANGE, lw=1.0, label=r"Im $\Sigma^\mathrm{eff}_{nK}$")
        x.set_xlim(-1, 1); x.set_ylim(-0.35, 0.35); x.set_xlabel(r"Énergie $\varepsilon - E_D$ (eV)")
        x.set_title(f"({'ab'[iv]}) " + ("non aligné" if var == "tel_quel" else "Kumagai–Oba") + r", $c_i$ = 1\,\%, $k = K$", fontsize=9)
    axs[0].set_ylabel(r"$\Sigma$ (eV)"); axs[0].legend(fontsize=6, loc="lower left")
    fig.tight_layout(); save(fig, "fig_kaasbjerg_sigma_K", a.outdir); plt.close(fig)


FIGURES = {"dos_c": fig_dos_c, "spectral_GKM": fig_spectral_GKM, "sensibilites": fig_sensibilites, "superposition": fig_superposition,
           "sigma_K": fig_sigma_K}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default=None, help="noms séparés par des virgules : " + ", ".join(FIGURES))
    ap.add_argument("--outdir", default=str(ROOT / "figures" / "electron_defect"))
    a = ap.parse_args()
    use_style()
    names = a.only.split(",") if a.only else list(FIGURES)
    for n in names:
        FIGURES[n](a); print(f"écrit {a.outdir}/fig_kaasbjerg_{n}.pdf/.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
