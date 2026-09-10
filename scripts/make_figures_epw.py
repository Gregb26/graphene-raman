#!/usr/bin/env python
"""
Figures du chapitre 5 (couplage électron-phonon, EPW), style figures/memoire.mplstyle, français, données lues dans results/epw/*.npz
et results/M/resonance_9x9.npz (Γ^ed). Trois figures :
  fig_epw_validation : (a) bandes DFT (bands.x) vs EPW/Wannier, (b) phonons matdyn vs EPW, chemin Γ–K–M–Γ    [validation_24k24q.npz]
  fig_epw_gamma      : (a) Γ^ep(ε) convergence (degaussw à 120², 240²), (b) Γ^ep(ε) à 300 K et 10 K (production)  [selfen_*.npz]
  fig_epw_vs_ed      : Γ^ep(300 K) et Γ^ed(c = 1 %) = c × Γ N_cells (matrice T, 9×9, R_cut 3, η 0.02) sur le même axe
Optionnel (--control) : fig_epw_g_control, |g| EPW vs DFPT aux q déjà appariés (sommes gauge-invariantes, q = M exclu : double comptage).
Usage : python scripts/make_figures_epw.py [--outdir figures] [--control] [--prod-tag prod]
"""
import argparse, os, glob, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.style.use("figures/memoire.mplstyle")
ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="figures"); ap.add_argument("--control", action="store_true"); ap.add_argument("--prod-tag", default="prod"); a = ap.parse_args()
C_DFT, C_EPW, C_T, C_10K, MUTED = "#52514e", "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
CONV_COL = {"120_dg0.01": "#eda100", "120_dg0.02": "#2a78d6", "120_dg0.05": "#e87ba4", "240_dg0.02": "#4a3aa7"}
LBL_E = r"Énergie $\varepsilon - E_D$ (eV)"
# chemin Γ–K–M–Γ : longueurs |ΓK| = 2/3, |KM| = 1/3, |MΓ| = 1/√3 (unités 2π/a)
L = np.array([2 / 3, 1 / 3, 1 / np.sqrt(3)]); TICKS = np.concatenate([[0], np.cumsum(L)]) / L.sum(); TLAB = [r"$\Gamma$", "K", "M", r"$\Gamma$"]
def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    w, h = fig.get_size_inches(); plt.close(fig); print(f"écrit {a.outdir}/{name}.pdf/.png ({w:.2f} × {h:.2f} po)")
def path_axis(ax):
    ax.set_xticks(TICKS); ax.set_xticklabels(TLAB); ax.set_xlim(0, 1)
    for t in TICKS[1:-1]: ax.axvline(t, color=MUTED, lw=0.6)

# ---------------- 1. validation bandes + phonons
V = np.load("results/epw/validation_24k24q.npz", allow_pickle=True)
ED = float(V["bands_ED_dft"]); s = V["bands_s"]; Ed = V["bands_E_dft"] - ED; Ee = V["bands_E_epw_on_dft"] - ED
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.4))
for n in range(Ed.shape[1]): a1.plot(s, Ed[:, n], color=C_DFT, lw=1.0, label="DFT (bands.x)" if n == 0 else None)
for n in range(Ee.shape[1]): a1.plot(s, Ee[:, n], color=C_EPW, lw=1.2, ls="--", label="EPW (Wannier, 24×24)" if n == 0 else None)
a1.axhline(0, color=MUTED, lw=0.8, ls=":"); a1.axhline(2.5, color=C_T, lw=0.8, ls=":"); a1.text(0.99, 2.6, r"fenêtre gelée $E_D+2{,}5$ eV", ha="right", va="bottom", fontsize=7, color=C_T)
a1.set_ylim(-21, 9); a1.set_ylabel(LBL_E); a1.legend(fontsize=7, loc="lower left"); path_axis(a1); panel(a1, "a")
CM2MEV = 1.0 / 8.06554; sq = V["ph_s"]; Fm = V["ph_F_matdyn"] * CM2MEV; Fe = V["ph_F_epw_on_matdyn"] * CM2MEV
for m in range(6): a2.plot(sq, Fm[:, m], color=C_DFT, lw=1.0, label="DFPT (matdyn)" if m == 0 else None); a2.plot(sq, Fe[:, m], color=C_EPW, lw=1.2, ls="--", label="EPW (24×24 q)" if m == 0 else None)
a2.set_ylabel(r"$\hbar\omega_{\mathbf{q}\nu}$ (meV)"); a2.set_ylim(0, 210); a2.legend(fontsize=7, loc="lower right"); path_axis(a2); panel(a2, "b")
fig.tight_layout(); save(fig, "fig_epw_validation")

# ---------------- 2. Γ^ep(ε) : convergence + température
def load_sel(f): R = np.load(f, allow_pickle=True); return R["eg"] - float(R["E_D"]), R["Gamma_e"] * 1e3, R
conv = {k: f"results/epw/selfen_{k}_T300.npz" for k in CONV_COL}; conv = {k: f for k, f in conv.items() if os.path.exists(f)}
prod = {T: f"results/epw/selfen_{a.prod_tag}_T{T}.npz" for T in (300, 10)}; prod = {T: f for T, f in prod.items() if os.path.exists(f)}
if conv:
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(6.5, 3.4), sharey=True)
    for k, f in conv.items():
        x, G, R = load_sel(f); b1.plot(x, G, color=CONV_COL[k], lw=1.2, label=rf"{int(R['n_mesh'])}$^2$, $\sigma$ = {float(R['degaussw']):.2f} eV")
    b1.set_xlabel(LBL_E); b1.set_ylabel(r"$\Gamma^{ep}(\varepsilon)$ (meV)"); b1.set_title(r"Convergence, $T$ = 300 K", loc="left", fontsize=9); b1.legend(fontsize=7); b1.axvline(0, color=MUTED, lw=0.8, ls=":"); panel(b1, "a")
    for T, f in prod.items():
        x, G, R = load_sel(f); b2.plot(x, G, color={300: C_T, 10: C_10K}[T], lw=1.2, label=rf"$T$ = {T} K")
    b2.set_xlabel(LBL_E); b2.set_title("Production (paramètres convergés)", loc="left", fontsize=9); b2.axvline(0, color=MUTED, lw=0.8, ls=":"); panel(b2, "b")
    if prod: b2.legend(fontsize=8)
    fig.tight_layout(); save(fig, "fig_epw_gamma")

# ---------------- 3. Γ^ep(300 K) vs Γ^ed(c = 1 %)
if 300 in prod and os.path.exists("results/M/resonance_9x9.npz"):
    x, G, R = load_sel(prod[300]); M = np.load("results/M/resonance_9x9.npz", allow_pickle=True); c = float(M["conc"]); xe = M["eg"] - float(M["E_D"])
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.semilogy(xe, c * M["Gamma_T"] * 1e3, color=C_EPW, label=rf"$\Gamma^{{ed}}$, lacune, $c$ = {c*100:.0f}\,\% (matrice $T$, 9$\times$9, $\eta$ = {float(M['eta'])} eV)")
    ax.semilogy(x, G, color=C_T, label=rf"$\Gamma^{{ep}}$, $T$ = 300 K (EPW, {int(R['n_mesh'])}$^2$, $\sigma$ = {float(R['degaussw']):.2f} eV)")
    ax.set_xlabel(LBL_E); ax.set_ylabel(r"$\Gamma$ (meV)"); ax.axvline(0, color=MUTED, lw=0.8, ls=":"); ax.legend(fontsize=8); ax.set_xlim(-3, 3)
    fig.tight_layout(); save(fig, "fig_epw_vs_ed")

# ---------------- optionnel : contrôle |g|
if a.control:
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.2))
    for ax, kname, letter in ((axs[0], "G", "a"), (axs[1], "K", "b")):
        g = V[f"g_{kname}"]; keep = (g[:, 4] > 20) & ~np.isclose(g[:, 0], 0.634, atol=0.003)     # q = M exclu (double comptage des groupes de modes)
        d, e = g[keep, 4], g[keep, 5]; rel = np.abs(e - d) / d
        ax.loglog(d, e, "o", ms=3, color=C_EPW, label=rf"{keep.sum()} sommes, $q \in \{{{', '.join(f'{v:.3f}' for v in np.unique(np.round(g[keep,0],3)))}\}}$")
        lim = [min(d.min(), e.min()) * 0.8, max(d.max(), e.max()) * 1.2]; ax.plot(lim, lim, color=MUTED, lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(r"$G$ DFPT (meV)"); ax.set_ylabel(r"$G$ EPW (meV)"); klab = r"$\Gamma$" if kname == "G" else "K"; ax.set_title(rf"$k$ = {klab} : médiane {np.median(rel)*100:.1f}\,\%, max {rel.max()*100:.1f}\,\%", loc="left", fontsize=9); ax.legend(fontsize=6, loc="lower right"); panel(ax, letter)
        print(f"[contrôle |g| k={kname}] {keep.sum()} sommes (G>20 meV, q=M exclu) : rel. médiane {np.median(rel):.2%}, max {rel.max():.2%}, |ΔG| max {np.abs(e-d).max():.2f} meV")
    fig.tight_layout(); save(fig, "fig_epw_g_control")
