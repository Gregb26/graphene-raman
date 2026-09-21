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
from matplotlib.ticker import FuncFormatter, ScalarFormatter, LogFormatterSciNotation
def fr(x, nd=2):
    """nombre formaté (POINT décimal : convention du mémoire depuis P18, 2026-09-20 ; l'ancienne virgule() n'est plus appliquée)."""
    return f"{x:.{nd}f}"
def virgule(ax):
    """conservée pour compatibilité : ne fait plus rien (point décimal sur les ticks, formateurs matplotlib par défaut)."""
    return None
ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="figures"); ap.add_argument("--control", action="store_true"); ap.add_argument("--prod-tag", default="prod"); ap.add_argument("--phself-tag", default="path_1200_dg0.02")
ap.add_argument("--val-tag", default="24k24q", help="validation_<tag>.npz (bandes, phonons, décroissance, |g|)"); ap.add_argument("--sel-suffix", default="", help="suffixe des selfen de convergence : selfen_<k><suffixe>_T300.npz")
ap.add_argument("--kohn-val-tags", default="24k24q,24k24q_mv0.02", help="tags des validation npz (matdyn) pour fig_epw_kohn_degauss, dans l'ordre"); ap.add_argument("--kohn-labels", default="degauss 0.002 Ry,degauss 0.02 Ry"); ap.add_argument("--dfpt-tag", default="24k24q", help="dfpt_path_freq_<tag>.npz (DFPT direct 16×16, 31 q)"); a = ap.parse_args()
C_DFT, C_EPW, C_T, C_10K, MUTED = "#52514e", "#2a78d6", "#eb6834", "#1baf7a", "#8a8984"
CONV_COL = {"120_dg0.01": "#eda100", "120_dg0.02": "#2a78d6", "120_dg0.05": "#e87ba4", "240_dg0.02": "#4a3aa7"}
LBL_E = r"Énergie $\varepsilon - E_D$ (eV)"
# chemin Γ–K–M–Γ : longueurs |ΓK| = 2/3, |KM| = 1/3, |MΓ| = 1/√3 (unités 2π/a)
L = np.array([2 / 3, 1 / 3, 1 / np.sqrt(3)]); TICKS = np.concatenate([[0], np.cumsum(L)]) / L.sum(); TLAB = [r"$\Gamma$", "K", "M", r"$\Gamma$"]
def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())
def save(fig, name):
    for ax in fig.axes: virgule(ax)
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    w, h = fig.get_size_inches(); plt.close(fig); print(f"écrit {a.outdir}/{name}.pdf/.png ({w:.2f} × {h:.2f} po)")
def path_axis(ax):
    ax.set_xticks(TICKS); ax.set_xticklabels(TLAB); ax.set_xlim(0, 1); ax.xaxis._etiquettes_fixes = True
    for t in TICKS[1:-1]: ax.axvline(t, color=MUTED, lw=0.6)

# ---------------- 1. validation bandes + phonons
V = np.load(f"results/epw/validation_{a.val_tag}.npz", allow_pickle=True)
ED = float(V["bands_ED_dft"]); s = V["bands_s"]; Ed = V["bands_E_dft"] - ED; Ee = V["bands_E_epw_on_dft"] - ED
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.5, 3.4))
for n in range(Ed.shape[1]): a1.plot(s, Ed[:, n], color=C_DFT, lw=1.0, label="DFT (bands.x)" if n == 0 else None)
for n in range(Ee.shape[1]): a1.plot(s, Ee[:, n], color=C_EPW, lw=1.2, ls="--", label="EPW (Wannier, 24×24)" if n == 0 else None)
a1.axhline(0, color=MUTED, lw=0.8, ls=":"); a1.axhline(2.5, color=C_T, lw=0.8, ls=":"); a1.text(0.99, 2.6, r"fenêtre gelée $E_D+2.5$ eV", ha="right", va="bottom", fontsize=7, color=C_T)
a1.set_ylim(-21, 9); a1.set_ylabel(LBL_E); a1.legend(fontsize=7, loc="lower left"); path_axis(a1); panel(a1, "a")
CM2MEV = 1.0 / 8.06554; sq = V["ph_s"]; Fm = V["ph_F_matdyn"] * CM2MEV; Fe = V["ph_F_epw_on_matdyn"] * CM2MEV
for m in range(6): a2.plot(sq, Fm[:, m], color=C_DFT, lw=1.0, label="DFPT (matdyn)" if m == 0 else None); a2.plot(sq, Fe[:, m], color=C_EPW, lw=1.2, ls="--", label="EPW (24×24 q)" if m == 0 else None)
a2.set_ylabel(r"$\hbar\omega_{\mathbf{q}\nu}$ (meV)"); a2.set_ylim(0, 210); a2.legend(fontsize=7, loc="lower right"); path_axis(a2); panel(a2, "b")
fig.tight_layout(); save(fig, "fig_epw_validation")

# ---------------- 2. Γ^ep(ε) : convergence + température
def load_sel(f): R = np.load(f, allow_pickle=True); return R["eg"] - float(R["E_D"]), R["Gamma_e"] * 1e3, R
conv = {k: f"results/epw/selfen_{k}{a.sel_suffix}_T300.npz" for k in CONV_COL}; conv = {k: f for k, f in conv.items() if os.path.exists(f)}
prod = {T: f"results/epw/selfen_{a.prod_tag}_T{T}.npz" for T in (300, 10)}; prod = {T: f for T, f in prod.items() if os.path.exists(f)}
if conv:
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(6.5, 3.4), sharey=True)
    for k, f in conv.items():
        x, G, R = load_sel(f); b1.plot(x, G, color=CONV_COL[k], lw=1.2, label=rf"{int(R['n_mesh'])}$^2$, $\sigma$ = {fr(float(R['degaussw']))} eV")
    b1.set_xlabel(LBL_E); b1.set_ylabel(r"$\Gamma^{ep}(\varepsilon)$ (meV)"); b1.set_title(r"Convergence, $T$ = 300 K", loc="left", fontsize=9); b1.legend(fontsize=7); b1.axvline(0, color=MUTED, lw=0.8, ls=":"); panel(b1, "a")
    for T, f in prod.items():
        x, G, R = load_sel(f); b2.plot(x, G, color={300: C_T, 10: C_10K}[T], lw=1.2, label=rf"$T$ = {T} K")
    b2.set_xlabel(LBL_E); b2.set_title(r"Production, 240$^2$, $\sigma$ = 0.02 eV", loc="left", fontsize=9); b2.axvline(0, color=MUTED, lw=0.8, ls=":"); panel(b2, "b")
    if prod: b2.legend(fontsize=8)
    fig.tight_layout(); save(fig, "fig_epw_gamma")

# ---------------- 3. Γ^ep(300 K) vs Γ^ed(c = 1 %)
if 300 in prod and os.path.exists("results/M/resonance_9x9.npz"):
    x, G, R = load_sel(prod[300]); M = np.load("results/M/resonance_9x9.npz", allow_pickle=True); c = float(M["conc"]); xe = M["eg"] - float(M["E_D"])
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    ax.semilogy(xe, c * M["Gamma_T"] * 1e3, color=C_EPW, label=rf"$\Gamma^{{ed}}$, lacune, $c$ = {c*100:.0f}\,\% (matrice $T$, 9$\times$9, $\eta$ = {fr(float(M['eta']))} eV)")
    ax.semilogy(x, G, color=C_T, label=rf"$\Gamma^{{ep}}$, $T$ = 300 K (EPW, {int(R['n_mesh'])}$^2$, $\sigma$ = {fr(float(R['degaussw']))} eV)")
    ax.set_xlabel(LBL_E); ax.set_ylabel(r"$\Gamma$ (meV)"); ax.axvline(0, color=MUTED, lw=0.8, ls=":"); ax.legend(fontsize=8); ax.set_xlim(-3, 3)
    fig.tight_layout(); save(fig, "fig_epw_vs_ed")

# ---------------- optionnel : contrôle |g|
if a.control:
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.2))
    for ax, kname, letter in ((axs[0], "G", "a"), (axs[1], "K", "b")):
        g = V[f"g_{kname}"]; keep = (g[:, 4] > 20) & ~np.isclose(g[:, 0], 0.634, atol=0.003)     # q = M exclu (double comptage des groupes de modes)
        d, e = g[keep, 4], g[keep, 5]; rel = np.abs(e - d) / d
        ax.loglog(d, e, "o", ms=3, color=C_EPW, label=rf"{keep.sum()} sommes, $q \in \{{{', '.join(fr(v, 3) for v in np.unique(np.round(g[keep,0],3)))}\}}$")
        lim = [min(d.min(), e.min()) * 0.8, max(d.max(), e.max()) * 1.2]; ax.plot(lim, lim, color=MUTED, lw=0.8); ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel(r"$G$ DFPT (meV)"); ax.set_ylabel(r"$G$ EPW (meV)"); klab = r"$\Gamma$" if kname == "G" else "K"; ax.set_title(rf"$k$ = {klab} : médiane {fr(np.median(rel)*100, 1)}\,\%, max {fr(rel.max()*100, 1)}\,\%", loc="left", fontsize=9); ax.legend(fontsize=6, loc="lower right"); panel(ax, letter)
        print(f"[contrôle |g| k={kname}] {keep.sum()} sommes (G>20 meV, q=M exclu) : rel. médiane {np.median(rel):.2%}, max {rel.max():.2%}, |ΔG| max {np.abs(e-d).max():.2f} meV")
    fig.tight_layout(); save(fig, "fig_epw_g_control")

# ---------------- 4. fig_epw_phonselfen : γ_qν le long de Γ–K–M–Γ (300 K) + valeurs clés aux deux T
PH = f"results/epw/phself_{a.phself_tag}.npz"
if os.path.exists(PH):
    import electron_defect_interaction.electron_phonon.phself as _ph
    P = np.load(PH, allow_pickle=True); T = P["T"]; s = P["s"]; om = P["omega"]; gam = P["gamma_fwhm"]; i300 = int(np.argmin(np.abs(T - 300))); i10 = int(np.argmin(np.abs(T - 10)))
    fig, (c1, c2) = plt.subplots(1, 2, figsize=(6.5, 3.4), gridspec_kw={"width_ratios": [2.2, 1]})
    BR_COL = ["#b5b4b0", "#8a8984", "#52514e", "#1baf7a", "#eda100", "#eb6834"]
    for m in range(6):
        g = gam[i300, :, m].copy()
        if m < 3: g[om[:, m] < 5.0] = np.nan                 # branches acoustiques : divergence q→0 sous smearing, masquées près de Γ
        c1.plot(s, g, color=BR_COL[m], lw=1.1, label=rf"$\nu$ = {m+1}")
    sp = _ph.special_points(P["q"]); iG, iK = sp["G"][0], sp["K"][0]
    Rp = dict(T=T, omega=om, gamma_epw=P["gamma_epw"]); mE = list(_ph.modes_E2g(Rp, iG)); mA = _ph.mode_A1p(Rp, iK)     # par fréquence / par caractère, jamais par index
    print(f"[phonselfen] E2g(Gamma) = modes {mE[0]+1}+{mE[1]+1} ({om[iG, mE].mean():.2f} meV), A1'(K) = mode {mA+1} ({om[iK, mA]:.2f} meV, plus grand gamma à K)")
    c1.annotate(r"E$_{2g}$", (s[iG], gam[i300, iG, mE[0]]), xytext=(6, 4), textcoords="offset points", fontsize=8)
    c1.annotate(r"A$_1'$", (s[iK], gam[i300, iK, mA]), xytext=(6, -2), textcoords="offset points", fontsize=8)
    c1.set_ylabel(r"$\gamma_{\mathbf{q}\nu}$ (meV, largeur totale)"); c1.set_title(rf"$T$ = 300 K, {int(P['nkf'])}$^2$ $k$, $\sigma$ = {fr(float(P['degaussw']))} eV", loc="left", fontsize=9)
    c1.legend(fontsize=7, loc="upper right", title="branches (tri en fréquence)", title_fontsize=7, ncol=2); path_axis(c1); panel(c1, "a")
    vals = {r"E$_{2g}$ ($\Gamma$)": (gam[i10, iG, mE].mean(), gam[i300, iG, mE].mean()), r"A$_1'$ (K)": (gam[i10, iK, mA], gam[i300, iK, mA])}
    x = np.arange(len(vals)); w = 0.36
    c2.bar(x - w / 2, [v[0] for v in vals.values()], w, color=C_10K, label=r"$T$ = 10 K"); c2.bar(x + w / 2, [v[1] for v in vals.values()], w, color=C_T, label=r"$T$ = 300 K")
    for xi, v in zip(x, vals.values()):
        for dx, val in ((-w / 2, v[0]), (w / 2, v[1])): c2.text(xi + dx, val, fr(val), ha="center", va="bottom", fontsize=7)
    c2.set_xticks(x); c2.set_xticklabels(list(vals)); c2.xaxis._etiquettes_fixes = True; c2.set_ylabel(r"$\gamma$ (meV, largeur totale)"); c2.legend(fontsize=7, loc="upper left"); c2.set_ylim(0, max(max(v) for v in vals.values()) * 1.35); panel(c2, "b")
    fig.tight_layout(); save(fig, "fig_epw_phonselfen")

# ---------------- 5. fig_epw_decay : décroissance de H, D et g en représentation de Wannier (chaîne 24k-24q, epw1/decay.*)
# Données : validation_24k24q.npz (decay_<q>_r en Å, decay_<q>_v en Ry tel qu'écrit par EPW), converties en eV (× RY2EV) sur l'axe y.
# Étiquettes en convention du mémoire (ch. 2) : atome α, direction μ ; en espace réel C = constantes de force (EPW « dynmat », masses identiques :
# facteur global) et pas d'indice de branche ν (il n'apparaît qu'après recombinaison par les vecteurs propres) ; g_{mn}, m final, n initial.
RY2EV = 13.605693122994
if all(f"decay_{k}_r" in V.files for k in ("H", "dynmat", "epmate", "epmatp")):
    a_A = 2.4659; ws_in = 24 * a_A / 2; ws_out = 24 * a_A / np.sqrt(3)          # demi-largeur (apothème) et rayon (sommet) de la cellule WS de la supercellule 24×24
    spec = [("H", r"$|R_e|$ (\AA)", r"max$_{nm}\,|H_{nm}(R_e)|$ (eV)", r"$H$ : hamiltonien", "#2a78d6"),
            ("dynmat", r"$|R_p|$ (\AA)", r"max$_{\alpha\mu,\alpha'\mu'}\,|C_{\alpha\mu,\alpha'\mu'}(R_p)|$ (eV)", r"$C$ : constantes de force", "#1baf7a"),
            ("epmate", r"$|R_e|$ (\AA)", r"max$_{mn,\alpha\mu}\,|g_{mn}^{(\alpha\mu)}(R_e,\,R_p)|$ (eV)", r"$g$ : côté électron ($R_e$)", "#eb6834"),
            ("epmatp", r"$|R_p|$ (\AA)", r"max$_{mn,\alpha\mu}\,|g_{mn}^{(\alpha\mu)}(R_e,\,R_p)|$ (eV)", r"$g$ : côté phonon ($R_p$)", "#eda100")]
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.6), sharex=True, sharey=False); axs = axs.ravel()
    stats = {}
    for i, (key, xl, yl, title, col) in enumerate(spec):
        r = V[f"decay_{key}_r"]; v = V[f"decay_{key}_v"] * RY2EV; ax = axs[i]
        ax.semilogy(r, v, "o", color=col, ms=3.2, zorder=2)
        ax.set_title(title, loc="left"); ax.set_ylabel(yl); ax.set_xlim(0, ws_out * 1.03); ax.grid(True, alpha=0.3); panel(ax, "abcd"[i])
        if i >= 2: ax.set_xlabel(xl)
        # chiffres pour le texte : maximum, plancher (médiane des R ≥ 0,9 R_max), distance où l'enveloppe (max glissant décroissant) passe sous 10 × plancher
        v0 = v[np.isclose(r, r.min())].max(); floor = np.median(v[r >= 0.9 * r.max()]); order = np.argsort(r); env = np.maximum.accumulate(v[order][::-1])[::-1]
        rc = r[order][np.argmax(env < 10 * floor)] if (env < 10 * floor).any() else np.nan
        stats[key] = (v0, floor, rc, np.log10(v0 / floor))
        print(f"[decay {key:7s}] max {v0:.3e} eV à R = {r.min():.2f} Å ; plancher {floor:.2e} eV (médiane R ≥ 0,9 R_max) ; chute de {np.log10(v0/floor):.1f} ordres ; enveloppe < 10 × plancher dès R = {rc:.1f} Å")
    print(f"[decay] cellule de Wigner-Seitz 24×24 : apothème (demi-largeur) {ws_in:.2f} Å, rayon (sommet) {ws_out:.2f} Å = R_max des fichiers")
    fig.tight_layout(); save(fig, "fig_epw_decay")

# ---------------- 6. fig_epw_kohn_degauss : deux branches optiques les plus hautes sur Γ–K–M–Γ, matdyn (degauss 0.002 / 0.02 Ry) et DFPT direct 16×16 (31 q)
KV = [(t, l) for t, l in zip(a.kohn_val_tags.split(","), a.kohn_labels.split(",")) if os.path.exists(f"results/epw/validation_{t}.npz")]
DF = f"results/epw/dfpt_path_freq_{a.dfpt_tag}.npz"
if KV and os.path.exists(DF):
    fig, ax = plt.subplots(figsize=(6.5, 3.6)); KCOL = ["#52514e", "#eb6834", "#2a78d6"]
    for i, (t, l) in enumerate(KV):
        Vt = np.load(f"results/epw/validation_{t}.npz", allow_pickle=True); Ft = Vt["ph_F_matdyn"] * CM2MEV
        for m in (4, 5): ax.plot(Vt["ph_s"], Ft[:, m], color=KCOL[i], lw=1.2, label=rf"matdyn 24$\times$24 q, {l}" if m == 4 else None)
    D = np.load(DF, allow_pickle=True); Fd = D["freq"] * CM2MEV
    for m in (4, 5): ax.plot(D["s"], Fd[:, m], "o", color=KCOL[2], ms=3.2, label=rf"DFPT direct, 16$\times$16 $k$ ({len(D['s'])} $q$)" if m == 4 else None)
    ax.set_ylabel(r"$\hbar\omega_{\mathbf{q}\nu}$ (meV)"); ax.set_title(r"Branches optiques $\nu$ = 5, 6 : anomalies de Kohn selon degauss", loc="left", fontsize=9); ax.legend(fontsize=7, loc="lower left"); path_axis(ax)
    for t, l in KV:
        Vt = np.load(f"results/epw/validation_{t}.npz", allow_pickle=True); sq = Vt["ph_s"]; F = Vt["ph_F_matdyn"]; iG = int(np.argmin(sq)); iK = int(np.argmin(np.abs(sq - TICKS[1]))); iM = int(np.argmin(np.abs(sq - TICKS[2])))
        print(f"[kohn {t} ({l})] matdyn omega(Gamma) = {F[iG, 4]:.2f}/{F[iG, 5]:.2f} cm^-1 ; omega(K), 6 modes = {np.round(F[iK], 2).tolist()} ; omega(M) LO/TO = {F[iM, 4]:.2f}/{F[iM, 5]:.2f} cm^-1")
    iK = int(np.argmin(np.abs(D["s"] - TICKS[1]))); iM = int(np.argmin(np.abs(D["s"] - TICKS[2])))
    print(f"[kohn dfpt {a.dfpt_tag}] omega(Gamma) = {D['freq'][0, 4]:.2f}/{D['freq'][0, 5]:.2f} ; omega(K), 6 modes = {np.round(D['freq'][iK], 2).tolist()} ; omega(M) = {D['freq'][iM, 4]:.2f}/{D['freq'][iM, 5]:.2f} cm^-1")
    fig.tight_layout(); save(fig, "fig_epw_kohn_degauss")
