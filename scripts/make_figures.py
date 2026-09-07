#!/usr/bin/env python
"""
Figures du mémoire (matplotlib, style figures/memoire.mplstyle, texte en français, données lues uniquement dans results/M/*.npz).
  fig_rcut          Γ N_cells vs R_cut, quatre tailles (grille et η gelés)
  fig_plateau       carte (grille × η) à R_cut gelé, taille de référence
  fig_level2        Γ N_cells vs N (niveau 2)
  fig_locality      sur-site pz–pz et décroissance ‖M_wR(R,R0)‖ : grille N×N (aliasée) vs dense
  fig_spectral      Γ(ε) sur-couche (Born vs T), Γ/ρ0, ρ0 vs ρ_dis, T̄(K), critère det/valeur propre, Γ à c = 0,1 %
  fig_M_map         |M_nn(k', k=K)| sur la zone de Brillouin (Ṽ = A_cell |M|, eV Å²), π et π*
  fig_Ved_boundary  V_ed^L du site de la lacune à la frontière de la super-cellule, quatre N
  fig_M_scaling     max|M| vs N et max|M| N_cells vs N (convention intensive)
Écrit aussi results/M/level1_summary.csv et level2_summary.csv.
Usage : python scripts/make_figures.py [--tag prod] [--sizes 5x5,7x7,8x8,9x9] [--outdir figures]
"""
import argparse, csv, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from electron_defect_interaction.config import load_production

plt.style.use("figures/memoire.mplstyle")
# palette catégorielle fixe (ordre du cycler) : taille -> teinte, jamais recyclée
COL = {"5x5": "#2a78d6", "7x7": "#eb6834", "8x8": "#1baf7a", "9x9": "#eda100"}
C_T, C_BORN, C_REF, C_DIS = "#2a78d6", "#eb6834", "#52514e", "#1baf7a"
INK = "#0b0b0b"; MUTED = "#8a8984"
LBL_E = r"Énergie $\varepsilon - E_D$ (eV)"
LBL_G = r"Taux d'amortissement $\Gamma\,N_\mathrm{cells}$ (meV)"

ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="prod"); ap.add_argument("--sizes", default="5x5,7x7,8x8,9x9"); ap.add_argument("--outdir", default="figures")
ap.add_argument("--resonance", default="results/M/resonance_9x9.npz"); ap.add_argument("--locality", default="results/M/mwr_locality.npz")
ap.add_argument("--analysis", default="results/M/M_analysis.npz")
a = ap.parse_args(); cfg = load_production(); sizes = a.sizes.split(","); os.makedirs(a.outdir, exist_ok=True)
RC, GRID, ETA = cfg["R_cut"], cfg["grid"], cfg["eta_eV"]
def lab(S): return S.replace("x", r"$\times$")
def panel(ax, letter):
    t = ax.get_title(loc="left"); ax.set_title(f"({letter}) {t}" if t else f"({letter})", loc="left", fontsize=ax.title.get_fontsize())
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig); print(f"écrit {a.outdir}/{name}.pdf/.png")

def load_map(S):
    f = f"results/M/specwd_{S}_{a.tag}.npz"
    if not os.path.exists(f): return None
    r = np.load(f)["results"]; return {(int(x[0]), int(x[1]), round(float(x[2]), 4)): (float(x[3]), float(x[4])) for x in r}
maps = {S: m for S in sizes if (m := load_map(S))}

if maps:
    rcs = sorted({k[0] for m in maps.values() for k in m}); grids = sorted({k[1] for m in maps.values() for k in m}); etas = sorted({k[2] for m in maps.values() for k in m}, reverse=True)
    with open("results/M/level1_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["size", "R_cut", "grid", "eta_eV", "median_Gamma_Ncells_meV", "argmax_E_minus_ED_eV"])
        for S, m in maps.items():
            for (rc, N, e), (med, er) in sorted(m.items()): w.writerow([S, rc, N, e, f"{med:.4f}", f"{er:.4f}"])
    with open("results/M/level2_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["size", "N", "N_cells", "R_cut", "grid", "eta_eV", "median_Gamma_Ncells_meV"])
        for S, m in maps.items():
            N = int(S.split("x")[0]); w.writerow([S, N, N * N, RC, GRID, ETA, f"{m[(RC, GRID, ETA)][0]:.4f}"])
    print("écrit results/M/level1_summary.csv, level2_summary.csv")

    # ---- fig_rcut
    fig, ax = plt.subplots()
    for S, m in maps.items():
        ax.plot(rcs, [m[(rc, GRID, ETA)][0] for rc in rcs], "-o", color=COL[S], label=lab(S))
    ax.set_xlabel(r"Rayon de coupure $R_\mathrm{cut}$ (mailles)"); ax.set_ylabel(LBL_G); ax.set_xticks(rcs)
    ax.set_title(rf"Grille $k$ {GRID}$\times${GRID}, $\eta$ = {ETA} eV", loc="left"); ax.legend(title="Super-cellule", ncol=2)
    save(fig, "fig_rcut")

    # ---- fig_plateau (séquentielle, une seule teinte)
    S = cfg["reference_size"]
    if S in maps:
        Z = np.array([[maps[S][(RC, N, e)][0] for e in etas] for N in grids]); ref = Z[-1, etas.index(ETA)]; D = (Z / ref - 1) * 100; vmax = max(5, np.abs(D).max())
        fig, ax = plt.subplots(figsize=(6.0, 4.4)); im = ax.imshow(np.abs(D), cmap="Blues", vmin=0, vmax=vmax, origin="lower")
        ax.set_xticks(range(len(etas))); ax.set_xticklabels([f"{e:g}" for e in etas]); ax.set_yticks(range(len(grids))); ax.set_yticklabels([rf"{N}$^2$" for N in grids]); ax.grid(False)
        for i in range(len(grids)):
            for j in range(len(etas)):
                ax.text(j, i, f"{Z[i,j]:.0f}\n({D[i,j]:+.1f}\\,\\%)", ha="center", va="center", fontsize=9, color=INK if abs(D[i, j]) < 0.6 * vmax else "white")
        ax.set_xlabel(r"Élargissement $\eta$ (eV)"); ax.set_ylabel(r"Grille $k$ de sortie")
        ax.set_title(rf"{lab(S)}, $R_\mathrm{{cut}}$ = {RC} : médiane de $\Gamma N_\mathrm{{cells}}$ (meV)", loc="left", fontsize=8)
        cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label(rf"Écart à ({GRID}$^2$, {ETA} eV) (\%)"); save(fig, "fig_plateau")

    # ---- fig_level2
    fig, ax = plt.subplots(); Ns = [int(S.split("x")[0]) for S in maps]; ys = [maps[S][(RC, GRID, ETA)][0] for S in maps]
    for S, N, y in zip(maps, Ns, ys): ax.plot([N], [y], "o", color=COL[S], label=lab(S))
    ax.plot(Ns, ys, "-", color=MUTED, lw=1, zorder=0); ax.axvspan(0, cfg["N_min"] - 0.5, color="#e6e6e3", alpha=0.5, lw=0)
    ax.set_xlabel(r"Taille de la super-cellule $N$ ($N\times N$)"); ax.set_ylabel(LBL_G); ax.set_xticks(Ns); ax.set_xlim(min(Ns) - 1, max(Ns) + 1)
    ax.set_title(rf"$R_\mathrm{{cut}}$ = {RC}, grille {GRID}$^2$, $\eta$ = {ETA} eV ; zone grise : $N<N_\mathrm{{min}}$", loc="left", fontsize=8); ax.legend(ncol=2); save(fig, "fig_level2")

# ---- fig_locality : dense (quatre N) vs grille N×N aliasée (trois N), sur-site pz–pz en légende ; 2×2 panneaux
if os.path.exists(a.locality):
    L = np.load(a.locality); have = [S for S in sizes if f"{S}_dense_w" in L]
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 5.6), sharey=True, sharex=True); axs = axs.ravel()
    for i, (ax, S) in enumerate(zip(axs, have)):
        ax.semilogy(L[f"{S}_dense_dist"], L[f"{S}_dense_w"], "o", color=COL[S], ms=3.2, zorder=2,
                    label=f"dense (zero-padding)\npz–pz = {float(L[f'{S}_dense_onsite_pzA']):.2f} eV")
        if f"{S}_coarse_w" in L:
            ax.semilogy(L[f"{S}_coarse_dist"], L[f"{S}_coarse_w"], "x", color=C_REF, ms=4, zorder=3,
                        label=f"grille {lab(S)} (aliasée)\npz–pz = {float(L[f'{S}_coarse_onsite_pzA']):.2f} eV")
        ax.set_title(lab(S), loc="left"); ax.legend(loc="upper right", fontsize=6.5, handletextpad=0.4, labelspacing=0.8); ax.set_ylim(2e-5, 60); panel(ax, "abcd"[i])
    for ax in axs[2:]: ax.set_xlabel(r"Distance $|R-R_0|$ ($a$)")
    for ax in axs[::2]: ax.set_ylabel(r"$\|M_{wR}(R,R_0)\|$ (eV)")
    fig.tight_layout(); save(fig, "fig_locality")

# ---- fig_spectral
if os.path.exists(a.resonance):
    R = np.load(a.resonance); ED = float(R["E_D"]); x = R["eg"] - ED; xg = R["egrid"] - ED; S = str(R["size"]); c = float(R["conc"])
    crit = a.resonance.replace("resonance_", "resonance_criteria_"); C = np.load(crit) if os.path.exists(crit) else None
    fig, axs = plt.subplots(3, 2, figsize=(6.5, 9.0)); (a1, a2), (a3, a4), (a5, a6) = axs
    a1.plot(x, R["Gamma_T"] * 1e3, color=C_T, label=r"matrice $T$ (exacte)"); a1.plot(x, R["Gamma_Born"] * 1e3, color=C_BORN, label="approximation de Born (2$^\\mathrm{e}$ ordre)")
    a1.set_ylabel(LBL_G); a1.set_title(rf"{lab(S)} : taux sur-couche, $\eta$ = {float(R['eta'])} eV", loc="left"); a1.legend()
    a2.plot(x, R["ratio"] / np.nanmax(R["ratio"]), color=C_T); a2.set_ylabel(r"$\Gamma/\rho_0$ (normalisé, $\propto|T|^2$)"); a2.set_title(r"$\Gamma(\varepsilon)/\rho_0(\varepsilon)$", loc="left")
    a3.plot(x, R["rho0"], color=C_REF, label=r"$\rho_0$ (cristal parfait, grille fine)"); a3.plot(x, R["rho_dis"], color=C_DIS, label=rf"$\rho_0 + c\,\delta\rho$, $c$ = {c*100:.0f}\,\%")
    a3.set_ylabel("Densité d'états (états/eV/cellule)"); a3.set_title(r"$\delta\rho = \rho_\mathrm{dis}-\rho_0$", loc="left"); a3.legend()
    tr = R["Tbar_tr"]; a4.plot(xg, tr.real, color=C_T, label=r"Re $\bar T_{\pi\pi}(K,K;\varepsilon)$"); a4.plot(xg, tr.imag, color=C_BORN, label="Im"); a4.axhline(0, color=MUTED, lw=0.8)
    a4.set_ylabel(r"$\bar T$ (eV, par défaut)"); a4.set_title(r"$\bar T$ en $K$, paire $\pi$ (trace/2)", loc="left"); a4.legend()
    if C is not None:
        xc = C["eg"] - float(C["E_D"])
        a5.semilogy(xc, np.exp(C["logdet_rel"]), color=C_T, label=r"$|\det[1-Vg_0]|$ / max"); a5.semilogy(xc, C["minlam"], color=C_BORN, label=r"$\min_i|\lambda_i(1-Vg_0)|$")
        a5.set_ylabel("Sans dimension"); a5.set_title("Critère de résonance", loc="left"); a5.legend(fontsize=7)
        cc = float(C["c_compare"]); a6.plot(C["x_c"], C["Gamma_c"] * 1e3, color=C_T)
        a6.set_ylabel(rf"$\Gamma$ (meV) à $c$ = {cc*100:.1f}\,\%"); a6.set_xlim(-1, 1)
        a6.set_title("Lacune, $c$ = %.1f\\,\\%% (vs Kaasbjerg fig. 17 : N substitutionnel,\nordre de grandeur seulement)" % (cc * 100), loc="left", fontsize=7.5)
    for i, ax in enumerate(axs.ravel()): ax.set_xlabel(LBL_E); ax.axvline(0, color=MUTED, lw=0.8, ls=":"); panel(ax, "abcdef"[i])
    fig.tight_layout(); save(fig, "fig_spectral")

# ---- figures de la matrice M (§4.1.5)
if os.path.exists(a.analysis):
    Z = np.load(a.analysis, allow_pickle=True)
    from scipy.spatial import Voronoi
    def bz_vertices(B):
        """Première zone de Brillouin (cellule de Wigner–Seitz réciproque) : sommets de la région de Voronoï de Γ."""
        pts = np.array([i * B[:2, 0] + j * B[:2, 1] for i in range(-2, 3) for j in range(-2, 3)]); vor = Voronoi(pts)
        reg = vor.regions[vor.point_region[np.argmin(np.linalg.norm(pts, axis=1))]]; v = vor.vertices[reg]
        return v[np.argsort(np.arctan2(v[:, 1], v[:, 0]))]
    # fig_M_map : (a) π, (b) π* : Ṽ = A_cell|M| ; (c),(d) composantes de M^L et M^NL le long de M total (ligne π), même échelle
    fig, axs = plt.subplots(2, 2, figsize=(6.5, 6.0), sharex=True, sharey=True, layout="constrained"); axs = axs.ravel()
    vmax = max(Z["map_Vpi"].max(), Z["map_Vpistar"].max()); hexa = bz_vertices(Z["map_B"])
    lo = min(0.0, Z["map_Lpar"].min(), Z["map_Npar"].min()); hi = max(Z["map_Lpar"].max(), Z["map_Npar"].max())
    specs = ((0, "map_Vpi", r"$n' = \pi(k')$ : $A_\mathrm{cell}|M|$", 0, vmax), (1, "map_Vpistar", r"$n' = \pi^*(k')$ : $A_\mathrm{cell}|M|$", 0, vmax),
             (2, "map_Lpar", r"$\pi$ : $M^L$ le long de $M$", lo, hi), (3, "map_Npar", r"$\pi$ : $M^{NL}$ le long de $M$", lo, hi))
    for i, key, title, vlo, vhi in specs:
        ax = axs[i]; sc = ax.scatter(Z["map_kx"], Z["map_ky"], c=Z[key], cmap="Blues", vmin=vlo, vmax=vhi, s=22, marker="h", linewidths=0)
        ax.add_patch(Polygon(hexa, closed=True, fill=False, ec=MUTED, lw=0.8)); ax.plot(*Z["map_K"], "o", mfc="none", mec=C_BORN, mew=1.2, ms=8)
        ax.annotate("$K$", Z["map_K"], xytext=(6, 4), textcoords="offset points", color=C_BORN, fontsize=10)
        ax.set_aspect("equal"); ax.set_title(title, loc="left", fontsize=10); panel(ax, "abcd"[i])
        if i >= 2: ax.set_xlabel(r"$k'_x$ (Å$^{-1}$)")
        if i % 2 == 0: ax.set_ylabel(r"$k'_y$ (Å$^{-1}$)")
        if i == 1: cb = fig.colorbar(sc, ax=axs[:2], shrink=0.85, pad=0.02); cb.set_label(r"$\tilde V$ (eV\,Å$^2$)")
        if i == 3: cb = fig.colorbar(sc, ax=axs[2:], shrink=0.85, pad=0.02); cb.set_label(r"$A_\mathrm{cell}\,\mathrm{Re}[\hat M^\dagger M^{X}]$ (eV\,Å$^2$)")
    save(fig, "fig_M_map")
    # fig_Ved_boundary : quatre N, axe normalisé
    fig, ax = plt.subplots()
    for S in sizes:
        if f"ved_{S}_line" in Z: ax.plot(Z[f"ved_{S}_x"], Z[f"ved_{S}_line"], color=COL[S], label=lab(S))
    ax.axhline(0, color=MUTED, lw=0.8); ax.set_xlabel(r"Distance au site / demi-largeur de la boîte (le long de $\mathbf a_1$)"); ax.set_ylabel(r"$V_\mathrm{ed}^L$ (eV)")
    ax.set_title(r"Potentiel de défaut local, du site de la lacune à la frontière", loc="left", fontsize=8); ax.legend(title="Super-cellule", ncol=2)
    ax.set_yscale("symlog", linthresh=0.01, linscale=0.4); ax.set_yticks([100, 10, 1, 0.1, 0.01, 0, -0.01, -0.1, -1]); ax.set_yticklabels(["$10^{2}$", "$10^{1}$", "$1$", "$10^{-1}$", "$10^{-2}$", "$0$", "$-10^{-2}$", "$-10^{-1}$", "$-1$"])
    ax.set_ylim(-2, 200); save(fig, "fig_Ved_boundary")
    # fig_M_scaling : deux panneaux côte à côte
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.2))
    Ns = [int(S[0]) for S in sizes if f"scale_{S}_dense" in Z]; dd = [float(Z[f"scale_{S}_dense"]) for S in sizes if f"scale_{S}_dense" in Z]; cc_ = [float(Z[f"scale_{S}_coarse"]) for S in sizes if f"scale_{S}_dense" in Z]
    nc = np.array(Ns) ** 2
    axs[0].plot(Ns, dd, "o-", color=C_T, label="dense (zero-padding)"); axs[0].plot(Ns, cc_, "x--", color=C_REF, label=r"grille $N\times N$")
    axs[0].set_ylabel(r"$\max|M|$ (eV)"); axs[0].set_ylim(0, max(dd + cc_) * 1.25)
    axs[1].plot(Ns, np.array(dd) * nc, "o-", color=C_T, label="dense (zero-padding)"); axs[1].plot(Ns, np.array(cc_) * nc, "x--", color=C_REF, label=r"grille $N\times N$")
    axs[1].set_ylabel(r"$\max|M|\cdot N_\mathrm{cells}$ (eV)")
    axs[0].set_title("Convention cellule unitaire : intensif", loc="left", fontsize=8); axs[1].set_title(r"Convention super-cellule $\times N_\mathrm{cells}$ : extensif", loc="left", fontsize=8)
    for i, ax in enumerate(axs): ax.set_xlabel(r"Taille de la super-cellule $N$"); ax.set_xticks(Ns); ax.legend(fontsize=7); panel(ax, "ab"[i])
    fig.tight_layout(); save(fig, "fig_M_scaling")
