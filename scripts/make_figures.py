#!/usr/bin/env python
"""
Thesis figures (matplotlib, uniform style, data from results/M/*.npz only).
  1 fig_rcut.pdf        Gamma*N_cells vs R_cut, four sizes (frozen grid/eta)
  2 fig_plateau.pdf     (grid x eta) map at R_cut 3, reference size
  3 fig_level2.pdf      Gamma*N_cells vs N (frozen R_cut/grid/eta)
  4 fig_locality.pdf    on-site + ||Mwr(R,R0)|| decay: coarse NxN (aliased) vs dense
  5 fig_spectral.pdf    Gamma(eps) on-shell, Born vs T; Gamma/rho0 (∝|T|^2); rho0 vs rho_dis (c=1%); Tbar(K)
Also writes results/M/level1_summary.csv and level2_summary.csv.
Usage: python scripts/make_figures.py [--tag prod] [--sizes 5x5,7x7,8x8,9x9] [--outdir figures]
"""
import argparse, csv, os, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from electron_defect_interaction.config import load_production

# fixed categorical order (validated reference palette): size -> hue, never cycled
COL = {"5x5": "#2a78d6", "7x7": "#eb6834", "8x8": "#1baf7a", "9x9": "#eda100"}
C_T, C_BORN, C_REF, C_DIS = "#2a78d6", "#eb6834", "#52514e", "#1baf7a"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e3", "grid.linewidth": 0.6,
                     "lines.linewidth": 1.6, "lines.markersize": 5, "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
                     "legend.frameon": False, "axes.edgecolor": "#8a8984", "xtick.color": "#52514e", "ytick.color": "#52514e"})
ap = argparse.ArgumentParser(); ap.add_argument("--tag", default="prod"); ap.add_argument("--sizes", default="5x5,7x7,8x8,9x9"); ap.add_argument("--outdir", default="figures")
ap.add_argument("--resonance", default="results/M/resonance_9x9.npz"); ap.add_argument("--locality", default="results/M/mwr_locality.npz")
a = ap.parse_args(); cfg = load_production(); sizes = a.sizes.split(","); os.makedirs(a.outdir, exist_ok=True)
RC, GRID, ETA = cfg["R_cut"], cfg["grid"], cfg["eta_eV"]

def load_map(S):
    f = f"results/M/specwd_{S}_{a.tag}.npz"
    if not os.path.exists(f): return None
    r = np.load(f)["results"]; return {(int(x[0]), int(x[1]), round(float(x[2]), 4)): (float(x[3]), float(x[4])) for x in r}
maps = {S: load_map(S) for S in sizes}; maps = {S: m for S, m in maps.items() if m}
def save(fig, name):
    for ext in ("pdf", "png"): fig.savefig(f"{a.outdir}/{name}.{ext}")
    plt.close(fig); print(f"wrote {a.outdir}/{name}.pdf/.png")

if maps:
    rcs = sorted({k[0] for m in maps.values() for k in m}); grids = sorted({k[1] for m in maps.values() for k in m}); etas = sorted({k[2] for m in maps.values() for k in m}, reverse=True)
    # --- CSV summaries
    with open("results/M/level1_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["size", "R_cut", "grid", "eta_eV", "median_Gamma_Ncells_meV", "argmax_E_minus_ED_eV"])
        for S, m in maps.items():
            for (rc, N, e), (med, er) in sorted(m.items()): w.writerow([S, rc, N, e, f"{med:.4f}", f"{er:.4f}"])
    with open("results/M/level2_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["size", "N", "N_cells", "R_cut", "grid", "eta_eV", "median_Gamma_Ncells_meV"])
        for S, m in maps.items():
            N = int(S.split("x")[0]); w.writerow([S, N, N * N, RC, GRID, ETA, f"{m[(RC, GRID, ETA)][0]:.4f}"])
    print("wrote results/M/level1_summary.csv, level2_summary.csv")
    # --- 1. Gamma*N_cells vs R_cut
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    for S, m in maps.items():
        y = [m[(rc, GRID, ETA)][0] for rc in rcs]; ax.plot(rcs, y, "-o", color=COL[S], label=S)
    ax.set_xlabel(r"$R_\mathrm{cut}$ (lattice shells)"); ax.set_ylabel(r"median $\Gamma\,N_\mathrm{cells}$ (meV)"); ax.set_xticks(rcs)
    ax.set_title(rf"grid {GRID}$\times${GRID}, $\eta$ = {ETA} eV", loc="left"); ax.legend(title="supercell", ncol=2); ax.set_xlim(rcs[0] - 0.2, rcs[-1] + 0.2)
    save(fig, "fig_rcut")
    # --- 2. plateau map at R_cut, reference size (sequential single hue)
    S = cfg["reference_size"]
    if S in maps:
        Z = np.array([[maps[S][(RC, N, e)][0] for e in etas] for N in grids]); ref = Z[-1, etas.index(ETA)] if ETA in etas else Z[-1, -1]
        D = (Z / ref - 1) * 100
        fig, ax = plt.subplots(figsize=(3.4, 2.6)); im = ax.imshow(np.abs(D), cmap="Blues", vmin=0, vmax=max(5, np.abs(D).max()), origin="lower")
        ax.set_xticks(range(len(etas))); ax.set_xticklabels([f"{e:g}" for e in etas]); ax.set_yticks(range(len(grids))); ax.set_yticklabels([f"{N}$^2$" for N in grids]); ax.grid(False)
        for i in range(len(grids)):
            for j in range(len(etas)): ax.text(j, i, f"{Z[i,j]:.2f}\n({D[i,j]:+.1f}%)", ha="center", va="center", fontsize=7, color="#0b0b0b" if abs(D[i, j]) < 0.6 * max(5, np.abs(D).max()) else "white")
        ax.set_xlabel(r"$\eta$ (eV)"); ax.set_ylabel("output k-grid"); ax.set_title(rf"{S}, $R_\mathrm{{cut}}$ = {RC}: median $\Gamma N_\mathrm{{cells}}$ (meV), rel. to ({GRID}$^2$, {ETA})", loc="left", fontsize=8)
        cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("|deviation| (%)"); save(fig, "fig_plateau")
    # --- 3. Level 2: Gamma*N_cells vs N
    fig, ax = plt.subplots(figsize=(3.4, 2.6)); Ns = [int(S.split("x")[0]) for S in maps]; ys = [maps[S][(RC, GRID, ETA)][0] for S in maps]
    for S, N, y in zip(maps, Ns, ys): ax.plot([N], [y], "o", color=COL[S], label=S)
    ax.plot(Ns, ys, "-", color="#8a8984", lw=1, zorder=0); ax.axvspan(0, cfg["N_min"] - 0.5, color="#e6e6e3", alpha=0.5, lw=0)
    ax.set_xlabel(r"supercell $N$ ($N\times N$)"); ax.set_ylabel(r"median $\Gamma\,N_\mathrm{cells}$ (meV)"); ax.set_xticks(Ns); ax.set_xlim(min(Ns) - 1, max(Ns) + 1)
    ax.set_title(rf"$R_\mathrm{{cut}}$ = {RC}, grid {GRID}$^2$, $\eta$ = {ETA} eV; shaded: $N<N_\min$", loc="left", fontsize=8); ax.legend(ncol=2); save(fig, "fig_level2")

# --- 4. locality: coarse (aliased) vs dense
if os.path.exists(a.locality):
    L = np.load(a.locality); have = [S for S in sizes if f"{S}_dense_w" in L]
    fig, axs = plt.subplots(1, len(have), figsize=(1.75 * len(have) + 0.6, 2.8), sharey=True)
    axs = np.atleast_1d(axs)
    for ax, S in zip(axs, have):
        ax.semilogy(L[f"{S}_dense_dist"], L[f"{S}_dense_w"], "o", color=COL[S], ms=3.5, label="dense (zero-padded)", zorder=2)
        txt = f"on-site pz(A)\ndense {float(L[f'{S}_dense_onsite_pzA']):.2f} eV"
        if f"{S}_coarse_w" in L:
            ax.semilogy(L[f"{S}_coarse_dist"], L[f"{S}_coarse_w"], "x", color="#52514e", ms=4, label=f"{S} grid (aliased)", zorder=3)
            txt += f"\n{S} grid {float(L[f'{S}_coarse_onsite_pzA']):.2f} eV"
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, ha="right", va="top", fontsize=6.5, color="#0b0b0b")
        ax.set_title(S, loc="left"); ax.set_xlabel(r"$|R-R_0|$ ($a$)"); ax.legend(loc="center right", fontsize=6, handletextpad=0.3)
        ax.set_ylim(2e-5, 40)
    axs[0].set_ylabel(r"$\|M_{wR}(R,R_0)\|$ (eV)"); save(fig, "fig_locality")

# --- 5. spectral: Born vs T, |T|^2 proxy, DOS change, Tbar(K), det/eigenvalue criterion, Gamma at c=0.1%
if os.path.exists(a.resonance):
    R = np.load(a.resonance); ED = float(R["E_D"]); x = R["eg"] - ED; xg = R["egrid"] - ED; S = str(R["size"]); c = float(R["conc"])
    crit = a.resonance.replace("resonance_", "resonance_criteria_"); C = np.load(crit) if os.path.exists(crit) else None
    fig, axs = plt.subplots(3, 2, figsize=(6.8, 7.4)); (a1, a2), (a3, a4), (a5, a6) = axs
    a1.plot(x, R["Gamma_T"] * 1e3, color=C_T, label="T-matrix (exact)"); a1.plot(x, R["Gamma_Born"] * 1e3, color=C_BORN, label="Born (2nd order)")
    a1.set_ylabel(r"$\Gamma\,N_\mathrm{cells}$ (meV)"); a1.set_title(f"{S}: on-shell rate, $\\eta$ = {float(R['eta'])} eV", loc="left"); a1.legend()
    a2.plot(x, R["ratio"] / np.nanmax(R["ratio"]), color=C_T); a2.set_ylabel(r"$\Gamma/\rho_0$ (norm., $\propto|T|^2$)"); a2.set_title(r"(a) $\Gamma(\varepsilon)/\rho_0(\varepsilon)$", loc="left")
    a3.plot(x, R["rho0"], color=C_REF, label=r"$\rho_0$ (pristine, fine grid)"); a3.plot(x, R["rho_dis"], color=C_DIS, label=rf"$\rho_0 + c\,\delta\rho$, $c$ = {c*100:.0f}%")
    a3.set_ylabel("DOS (states / eV / cell / spin)"); a3.set_title(r"(b) $\delta\rho = \rho_\mathrm{dis}-\rho_0$", loc="left"); a3.legend()
    tr = R["Tbar_tr"]; a4.plot(xg, tr.real, color=C_T, label=r"Re $\bar T_{\pi\pi}(K,K;\varepsilon)$"); a4.plot(xg, tr.imag, color=C_BORN, label="Im"); a4.axhline(0, color="#8a8984", lw=0.8)
    a4.set_ylabel(r"$\bar T$ (eV, per defect)"); a4.set_title(r"(c) $\bar T$ at $K$, $\pi$ pair (trace/2)", loc="left"); a4.legend()
    if C is not None:
        xc = C["eg"] - float(C["E_D"])
        a5.semilogy(xc, np.exp(C["logdet_rel"]), color=C_T, label=r"$|\det[1-Vg_0]|$ / max"); a5.semilogy(xc, C["minlam"], color=C_BORN, label=r"$\min_i|\lambda_i(1-Vg_0)|$")
        a5.set_ylabel("dimensionless"); a5.set_title(r"(d) resonance criterion", loc="left"); a5.legend(fontsize=7)
        cc = float(C["c_compare"]); a6.plot(C["x_c"], C["Gamma_c"] * 1e3, color=C_T)
        a6.set_ylabel(rf"$\Gamma$ (meV) at $c$ = {cc*100:.1f}%"); a6.set_xlim(-1, 1)
        a6.set_title("(e) vacancy, $c$ = %.1f%% (vs Kaasbjerg Fig. 17: subst. N,\norder of magnitude only)" % (cc * 100), loc="left", fontsize=7.5)
    for ax in axs.ravel(): ax.set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax.axvline(0, color="#8a8984", lw=0.8, ls=":")
    fig.tight_layout(); save(fig, "fig_spectral")
