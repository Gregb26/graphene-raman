#!/usr/bin/env python
"""
compute_convergence.py
    Capstone synthesis: how the single-defect observables converge towards the dilute
    (isolated-defect) limit as the supercell grows. The supercell size N sets the defect
    concentration c = 1/N^2 (one defect per N x N unit cells), so N -> infinity is dilute.

    Collects the per-size outputs produced by compute_tmatrix.py (dos_<N>.npz) and
    compute_spectral.py (gamma_<N>.npz) and plots:
        (1) the defect-induced Delta-DOS(E) for every size (spectral signature),
        (2) scalar convergence metrics vs N:  int |Delta-DOS| dE  and  median Gamma.

    Usage: python scripts/compute_convergence.py [--plot]
"""
import argparse
import glob
import os
import re

import numpy as np

SIZE_RE = re.compile(r"_(\d+)x(\d+)\.npz$")


def collect(prefix):
    out = {}
    for f in sorted(glob.glob(f"results/M/{prefix}_*x*.npz")):
        m = SIZE_RE.search(f)
        if not m:
            continue
        N = int(m.group(1))
        out[N] = np.load(f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    dos = collect("dos")
    gam = collect("gamma")
    Ns = sorted(dos)
    if not Ns:
        raise SystemExit("no dos_*.npz found; run compute_tmatrix.py first")

    # (1) convergence metric from Delta-DOS: total rearranged spectral weight
    weight = {N: float(np.trapz(np.abs(dos[N]["ddos"]), dos[N]["eps"])) for N in Ns}
    # (2) scattering-rate summary (median over states), meV
    Ngam = sorted(gam)
    gmed = {N: float(np.median(np.abs(gam[N]["gamma"]))) * 1e3 for N in Ngam}

    print(f"{'N':>4} {'conc 1/N^2':>12} {'int|dDOS|':>12} {'median Gamma(meV)':>18}")
    for N in Ns:
        g = f"{gmed[N]:.3f}" if N in gmed else "-"
        print(f"{N:>4} {1.0/N**2:>12.4e} {weight[N]:>12.4f} {g:>18}")

    np.savez("results/M/convergence.npz",
             N=np.array(Ns), conc=np.array([1.0 / N**2 for N in Ns]),
             ddos_weight=np.array([weight[N] for N in Ns]),
             gamma_N=np.array(Ngam), gamma_median_meV=np.array([gmed[N] for N in Ngam]))
    print("saved results/M/convergence.npz")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("viridis")
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))

        for i, N in enumerate(Ns):
            c = cmap(i / max(1, len(Ns) - 1))
            ax[0].plot(dos[N]["eps"], dos[N]["ddos"], color=c, lw=1.2, label=f"{N}x{N}")
        ax[0].axvline(0, ls=":", c="grey"); ax[0].axhline(0, ls=":", c="grey")
        ax[0].set_xlabel(r"$E - E_F$ (eV)"); ax[0].set_ylabel(r"$\Delta$DOS (1/eV)")
        ax[0].set_title("Defect-induced DOS vs supercell size"); ax[0].legend(fontsize=8, ncol=2)

        axb = ax[1]
        axb.plot(Ns, [weight[N] for N in Ns], "o-", color="tab:blue", label=r"$\int|\Delta$DOS$|\,dE$")
        axb.set_xlabel("supercell size N"); axb.set_ylabel(r"$\int|\Delta$DOS$|\,dE$", color="tab:blue")
        axb.tick_params(axis="y", labelcolor="tab:blue")
        if Ngam:
            axr = axb.twinx()
            axr.plot(Ngam, [gmed[N] for N in Ngam], "s--", color="tab:red", label=r"median $\Gamma$")
            axr.set_ylabel(r"median $\Gamma$ (meV)", color="tab:red")
            axr.tick_params(axis="y", labelcolor="tab:red")
        axb.set_title("Convergence to the dilute limit")
        fig.tight_layout()
        fig.savefig("results/M/convergence.png", dpi=150)
        print("saved results/M/convergence.png")


if __name__ == "__main__":
    main()
