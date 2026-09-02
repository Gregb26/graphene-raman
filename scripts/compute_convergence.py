#!/usr/bin/env python
"""
compute_convergence.py
    Capstone synthesis in the SUPERCELL Bloch-normalization convention. The size-independent
    (per-defect) observables are tracked vs supercell size N (defect concentration 1/N_cells):

        * median per-defect scattering rate  Gamma*N_cells  (meV),
        * vacancy-resonance position near the Dirac point  (argmax |Delta-DOS| for |E-E_F|<E_win).

    Interpretation criteria (guardrail 3):
        (a) Gamma*N_cells plateaus between 7x7 and 8x8 (relative change <= ~5-10%),
        (b) the resonance position is stable between those two sizes.
    If (a) fails but (b) holds, eta is likely still too large for the big cells -> rerun the two
    largest sizes at 2-3 eta values before concluding (this script only flags it).

    Usage: python scripts/compute_convergence.py [--plot --ewin 1.5]
"""
import argparse
import glob
import re

import numpy as np

SIZE_RE = re.compile(r"_(\d+)x(\d+)\.npz$")


def collect(prefix):
    out = {}
    for f in sorted(glob.glob(f"results/M/{prefix}_*x*.npz")):
        m = SIZE_RE.search(f)
        if m:
            out[int(m.group(1))] = np.load(f)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ewin", type=float, default=1.5, help="|E-E_F| window (eV) to locate the resonance")
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    dos = collect("dos")
    gam = collect("gamma")
    Ns = sorted(set(dos) & set(gam))
    if not Ns:
        raise SystemExit("need matching dos_*.npz and gamma_*.npz (run compute_tmatrix + compute_spectral)")

    rows = {}
    for N in Ns:
        g = gam[N]
        gperdef = g["gamma_perdef"] if "gamma_perdef" in g else g["gamma"] * (N * N)
        med = float(np.median(np.abs(gperdef))) * 1e3                     # meV
        d = dos[N]
        eps, ddos = d["eps"], d["ddos"]
        win = np.abs(eps) <= args.ewin
        e_res = float(eps[win][np.argmax(np.abs(ddos[win]))]) if win.any() else np.nan
        eta = float(g["eta"]) if "eta" in g else np.nan
        rows[N] = dict(gperdef=med, e_res=e_res, eta=eta)

    print(f"{'N':>4} {'eta(eV)':>9} {'median G*Ncells(meV)':>22} {'resonance E-EF(eV)':>20}")
    for N in Ns:
        r = rows[N]
        print(f"{N:>4} {r['eta']:>9.4f} {r['gperdef']:>22.2f} {r['e_res']:>20.3f}")

    # guardrail-3 interpretation on 7x7 vs 8x8 if available
    if 7 in rows and 8 in rows:
        g7, g8 = rows[7]["gperdef"], rows[8]["gperdef"]
        plateau = abs(g8 - g7) / max(1e-30, abs(g7))
        dres = abs(rows[8]["e_res"] - rows[7]["e_res"])
        print(f"\n[7x7 -> 8x8] G*Ncells change = {plateau*100:.1f}%  |  resonance shift = {dres*1000:.0f} meV")
        oka = plateau <= 0.10
        okb = dres <= 0.05
        print(f"  (a) plateau (<=10%): {'PASS' if oka else 'FAIL'}   (b) resonance stable (<=50 meV): {'PASS' if okb else 'FAIL'}")
        if not oka and okb:
            print("  -> (a) fails while (b) holds: eta likely still too large for large cells. "
                  "Rerun 7x7/8x8 at 2-3 eta values (compute_spectral.py --eta ...) before concluding.")
    else:
        print("\n(need both 7x7 and 8x8 for the plateau/resonance interpretation)")

    np.savez("results/M/convergence.npz",
             N=np.array(Ns),
             gperdef_meV=np.array([rows[N]["gperdef"] for N in Ns]),
             e_res=np.array([rows[N]["e_res"] for N in Ns]),
             eta=np.array([rows[N]["eta"] for N in Ns]))
    print("saved results/M/convergence.npz")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cmap = plt.get_cmap("viridis")
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
        for i, N in enumerate(Ns):
            ax[0].plot(dos[N]["eps"], dos[N]["ddos"], color=cmap(i / max(1, len(Ns) - 1)), label=f"{N}x{N}")
        ax[0].axvline(0, ls=":", c="grey"); ax[0].set_xlim(-args.ewin*2, args.ewin*2)
        ax[0].set_xlabel(r"$E-E_F$ (eV)"); ax[0].set_ylabel(r"$\Delta$DOS (1/eV)")
        ax[0].set_title("Defect DOS near E_Dirac"); ax[0].legend(fontsize=8, ncol=2)
        ax[1].plot(Ns, [rows[N]["gperdef"] for N in Ns], "o-")
        ax[1].set_xlabel("N"); ax[1].set_ylabel(r"median $\Gamma\cdot N_{cells}$ (meV)")
        ax[1].set_title("Per-defect rate (should plateau)")
        ax[2].plot(Ns, [rows[N]["e_res"] for N in Ns], "s-", color="tab:red")
        ax[2].set_xlabel("N"); ax[2].set_ylabel(r"resonance $E-E_F$ (eV)")
        ax[2].set_title("Vacancy resonance position")
        fig.tight_layout(); fig.savefig("results/M/convergence.png", dpi=150)
        print("saved results/M/convergence.png")


if __name__ == "__main__":
    main()
