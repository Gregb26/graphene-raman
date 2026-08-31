#!/usr/bin/env python
"""
compute_spectral.py
    From the electron-defect scattering matrix M, evaluate the on-shell single-defect
    scattering rate (inverse lifetime) of each Bloch state:

        Gamma_{nk} = -2 Im T_{nk,nk}(eps = eps_{nk})

    where T(eps) = [1 - M G0(eps)]^{-1} M (single_defect.compute_T). This is the defect-
    induced level width; tau_{nk} = hbar / Gamma_{nk} is the lifetime (per single defect,
    up to the impurity concentration). Saves Gamma_{nk}, tau_{nk} and a plot vs E - E_F.

    Usage: python scripts/compute_spectral.py --size 5x5 [--eta 0.02 --plot]
"""
import argparse
import numpy as np

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.defects.many_body.single_defect import compute_T

HA2EV = 27.211386245988
HBAR_EVS = 6.582119569e-16  # hbar in eV*s


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--size", default="5x5")
    p.add_argument("--eta", type=float, default=0.02, help="broadening in eV")
    p.add_argument("--chunk", type=int, default=16, help="on-shell energies per batch")
    p.add_argument("--plot", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    uc = f"data/graphene/unit_cell/qe/defect_{args.size}.save"
    M = np.load(f"results/M/M_ed_{args.size}.npy")               # (nb, nk, nb, nk)
    eigs = qe_io.get_eigenvalues(uc, shift_Fermi=True) * HA2EV   # (nb, nk) eV rel. E_F
    nb, nk = eigs.shape
    assert M.shape == (nb, nk, nb, nk)
    N = nb * nk
    eflat = eigs.reshape(N)                                      # on-shell energies, one per state

    gamma = np.zeros(N)
    for s in range(0, N, args.chunk):
        e = eflat[s:s + args.chunk]
        T = compute_T(e, eigs, M, args.eta)                     # (ne, nb, nk, nb, nk)
        Tmat = T.reshape(len(e), N, N)
        # on-shell diagonal: state a = global index s+j evaluated at its own energy e[j]
        for j in range(len(e)):
            a = s + j
            gamma[a] = -2.0 * Tmat[j, a, a].imag
        print(f"  states {s}-{min(s+args.chunk, N)}/{N} done", flush=True)

    gamma = gamma.reshape(nb, nk)                               # eV
    gamma_pos = np.clip(gamma, 1e-12, None)
    tau = HBAR_EVS / gamma_pos                                  # s (per single defect)

    out = f"results/M/gamma_{args.size}.npz"
    np.savez(out, gamma=gamma, tau=tau, eigs=eigs, eta=args.eta, size=args.size)
    print(f"saved {out}  (nb={nb}, nk={nk}, eta={args.eta} eV)")
    print(f"  Gamma range: {gamma.min():.2e} .. {gamma.max():.2e} eV   "
          f"median tau = {np.median(tau):.2e} s")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6, 4))
        sc = ax.scatter(eigs.ravel(), gamma.ravel() * 1e3, s=12, c=gamma.ravel() * 1e3,
                        cmap="viridis")
        ax.axvline(0, ls=":", c="grey")
        ax.set_xlabel(r"$E_{nk} - E_F$ (eV)"); ax.set_ylabel(r"$\Gamma_{nk}$ (meV)")
        ax.set_title(f"Single-defect scattering rate, {args.size} supercell")
        fig.colorbar(sc, label=r"$\Gamma$ (meV)")
        fig.tight_layout()
        png = f"results/M/gamma_{args.size}.png"
        fig.savefig(png, dpi=150)
        print(f"saved {png}")


if __name__ == "__main__":
    main()
