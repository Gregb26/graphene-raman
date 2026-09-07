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

from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.defects.many_body.single_defect import compute_T

HA2EV = 27.211386245988
HBAR_EVS = 6.582119569e-16  # hbar in eV*s


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--size", default="5x5")
    p.add_argument("--eta", type=float, default=None, help="broadening in eV; default = max(2x level spacing near E_F, 0.05)")
    p.add_argument("--chunk", type=int, default=16, help="on-shell energies per batch")
    p.add_argument("--plot", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    uc = f"data/graphene/unit_cell/qe/defect_{args.size}.save"
    # require the supercell-normalized M (matrix_io refuses an un-normalized/untagged one).
    M = matrix_io.load_M_checked(f"results/M/M_ed_{args.size}_norm.npy", units=matrix_io.EV)   # (nb, nk, nb, nk), eV
    # Structural k-pairing: eps aligned to M's k-grid (asserts count/order match, see qe_io).
    eigs = qe_io.aligned_eigenvalues(uc, nk_expected=M.shape[1], shift_Fermi=True) * HA2EV
    nb, nk = eigs.shape
    N = nb * nk
    N_cells = nk                                                 # one unit cell per k-point

    # eta must exceed ~2x the mean level spacing of the finite k-grid near E_F (else the on-shell
    # rate is a broadening artefact, not a proper continuum average).
    win = (np.abs(eigs) <= 3.0)
    level_spacing = 6.0 / max(1, int(win.sum()))                 # ~ 6 eV window / #levels
    eta = args.eta if args.eta is not None else max(2.0 * level_spacing, 0.05)
    if eta < 2.0 * level_spacing:
        print(f"WARNING: eta={eta:.3f} eV < 2x level spacing ({2*level_spacing:.3f} eV) for nk={nk}: "
              f"on-shell rate undersampled -- use a denser k-grid (bigger supercell) or larger eta.", flush=True)
    print(f"[{args.size}] nk={nk}, eta={eta:.4f} eV, level_spacing~{level_spacing:.4f} eV", flush=True)

    eflat = eigs.reshape(N)                                      # on-shell energies, one per state
    gamma = np.zeros(N)
    for s in range(0, N, args.chunk):
        e = eflat[s:s + args.chunk]
        T = compute_T(e, eigs, M, eta)                          # (ne, nb, nk, nb, nk)
        Tmat = T.reshape(len(e), N, N)
        # on-shell diagonal: state a = global index s+j evaluated at its own energy e[j]
        for j in range(len(e)):
            a = s + j
            gamma[a] = -2.0 * Tmat[j, a, a].imag
        print(f"  states {s}-{min(s+args.chunk, N)}/{N} done", flush=True)

    gamma = gamma.reshape(nb, nk)                               # eV, at the array concentration 1/N_cells
    gamma_perdef = gamma * N_cells                              # per-defect (concentration-normalised)
    gamma_pos = np.clip(gamma, 1e-12, None)
    tau = HBAR_EVS / gamma_pos                                  # s (at the array concentration)

    out = f"results/M/gamma_{args.size}.npz"
    np.savez(out, gamma=gamma, gamma_perdef=gamma_perdef, tau=tau,
             eigs=eigs, eta=eta, N_cells=N_cells, size=args.size)
    print(f"saved {out}  (nb={nb}, nk={nk}, eta={eta:.4f} eV)")
    print(f"  Gamma(array): median {np.median(np.abs(gamma))*1e3:.3f} meV | "
          f"per-defect Gamma*N_cells: median {np.median(np.abs(gamma_perdef))*1e3:.1f} meV", flush=True)

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
