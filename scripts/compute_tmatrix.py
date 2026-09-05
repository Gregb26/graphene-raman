#!/usr/bin/env python
"""
compute_tmatrix.py
    From a computed electron-defect scattering matrix M = M^L + M^NL, evaluate the
    single-defect T-matrix and interacting Green's function (single_defect.py) and
    extract the density of states:

        A0(e) = -1/pi sum_{n,k} Im G0_{nk}(e)         (pristine)
        A(e)  = -1/pi sum_{n,k} Im G_{nk,nk}(e)       (with the defect, via T-matrix)

    The difference A - A0 is the defect-induced change in the DOS (bound states /
    resonances). Energies are relative to the Fermi level (eV).

    Usage:
        python scripts/compute_tmatrix.py --size 7x7 [--ne 400 --eta 0.02 --plot]
"""
import argparse
import numpy as np

from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.defects.many_body.single_defect import compute_G0, compute_G

HA2EV = 27.211386245988


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--size", default="7x7", help="supercell size, e.g. 7x7 (needs results/M/M_ed_<size>.npy)")
    p.add_argument("--ne", type=int, default=600, help="minimum number of energy points (raised if needed for dE<=eta/4)")
    p.add_argument("--eta", type=float, default=None, help="broadening in eV; default = max(2x level spacing, 0.05)")
    p.add_argument("--emin", type=float, default=-3.0, help="min energy in eV rel. E_F (tight window near E_F)")
    p.add_argument("--emax", type=float, default=3.0, help="max energy in eV rel. E_F")
    p.add_argument("--chunk", type=int, default=10, help="energies per batch (memory)")
    p.add_argument("--plot", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    uc = f"data/graphene/unit_cell/qe/defect_{args.size}.save"
    # require the supercell-normalized M (matrix_io refuses an un-normalized/untagged one).
    M = matrix_io.load_M_checked(f"results/M/M_ed_{args.size}_norm.npy") * HA2EV   # (nb, nk, nb, nk); Ha -> eV to match eigs
    # Structural k-pairing: eps aligned to M's k-grid; asserts the k-count/order match (see qe_io).
    eigs = qe_io.aligned_eigenvalues(uc, nk_expected=M.shape[1], shift_Fermi=True) * HA2EV
    nb, nk = eigs.shape

    emin, emax = args.emin, args.emax
    # mean level spacing of the finite k-grid inside the window: eta must exceed it or the DOS is
    # just a comb of delta-like spikes rather than a smooth spectrum.
    in_win = (eigs >= emin) & (eigs <= emax)
    n_lev = max(1, int(in_win.sum()))
    level_spacing = (emax - emin) / n_lev
    eta = args.eta if args.eta is not None else max(2.0 * level_spacing, 0.05)

    # energy grid fine enough to resolve the Lorentzians: dE <= eta/4
    ne = max(args.ne, int(np.ceil(4.0 * (emax - emin) / eta)))
    eps = np.linspace(emin, emax, ne)
    dE = eps[1] - eps[0]
    if eta < 2.0 * level_spacing:
        print(f"WARNING: eta={eta:.3f} eV < 2x mean level spacing ({2*level_spacing:.3f} eV) for "
              f"nk={nk}: DOS undersampled -- use a larger supercell (denser k-grid) or larger eta.", flush=True)
    if eta < dE:
        print(f"WARNING: eta={eta:.3f} eV < energy grid spacing dE={dE:.3f} eV.", flush=True)
    print(f"[{args.size}] window [{emin},{emax}] eV, ne={ne}, dE={dE:.4f}, eta={eta:.4f}, "
          f"level_spacing={level_spacing:.4f} eV", flush=True)

    dos0 = np.zeros(ne)
    dos = np.zeros(ne)
    for s in range(0, ne, args.chunk):
        e = eps[s:s + args.chunk]
        G0 = compute_G0(e, eigs, eta)                             # (ne, nb, nk)
        dos0[s:s + args.chunk] = -1.0 / np.pi * np.sum(G0.imag, axis=(1, 2))
        G = compute_G(e, eigs, M, eta)                           # (ne, nb, nk, nb, nk)
        Gdiag = np.einsum("eabab->eab", G)                       # (ne, nb, nk)
        dos[s:s + args.chunk] = -1.0 / np.pi * np.sum(Gdiag.imag, axis=(1, 2))
        print(f"  energies {s}-{min(s+args.chunk, ne)}/{ne} done", flush=True)

    out = f"results/M/dos_{args.size}.npz"
    np.savez(out, eps=eps, dos0=dos0, dos=dos, ddos=dos - dos0, eta=eta, size=args.size)
    print(f"saved {out}  (size={args.size}, nb={nb}, nk={nk}, eta={eta} eV)")

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
        ax[0].plot(eps, dos0, label="pristine A0", color="k")
        ax[0].plot(eps, dos, label="with defect A", color="r")
        ax[0].axvline(0, ls=":", c="grey"); ax[0].set_ylabel("DOS (1/eV)"); ax[0].legend()
        ax[1].plot(eps, dos - dos0, color="b"); ax[1].axhline(0, ls=":", c="grey")
        ax[1].axvline(0, ls=":", c="grey")
        ax[1].set_ylabel(r"$\Delta$DOS"); ax[1].set_xlabel(r"$E - E_F$ (eV)")
        fig.suptitle(f"Single-defect DOS, {args.size} supercell")
        fig.tight_layout()
        png = f"results/M/dos_{args.size}.png"
        fig.savefig(png, dpi=150)
        print(f"saved {png}")


if __name__ == "__main__":
    main()
