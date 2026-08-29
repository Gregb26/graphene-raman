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

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.defects.many_body.single_defect import compute_G0, compute_G

HA2EV = 27.211386245988


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--size", default="7x7", help="supercell size, e.g. 7x7 (needs results/M/M_ed_<size>.npy)")
    p.add_argument("--ne", type=int, default=400, help="number of energy points")
    p.add_argument("--eta", type=float, default=0.02, help="broadening in eV")
    p.add_argument("--emin", type=float, default=None, help="min energy in eV (rel. to E_F); default from bands")
    p.add_argument("--emax", type=float, default=None, help="max energy in eV; default from bands")
    p.add_argument("--chunk", type=int, default=10, help="energies per batch (memory)")
    p.add_argument("--plot", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    uc = f"data/graphene/unit_cell/qe/defect_{args.size}.save"
    M = np.load(f"results/M/M_ed_{args.size}.npy")                 # (nb, nk, nb, nk)
    eigs_ha = qe_io.get_eigenvalues(uc, shift_Fermi=True)          # (nb, nk), Hartree rel. E_F
    eigs = eigs_ha * HA2EV                                         # eV
    nb, nk = eigs.shape
    assert M.shape == (nb, nk, nb, nk), f"M {M.shape} vs eigs {eigs.shape}"

    eta = args.eta
    emin = args.emin if args.emin is not None else float(eigs.min()) - 1.0
    emax = args.emax if args.emax is not None else float(eigs.max()) + 1.0
    eps = np.linspace(emin, emax, args.ne)

    dos0 = np.zeros(args.ne)
    dos = np.zeros(args.ne)
    for s in range(0, args.ne, args.chunk):
        e = eps[s:s + args.chunk]
        G0 = compute_G0(e, eigs, eta)                             # (ne, nb, nk)
        dos0[s:s + args.chunk] = -1.0 / np.pi * np.sum(G0.imag, axis=(1, 2))
        G = compute_G(e, eigs, M, eta)                           # (ne, nb, nk, nb, nk)
        Gdiag = np.einsum("eabab->eab", G)                       # (ne, nb, nk)
        dos[s:s + args.chunk] = -1.0 / np.pi * np.sum(Gdiag.imag, axis=(1, 2))
        print(f"  energies {s}-{min(s+args.chunk, args.ne)}/{args.ne} done", flush=True)

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
