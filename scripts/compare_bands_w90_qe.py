"""
compare_bands_w90_qe.py
    Overlay the Wannier90-plotted bands (`wannier_band.dat`, from bands_plot=.true.) on top of the
    Quantum ESPRESSO DFT bands (`bands.dat`, from bands.x) along the same high-symmetry path.

    Both files are read directly; no Hamiltonian interpolation here (unlike compare_bands_qe.py, which
    re-interpolates H(R) onto the QE path). This just compares the two band files as produced.

    Two reference issues are handled so the curves actually lie on top of each other:
      * Energy:  the separate bands.x run and the nscf run that seeded Wannier90 differ by a rigid
                 energy-reference offset (~2.4 eV for graphene). Each band structure is shifted so its
                 OWN Dirac point sits at 0. The Dirac energy is found WITHOUT needing the K coordinate:
                 it is the midpoint of the smallest gap between the top-valence and bottom-conduction
                 bands along the path (graphene: bands 4/5, i.e. 0-indexed 3/4).
      * x-axis:  wannier_band.dat distances are in Ang^-1, bands.dat k is in 2*pi/alat. Both path
                 parametrizations are rescaled to [0, 1]; this lines up the high-symmetry corners
                 provided the two runs sample the SAME path with proportional segment sampling.

    File formats:
      bands.dat (bands.x): header '&plot nbnd=.., nks=.. /', then per k-point a line with cartesian k
                           (2*pi/alat) followed by nbnd eigenvalues in eV.
      wannier_band.dat:    blocks separated by blank lines, one block per band, each line "dist energy"
                           with energy in eV.
"""
import re
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import _paths
from _bands import SYM_POINTS, path_corners, label_corner, dirac_energy
from electron_defect_interaction.io import qe_io

HA_TO_EV = 27.211386245988

# Defaults (override on the command line). 11x11 ships both bands.dat and wannier_band.dat locally.
DATA = _paths.uc("11x11")
QE_BANDS = f"{DATA}/bands.dat"
W90_BANDS = f"{DATA}/wannier_band.dat"
OUT = "results/compare_bands_w90_qe.png"

# graphene: 8 valence electrons -> 4 occupied bands; the pi/pi* Dirac cone touches between the
# 4th and 5th bands (0-indexed 3 and 4).
NOCC = 4


def read_w90_bands(path):
    """Return x (nks,) cumulative distance and eig (nks,nbnd) in eV from a wannier_band.dat file."""
    blocks = re.split(r"\n\s*\n", open(path).read().strip())
    cols = [np.array(b.split(), dtype=float).reshape(-1, 2) for b in blocks if b.strip()]
    x = cols[0][:, 0]
    eig = np.column_stack([c[:, 1] for c in cols])           # (nks, nbnd)
    return x, eig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qe", default=QE_BANDS, help="QE bands.x bands.dat")
    ap.add_argument("--w90", default=W90_BANDS, help="Wannier90 wannier_band.dat")
    ap.add_argument("--save", default=DATA, help=".save dir (for cart->reduced k conversion / corner labels)")
    ap.add_argument("--out", default=OUT, help="output figure path")
    ap.add_argument("--nocc", type=int, default=NOCC, help="number of occupied bands (Dirac alignment)")
    args = ap.parse_args()

    # --- read both band files ---
    kcart, qe = qe_io.get_qe_bands(args.qe)                  # (nks_q,3), (nks_q, nbnd) eV
    xw, w90 = read_w90_bands(args.w90)                       # (nks_w,),  (nks_w, nw)  eV

    # --- align each on its own Dirac point (removes the rigid energy-reference offset) ---
    Ed_qe = dirac_energy(qe, args.nocc)
    Ed_w90 = dirac_energy(w90, args.nocc)
    qe = qe - Ed_qe
    w90 = w90 - Ed_w90
    print(f"QE : {qe.shape[0]} k, {qe.shape[1]} bands, Dirac at {Ed_qe:.3f} eV")
    print(f"W90: {w90.shape[0]} k, {w90.shape[1]} bands, Dirac at {Ed_w90:.3f} eV")
    print(f"rigid energy-reference offset (QE - W90) = {Ed_qe - Ed_w90:.3f} eV")

    # --- common x-axis: rescale both path parametrizations to [0, 1] ---
    dq = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(kcart, axis=0), axis=1))])
    xq = dq / dq[-1]
    xw = (xw - xw[0]) / (xw[-1] - xw[0])

    # --- high-symmetry ticks from the QE path (in the same normalized units) ---
    kred = qe_io.kcart_to_kred(kcart, args.save)
    ticks, labels = [], []
    for c in path_corners(kcart):
        if not ticks or abs(xq[c] - ticks[-1]) > 1e-9:
            ticks.append(xq[c])
            labels.append(label_corner(kred[c]))
    print("path corners:", " -> ".join(l or "?" for l in labels))

    # --- plot: QE DFT black solid, Wannier90 red dashed ---
    fig, ax = plt.subplots(figsize=(6, 5))
    for ib in range(qe.shape[1]):
        ax.plot(xq, qe[:, ib], "k-", lw=1.0, label="DFT (QE)" if ib == 0 else None)
    for iw in range(w90.shape[1]):
        ax.plot(xw, w90[:, iw], "r--", lw=1.3, label="Wannier90" if iw == 0 else None)
    for x in ticks:
        ax.axvline(x, color="k", lw=0.5, alpha=0.4)
    ax.axhline(0.0, color="b", lw=0.6, ls=":", alpha=0.7, label="Dirac point")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlim(0, 1)
    ax.set_ylabel(r"$E - E_\mathrm{Dirac}$ (eV)")
    ax.set_title("QE DFT vs Wannier90 bands (graphene)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"saved {args.out}")


if __name__ == "__main__":
    main()
