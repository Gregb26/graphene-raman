"""
compare_bands_qe.py
    Compare Quantum ESPRESSO DFT bands (bands.x 'bands.dat') with the Wannier90-interpolated bands along
    the same k-path. This is the QE analogue of scripts/compare_bands.py (which is ABINIT/EIG.nc only and,
    importantly, ignores the Wigner-Seitz degeneracies). Here we use read_w90_HR WITH ndegen and Hwr_to_Hwk.

    bands.dat (bands.x raw format): header line '&plot nbnd=.., nks=.. /', then for each k-point a line
    with the cartesian k (units 2*pi/alat) followed by nbnd eigenvalues in eV. We convert k to reduced
    coordinates (mirroring qe_io.get_k_red), Fourier-interpolate H(R) -> H(k) -> eigenvalues, and align
    both DFT and Wannier to the QE Fermi level.

    Unlike the coarse-grid check in validate_wannier_bands.py, this compares the interpolation BETWEEN the
    coarse k-points (along a continuous path), which is the actual test of Wannier interpolation quality.
"""
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.wannier_io import read_w90_HR
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk

HA_TO_EV = 27.211386245988

DATA = "data/graphene/unit_cell/qe/defect_5x5.save"
BANDS_DAT = f"{DATA}/bands.dat"
TB_PATH = f"{DATA}/wannier_tb.dat"

# hexagonal high-symmetry points (reduced coords) + symmetry images, to label the path corners
SYM_POINTS = {
    r"$\Gamma$": [(0, 0, 0)],
    "K": [(2/3, 1/3, 0), (1/3, 2/3, 0), (-1/3, 1/3, 0), (1/3, -1/3, 0), (-1/3, -2/3, 0), (2/3, -1/3, 0)],
    "M": [(1/2, 0, 0), (0, 1/2, 0), (1/2, 1/2, 0), (-1/2, 1/2, 0), (1/2, -1/2, 0), (0, -1/2, 0)],
}


def read_qe_bands(path):
    """Return kcart (nks,3) in 2*pi/alat units and eig (nks,nbnd) in eV from a bands.x bands.dat file."""
    with open(path) as f:
        header = f.readline()
        nbnd = int(re.search(r"nbnd=\s*(\d+)", header).group(1))
        nks = int(re.search(r"nks=\s*(\d+)", header).group(1))
        data = np.array(f.read().split(), dtype=float).reshape(nks, 3 + nbnd)
    return data[:, :3], data[:, 3:]


def cart_to_red(kcart, save_dir):
    """QE cartesian k (2*pi/alat) -> reduced coords; mirrors qe_io.get_k_red (k_cart = B @ k_red)."""
    B, _ = qe_io.get_B_volume(save_dir)
    alat = float(ET.parse(f"{save_dir}/data-file-schema.xml").getroot()
                 .find(".//output/atomic_structure").get("alat"))
    return (kcart * (2 * np.pi / alat)) @ np.linalg.inv(B).T


def fermi_eV(save_dir):
    root = ET.parse(f"{save_dir}/data-file-schema.xml").getroot()
    return float(root.find(".//output/band_structure/fermi_energy").text) * HA_TO_EV


def path_corners(kcart):
    """Indices of the path endpoints and direction-change vertices (high-symmetry corners)."""
    seg = np.diff(kcart, axis=0)
    n = np.linalg.norm(seg, axis=1, keepdims=True)
    t = seg / np.where(n > 0, n, 1.0)
    corners = [0]
    for i in range(1, len(t)):
        if np.dot(t[i], t[i - 1]) < 0.999:        # tangent direction changes -> a corner
            corners.append(i)
    corners.append(len(kcart) - 1)
    return corners


def label_corner(kred, tol=2e-2):
    for name, imgs in SYM_POINTS.items():
        for s in imgs:
            d = kred - np.array(s)
            if np.linalg.norm(d - np.round(d)) < tol:   # match modulo a reciprocal lattice vector
                return name
    return ""


def nearest_index(kred, targets):
    """Index of the path point closest (modulo a reciprocal lattice vector) to any of `targets`."""
    best_i, best_d = 0, np.inf
    for i, kr in enumerate(kred):
        for s in targets:
            d = kr - np.array(s)
            dd = np.linalg.norm(d - np.round(d))
            if dd < best_d:
                best_d, best_i = dd, i
    return best_i


def main():
    # --- read DFT bands and interpolate Wannier bands on the same path ---
    kcart, dft = read_qe_bands(BANDS_DAT)                 # (nks,3), (nks,nbnd) eV (absolute)
    kred = cart_to_red(kcart, DATA)

    Hwr, Rw, ndegen = read_w90_HR(TB_PATH)
    wann = Hwr_to_Hwk(Hwr, Rw, kred, ndegen=ndegen)[1]    # (nks, nw) eV (absolute)
    nks, nbnd = dft.shape
    nw = wann.shape[1]

    # Reference both band structures to their OWN Dirac point at K. The separate bands.x run and the
    # nscf run that seeded Wannier90 differ by a rigid energy-reference offset (~2.4 eV here, verified
    # constant across the low bands); aligning on the graphene Dirac point removes it physically.
    # Graphene: 8 valence e- -> 4 occupied bands, the pi/pi* cone touches between band 3 and 4 (0-indexed).
    iK = nearest_index(kred, SYM_POINTS["K"])
    nocc = 4
    Ef_dft = np.sort(dft[iK])[nocc - 1:nocc + 1].mean()
    Ef_wan = np.sort(wann[iK])[nocc - 1:nocc + 1].mean()
    dft = dft - Ef_dft
    wann = wann - Ef_wan
    print(f"path: {nks} k-points, nbnd_DFT={nbnd}, nw={nw}")
    print(f"Dirac alignment at K: E_F(DFT)={Ef_dft:.3f} eV, E_F(Wannier)={Ef_wan:.3f} eV, "
          f"rigid offset={Ef_dft - Ef_wan:.3f} eV")

    # --- quantitative agreement along the path: each Wannier band vs the nearest DFT band ---
    # (robust to band crossings; outside the frozen window the disentangled subspace can drift off
    #  individual DFT bands, so we report the median and the fraction within 50 meV, not a raw RMS)
    dmin = np.array([[np.min(np.abs(dft[ik] - e)) for e in wann[ik]] for ik in range(nks)])
    print(f"Wannier vs DFT along the path (nearest band): median={np.median(dmin)*1e3:.1f} meV, "
          f"90th pct={np.percentile(dmin, 90)*1e3:.1f} meV, frac<50meV={np.mean(dmin < 0.05):.2f}")

    # --- curvilinear distance + high-symmetry ticks ---
    dist = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(kcart, axis=0), axis=1))])
    ticks, labels = [], []
    for c in path_corners(kcart):
        if not ticks or abs(dist[c] - ticks[-1]) > 1e-9:
            ticks.append(dist[c])
            labels.append(label_corner(kred[c]))
    print("path corners:", " -> ".join(l or "?" for l in labels))

    # --- plot (DFT black solid, Wannier red dashed) ---
    fig, ax = plt.subplots(figsize=(6, 5))
    for ib in range(nbnd):
        ax.plot(dist, dft[:, ib], "k-", lw=1.0, label="DFT (QE)" if ib == 0 else None)
    for iw in range(nw):
        ax.plot(dist, wann[:, iw], "r--", lw=1.3, label="Wannier" if iw == 0 else None)
    for x in ticks:
        ax.axvline(x, color="k", lw=0.5, alpha=0.4)
    ax.axhline(0.0, color="b", lw=0.6, ls=":", alpha=0.7, label=r"$E_F$")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels)
    ax.set_xlim(dist[0], dist[-1])
    ax.set_ylim(dft.min() - 0.5, dft.max() + 0.5)
    ax.set_ylabel(r"$E - E_F$ (eV)")
    ax.set_title("QE DFT vs Wannier-interpolated bands (graphene)")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    fig.tight_layout()
    out = "results/compare_bands_qe.png"
    fig.savefig(out, dpi=200)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
