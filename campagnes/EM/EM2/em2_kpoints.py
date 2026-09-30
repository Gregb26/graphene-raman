"""
EM2, step A.1: k list of the direct DFT check (QE `K_POINTS crystal`).

- 72 ring points from ring_kpoints_crystal (F17): {2.33: 48, 1.96: 12, 2.54: 12} eV around K, mu = E_D;
- 12 points on a circle of radius Q_DIRAC = 0.005 1/Angstrom around K (unit check independent of Wannier:
  hbar v_F^DFT = <Delta eps / (2 q)>_theta from the QE eigenvalues).

Checks first that CELL_PARAMETERS of the reference nscf.in is unit_cell_cart of wannier.win (same vectors,
same order), otherwise the crystal coordinates would not mean the same thing, and that tb.lattice is that
cell in Angstrom.

Writes em2_kpoints_crystal.txt (K_POINTS block, weight 1) and em2_kpoints_table.txt (index, set, hw, theta,
k_red, k_cart, |k - K|). Run: python em2_kpoints.py (repo from $GRAPHENE_RAMAN).
"""

import hashlib
import os
import re
from pathlib import Path

import numpy as np
from electron_defect_interaction.electron_photon import make_grid_tb, make_wannier_tb, ring_kpoints_crystal

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
REF = REPO.parent / "graphene" / "qe" / "defects" / "unit_cell" / "27x27"   # read only
TB = REPO / "wannier" / "27x27" / "wannier_tb.dat"
E_D = -4.238895                                   # eV, Dirac point (EM.md section 2)
NPOINTS = {2.33: 48, 1.96: 12, 2.54: 12}          # eV -> points on the ring
FIRST = np.array([0.73959906, 0.40626573, 0.0])   # first point at 2.33 eV (EM.md, F17)
Q_DIRAC, N_DIRAC = 0.005, 12                      # 1/Angstrom, points on the Dirac circle
BOHR_A = 0.529177210903                           # Angstrom (CODATA 2018, QE 7.5)


def read_cell(text, start, stop):
    """Three rows of floats between the line matching `start` and the line matching `stop`."""
    lines = text.splitlines()
    i0 = next(i for i, l in enumerate(lines) if re.match(start, l.strip(), re.I))
    rows = []
    for l in lines[i0 + 1:]:
        if re.match(stop, l.strip(), re.I):
            break
        parts = l.split()
        if len(parts) == 3:
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                pass    # 'bohr' unit line
    assert len(rows) == 3, f"expected 3 rows after {start}, got {len(rows)}"
    return np.array(rows)


# --- cell check (step 0.3, redone here so that the list and the check travel together) ---
nscf = (REF / "nscf.in").read_text()
win = (REF / "wannier.win").read_text()
assert re.search(r"CELL_PARAMETERS\s+bohr", nscf), "nscf.in CELL_PARAMETERS not in bohr"
cell_qe = read_cell(nscf, r"CELL_PARAMETERS", r"ATOMIC_POSITIONS")
cell_w90 = read_cell(win, r"begin unit_cell_cart", r"end unit_cell_cart")
assert re.search(r"begin unit_cell_cart\s*\n\s*bohr", win, re.I), "unit_cell_cart not in bohr"
assert np.array_equal(cell_qe, cell_w90), f"cells differ:\n{cell_qe}\n{cell_w90}"

tb = make_wannier_tb(TB)
lat_err = np.abs(tb.lattice.T - cell_qe * BOHR_A).max()
assert lat_err < 1e-6, f"tb.lattice is not the nscf cell: {lat_err}"
K = make_grid_tb(tb, 3).K

# --- ring points (F17) and Dirac circle ---
rings = ring_kpoints_crystal(tb, K, E_D, NPOINTS)
assert np.abs(rings[2.33][0][0] - FIRST).max() < 5e-9, rings[2.33][0][0]

theta_D = np.linspace(0, 2 * np.pi, N_DIRAC, endpoint=False)
k_D = K + Q_DIRAC * np.column_stack((np.cos(theta_D), np.sin(theta_D), np.zeros(N_DIRAC)))
k_red_D = k_D @ tb.lattice / (2 * np.pi)

sets = [("ring", hw, *rings[hw]) for hw in NPOINTS] + [("dirac", 0.0, k_red_D, theta_D)]
k_red = np.vstack([s[2] for s in sets])
B = 2 * np.pi * np.linalg.inv(tb.lattice).T        # columns b_j (a_i . b_j = 2 pi delta_ij)
k_cart = k_red @ B.T
K_red = K @ tb.lattice / (2 * np.pi)
assert np.abs(K_red - [2 / 3, 1 / 3, 0]).max() < 1e-12, K_red
assert np.abs(k_red[:, 2]).max() < 1e-14

# --- files ---
nk = k_red.shape[0]
with open(HERE / "em2_kpoints_crystal.txt", "w") as f:
    f.write(f"K_POINTS crystal\n{nk}\n")
    for k in k_red:
        f.write(f"  {k[0]:18.14f} {k[1]:18.14f} {k[2]:18.14f}  1.0\n")

with open(HERE / "em2_kpoints_table.txt", "w") as f:
    f.write(f"# EM2 k list: {nk} points, K_POINTS crystal order. E_D = {E_D} eV, K = {K} 1/A (crystal {K_red})\n")
    f.write(f"# _tb.dat {TB.relative_to(REPO)} sha256 {hashlib.sha256(TB.read_bytes()).hexdigest()[:16]}\n")
    f.write("# ik set hw(eV) theta(deg) k1 k2 k3 (crystal) kx ky kz (1/A) |k-K|(1/A)\n")
    ik = 0
    for name, hw, kr, th in sets:
        for j in range(kr.shape[0]):
            kc = kr[j] @ B.T
            f.write(f"{ik + 1:3d} {name:5s} {hw:5.2f} {np.degrees(th[j]):9.4f} "
                    f"{kr[j, 0]:17.14f} {kr[j, 1]:17.14f} {kr[j, 2]:17.14f} "
                    f"{kc[0]:15.12f} {kc[1]:15.12f} {kc[2]:15.12f} {np.linalg.norm(kc - K):14.12f}\n")
            ik += 1

print(f"cells identical (bohr), tb.lattice = cell x {BOHR_A} to {lat_err:.1e} A")
print(f"K = {K} 1/A, crystal {K_red}")
for name, hw, kr, th in sets:
    q = np.linalg.norm(kr @ B.T - K, axis=1)
    centre = (kr @ B.T).mean(axis=0) - K
    print(f"{name:5s} {hw:4.2f} eV: {kr.shape[0]:2d} k, |k-K| {q.min():.6f}-{q.max():.6f} 1/A, centre - K {np.abs(centre).max():.1e}")
print(f"first point at 2.33 eV: {rings[2.33][0][0]}")
print(f"wrote {nk} k points")
