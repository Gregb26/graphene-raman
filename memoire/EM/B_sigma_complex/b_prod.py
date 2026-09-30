"""
Perspective B (memoire/EM/EM.md section 11): complex sigma(omega)/sigma_0 of the 27 x 27 data, full velocity,
gaussian_complex kernel, eta = 0.04 eV, N = 1200 (as M4), hw = 0.20 ... 6.00 eV, four cases sharing each
diagonalization: undoped at T = 0 and 300 K, n and p doping (mu = E_D +/- 0.3 eV) at 300 K. The block loop is
that of sigma_on_grid (compute_velocity -> kubo_accumulate -> kubo_normalize), checked against it on a small grid.

Run from anywhere:  .venv/bin/python memoire/EM/B_sigma_complex/b_prod.py   (outputs next to this script)
"""

import hashlib
import subprocess
import time
from pathlib import Path

import numpy as np
from electron_defect_interaction.electron_photon import (compute_velocity, gaussian_complex, kubo_accumulate,
                                                         kubo_normalize, make_grid_tb, make_wannier_tb, sigma_on_grid)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
E_D = -4.238895                                  # eV (EM.md section 2)
N, ETA, CHUNK = 1200, 0.04, 10000                # grid, eV, k per block (complex kernel: ~0.7 GB per block)
HW = np.arange(20, 601) / 100                    # eV
CASES = {"undoped_T0": (E_D, 0.0), "undoped_300K": (E_D, 0.025), "n_300K": (E_D + 0.3, 0.025), "p_300K": (E_D - 0.3, 0.025)}

tb = make_wannier_tb(REPO / "wannier" / "27x27" / "wannier_tb.dat")
A_cell = np.linalg.norm(np.cross(tb.lattice[:, 0], tb.lattice[:, 1]))


def run(n, hw, chunk):
    grid = make_grid_tb(tb, n)
    S = {c: 0 for c in CASES}
    for b in range(0, grid.nk, chunk):
        _, eps, _, hv = compute_velocity(tb, grid.k_cart[b:b + chunk], "berry")
        for c, (mu, kT) in CASES.items():
            S[c] = S[c] + kubo_accumulate(eps, hv, hw, mu, ETA, kT=kT, kernel=gaussian_complex)
    return {c: kubo_normalize(S[c], hw, grid.nk, A_cell, kernel=gaussian_complex) for c in CASES}


# the shared loop is sigma_on_grid, case by case
small = run(120, HW[::50], 3000)
for c, (mu, kT) in CASES.items():
    ref = sigma_on_grid(tb, 120, HW[::50], mu=mu, eta=ETA, kT=kT, kernel=gaussian_complex)
    assert np.allclose(small[c], ref, rtol=0, atol=1e-13), c
print("shared loop == sigma_on_grid (N = 120)", flush=True)

t0 = time.time()
sig = run(N, HW, CHUNK)
src = REPO / "src" / "electron_defect_interaction" / "electron_photon" / "kubo.py"
git = lambda *a: subprocess.run(["git", "-C", str(REPO), *a], capture_output=True, text=True).stdout.strip()
np.savez(HERE / "b_sigma_complex_N1200_eta0.04.npz", hw=HW, cases=np.array(list(CASES)),
         mu=np.array([m for m, _ in CASES.values()]), kT=np.array([k for _, k in CASES.values()]),
         sigma=np.stack([sig[c] for c in CASES]), N=N, eta=ETA, E_D=E_D, kernel="gaussian_complex",
         commit=git("rev-parse", "--short", "HEAD"), src_dirty=bool(git("status", "--porcelain", "src")),
         kubo_sha256=hashlib.sha256(src.read_bytes()).hexdigest()[:16])
iL = [int(np.flatnonzero(np.isclose(HW, e))[0]) for e in (1.96, 2.33, 2.54, 4.05)]
print(f"done in {time.time() - t0:.0f} s")
for c in CASES:
    print(f"{c:13s} xx at 1.96 2.33 2.54 4.05:", " ".join(f"{z.real:.5f}{z.imag:+.5f}i" for z in sig[c][iL, 0, 0]),
          f"| max|yy-xx| {np.abs(sig[c][:, 1, 1] - sig[c][:, 0, 0]).max():.1e}")
