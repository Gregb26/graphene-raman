"""
M4 pilot, convergence sweep of sigma(omega)/sigma_0 on the 27 x 27 data (EM.md, section 8).
Same loop as sigma_on_grid (make_grid_tb -> blocks -> compute_velocity -> kubo_accumulate -> kubo_normalize),
but the three eta share each diagonalization. Read-only on the repo; outputs next to this script.
"""
import time, sys
from pathlib import Path
import numpy as np
from electron_defect_interaction.electron_photon import (make_wannier_tb, centres_only, make_grid_tb,
                                                         compute_velocity, kubo_accumulate, kubo_normalize)

REPO = Path(__file__).resolve().parents[4]    # repo root (memoire/EM/M4_sigma/pilote/)
OUT = Path(__file__).resolve().parent
E_D = -4.238895                          # eV, Dirac point (W90_REF)
NS = [900, 1200, 1800]
ETAS = [0.02, 0.04, 0.08]                # eV, Gaussian standard deviation
HW = np.arange(20, 601) / 100            # eV, 0.20 ... 6.00, exact decimals (1.96, 2.33, 2.54 on the grid)
CHUNK = 50_000

tb = make_wannier_tb(REPO / "wannier/27x27/wannier_tb.dat")
variants = {"full": (tb, "berry"), "centres_only": (centres_only(tb), "berry"), "no_berry": (tb, "no_berry")}
A_cell = np.linalg.norm(np.cross(tb.lattice[:, 0], tb.lattice[:, 1])) # as in sigma_on_grid

for N in NS:
    grid = make_grid_tb(tb, N)
    for name, (model, mode) in variants.items():
        t0 = time.time()
        S = np.zeros((len(ETAS), len(HW), 3, 3))
        for b in range(0, grid.nk, CHUNK):
            _, eps, _, hv = compute_velocity(model, grid.k_cart[b:b+CHUNK], mode)
            for i, eta in enumerate(ETAS):
                S[i] += kubo_accumulate(eps, hv, HW, mu=E_D, eta=eta)
        sigma = np.stack([kubo_normalize(S[i], HW, grid.nk, A_cell) for i in range(len(ETAS))])
        np.savez(OUT / f"em_sigma_{name}_N{N}.npz", hw=HW, eta=np.array(ETAS), sigma=sigma, N=N, mu=E_D,
                 variant=name, seconds=time.time() - t0)
        iL = [np.flatnonzero(HW == e)[0] for e in (1.96, 2.33, 2.54)]
        print(f"N={N:5d} {name:13s} {time.time()-t0:6.1f} s | sigma_xx(1.96, 2.33, 2.54), eta=0.04: "
              f"{np.round(sigma[1, iL, 0, 0], 5)}", flush=True)
print("DONE", flush=True)
