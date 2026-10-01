# the sweep loop (three eta per diagonalization) reproduces sigma_on_grid at small N
from pathlib import Path
import numpy as np
from graphene_raman.electron_photon import make_wannier_tb, make_grid_tb, compute_velocity, kubo_accumulate, kubo_normalize, sigma_on_grid
tb = make_wannier_tb(Path(__file__).resolve().parents[4] / "results/wannier/27x27/wannier_tb.dat"); E_D = -4.238895
HW = np.arange(20, 601) / 100; N = 90; grid = make_grid_tb(tb, N)
A = np.linalg.norm(np.cross(tb.lattice[:, 0], tb.lattice[:, 1]))
S = sum(kubo_accumulate(*compute_velocity(tb, grid.k_cart[b:b+2000], 'berry')[1:4:2], HW, mu=E_D, eta=0.04) for b in range(0, grid.nk, 2000))
ref = sigma_on_grid(tb, N, HW, mu=E_D, eta=0.04, mode='berry')
print("max |diff|", np.abs(kubo_normalize(S, HW, grid.nk, A) - ref).max())
