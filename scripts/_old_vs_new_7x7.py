import numpy as np
from electron_defect_interaction.io import matrix_io
old = np.load("results/M/obsolete_grid_7x7/M_L_dense_7x7.npy", mmap_mode="r"); new = np.load("results/M/M_L_dense_7x7.npy", mmap_mode="r")
oM = np.load("results/M/obsolete_grid_7x7/M_dense_7x7.npy", mmap_mode="r"); nM = np.load("results/M/M_dense_7x7.npy", mmap_mode="r")
rng = np.random.default_rng(0); ks = rng.choice(784, 40, replace=False); dL = dM = 0.0; svL = svM = 0.0
for i in ks:
    for j in ks:
        a, b = np.array(old[:16, i, :16, j]), np.array(new[:16, i, :16, j]); dL = max(dL, np.abs(a - b).max()); svL = max(svL, np.abs(np.linalg.svd(a, compute_uv=False) - np.linalg.svd(b, compute_uv=False)).max() / np.linalg.svd(b, compute_uv=False)[0])
        a, b = np.array(oM[:16, i, :16, j]), np.array(nM[:16, i, :16, j]); dM = max(dM, np.abs(a - b).max()); svM = max(svM, np.abs(np.linalg.svd(a, compute_uv=False) - np.linalg.svd(b, compute_uv=False)).max() / np.linalg.svd(b, compute_uv=False)[0])
mL = max(np.abs(np.array(new[:16, i, :16, j])).max() for i in ks[:5] for j in ks[:5]); mM = np.abs(np.array(nM[:16, ks[0], :16, ks[0]])).max()
dg = lambda X: np.array([X[n, k, n, k] for n in range(20) for k in range(0, 784, 8)])
print(f"[7x7 dense] old (ix%30 on 216) vs new (resampled 217), 40x40 k pairs, bands 1-16: M^L max|diff| = {dL*27.2114:.4f} eV (max|M^L| ~ {mL*27.2114:.3f} eV), SV mismatch {svL:.3e}; "
      f"full M: max|diff| = {dM*27.2114:.4f} eV, SV mismatch {svM:.3e}; mean diag M^L old {dg(old).real.mean()*27.2114:+.4f} eV new {dg(new).real.mean()*27.2114:+.4f} eV")
