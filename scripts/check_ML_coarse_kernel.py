"""Plumbing check: M_L from the node-shared kernel on the COARSE grid must equal the reference M_L_<size>.npy
(produced by compute_ML_R_mpi) to machine precision -- same wavefunctions, same potential, same convention."""
import sys, numpy as np
from electron_defect_interaction.config import load_production, results_dir
RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)
size = sys.argv[1]
from electron_defect_interaction.io import matrix_io
for _f in (f"{RES}/M_L_{size}.npy", f"{RES}/M_L_dense_{size}_coarsecheck.npy"): matrix_io.check_manifest(_f, require_normalization=matrix_io.M_NORM_V2)
A = np.load(f"{RES}/M_L_{size}.npy"); B = np.load(f"{RES}/M_L_dense_{size}_coarsecheck.npy")
print("shapes", A.shape, B.shape)
d = np.max(np.abs(A - B)); m = np.max(np.abs(A))
print(f"max|A-B| = {d:.3e}   max|A| = {m:.3e}   rel = {d/m:.3e}")
H = np.max(np.abs(B.reshape(A.shape[0]*A.shape[1], -1) - B.reshape(A.shape[0]*A.shape[1], -1).conj().T))
print(f"Hermiticity |M - M^dag| max = {H:.3e}")
print("RESULT:", "PASS" if d/m < 1e-10 else "FAIL")
