"""Plumbing check: M_L from the node-shared kernel on the COARSE grid must equal the reference M_L_<size>.npy
(produced by compute_ML_R_mpi) to machine precision -- same wavefunctions, same potential, same convention."""
import sys, numpy as np
size = sys.argv[1]
A = np.load(f"results/M/M_L_{size}.npy"); B = np.load(f"results/M/M_L_dense_{size}_coarsecheck.npy")
print("shapes", A.shape, B.shape)
d = np.max(np.abs(A - B)); m = np.max(np.abs(A))
print(f"max|A-B| = {d:.3e}   max|A| = {m:.3e}   rel = {d/m:.3e}")
H = np.max(np.abs(B.reshape(A.shape[0]*A.shape[1], -1) - B.reshape(A.shape[0]*A.shape[1], -1).conj().T))
print(f"Hermiticity |M - M^dag| max = {H:.3e}")
print("RESULT:", "PASS" if d/m < 1e-10 else "FAIL")
