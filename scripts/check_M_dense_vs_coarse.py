"""Physics check of the dense (zero-padded) M^L: at the k-points the dense (pN)x(pN) grid shares with the
coarse NxN grid, the nb x nb blocks M[:,k',:,k] must agree with the coarse M_L up to the band gauge
(different nscf runs) -> compare their SINGULAR VALUES (gauge invariant) and the gauge-invariant
Frobenius norms. Usage: python scripts/check_M_dense_vs_coarse.py 5x5"""
import sys, json, numpy as np
from electron_defect_interaction.io import qe_io
size = sys.argv[1]; N = int(size.split("x")[0])
D = {"5x5": 25, "6x6": 24, "7x7": 28, "8x8": 32, "9x9": 27, "12x12": 24}[size]
uc_c = f"data/graphene/unit_cell/qe/defect_{size}.save"
uc_d = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"
kc = qe_io.get_k_red(uc_c); kd = qe_io.get_k_red(uc_d)
Mc = np.load(f"results/M/M_L_{size}.npy") if len(sys.argv) < 3 else np.load(sys.argv[2])
Md = np.load(f"results/M/M_L_dense_{size}.npy") if len(sys.argv) < 4 else np.load(sys.argv[3])
key = lambda k: tuple(np.round(np.mod(k + 1e-9, 1.0), 6))
idx_d = {key(k): i for i, k in enumerate(kd)}
pairs = [(ic, idx_d[key(k)]) for ic, k in enumerate(kc) if key(k) in idx_d]
print(f"coarse nk={len(kc)} dense nk={len(kd)}  coincident k: {len(pairs)} (expected {N*N})")
worst = 0.0; worst_n = 0.0
for ic, id_ in pairs:
    for jc, jd in pairs:
        nbc = min(Mc.shape[0], Md.shape[0]); A = Mc[:nbc, ic, :nbc, jc]; B = Md[:nbc, id_, :nbc, jd]   # coarse (16 bands) vs dense (20)
        sa = np.linalg.svd(A, compute_uv=False); sb = np.linalg.svd(B, compute_uv=False)
        worst = max(worst, np.max(np.abs(sa - sb)) / max(1e-30, sa[0]))
        worst_n = max(worst_n, abs(np.linalg.norm(A) - np.linalg.norm(B)) / max(1e-30, np.linalg.norm(A)))
print(f"max rel. singular-value mismatch over coincident blocks: {worst:.3e}")
print(f"max rel. Frobenius-norm mismatch: {worst_n:.3e}")
print(f"max|M_L dense| = {np.max(np.abs(Md)):.4e}  max|M_L coarse| = {np.max(np.abs(Mc)):.4e}  (both unit_cell-normalised: same order)")
print("RESULT:", "PASS" if worst < 1e-4 else "CHECK (nscf gauge/convergence tolerance exceeded)")
