"""Q1: on-site pz-pz of V_loc (dense) in eV for the four sizes; Q2: M^L vs M^NL at K for the pi pair (9x9, 27x27 grid contains K);
plus max|M_L|, max|M_NL|, max|M| per size to explain 0.26 (5x5) vs 8.8e-3 (9x9)."""
import numpy as np, json
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
HA = 27.211386245988
PF = {"5x5": 25, "7x7": 28, "8x8": 32, "9x9": 27}
PZ = (3, 4)   # projections: 3 sp2 on atom A (idx 0..2), pz on A (3), pz on B (4)

print("=== Q2/magnitudes: max|M_L|, max|M_NL|, max|M| (Ha) and mean diag of M_L (mean-potential shift)")
for S, D in PF.items():
    ML = np.load(f"results/M/M_L_dense_{S}.npy", mmap_mode="r"); MN = np.load(f"results/M/M_NL_dense_{S}.npy", mmap_mode="r"); M = np.load(f"results/M/M_dense_{S}.npy", mmap_mode="r")
    nb, nk = ML.shape[:2]
    dL = np.array([ML[n, k, n, k] for n in range(nb) for k in range(nk)])
    print(f"[{S}] max|M_L|={np.abs(ML).max():.4e}  max|M_NL|={np.abs(MN).max():.4e}  max|M|={np.abs(M).max():.4e}  "
          f"mean diag M_L={dL.mean().real:.4e} Ha ({dL.mean().real*HA:.3f} eV)  max|offdiag-ish M_L|={np.abs(ML[:, 0, :, 1:]).max():.4e}", flush=True)
    del ML, MN, M

print("\n=== Q2: M^L vs M^NL at K, pi pair (9x9 dense, 27x27)")
S, D = "9x9", 27
uc = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"
k = qe_io.get_k_red(uc); kk, eps = qe_io.get_k_eigenvalues(uc, False); eps = np.asarray(eps)
if eps.shape[0] != len(kk): eps = eps.T
for Kp, lab in (((2/3, 1/3, 0), "K"), ((1/3, 2/3, 0), "K'")):   # Dirac points of the 60-deg QE cell
    d = np.linalg.norm(np.mod(k - np.array(Kp) + 0.5, 1.0) - 0.5, axis=1); iK = int(np.argmin(d))
    print(f"[{lab}] k_red={k[iK]} (dist {d[iK]:.2e}); eps bands 3..5 (eV): {np.round(eps[iK, 2:6]*HA, 4)}")
    ML = np.load(f"results/M/M_L_dense_{S}.npy", mmap_mode="r"); MN = np.load(f"results/M/M_NL_dense_{S}.npy", mmap_mode="r")
    bL = np.array(ML[3:5, iK, 3:5, iK]); bN = np.array(MN[3:5, iK, 3:5, iK])
    dL = np.abs(np.diag(bL)); dN = np.abs(np.diag(bN))
    print(f"[{lab}] |M_L[n,K,n,K]| eV: {np.round(dL*HA, 4)}  mean {dL.mean()*HA:.4f} eV   (Re: {np.round(np.diag(bL).real*HA, 4)})")
    print(f"[{lab}] |M_NL[n,K,n,K]| eV: {np.round(dN*HA, 4)}  mean {dN.mean()*HA:.4f} eV   (Re: {np.round(np.diag(bN).real*HA, 4)})")
    print(f"[{lab}] ||2x2 pi block||_F: M_L {np.linalg.norm(bL)*HA:.4f} eV  M_NL {np.linalg.norm(bN)*HA:.4f} eV  ratio NL/L {np.linalg.norm(bN)/np.linalg.norm(bL):.2f}", flush=True)
    del ML, MN

print("\n=== Q1: on-site V_loc (dense) pz-pz in eV")
for S, D in PF.items():
    uc = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"; W = f"wannier/{D}x{D}"
    M = matrix_io.load_M_checked(f"results/M/M_dense_{S}.npy", require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    k = qe_io.get_k_red(uc); MP = _infer_mp_grid(k)
    U, kU = read_w90_mat(f"{W}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(f"{W}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
    Mwk = Mbk_to_Mwk(M, U, Ud); del M
    Mwr, R = Mwk_to_Mwr(Mwk, k, MP); del Mwk
    Rn, Rd = lt.recenter_mwr(Mwr, R, MP)
    i0 = int(np.argmin(np.abs(Rn).sum(1))); on = Mwr[:, i0, :, i0]
    a, b = PZ
    print(f"[{S}] R_d={Rd.tolist()} ||on-site 5x5||={np.linalg.norm(on):.4f} Ha = {np.linalg.norm(on)*HA:.3f} eV")
    print(f"[{S}] pz(A)-pz(A) = {on[a,a].real*HA:+.4f} eV   pz(B)-pz(B) = {on[b,b].real*HA:+.4f} eV   pz(A)-pz(B) = {on[a,b].real*HA:+.4f}{on[a,b].imag*HA:+.4f}i eV   |Im diag| max {np.abs(np.diag(on).imag).max()*HA:.1e}")
    print(f"[{S}] sigma diag (eV): {np.round(np.diag(on)[:3].real*HA, 4)}   full on-site block |.| (eV):\n{np.round(np.abs(on)*HA, 3)}", flush=True)
