"""Figure-4 data: on-site pz-pz and ||Mwr(R,R0)|| decay, coarse NxN (aliased) vs dense zero-padded M, per size. eV.
R10 : Mwr aligned (M_W(R,R) - C_N on the N x N box, approximation (i), C_N from the config block "alignment"); abscissa = true distance of the
Wigner-Seitz image to the defect (ws_images, P-b2)."""
import numpy as np
from graphene_raman.io import qe_io, matrix_io
from graphene_raman.io.wannier_io import read_w90_mat
from graphene_raman.wannier.wannier_interpolation import _infer_mp_grid, _match_kpoint_order, ws_images
from graphene_raman.defects.many_body import cluster_tmatrix as ct
from graphene_raman.config import load_production, dense_paths, HA2EV, results_dir, matrices_dir, alignment_C, wannier_dir
RES = results_dir(load_production(verbose=False))          # R10 : produits (results/M2_plateau)
MAT = matrices_dir(load_production(verbose=False))         # R10 : matrices M2 brutes (results/M2, lecture seule)
import os
cfg = load_production(); out = {}
def run(uc, mfile, wdir, tag, S):
    M = matrix_io.load_M_checked(mfile, require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    k = qe_io.get_k_red(uc); MP = _infer_mp_grid(k)
    U, kU = read_w90_mat(f"{wdir}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(f"{wdir}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
    C_N = alignment_C(cfg, S); d = ct.defect_mwr(M, U, Ud, k, MP, n_box=int(S.split("x")[0]), C_N=C_N); Mwr, Rn = d["Mwr"], d["Rn"]   # R10 : approximation (i)
    dist, wt = ct.mwr_locality(Mwr, Rn); i0 = int(np.argmin(np.abs(Rn).sum(1))); on = Mwr[:, i0, :, i0]
    # cartesian distance in units of a (hexagonal lattice: |R| = a sqrt(i^2 + j^2 - i j) for a 60-deg cell? use metric from the cell)
    A = np.array([[4.0354919061, -2.3298923383], [4.0354919061, 2.3298923383]]) / 4.6597846766   # rows a1,a2 in units of a
    A_cols = np.eye(3); A_cols[:2, :2] = A.T                                                          # columns a1, a2 (units of a), a3 unused (R_z = 0)
    dcart = ws_images(Rn, np.zeros(3, int), MP, A_cols)["dist"]                                          # R10 (P-b2) : Wigner-Seitz image, was |Rn @ A|
    order = np.argsort(dcart); wR = np.array([np.linalg.norm(Mwr[:, i, :, i0]) for i in range(len(Rn))])
    out[f"{tag}_dist"] = dcart[order]; out[f"{tag}_w"] = wR[order]; out[f"{tag}_onsite_pzA"] = on[3, 3].real; out[f"{tag}_onsite_pzB"] = on[4, 4].real
    out[f"{tag}_onsite_norm"] = np.linalg.norm(on); out[f"{tag}_onsite_pzvac"] = max(on[3, 3].real, on[4, 4].real); out[f"{tag}_vac_sublattice"] = "A" if on[3, 3].real >= on[4, 4].real else "B"
    print(f"[{tag}] C_N = {C_N*1e3:+.4f} meV; on-site pz(A)={on[3,3].real:.4f} pz(B)={on[4,4].real:.4f} ||on||={np.linalg.norm(on):.4f} eV; first shells:",
          [(round(float(d), 2), round(float(w), 4)) for d, w in zip(dcart[order][:8], wR[order][:8])], flush=True)
for S in ("5x5", "6x6", "7x7", "8x8", "9x9", "12x12"):
    dp = dense_paths(cfg, S)
    if not os.path.exists(dp["mfile"]) or not os.path.isdir(dp["wdir"]): print(f"[{S}] dense M or wannier missing; skipped"); continue
    run(dp["uc"], dp["mfile"], dp["wdir"], f"{S}_dense", S)
    if os.path.isdir(wannier_dir(cfg, S)):
        run(f"data/graphene/unit_cell/qe/defect_{S}.save", f"{MAT}/M_ed_{S}.npy", wannier_dir(cfg, S), f"{S}_coarse", S)
    else:
        print(f"[{S}] no coarse wannierization ({wannier_dir(cfg, S)}); dense only")
np.savez(f"{RES}/mwr_locality.npz", **out); print("saved <results_dir>/mwr_locality.npz")
