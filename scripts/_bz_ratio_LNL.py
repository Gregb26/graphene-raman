import sys, numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat
from electron_defect_interaction.wannier.wannier_interpolation import _match_kpoint_order
from electron_defect_interaction.config import load_production, dense_paths, HA2EV
cfg = load_production(verbose=False)
for S in sys.argv[1:]:
    dp = dense_paths(cfg, S); kd = qe_io.get_k_red(dp["uc"]); kk, eps = qe_io.get_k_eigenvalues(dp["uc"], False); eps = np.asarray(eps)
    if eps.shape[0] != len(kd): eps = eps.T
    U, kU = read_w90_mat(f"{dp['wdir']}/wannier_u.mat"); U = U[_match_kpoint_order(kU, kd)]; Ud, kUd = read_w90_mat(f"{dp['wdir']}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, kd)]
    V = np.einsum("kbw,kwv->kbv", Ud, U); w = np.abs(V[:, :, 3]) ** 2 + np.abs(V[:, :, 4]) ** 2                     # pz weight per (k, band)
    pi_idx = np.array([sorted(np.argsort(-w[k])[:2], key=lambda n: eps[k, n]) for k in range(len(kd))])
    ML = np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r"); MN = np.load(dp["mfile"].replace("M_dense_", "M_NL_dense_"), mmap_mode="r"); nk = len(kd)
    fL = np.zeros((nk, nk)); fN = np.zeros((nk, nk))
    for ik in range(nk):
        rowL = np.array(ML[:, :, :, ik]) * HA2EV; rowN = np.array(MN[:, :, :, ik]) * HA2EV
        for jk in range(nk):
            fL[jk, ik] = np.linalg.norm(rowL[pi_idx[jk]][:, pi_idx[ik]]); fN[jk, ik] = np.linalg.norm(rowN[pi_idx[jk]][:, pi_idx[ik]])
    print(f"[{S}] pi/pi* subspace, all (k',k) on {int(np.sqrt(nk))}x{int(np.sqrt(nk))}: <|M^NL|_F>/<|M^L|_F> = {fN.mean()/fL.mean():.2f}; <|M^NL|_F> = {fN.mean():.4f} eV, <|M^L|_F> = {fL.mean():.4f} eV; diag k'=k: {np.diag(fN).mean()/np.diag(fL).mean():.2f}; ratio min {(fN/fL).min():.2f} max {(fN/fL).max():.2f}", flush=True)
