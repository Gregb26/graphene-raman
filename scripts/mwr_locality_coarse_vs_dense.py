"""Figure-4 data: on-site pz-pz and ||Mwr(R,R0)|| decay, coarse NxN (aliased) vs dense zero-padded M, per size. eV."""
import numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from electron_defect_interaction.config import load_production, dense_paths, HA2EV
import os
cfg = load_production(); out = {}
def run(uc, mfile, wdir, tag):
    M = matrix_io.load_M_checked(mfile, require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)
    k = qe_io.get_k_red(uc); MP = _infer_mp_grid(k)
    U, kU = read_w90_mat(f"{wdir}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(f"{wdir}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
    Mwr, R = Mwk_to_Mwr(Mbk_to_Mwk(M, U, Ud), k, MP); Rn, Rd = lt.recenter_mwr(Mwr, R, MP)
    dist, wt = lt.mwr_locality(Mwr, Rn); i0 = int(np.argmin(np.abs(Rn).sum(1))); on = Mwr[:, i0, :, i0]
    # cartesian distance in units of a (hexagonal lattice: |R| = a sqrt(i^2 + j^2 - i j) for a 60-deg cell? use metric from the cell)
    A = np.array([[4.0354919061, -2.3298923383], [4.0354919061, 2.3298923383]]) / 4.6597846766   # rows a1,a2 in units of a
    dcart = np.linalg.norm(Rn[:, :2] @ A, axis=1)
    order = np.argsort(dcart); wR = np.array([np.linalg.norm(Mwr[:, i, :, i0]) for i in range(len(Rn))])
    out[f"{tag}_dist"] = dcart[order]; out[f"{tag}_w"] = wR[order]; out[f"{tag}_onsite_pzA"] = on[3, 3].real; out[f"{tag}_onsite_pzB"] = on[4, 4].real
    out[f"{tag}_onsite_norm"] = np.linalg.norm(on); out[f"{tag}_onsite_pzvac"] = max(on[3, 3].real, on[4, 4].real); out[f"{tag}_vac_sublattice"] = "A" if on[3, 3].real >= on[4, 4].real else "B"
    print(f"[{tag}] on-site pz(A)={on[3,3].real:.4f} pz(B)={on[4,4].real:.4f} ||on||={np.linalg.norm(on):.4f} eV; first shells:",
          [(round(float(d), 2), round(float(w), 4)) for d, w in zip(dcart[order][:8], wR[order][:8])], flush=True)
for S in ("5x5", "6x6", "7x7", "8x8", "9x9", "12x12"):
    dp = dense_paths(cfg, S)
    if not os.path.exists(dp["mfile"]) or not os.path.isdir(dp["wdir"]): print(f"[{S}] dense M or wannier missing; skipped"); continue
    run(dp["uc"], dp["mfile"], dp["wdir"], f"{S}_dense")
    if os.path.isdir(f"wannier/{S}"):
        run(f"data/graphene/unit_cell/qe/defect_{S}.save", f"results/M/M_ed_{S}.npy", f"wannier/{S}", f"{S}_coarse")
    else:
        print(f"[{S}] no coarse wannierization (wannier/{S}); dense only")
np.savez("results/M/mwr_locality.npz", **out); print("saved results/M/mwr_locality.npz")
