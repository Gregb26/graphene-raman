"""
REAL golden test (blocking): on the coarse 5x5 grid, the local Wannier t-matrix rate must equal the
dense compute_T rate when both are posed in the SAME 5-WF subspace: M is projected to the Wannier
smooth-Bloch gauge (5 bands, Wannier-interpolated energies), dense uses M/N_cells (supercell Bloch
norm), local uses the intensive Mwr from M_raw with R_local = full MP-dual grid (exact) and the
same coarse internal grid. Also reports the on-site V_loc and positivity. Usage: python ... [5x5]
"""
import sys
import numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import (
    Mbk_to_Mwk, Mwk_to_Mwr, Mwk_to_Mbk, _infer_mp_grid, _match_kpoint_order)
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.defects.many_body.single_defect import compute_T
from electron_defect_interaction.defects.many_body import local_tmatrix as lt

N = sys.argv[1] if len(sys.argv) > 1 else "5x5"
dense = len(sys.argv) > 2 and sys.argv[2] == "--dense"
eta = 0.10
PF = {"5x5": (5, 25), "7x7": (4, 28), "8x8": (4, 32), "9x9": (3, 27)}
if dense:
    D = PF[N][1]; W = f"wannier/{D}x{D}"
    uc = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"
    MFILE = f"results/M/M_dense_{N}.npy"
else:
    W = f"wannier/{N}"; uc = f"data/graphene/unit_cell/qe/defect_{N}.save"; MFILE = f"results/M/M_ed_{N}.npy"
k = qe_io.get_k_red(uc); Nc = len(k); MP = _infer_mp_grid(k)
U, kU = read_w90_mat(f"{W}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
Ud, kUd = read_w90_mat(f"{W}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
Hwr, Rw, nd = read_w90_HR(f"{W}/wannier_tb.dat")
M_raw = matrix_io.load_M_checked(MFILE, require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)   # eV (Wannier H is in eV)

# Wannier-gauge M (intensive), real space, recentered; 5-band smooth-Bloch projection on the coarse grid
Mwk = Mbk_to_Mwk(M_raw, U, Ud)
Mwr, R = Mwk_to_Mwr(Mwk, k, MP)
Rn, Rd = lt.recenter_mwr(Mwr, R, MP)
Mbk5 = Mwk_to_Mbk(Mwk, Hwr, Rw, k, ndegen=nd)              # (5, Nc, 5, Nc), intensive
_, Ew, _ = Hwr_to_Hwk(Hwr, Rw, k, ndegen=nd)               # (Nc, 5) eV
eigs = Ew.T                                                # (5, Nc)

# dense (supercell Bloch norm): compute_T on Mbk5/N_cells, on-shell, x N_cells -> per-defect
nw = 5
Gd = np.zeros((nw, Nc))
for ik in range(Nc):
    for n in range(nw):
        T = compute_T(np.array([Ew[ik, n]]), eigs, Mbk5 / Nc, eta)
        Gd[n, ik] = -2.0 * T[0, n, ik, n, ik].imag
Gd *= Nc

# local (intensive V_loc from M_raw), full R_cut, same coarse internal grid
V_loc, res = lt.extract_V_loc(Mwr, Rn, Rn)                 # all R -> exact
Gl = lt.scattering_rate(Hwr, Rw, nd, V_loc, Rn, k, eta, k_int=k)

rel = np.max(np.abs(Gl - Gd)) / max(1e-30, np.max(np.abs(Gd)))
i0 = int(np.argmin(np.abs(Rn).sum(1))); on = Mwr[:, i0, :, i0]
print(f"[{N}] R_d={Rd.tolist()}  on-site ||V_loc(0,0)||={np.linalg.norm(on):.4f} eV  herm_res={res:.1e}")
print(f"[{N}] Gamma range dense*Nc: {Gd.min():.3e}..{Gd.max():.3e}  local: {Gl.min():.3e}..{Gl.max():.3e}")
print(f"[{N}] REAL GOLDEN: max|local - dense*N_cells| rel = {rel:.3e}   (must be ~1e-10)")
print(f"[{N}] positivity min Gamma(local) = {Gl.min():.3e}")
print("RESULT:", "PASS" if (rel < 1e-8 and Gl.min() >= -1e-8) else "FAIL")
