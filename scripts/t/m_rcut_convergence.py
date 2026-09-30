#!/usr/bin/env python
"""
Convergence of the R-space truncation at the level of M itself (§4.1.5). For a size S: Mwr (dense M, recentred, U of the
dense wannierization) -> M(k'_f, k_f) on a fine NfxNf grid by the inverse Fourier transform, with the sums over R, R'
restricted to |R|,|R'| <= R_cut (reduced-coordinate norm, as in production) and with the full R set (reference).
R10 : Mwr aligned (M_W(R,R) - C_N on the N x N box, approximation (i), C_N from the config block "alignment"); off-grid phases through the
Wigner-Seitz images of the recentred labels (ws_images / ws_phase, P-c4), the R_cut masks stay on the recentred labels.
Reported per R_cut: max|dM|/max|M_ref| over the pi/pi* blocks of all (k',k); max relative singular-value mismatch of the
2x2 pi/pi* blocks (gauge invariant); same two measures on the diagonal elements M_nn(k,k) only. eV. Writes
<results_dir>/m_rcut_convergence.csv (appends one row per (size, R_cut)).
Usage: python scripts/t/m_rcut_convergence.py --size 9x9 [--nf 60] [--rcuts 0,1,2,3]
"""
import argparse, csv, os, numpy as np
from graphene_raman.io import qe_io, matrix_io
from graphene_raman.io.wannier_io import read_w90_mat, read_w90_tb
from graphene_raman.wannier.wannier_interpolation import _infer_mp_grid, _match_kpoint_order, ws_images, ws_phase
from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk
from graphene_raman.defects.many_body import local_tmatrix as lt
from graphene_raman.config import load_production, dense_paths, results_dir, alignment_C
RES = results_dir(load_production(verbose=False))          # R10 : produits (results/M2_plateau) ; matrices par matrices_dir
ap = argparse.ArgumentParser(); ap.add_argument("--size", default="9x9"); ap.add_argument("--nf", type=int, default=60); ap.add_argument("--rcuts", default="0,1,2,3"); a = ap.parse_args()
cfg = load_production(); dp = dense_paths(cfg, a.size); S = a.size
M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k); k = np.round(k * np.asarray(MP)) / np.asarray(MP)   # exact MP k (XML rounding)
U, kU = read_w90_mat(f"{dp['wdir']}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]; Ud, kUd = read_w90_mat(f"{dp['wdir']}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
Hwr, Rw, nd, _, _ = read_w90_tb(f"{dp['wdir']}/wannier_tb.dat")
C_N = alignment_C(cfg, S); d = lt.defect_mwr(M, U, Ud, k, MP, n_box=int(S.split("x")[0]), C_N=C_N); del M    # R10 : approximation (i)
Mwr, Rn, Rd = d["Mwr"], d["Rn"], d["R_d"]; del d; lt.mwr_locality(Mwr, Rn)
nw, nR = Mwr.shape[0], Mwr.shape[1]; print(f"[{S}] Mwr {Mwr.shape}, R grid {MP}, defect at R_d={Rd.tolist()} (recentred); C_N = {C_N*1e3:+.4f} meV (config alignment)", flush=True)
A_uc, _ = qe_io.get_A_volume(dp["uc"]); ws = ws_images(Rn, np.zeros(3, int), MP, A_uc)                # R10 (P-c4) : images de Wigner-Seitz, cellule en bohr
print(f"[{S}] Wigner-Seitz: {int((ws['n_tie'] > 1).sum())} labels with tied images, {len(ws['w'])} images for {nR} labels", flush=True)
kf = lt.mp_grid(a.nf, a.nf, 1); nk = len(kf)
# pi/pi* eigenvectors of H_W(k) on the fine grid (same gauge for reference and truncations)
_, Ef, Uf = Hwr_to_Hwk(Hwr, Rw, kf, ndegen=nd)                         # Uf (nk, nw, nw), bands sorted by energy
gap = Ef[:, 4] - Ef[:, 3]; iD = int(np.argmin(gap)); ED = 0.5 * (Ef[iD, 3] + Ef[iD, 4]); print(f"[{S}] E_D (fine grid) = {ED:.4f} eV", flush=True)
Upi = np.ascontiguousarray(Uf[:, :, 3:5])                             # (nk, nw, 2)
# Fourier matrices: Mwk[w k', w' k] = sum_{R,R'} e^{-2pi i k'.R} Mwr[w R, w' R'] e^{+2pi i k.R'}   (Mwr_to_Mwk convention)
Pm = ws_phase(kf, ws, nR, -1)                                          # (nk, nR): conj side for k' (Wigner-Seitz images, R10)
Pp = ws_phase(kf, ws, nR, +1)                                          # (nk, nR) for k
def fine_pi(mask):
    """pi/pi* blocks (2, nk, 2, nk) of M on the fine grid from Mwr restricted to R,R' in mask."""
    Wm = Mwr.copy(); keep = np.asarray(mask, bool); Wm[:, ~keep, :, :] = 0; Wm[:, :, :, ~keep] = 0
    A = Wm.transpose(1, 0, 3, 2).reshape(nR * nw, nR * nw)             # index R*nw+w
    # right multiply by P_k (R' -> k), then left by P_k'^* (R -> k'), band-rotate to pi/pi* per k
    X = A.reshape(nR * nw, nR, nw) ; Y = np.einsum("aRw,kR->akw", X, Pp, optimize=True).reshape(nR, nw, nk, nw)   # (R, w, k, w')
    Y = np.einsum("Rwkv,kvn->Rwkn", Y, Upi, optimize=True)               # rotate ket to pi/pi*: (R, w, k, n)
    Z = np.einsum("Rwkn,KR->Kwkn", Y, Pm, optimize=True)                 # (k', w, k, n)
    return np.einsum("Kwkn,Kwm->mKnk", Z, Upi.conj(), optimize=True)     # bra rotation: (m, k', n, k)
dist = np.linalg.norm(Rn, axis=1)
ref = fine_pi(np.ones(nR, bool)); mref = np.abs(ref).max(); d_ref = np.array([[ref[n, i, n, i] for i in range(nk)] for n in range(2)]); mdref = np.abs(d_ref).max()
sv_ref = np.linalg.svd(ref.transpose(1, 3, 0, 2), compute_uv=False)     # (k', k, 2)
print(f"[{S}] reference (all {nR} R): max|M_pi| = {mref:.4f} eV, max|M_nn(k,k)| = {mdref:.4f} eV, mean sigma_max = {sv_ref[:, :, 0].mean():.4f} eV", flush=True)
rows = []
for rc in [int(x) for x in a.rcuts.split(",")]:
    mask = dist <= rc + 1e-9; Mt = fine_pi(mask)
    e_max = np.abs(Mt - ref).max() / mref
    sv = np.linalg.svd(Mt.transpose(1, 3, 0, 2), compute_uv=False); e_sv = (np.abs(sv - sv_ref).max(-1) / np.maximum(sv_ref[:, :, 0], 1e-12)).max()
    d_t = np.array([[Mt[n, i, n, i] for i in range(nk)] for n in range(2)]); e_diag = np.abs(d_t - d_ref).max() / mdref
    e_diag_sv = np.abs(np.abs(d_t) - np.abs(d_ref)).max() / mdref
    print(f"[{S}] R_cut={rc}: {mask.sum():3d} sites; max|dM|/max|M| = {e_max:.3e}; SV mismatch (pi/pi* 2x2, all k',k) = {e_sv:.3e}; diagonal M_nn(k,k): max|d|/max = {e_diag:.3e} (|.|: {e_diag_sv:.3e})", flush=True)
    rows.append([S, rc, int(mask.sum()), a.nf, f"{e_max:.4e}", f"{e_sv:.4e}", f"{e_diag:.4e}", f"{e_diag_sv:.4e}", f"{mref:.4f}", f"{mdref:.4f}"])
f = f"{RES}/m_rcut_convergence.csv"; new = not os.path.exists(f)
with open(f, "a", newline="") as fh:
    w = csv.writer(fh)
    if new: w.writerow(["size", "R_cut", "n_sites", "fine_grid", "max_dM_over_maxM", "sv_mismatch_pi_blocks", "diag_max_dM_over_max", "diag_abs_mismatch", "maxM_pi_eV", "maxMdiag_eV"])
    w.writerows(rows)
print(f"appended {len(rows)} rows to {f}")
