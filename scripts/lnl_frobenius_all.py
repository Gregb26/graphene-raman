#!/usr/bin/env python
"""
lnl_frobenius_all.py -- <|M^NL|_F>/<|M^L|_F> on the full pi/pi* subspace (2x2 blocks, ALL (k',k) pairs of the dense
grid) for every size with dense M^L / M^NL files (tab:L_NL). Same metric as the "[BZ avg]" block of analyze_M.py
(which only does the reference size). No new physics run; reads the dense Hartree M files (mmap) and the dense
wannierization U, U_dis to select the pi/pi* pair per k by pz weight. Output: results/M/lnl_frobenius.csv.
"""
import csv, os, sys, time, numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat
from electron_defect_interaction.wannier.wannier_interpolation import _match_kpoint_order
from electron_defect_interaction.config import load_production, dense_paths, HA2EV

cfg = load_production()
SIZES = sys.argv[1].split(",") if len(sys.argv) > 1 else ["5x5", "6x6", "7x7", "8x8", "9x9", "12x12"]
OUT = "results/M/lnl_frobenius.csv"

def mmap_M(path):
    meta = matrix_io.read_manifest(path); assert meta and meta.get("units") == matrix_io.HARTREE, f"{path}: untagged/non-Hartree sidecar"
    return np.load(path, mmap_mode="r")
def wannier_V(wdir, k):
    U, kU = read_w90_mat(f"{wdir}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(f"{wdir}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
    return np.einsum("kbw,kwv->kbv", Ud, U)                                 # V[k, band, wf]
def pi_pair(V, eps, ik):
    w = np.abs(V[ik][:, 3]) ** 2 + np.abs(V[ik][:, 4]) ** 2
    top = np.argsort(-w)[:2]; return tuple(sorted(top, key=lambda n: eps[ik, n])), w[top]

rows = []
for S in SIZES:
    dp = dense_paths(cfg, S); fL_path = dp["mfile"].replace("M_dense_", "M_L_dense_"); fN_path = dp["mfile"].replace("M_dense_", "M_NL_dense_")
    if not (os.path.exists(fL_path) and os.path.exists(fN_path)):
        print(f"[skip] {S}: no dense L/NL files"); continue
    t0 = time.time()
    kd = qe_io.get_k_red(dp["uc"]); _, eps = qe_io.get_k_eigenvalues(dp["uc"], False); eps = np.asarray(eps)
    if eps.shape[0] != len(kd): eps = eps.T
    eps = eps * HA2EV
    V = wannier_V(dp["wdir"], kd); nk = len(kd)
    pi_idx = np.zeros((nk, 2), int); wmin = np.zeros(nk)
    for ik in range(nk):
        (a, b), w = pi_pair(V, eps, ik); pi_idx[ik] = (a, b); wmin[ik] = w.min()
    ML = mmap_M(fL_path); MN = mmap_M(fN_path); nb = ML.shape[0]
    assert ML.shape == MN.shape == (nb, nk, nb, nk), (ML.shape, nk)
    fL = np.zeros((nk, nk)); fN = np.zeros((nk, nk)); jj = np.arange(nk)[:, None, None]
    for ik in range(nk):
        rowL = np.asarray(ML[:, :, :, ik]) * HA2EV; rowN = np.asarray(MN[:, :, :, ik]) * HA2EV      # (nb, nk', nb)
        bl = rowL[pi_idx[:, :, None], jj, pi_idx[ik][None, None, :]]                             # (nk', 2, 2)
        bn = rowN[pi_idx[:, :, None], jj, pi_idx[ik][None, None, :]]
        fL[:, ik] = np.linalg.norm(bl.reshape(nk, 4), axis=1); fN[:, ik] = np.linalg.norm(bn.reshape(nk, 4), axis=1)
    r = dict(size=S, N=int(S.split("x")[0]), D=dp["D"], nk=nk, nb=nb, min_pz_weight=float(wmin.min()),
             ratio_full=fN.mean() / fL.mean(), mean_fN_eV=fN.mean(), mean_fL_eV=fL.mean(),
             ratio_diag=np.diag(fN).mean() / np.diag(fL).mean(), ratio_min=(fN / fL).min(), ratio_max=(fN / fL).max(),
             ratio_median_pairs=float(np.median(fN / fL)))
    rows.append(r)
    print(f"[{S}] D={dp['D']} nk={nk} nb={nb}: <|M^NL|_F>/<|M^L|_F> = {r['ratio_full']:.3f} (<F_N> {r['mean_fN_eV']:.4f} eV, <F_L> {r['mean_fL_eV']:.4f} eV); "
          f"diag k'=k {r['ratio_diag']:.3f}; pairwise ratio min {r['ratio_min']:.3f} median {r['ratio_median_pairs']:.3f} max {r['ratio_max']:.3f}; "
          f"min pz weight {r['min_pz_weight']:.2f}; {time.time()-t0:.0f} s", flush=True)
    del ML, MN
if rows:
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT}")
