#!/usr/bin/env python
"""
resonance_metrics.py -- reference size (frozen config), +-e_window around E_D, exact local t-matrix:
  (a) Gamma_T(eps)/rho0(eps) on-shell (∝ |T|^2)          (b) delta_rho(eps) = (1/pi) Im Tr[t(eps) dg0/deps], rho_dis at c
  (c) Tbar_nn(K,K;eps) for the pi pair (Re/Im, pole crossing)
  Born vs T: Gamma_Born(eps) = -2 Im <V g0 V> (2nd-order term; the 1st-order <V> is real, contributes 0) vs Gamma_T(eps)
  Potential alignment: Gamma_T(eps) with M and with M - mean_diag(M^L) * 1 (G~=0 component removed)
All energies eV. Output npz: results/M/resonance_<size>.npz. Hard gauge gate + frozen config.
Usage: python scripts/resonance_metrics.py [--size 9x9] [--rho0-grid 900]
"""
import argparse, numpy as np
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from electron_defect_interaction.config import load_production, dense_paths, HA2EV

ap = argparse.ArgumentParser(); ap.add_argument("--size", default=None); ap.add_argument("--rho0-grid", type=int, default=900)
ap.add_argument("--out", default=None); a = ap.parse_args()
cfg = load_production(); S = a.size or cfg["reference_size"]; dp = dense_paths(cfg, S)
paths = wannier_provenance.load_wannier_checked(dp["manifest"]); print(f"[gauge] provenance OK: {dp['manifest']}", flush=True)
rc, N, eta, nk_int, ew, npe, conc = cfg["R_cut"], cfg["grid"], cfg["eta_eV"], cfg["nk_int"], cfg["e_window_eV"], cfg["ne_per_eta"], cfg["defect_concentration_for_dos"]

M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL) * HA2EV          # Ha -> eV, once
ML_diag_mean = float(np.mean([np.load(dp["mfile"].replace("M_dense_", "M_L_dense_"), mmap_mode="r")[n, k, n, k].real
                              for n in range(M.shape[0]) for k in range(0, M.shape[1], 7)])) * HA2EV
print(f"[align] mean diag of M^L (G~=0 component) = {ML_diag_mean*1e3:.2f} meV (unit-cell norm)", flush=True)
k_coarse = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k_coarse)
U, kU = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(kU, k_coarse)]
Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_mp := _match_kpoint_order(kUd, k_coarse)]
Hwr, Rw, nd = read_w90_HR(paths["tb"])

def vloc_from(Mb):
    Mwk = Mbk_to_Mwk(Mb, U, Ud); Mwr, R = Mwk_to_Mwr(Mwk, k_coarse, MP); Rn, Rd = lt.recenter_mwr(Mwr, R, MP)
    lt.mwr_locality(Mwr, Rn)
    Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]
    V, res = lt.extract_V_loc(Mwr, Rn, Rloc); return V, Rloc
V, Rloc = vloc_from(M)
eye = np.zeros_like(M); nb, nk = M.shape[:2]
for n in range(nb):
    for k in range(nk): eye[n, k, n, k] = 1.0
V_sub, _ = vloc_from(M - ML_diag_mean * eye); del eye
nL, nw = len(Rloc), Hwr.shape[1]; print(f"[cluster] R_cut={rc}: nL={nL} sites, dim={nL*nw}", flush=True)

# bands, Dirac point, output states
k_int = lt.mp_grid(nk_int, nk_int, 1); Hwk_int, _, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=nd)
_, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
k_out = lt.mp_grid(N, N, 1); _, E_out, U_out = lt.Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=nd)
ph_out = lt._phase(k_out, Rloc); phi = np.einsum("kL,kwn->knLw", ph_out, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)
sel = np.abs(E_out - E_D) <= ew
print(f"[dirac] E_D = {E_D:.4f} eV; {int(sel.sum())} on-shell states in +-{ew} eV on {N}x{N}", flush=True)

# energy grid, g0 and dg0/de (exact, batched)
de = eta / npe; egrid = np.arange(E_D - ew - eta, E_D + ew + eta + de, de); nE = len(egrid)
g0 = lt.local_green_batch(Hwk_int, k_int, Rloc, egrid, eta); g0p = lt.local_green_batch(Hwk_int, k_int, Rloc, egrid, eta, deriv=True)
print(f"[g0] {nE} energies, spacing {de*1e3:.2f} meV", flush=True)
I = np.eye(nL * nw)
t_T = np.array([V @ np.linalg.solve(I - g0[j] @ V, I) for j in range(nE)])
t_B = np.array([V + V @ g0[j] @ V for j in range(nE)])                          # Born, first two terms
t_S = np.array([V_sub @ np.linalg.solve(I - g0[j] @ V_sub, I) for j in range(nE)])
drho = np.array([(1.0 / np.pi) * np.trace(t_T[j] @ g0p[j]).imag for j in range(nE)])   # per defect, per spin, states/eV
print(f"[drho] integral over the window = {np.trapz(drho, egrid):+.4f} states (Friedel: -> 0 over the full band)", flush=True)

# per-state on-shell rates (nearest energy grid point, as in production)
def onshell(tc):
    g = np.full(E_out.shape, np.nan); j_all = np.rint((E_out - egrid[0]) / de).astype(int)
    for j in np.unique(j_all[sel]):
        m = sel & (j_all == j); P = phi[m]                      # (m, dim)
        g[m] = -2.0 * np.einsum("mi,mi->m", P.conj() @ tc[j], P).imag
    return g
G_T, G_B, G_S = onshell(t_T), onshell(t_B), onshell(t_S)
# energy-resolved (Lorentzian-weighted average of the on-shell states) and pristine DOS
eg = egrid[::2]
def lor(x): return (eta / np.pi) / (x * x + eta * eta)
Es, W = E_out[sel], None
def eres(g):
    w = lor(eg[:, None] - Es[None, :]); return (w * g[sel][None, :]).sum(1) / w.sum(1)
GT_e, GB_e, GS_e = eres(G_T), eres(G_B), eres(G_S)
rho0_240 = np.array([lor(e - E_out).sum() / len(k_out) for e in eg])
kf = lt.mp_grid(a.rho0_grid, a.rho0_grid, 1); _, E_f, _ = lt.Hwr_to_Hwk(Hwr, Rw, kf, ndegen=nd)
rho0 = np.array([lor(e - E_f).sum() / len(kf) for e in eg])                       # per unit cell, per spin
drho_e = np.interp(eg, egrid, drho); rho_dis = rho0 + conc * drho_e
ratio = GT_e / rho0
# (c) Tbar at K for the pi pair
iK = int(np.argmin(np.linalg.norm(np.mod(k_out - np.array([1/3, 1/3, 0]) + 0.5, 1) - 0.5, axis=1)))
PK = phi[iK, 3:5]                                                                  # (2, dim)
Tbar = np.array([PK.conj() @ t_T[j] @ PK.T for j in range(nE)])                    # (nE, 2, 2)
tr = 0.5 * np.trace(Tbar, axis1=1, axis2=2)
zc = np.where(np.diff(np.sign(tr.real)) != 0)[0]
print(f"[K] k_out[{iK}]={k_out[iK]}, E(pi,pi*)@K = {E_out[iK,3]:.4f}, {E_out[iK,4]:.4f} eV", flush=True)

def peak(x, y, lo=None):
    m = np.isfinite(y); i = int(np.argmax(y[m])); return float(x[m][i])
res = {"E_D": E_D, "peak_ratio": peak(eg, ratio) - E_D, "peak_drho": peak(egrid, drho) - E_D, "peak_rho_dis": peak(eg, rho_dis) - E_D,
       "peak_GT": peak(eg, GT_e) - E_D, "peak_GB": peak(eg, GB_e) - E_D,
       "Tbar_zero_crossings": [float(egrid[i] - E_D) for i in zc], "peak_absTbar": float(egrid[int(np.argmax(np.abs(tr)))] - E_D),
       "peak_ImTbar": float(egrid[int(np.argmax(-tr.imag))] - E_D)}
print("\n=== RESONANCE (eV rel. E_D; >0 above Dirac) ===")
for k, v in res.items(): print(f"  {k:22s} {v}")
m = np.abs(eg - E_D) <= ew
print(f"\n=== Born vs T (+-{ew} eV): median Gamma_T = {np.nanmedian(GT_e[m])*1e3:.2f} meV, median Gamma_Born = {np.nanmedian(GB_e[m])*1e3:.2f} meV, "
      f"ratio Born/T median {np.nanmedian(GB_e[m]/GT_e[m]):.3f}, min {np.nanmin(GB_e[m]/GT_e[m]):.3f}, max {np.nanmax(GB_e[m]/GT_e[m]):.3f}")
print(f"=== alignment: Gamma_T with vs without G~=0 shift ({ML_diag_mean*1e3:.1f} meV): max rel |diff| = {np.nanmax(np.abs(GS_e[m]-GT_e[m])/GT_e[m]):.3e}, "
      f"median rel = {np.nanmedian(np.abs(GS_e[m]-GT_e[m])/GT_e[m]):.3e}; per-state max rel = {np.nanmax(np.abs(G_S-G_T)/np.abs(G_T)):.3e}")
out = a.out or f"results/M/resonance_{S}.npz"
np.savez(out, size=S, E_D=E_D, eta=eta, R_cut=rc, grid=N, nk_int=nk_int, conc=conc, egrid=egrid, eg=eg,
         Gamma_T=GT_e, Gamma_Born=GB_e, Gamma_T_noshift=GS_e, rho0=rho0, rho0_240=rho0_240, ratio=ratio, drho=drho, rho_dis=rho_dis,
         Tbar=Tbar, Tbar_tr=tr, E_out=E_out, G_T=G_T, G_B=G_B, G_S=G_S, ML_diag_mean=ML_diag_mean, **{k: v for k, v in res.items() if k != "Tbar_zero_crossings"},
         Tbar_zero_crossings=np.array(res["Tbar_zero_crossings"]))
print(f"saved {out}")
