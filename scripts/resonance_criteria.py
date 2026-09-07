#!/usr/bin/env python
"""
resonance_criteria.py -- definitive resonance criterion + sum rule + concentration-scaled rate (reference size, frozen config).
  * |det[1 - V_loc g0(eps)]| and the eigenvalue of [1 - V_loc g0(eps)] closest to zero, on +-e_window around E_D:
    local minima (energy rel. E_D, depth) of both.
  * Friedel sum rule: int delta_rho d(eps) over the FULL Wannier bandwidth, delta_rho = (1/pi) Im Tr[t dg0/deps],
    cross-checked with Lloyd's formula delta_rho = -(1/pi) d/deps Im ln det[1 - g0 V]; where the +-e_window excess is compensated.
  * Gamma_T(eps) at c = 0.1 % on +-1 eV (from resonance_<size>.npz), for the order-of-magnitude comparison with
    Kaasbjerg Fig. 17 (vacancy here vs substitutional N there).
Output: results/M/resonance_criteria_<size>.npz
"""
import argparse, numpy as np
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from electron_defect_interaction.config import load_production, dense_paths

ap = argparse.ArgumentParser(); ap.add_argument("--size", default=None); ap.add_argument("--c-compare", type=float, default=1e-3)
ap.add_argument("--band-de", type=float, default=None, help="energy spacing for the full-band sum rule (default eta/4)"); a = ap.parse_args()
cfg = load_production(); S = a.size or cfg["reference_size"]; dp = dense_paths(cfg, S)
paths = wannier_provenance.load_wannier_checked(dp["manifest"]); print(f"[gauge] provenance OK: {dp['manifest']}", flush=True)
rc, eta, nk_int, ew, npe = cfg["R_cut"], cfg["eta_eV"], cfg["nk_int"], cfg["e_window_eV"], cfg["ne_per_eta"]

M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)
k_coarse = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k_coarse)
U, kU = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(kU, k_coarse)]
Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_kpoint_order(kUd, k_coarse)]
Hwr, Rw, nd = read_w90_HR(paths["tb"])
Mwr, R = Mwk_to_Mwr(Mbk_to_Mwk(M, U, Ud), k_coarse, MP); Rn, Rd = lt.recenter_mwr(Mwr, R, MP); lt.mwr_locality(Mwr, Rn); del M
Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]; V, _ = lt.extract_V_loc(Mwr, Rn, Rloc); dim = V.shape[0]
k_int = lt.mp_grid(nk_int, nk_int, 1); Hwk_int, E_int, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=nd)
_, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
print(f"[setup] {S}: R_cut={rc} dim={dim}, E_D={E_D:.4f} eV, Wannier bands span [{E_int.min():.2f}, {E_int.max():.2f}] eV", flush=True)
I = np.eye(dim)

# --- 1. det / eigenvalue criterion on +-ew
de = eta / npe; eg = np.arange(E_D - ew, E_D + ew + de / 2, de)
g0 = lt.local_green_batch(Hwk_int, k_int, Rloc, eg, eta)
A = I[None] - V[None] @ g0                                        # 1 - V g0
sign, logabs = np.linalg.slogdet(A); lam = np.linalg.eigvals(A)
minlam = np.abs(lam).min(1); ilam = np.abs(lam).argmin(1); lam_min = lam[np.arange(len(eg)), ilam]
logdet_rel = logabs - logabs.max()                                 # log|det| relative to its max on the window
def local_minima(y, x, n=6):
    idx = [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] <= y[i + 1]]
    idx = sorted(idx, key=lambda i: y[i])[:n]; return sorted(idx, key=lambda i: x[i])
print("\n=== criterion 1: local minima of |det[1 - V g0]| on +-%.0f eV (energy rel. E_D; depth = |det|/max|det| on the window)" % ew)
for i in local_minima(logdet_rel, eg):
    print(f"   eps-E_D = {eg[i]-E_D:+.3f} eV   |det|/max = {np.exp(logdet_rel[i]):.3e}   log10 = {logdet_rel[i]/np.log(10):+.2f}   min|lambda| there = {minlam[i]:.4f}")
print("=== criterion 2: local minima of min_i |lambda_i(1 - V g0)| on +-%.0f eV" % ew)
for i in local_minima(minlam, eg):
    print(f"   eps-E_D = {eg[i]-E_D:+.3f} eV   min|lambda| = {minlam[i]:.4f}   lambda = {lam_min[i].real:+.4f}{lam_min[i].imag:+.4f}i")
print(f"   global: min|lambda| = {minlam.min():.4f} at {eg[minlam.argmin()]-E_D:+.3f} eV; min |det|/max = {np.exp(logdet_rel.min()):.3e} at {eg[logdet_rel.argmin()]-E_D:+.3f} eV", flush=True)

# --- 2. Friedel sum rule over the full Wannier bandwidth (chunked; g0, dg0/de, t, Lloyd phase)
deb = a.band_de or eta / 4; lo, hi = E_int.min() - 1.5, E_int.max() + 1.5
eb = np.arange(lo, hi + deb / 2, deb); nb_ = len(eb); drho = np.zeros(nb_); phase = np.zeros(nb_)
print(f"[sum rule] {nb_} energies on [{lo:.2f}, {hi:.2f}] eV, spacing {deb*1e3:.1f} meV", flush=True)
for s0 in range(0, nb_, 400):
    ee = eb[s0:s0 + 400]; gb = lt.local_green_batch(Hwk_int, k_int, Rloc, ee, eta); gp = lt.local_green_batch(Hwk_int, k_int, Rloc, ee, eta, deriv=True)
    for j in range(len(ee)):
        Aj = I - gb[j] @ V; t = V @ np.linalg.solve(Aj, I)
        drho[s0 + j] = np.trace(t @ gp[j]).imag / np.pi
        sg, la = np.linalg.slogdet(Aj); phase[s0 + j] = np.angle(sg)
ph = np.unwrap(phase); drho_lloyd = -np.gradient(ph, eb) / np.pi
tot = np.trapz(drho, eb); tot_l = np.trapz(drho_lloyd, eb)
win = np.abs(eb - E_D) <= ew; inwin = np.trapz(drho[win], eb[win])
print(f"\n=== sum rule: int delta_rho = {tot:+.4f} states (Tr[t g0'])   |   {tot_l:+.4f} (Lloyd)   |   inside +-{ew:.0f} eV: {inwin:+.4f}")
print(f"    max |drho_Tr - drho_Lloyd| = {np.abs(drho - drho_lloyd).max():.3e} states/eV (agreement of the two formulas)")
# where is the compensation: contiguous intervals of sign(drho) with their integrals, sorted by |integral|
sgn = np.sign(drho); edges = np.where(np.diff(sgn) != 0)[0] + 1; bounds = np.concatenate([[0], edges, [nb_]])
segs = [(eb[b0] - E_D, eb[b1 - 1] - E_D, np.trapz(drho[b0:b1], eb[b0:b1])) for b0, b1 in zip(bounds[:-1], bounds[1:])]
segs = sorted(segs, key=lambda s: -abs(s[2]))[:10]
print("    largest contributions (interval rel. E_D, integral):")
for s in sorted(segs): print(f"      [{s[0]:+7.2f}, {s[1]:+7.2f}] eV : {s[2]:+.4f}")
cum = np.concatenate([[0], np.cumsum(0.5 * (drho[1:] + drho[:-1]) * np.diff(eb))])
print(f"    cumulative int at window edges: {np.interp(E_D - ew, eb, cum):+.4f} (at -{ew:.0f} eV), {np.interp(E_D + ew, eb, cum):+.4f} (at +{ew:.0f} eV), {cum[-1]:+.4f} (band top)", flush=True)

# --- 3. Gamma_T at c = c_compare on +-1 eV (from resonance_<size>.npz)
Rz = np.load(f"results/M/resonance_{S}.npz"); x = Rz["eg"] - float(Rz["E_D"]); m1 = np.abs(x) <= 1.0
Gc = Rz["Gamma_T"] * a.c_compare
print(f"\n=== Gamma_T at c = {a.c_compare*100:.1f}% on +-1 eV (VACANCY here vs SUBSTITUTIONAL N in Kaasbjerg Fig. 17: order of magnitude only):")
print(f"    min {np.nanmin(Gc[m1])*1e3:.2f} meV at {x[m1][np.nanargmin(Gc[m1])]:+.2f} eV, max {np.nanmax(Gc[m1])*1e3:.2f} meV at {x[m1][np.nanargmax(Gc[m1])]:+.2f} eV; "
      f"at E_D {np.interp(0, x, Gc)*1e3:.2f} meV, at -0.5 eV {np.interp(-0.5, x, Gc)*1e3:.2f} meV, at +0.5 eV {np.interp(0.5, x, Gc)*1e3:.2f} meV; hbar/Gamma at +-0.3 eV: "
      f"{658.2/ (np.interp(-0.3, x, Gc)*1e3):.0f} fs / {658.2/(np.interp(0.3, x, Gc)*1e3):.0f} fs")
np.savez(f"results/M/resonance_criteria_{S}.npz", size=S, E_D=E_D, eta=eta, R_cut=rc, eg=eg, logdet_rel=logdet_rel, minlam=minlam, lam_min=lam_min,
         eb=eb, drho_band=drho, drho_lloyd=drho_lloyd, cum_drho=cum, sumrule=tot, sumrule_lloyd=tot_l, sumrule_window=inwin,
         x_c=x[m1], Gamma_c=Gc[m1], c_compare=a.c_compare)
print(f"saved results/M/resonance_criteria_{S}.npz")
