#!/usr/bin/env python
"""
rcut_resigma.py -- complement to the R_cut truncation study (§4.1.5, reference size): on-shell self-energy
Sigma_nk = <nk| t(eps_nk) |nk> (COMPLEX) of the local Wannier t-matrix, per R_cut, on the frozen (grid, eta).
  Re Sigma_nk  : level shift (real part of the truncated tail)
  Gamma_nk     = -2 Im Sigma_nk  (same quantity as compute_spectral_wannier.py, per defect)
Reports medians over the +-e_window_eV window around E_D. Same loading/gauge/recentering chain and the same
on-shell nearest-energy-grid evaluation as compute_spectral_wannier.py. Output: npz tagged units='eV', E_D.
Usage: rcut_resigma.py --size 9x9 --rcut 0,1,2,3 [--grid 240 --eta 0.02] --out results/M/resigma_9x9_rc0123.npz
"""
import argparse, numpy as np
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from electron_defect_interaction.config import load_production, dense_paths

p = argparse.ArgumentParser(); p.add_argument("--size", required=True); p.add_argument("--rcut", required=True)
p.add_argument("--grid", type=int, default=None); p.add_argument("--eta", type=float, default=None); p.add_argument("--out", required=True)
a = p.parse_args(); cfg = load_production()
N = a.grid or int(cfg["grid"]); eta = a.eta or float(cfg["eta_eV"]); nk_int = int(cfg["nk_int"]); ew = float(cfg["e_window_eV"]); npe = int(cfg["ne_per_eta"])
dp = dense_paths(cfg, a.size); paths = wannier_provenance.load_wannier_checked(dp["manifest"]); print(f"[gauge] provenance OK: {dp['manifest']}", flush=True)
M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)
k_coarse = qe_io.get_k_red(dp["uc"])
U, k_U = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(k_U, k_coarse)]
U_dis, k_Ud = read_w90_mat(paths["u_dis"]); U_dis = U_dis[_match_kpoint_order(k_Ud, k_coarse)]
Hwr, Rw, ndegen = read_w90_HR(paths["tb"])
Mwk = Mbk_to_Mwk(M, U, U_dis); MP = _infer_mp_grid(k_coarse); Mwr, R_mwr = Mwk_to_Mwr(Mwk, k_coarse, MP)
R_mwr, R_d = lt.recenter_mwr(Mwr, R_mwr, MP); lt.mwr_locality(Mwr, R_mwr); del M, Mwk
print(f"[recenter] R_d={R_d.tolist()}; MP={MP}", flush=True)
_, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=ndegen); gap = E_ref[:, 4] - E_ref[:, 3]
iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4])); win = (E_D - ew, E_D + ew)
print(f"[dirac] E_D = {E_D:.4f} eV; window {win}; grid {N}, eta {eta}, nk_int {nk_int}", flush=True)

k_int = lt.mp_grid(nk_int, nk_int, 1); k_out = lt.mp_grid(N, N, 1)
Hwk_int, _, _ = lt.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=ndegen); _, E_out, U_out = lt.Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)
nw = Hwr.shape[1]; sel = (E_out >= win[0]) & (E_out <= win[1]); E_sel = E_out[sel]
de = eta / npe; egrid = np.arange(E_sel.min() - eta, E_sel.max() + eta + de, de)
jn = np.abs(egrid[None, None, :] - E_out[:, :, None]).argmin(-1)                        # nearest grid energy per (k, n)
out = dict(size=a.size, E_D=E_D, eta=eta, grid=N, nk_int=nk_int, e_window=ew, units="eV", E_out=E_out.T, egrid=egrid,
           note="Sigma[n,k] = <nk|t(eps_nk)|nk> per defect (intensive V_loc), Gamma = -2 Im Sigma; energies in eV, E_D = Wannier Dirac point")
print(f"\n{'Rcut':>5} {'nL':>4} {'dim':>5} {'med ReS(meV)':>13} {'med |ReS|(meV)':>15} {'med Gamma(meV)':>15} {'mean ReS':>10} {'mean Gamma':>11} {'|ReS|/Gamma med':>16} {'E_res-E_D':>10}")
for rc in [float(x) for x in a.rcut.split(",")]:
    Rloc = R_mwr[np.linalg.norm(R_mwr, axis=1) <= rc + 1e-9]; V_loc, _ = lt.extract_V_loc(Mwr, R_mwr, Rloc); nL = len(Rloc)
    ph = lt._phase(k_out, Rloc); phi = np.einsum("kL,kwn->knLw", ph, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)
    g0_all = lt.local_green_batch(Hwk_int, k_int, Rloc, egrid, eta)
    t_cache = [lt.local_t(V_loc, g0_all[j]) for j in range(len(egrid))]; del g0_all
    Sig = np.full((nw, len(k_out)), np.nan + 0j)
    for ik in range(len(k_out)):
        for n in range(nw):
            if sel[ik, n]:
                v = phi[ik, n]; Sig[n, ik] = v.conj() @ t_cache[jn[ik, n]] @ v
    G = -2 * Sig.imag; R = Sig.real; m = sel.T                                  # states inside the window only
    assert G[m].min() > -1e-8, f"positivity: min Gamma {G[m].min():.2e}"
    E = E_out.T; mm = m & (np.abs(E - E_D) <= 1.5); e_res = float(E[mm][np.argmax(G[mm])] - E_D)
    out[f"Sigma_rc{rc:.0f}"] = Sig; out[f"nL_rc{rc:.0f}"] = nL
    print(f"{rc:>5.0f} {nL:>4} {nL*nw:>5} {np.median(R[m])*1e3:>13.2f} {np.median(np.abs(R[m]))*1e3:>15.2f} {np.median(G[m])*1e3:>15.2f} "
          f"{R[m].mean()*1e3:>10.2f} {G[m].mean()*1e3:>11.2f} {np.median(np.abs(R[m])/G[m]):>16.3f} {e_res:>10.3f}", flush=True)
    del t_cache, phi
np.savez(a.out, **out); print(f"saved {a.out}")
