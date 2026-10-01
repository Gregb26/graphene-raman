#!/usr/bin/env python
"""
rcut_resigma.py -- complement to the R_cut truncation study (§4.1.5, reference size): on-shell self-energy
Sigma_nk = <nk| t(eps_nk) |nk> (COMPLEX) of the local Wannier t-matrix, per R_cut, on the frozen (grid, eta).
  Re Sigma_nk  : level shift (real part of the truncated tail)
  Gamma_nk     = -2 Im Sigma_nk  (same quantity as compute_spectral_wannier.py, per defect)
Reports medians over the +-e_window_eV window around E_D. Same loading/gauge/recentering chain and the same
on-shell nearest-energy-grid evaluation as compute_spectral_wannier.py. Output: npz tagged units='eV', E_D.
Usage: rcut_resigma.py --size 9x9 --rcut 0,1,2,3 [--grid 240 --eta 0.02] --out <results_dir>/resigma_9x9_rc0123.npz
"""
import argparse, numpy as np
from graphene_raman.io import qe_io, matrix_io, wannier_provenance
from graphene_raman.io.wannier_io import read_w90_mat, read_w90_tb
from graphene_raman.wannier.wannier_interpolation import _infer_mp_grid, _match_kpoint_order
from graphene_raman.wannier.wannier_hamiltonian import dirac_point
from graphene_raman.defects.many_body import cluster_tmatrix as ct
from graphene_raman.config import load_production, dense_paths, alignment_C

p = argparse.ArgumentParser(); p.add_argument("--size", required=True); p.add_argument("--rcut", required=True)
p.add_argument("--grid", type=int, default=None); p.add_argument("--eta", type=float, default=None); p.add_argument("--out", required=True); p.add_argument("--npe", type=int, default=None, help="ne_per_eta override (default: frozen config)")
p.add_argument("--nk-int", type=int, default=None, help="internal k-grid density N (NxN) for g0 override (default: frozen config); convergence test of N_k^int")
a = p.parse_args(); cfg = load_production()
N = a.grid or int(cfg["grid"]); eta = a.eta or float(cfg["eta_eV"]); nk_int = a.nk_int or int(cfg["nk_int"]); ew = float(cfg["e_window_eV"]); npe = a.npe or int(cfg["ne_per_eta"])
dp = dense_paths(cfg, a.size); paths = wannier_provenance.load_wannier_checked(dp["manifest"]); print(f"[gauge] provenance OK: {dp['manifest']}", flush=True)
M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
k_coarse = qe_io.get_k_red(dp["uc"])
U, k_U = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(k_U, k_coarse)]
U_dis, k_Ud = read_w90_mat(paths["u_dis"]); U_dis = U_dis[_match_kpoint_order(k_Ud, k_coarse)]
Hwr, Rw, ndegen, _, _ = read_w90_tb(paths["tb"])
MP = _infer_mp_grid(k_coarse); C_N = alignment_C(cfg, a.size)                           # R10 : M_W(R,R) - C_N sur la boîte (approximation (i))
d = ct.defect_mwr(M, U, U_dis, k_coarse, MP, n_box=int(a.size.split("x")[0]), C_N=C_N); Mwr, R_mwr, R_d = d["Mwr"], d["Rn"], d["R_d"]
ct.mwr_locality(Mwr, R_mwr); del M, d
print(f"[recenter] R_d={R_d.tolist()}; MP={MP}; C_N = {C_N*1e3:+.4f} meV (config alignment)", flush=True)
_, E_ref, _ = ct.Hwr_to_Hwk(Hwr, Rw, ct.mp_grid(90, 90, 1), ndegen=ndegen)
E_D, _ = dirac_point(E_ref); win = (E_D - ew, E_D + ew)
print(f"[dirac] E_D = {E_D:.4f} eV; window {win}; grid {N}, eta {eta}, nk_int {nk_int}, ne_per_eta {npe} (de = {eta/npe*1e3:.3f} meV)", flush=True)

k_int = ct.mp_grid(nk_int, nk_int, 1); k_out = ct.mp_grid(N, N, 1)
Hwk_int, _, _ = ct.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=ndegen); _, E_out, U_out = ct.Hwr_to_Hwk(Hwr, Rw, k_out, ndegen=ndegen)
nw = Hwr.shape[1]; sel = (E_out >= win[0]) & (E_out <= win[1]); E_sel = E_out[sel]
de = eta / npe; egrid = np.arange(E_sel.min() - eta, E_sel.max() + eta + de, de)
jn = np.abs(egrid[None, None, :] - E_out[:, :, None]).argmin(-1)                        # nearest grid energy per (k, n)
out = dict(size=a.size, E_D=E_D, eta=eta, grid=N, nk_int=nk_int, e_window=ew, units="eV", E_out=E_out.T, egrid=egrid,
           ne_per_eta=npe, note="Sigma[n,k] = <nk|t(eps_nk)|nk> per defect (intensive M_cluster), Gamma = -2 Im Sigma; energies in eV, E_D = Wannier Dirac point")
print(f"\n{'Rcut':>5} {'nL':>4} {'dim':>5} {'med ReS(meV)':>13} {'med |ReS|(meV)':>15} {'med Gamma(meV)':>15} {'mean ReS':>10} {'mean Gamma':>11} {'|ReS|/Gamma med':>16} {'E_res-E_D':>10}")
for rc in [float(x) for x in a.rcut.split(",")]:
    R_cluster = ct.cluster_cells(R_mwr, rc); M_cluster, _ = ct.cluster_potential(Mwr, R_mwr, R_cluster); nL = len(R_cluster)
    ph = ct._phase(k_out, R_cluster); phi = np.einsum("kL,kwn->knLw", ph, U_out, optimize=True).reshape(len(k_out), nw, nL * nw)
    g0_all = ct.cluster_green_batch(Hwk_int, k_int, R_cluster, egrid, eta)
    t_cache = [ct.cluster_t(M_cluster, g0_all[j]) for j in range(len(egrid))]; del g0_all
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
