#!/usr/bin/env python
"""
ks_reconstruction_all.py -- KS reconstruction test with the PRISTINE supercell potential (test C of
test_ks_reconstruction.py) on every unit-cell .save: coarse N×N (defect_N.save) and dense D×D
(defect_uc_dense_D.save), N = 5, 7, 8, 9.
    H_mn(k) = T_mn(k) + <psi_mk|V_p|psi_nk> + V^NL_mn(k)   must equal   eps_nk delta_mn (+ constant offset)
with the kinetic energy T = 1/2 |k+G|^2 added explicitly. k-chunked (psi never materialised for all k).
Reports per case: offset, max/mean |diag - eps - offset|, per-band max/median (meV), off-diagonal max
(absolute, and relative to max|eps|). Saves results/M/ks_reconstruction.npz (figure fig_ks_reconstruction).
"""
import numpy as np
from scipy.interpolate import CubicSpline
from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.pseudo_io import read_upf, fq_from_fr
from electron_defect_interaction.utils.planewaves import mask_invalid_G
from electron_defect_interaction.utils.lattice import red_to_cart
from electron_defect_interaction.defects.non_local import build_K_vectors, compute_phase, compute_angular_part
from electron_defect_interaction.wavefunctions.wfk import compute_psi_nk
from electron_defect_interaction.config import load_production, dense_paths, HA2EV
import sys; sys.path.insert(0, "scripts"); from test_ks_reconstruction import sc_pot_on_uc_grid

cfg = load_production(); DATA = "data/graphene"; out = {}
import argparse, os
ap = argparse.ArgumentParser(); ap.add_argument("--sizes", default="5x5,7x7,8x8,9x9"); ap.add_argument("--merge", action="store_true", help="merge into an existing results/M/ks_reconstruction.npz")
args = ap.parse_args()
if args.merge and os.path.exists("results/M/ks_reconstruction.npz"): out.update(dict(np.load("results/M/ks_reconstruction.npz")))

def sc_pot_on_uc_grid_fourier(uc_save, sc_save, pot_sc_file):
    """Non-commensurate FFT grids: restrict V_p to its unit-cell-periodic Fourier components G_sc = N.G_uc and
    resample on the unit-cell grid (exact for a unit-cell-periodic V_p). Diagnostic = norm of the dropped
    (non-periodic) components relative to the kept ones."""
    A_uc, _ = qe_io.get_A_volume(uc_save); A_sc, _ = qe_io.get_A_volume(sc_save)
    N = np.rint(np.diag(np.linalg.inv(A_uc) @ A_sc)).astype(int); n_uc = np.array(qe_io.get_ngfft(uc_save))
    V_sc, _ = qe_io.get_pot(pot_sc_file, subtract_mean=False); V_sc = V_sc.transpose(2, 1, 0); n_sc = np.array(V_sc.shape)
    Vg = np.fft.fftn(V_sc) / V_sc.size; Vuc_g = np.zeros(tuple(n_uc), complex); kept = 0.0
    for m1 in range(-(n_uc[0] // 2), (n_uc[0] - 1) // 2 + 1):
        for m2 in range(-(n_uc[1] // 2), (n_uc[1] - 1) // 2 + 1):
            for m3 in range(-(n_uc[2] // 2), (n_uc[2] - 1) // 2 + 1):
                i = (np.array([m1, m2, m3]) * N)
                if np.all(np.abs(i) <= (n_sc - 1) // 2):
                    Vuc_g[m1 % n_uc[0], m2 % n_uc[1], m3 % n_uc[2]] = Vg[i[0] % n_sc[0], i[1] % n_sc[1], i[2] % n_sc[2]]
    # dropped = everything not on the N-sublattice of G (non-periodic part), measured in real space
    fsigned = [np.where(np.arange(n) < (n + 1) // 2, np.arange(n), np.arange(n) - n) for n in n_sc]      # signed FFT frequencies
    mask = np.zeros(tuple(n_sc), bool); idx = np.ix_(*[f % Ni == 0 for f, Ni in zip(fsigned, N)]); mask[idx] = True
    dropped = np.fft.ifftn(np.where(mask, 0.0, Vg) * V_sc.size).real
    V_uc = np.fft.ifftn(Vuc_g).real * Vuc_g.size
    return V_uc, float(np.abs(dropped).max())

def reconstruct(uc_save, V_uc, upf_file, k_chunk=27):
    C_nkg, nG = qe_io.get_C_nk(uc_save); G_red = qe_io.get_G_red(uc_save); k_red = qe_io.get_k_red(uc_save)
    B_uc, _ = qe_io.get_B_volume(uc_save); A_uc, Om = qe_io.get_A_volume(uc_save); ecut = float(qe_io.get_ecut(uc_save))
    tau = red_to_cart(qe_io.get_x_red(uc_save), A_uc); eps = qe_io.get_eigenvalues(uc_save); ngfft = qe_io.get_ngfft(uc_save)
    nb, nk, _ = C_nkg.shape; keep = mask_invalid_G(nG); C = np.where(keep, C_nkg, 0.0)
    # kinetic
    K = red_to_cart(k_red[:, None, :] + G_red, B_uc); K2 = np.where(keep, np.sum(K ** 2, axis=2), 0.0)
    T = np.einsum("nkg,kg,mkg->knm", C.conj(), K2, C, optimize=True) * 0.5
    # non-local
    ekb_li, fr_li, rgrid, lmax, imax, _ = read_upf(upf_file)
    Kc, Kn, Kh = build_K_vectors(k_red, G_red, keep, B_uc)
    q = np.linspace(0, 2 * np.sqrt(2 * ecut), 2000); Fq = CubicSpline(q, fq_from_fr(rgrid, fr_li, q), axis=-1, extrapolate=False); F = Fq(Kn)
    phase = compute_phase(Kc, tau); Y = compute_angular_part(Kh, lmax); pref = 4 * np.pi / np.sqrt(Om)
    Bp = pref * np.einsum("nkg,likg,kglm,ksg->nkslim", np.conj(C), F, Y, phase, optimize=True)
    NL = np.einsum("li,pkslia,qkslia->kpq", ekb_li, Bp, np.conj(Bp), optimize=True)
    # local, k-chunked
    Nr = int(np.prod(ngfft)); dV = Om / Nr; Vr = np.asarray(V_uc).reshape(Nr); L = np.zeros((nk, nb, nb), complex)
    for s in range(0, nk, k_chunk):
        sl = slice(s, s + k_chunk)
        psi, _ = compute_psi_nk(C_nkg[:, sl], nG[sl], G_red[sl], k_red[sl], Om, ngfft=ngfft)
        for j in range(psi.shape[1]):
            P = psi[:, j].reshape(nb, Nr); L[s + j] = dV * (P.conj() * Vr) @ P.T
        del psi
    return T + L + NL, np.asarray(eps)          # (nk, nb, nb) Ha ; (nb, nk) Ha

def analyse(tag, H, eps):
    nk, nb, _ = H.shape; diag = np.real(np.einsum("kii->ki", H)).T                   # (nb, nk)
    off = H - np.einsum("ki,ij->kij", np.einsum("kii->ki", H), np.eye(nb)); offmax = float(np.abs(off).max())
    raw = diag - eps; shift = float(raw.mean()); dev = raw - shift
    per_band_max = np.abs(dev).max(1) * HA2EV * 1e3; per_band_med = np.median(np.abs(dev), 1) * HA2EV * 1e3
    print(f"[{tag}] nb={nb} nk={nk}: offset diag-eps = {shift*HA2EV:+.4f} eV; max|diag-eps| raw = {np.abs(raw).max()*HA2EV*1e3:.2f} meV; "
          f"after offset: max {np.abs(dev).max()*HA2EV*1e3:.3f} meV, mean {np.abs(dev).mean()*HA2EV*1e3:.3f} meV, median {np.median(np.abs(dev))*HA2EV*1e3:.3f} meV; "
          f"off-diag max {offmax*HA2EV*1e3:.3f} meV = {offmax/np.abs(eps).max():.2e} x max|eps|", flush=True)
    print(f"[{tag}] per band (meV) max: {np.round(per_band_max, 3).tolist()}\n[{tag}] per band (meV) median: {np.round(per_band_med, 3).tolist()}", flush=True)
    out.update(**{f"{tag}_offset_eV": shift * HA2EV, f"{tag}_max_meV": np.abs(dev).max() * HA2EV * 1e3, f"{tag}_mean_meV": np.abs(dev).mean() * HA2EV * 1e3,
                  f"{tag}_median_meV": np.median(np.abs(dev)) * HA2EV * 1e3, f"{tag}_raw_max_meV": np.abs(raw).max() * HA2EV * 1e3,
                  f"{tag}_offdiag_max_meV": offmax * HA2EV * 1e3, f"{tag}_offdiag_rel": offmax / np.abs(eps).max(),
                  f"{tag}_band_max_meV": per_band_max, f"{tag}_band_median_meV": per_band_med, f"{tag}_nb": nb, f"{tag}_nk": nk})

for S in args.sizes.split(","):
    sc_p = f"{DATA}/supercell/qe/defect_{S}_p.save"; pot_p = f"{sc_p}/Vks_{S}_p"; upf = f"{DATA}/unit_cell/qe/defect_{S}.save/C.upf"
    for kind, uc in (("coarse", f"{DATA}/unit_cell/qe/defect_{S}.save"), ("dense", dense_paths(cfg, S)["uc"])):
        tag = f"{S}_{kind}"
        try:
            try:
                V_uc, spread = sc_pot_on_uc_grid(uc, sc_p, pot_p); how = "restriction directe (grilles commensurables)"
            except ValueError as e:
                V_uc, spread = sc_pot_on_uc_grid_fourier(uc, sc_p, pot_p); how = f"restriction de Fourier (grilles non commensurables : {e})"
            print(f"[{tag}] {how}; V_p unit-cell periodicity spread / dropped non-periodic part = {spread*HA2EV*1e3:.3f} meV", flush=True)
            H, eps = reconstruct(uc, V_uc, upf); analyse(tag, H, eps); out[f"{tag}_spread_meV"] = spread * HA2EV * 1e3
        except Exception as e:
            print(f"[{tag}] FAILED: {e}", flush=True)
np.savez("results/M/ks_reconstruction.npz", **out); print("saved results/M/ks_reconstruction.npz")
