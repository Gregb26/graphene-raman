#!/usr/bin/env python
"""
analyze_M.py -- post-processing of the M matrix for §4.1.5 (no new physics runs). Produces results/M/M_analysis.npz
and results/M/M_tests_summary.csv, read by scripts/make_figures.py (fig_M_map, fig_Ved_boundary, fig_M_scaling).
  1. |M| map on the BZ for k=K (reference size, dense grid), Kaasbjerg Fig. 3 convention: Vt = A_cell |M| (eV A^2),
     unit-cell Bloch normalization; pi/pi* selected by pz Wannier weight; gauge-invariant sum over the degenerate pair at K.
  2. V_ed^L along a lattice-vector line from the vacancy to the supercell boundary, four N; boundary/site ratio.
  3. Re M^L and Re M^NL at (the grid point nearest to) K, pi-pair trace/2, four N.
  4. max|M| (bands 1-16) vs N and max|M| N_cells vs N, dense and coarse.
  5. tests table: hermiticity, padding (plumbing, coincident-k SV, nb20 vs nb16), Wannier closure, scaling.
"""
import csv, json, numpy as np
from scipy.ndimage import map_coordinates
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, Mwr_to_Mwk, Mwk_to_Mbk, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.config import load_production, dense_paths, HA2EV
BOHR = 0.529177210903
cfg = load_production(); SIZES = ["5x5", "7x7", "8x8", "9x9"]; REF = cfg["reference_size"]
DATA = "data/graphene"; out = {}; tests = []

def mmap_M(path):
    """mmap a Hartree-tagged M without loading it (gate: manifest must say units=hartree)."""
    meta = matrix_io.read_manifest(path); assert meta and meta.get("units") == matrix_io.HARTREE, f"{path}: untagged/non-Hartree sidecar"
    return np.load(path, mmap_mode="r")
def wannier_V(wdir, k):
    U, kU = read_w90_mat(f"{wdir}/wannier_u.mat"); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(f"{wdir}/wannier_u_dis.mat"); Ud = Ud[_match_kpoint_order(kUd, k)]
    return np.einsum("kbw,kwv->kbv", Ud, U), U, Ud                       # V[k, band, wf]
def pi_pair(V, eps, ik):
    w = (np.abs(V[ik][:, 3]) ** 2 + np.abs(V[ik][:, 4]) ** 2)          # pz(A)+pz(B) weight per DFT band
    top = np.argsort(-w)[:2]; return tuple(sorted(top, key=lambda n: eps[ik, n])), w[top]
def kdist(k, k0):
    d = np.mod(k - np.asarray(k0) + 0.5, 1.0) - 0.5; return np.linalg.norm(d, axis=1)

# ---------------- 1. |M| map at k=K, reference size, dense grid
dp = dense_paths(cfg, REF); kd = qe_io.get_k_red(dp["uc"]); kk, eps = qe_io.get_k_eigenvalues(dp["uc"], False); eps = np.asarray(eps)
if eps.shape[0] != len(kd): eps = eps.T
eps = eps * HA2EV
A_uc, _ = qe_io.get_A_volume(dp["uc"]); A_cell = np.linalg.norm(np.cross(A_uc[:, 0], A_uc[:, 1])) * BOHR ** 2      # A^2
B = 2 * np.pi * np.linalg.inv(A_uc).T / BOHR                                                                        # columns b_i, A^-1
iK = int(np.argmin(kdist(kd, tuple(cfg["K_red"])))); V9, U9, Ud9 = wannier_V(dp["wdir"], kd)
pairK, wK = pi_pair(V9, eps, iK)
M = mmap_M(dp["mfile"]); nb, nk = M.shape[:2]
MK = np.array(M[:, :, :, iK]) * HA2EV                                                                               # (nb, nk', nb_K) eV
Vpi = np.zeros(nk); Vps = np.zeros(nk); ipi = np.zeros(nk, int); ips = np.zeros(nk, int); wmin = np.zeros(nk)
for ik in range(nk):
    (a, b), w = pi_pair(V9, eps, ik); ipi[ik], ips[ik], wmin[ik] = a, b, w.min()
    Vpi[ik] = A_cell * np.sqrt(sum(abs(MK[a, ik, n]) ** 2 for n in pairK)); Vps[ik] = A_cell * np.sqrt(sum(abs(MK[b, ik, n]) ** 2 for n in pairK))
kc = (B[:2, :2] @ kd[:, :2].T).T
imgs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)]); Gc = (B[:2, :2] @ imgs.T).T
kfold = np.array([kc[i] + Gc[np.argmin(np.linalg.norm(kc[i] + Gc, axis=1))] for i in range(nk)])
corners = np.array([B[:2, :2] @ np.array(c) for c in [(1/3, 1/3), (1/3, -2/3), (-2/3, 1/3), (-1/3, -1/3), (-1/3, 2/3), (2/3, -1/3)]])
corners = corners[np.argsort(np.arctan2(corners[:, 1], corners[:, 0]))]
print(f"[map] {REF} dense {int(np.sqrt(nk))}x{int(np.sqrt(nk))}: K index {iK} k={kd[iK]}, pair at K = bands {pairK} (pz weight {np.round(wK,3)}), "
      f"A_cell = {A_cell:.3f} A^2; Vt_pi max {Vpi.max():.2f} eV A^2 (at K: {Vpi[iK]:.2f}), Vt_pi* max {Vps.max():.2f} (at K: {Vps[iK]:.2f}); min pz-weight of selected pair {wmin.min():.2f}", flush=True)
out.update(map_kx=kfold[:, 0], map_ky=kfold[:, 1], map_Vpi=Vpi, map_Vpistar=Vps, map_K=kfold[iK], map_corners=corners, map_A_cell=A_cell, map_B=B, map_eps_pair=eps[np.arange(nk), ipi], map_eps_pairstar=eps[np.arange(nk), ips])
# --- L / NL decomposition on the BZ at k=K (pi row), gauge-invariant: component of M^X along the total M
MLK = np.array(mmap_M(dp["mfile"].replace("M_dense_", "M_L_dense_"))[:, :, :, iK]) * HA2EV; MNK = np.array(mmap_M(dp["mfile"].replace("M_dense_", "M_NL_dense_"))[:, :, :, iK]) * HA2EV
Lpar = np.zeros(nk); Npar = np.zeros(nk); Labs = np.zeros(nk); Nabs = np.zeros(nk)
for ik in range(nk):
    v = np.array([MK[ipi[ik], ik, n] for n in pairK]); vL = np.array([MLK[ipi[ik], ik, n] for n in pairK]); vN = np.array([MNK[ipi[ik], ik, n] for n in pairK])
    nv = np.linalg.norm(v); Lpar[ik] = A_cell * (np.vdot(v, vL).real / nv); Npar[ik] = A_cell * (np.vdot(v, vN).real / nv)
    Labs[ik] = A_cell * np.linalg.norm(vL); Nabs[ik] = A_cell * np.linalg.norm(vN)
print(f"[map L/NL] pi row at k=K: <|M^NL|>/<|M^L|> over the BZ = {Nabs.mean()/Labs.mean():.2f}; components along M: <M^NL_par> = {Npar.mean():.2f}, <M^L_par> = {Lpar.mean():.2f} eV A^2 "
      f"(ratio {Npar.mean()/Lpar.mean():.2f}); min over k' of |M^NL|/|M^L| = {(Nabs/Labs).min():.2f}, max = {(Nabs/Labs).max():.2f}; M^L_par < 0 at {int((Lpar<0).sum())} of {nk} k'", flush=True)
# --- full pi-subspace double average over (k',k)
ML = mmap_M(dp["mfile"].replace("M_dense_", "M_L_dense_")); MN = mmap_M(dp["mfile"].replace("M_dense_", "M_NL_dense_"))
pi_idx = np.stack([ipi, ips], 1)                                                                      # (nk, 2) pi/pi* band index per k
fL = np.zeros((nk, nk)); fN = np.zeros((nk, nk))
for ik in range(nk):
    rowL = np.array(ML[:, :, :, ik]) * HA2EV; rowN = np.array(MN[:, :, :, ik]) * HA2EV                 # (nb, nk', nb)
    for jk in range(nk):
        bl = rowL[pi_idx[jk]][:, pi_idx[ik]]; bn = rowN[pi_idx[jk]][:, pi_idx[ik]]
        fL[jk, ik] = np.linalg.norm(bl); fN[jk, ik] = np.linalg.norm(bn)
print(f"[BZ avg] pi subspace (2x2 blocks pi/pi*, all (k',k)): <|M^NL|_F>/<|M^L|_F> = {fN.mean()/fL.mean():.2f}; <|M^NL|_F> = {fN.mean():.4f} eV, <|M^L|_F> = {fL.mean():.4f} eV; "
      f"diagonal k'=k only: {np.diag(fN).mean()/np.diag(fL).mean():.2f}; ratio min {(fN/fL).min():.2f} max {(fN/fL).max():.2f}", flush=True)
out.update(map_Lpar=Lpar, map_Npar=Npar, map_Labs=Labs, map_Nabs=Nabs, bz_ratio_row=Nabs.mean()/Labs.mean(), bz_ratio_par=Npar.mean()/Lpar.mean(),
           bz_ratio_full=fN.mean()/fL.mean(), bz_ratio_diag=np.diag(fN).mean()/np.diag(fL).mean(), bz_fL_mean=fL.mean(), bz_fN_mean=fN.mean())
del MK, MLK, MNK, ML, MN

# ---------------- 2. V_ed^L along a lattice-vector line, four N
for S in SIZES:
    n = int(S[0]); scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; scd = f"{DATA}/supercell/qe/defect_{S}_d.save"
    A_sc, _ = qe_io.get_A_volume(scd); xp = np.mod(qe_io.get_x_red(scp), 1.0); xd = np.mod(qe_io.get_x_red(scd), 1.0)
    dmin = np.array([np.min(np.linalg.norm(np.mod(xd - p + 0.5, 1) - 0.5, axis=1)) for p in xp]); s_vac = xp[int(np.argmax(dmin))]
    Vp, ng = qe_io.get_pot(f"{scp}/Vks_{S}_p", subtract_mean=False, to_hartree=True); Vd, _ = qe_io.get_pot(f"{scd}/Vks_{S}_d", subtract_mean=False, to_hartree=True)
    dV = (Vd - Vp).transpose(2, 1, 0) * HA2EV; nr = np.array(dV.shape)                                     # [ix, iy, iz], eV
    t = np.linspace(0, 0.5, 401); s = s_vac[None, :] + t[:, None] * np.array([1.0, 0, 0])[None, :]
    line = map_coordinates(dV, (s * nr[None, :]).T, order=1, mode="wrap")
    site = float(map_coordinates(dV, (s_vac * nr)[:, None], order=1, mode="wrap")[0])
    i1 = int(np.round((s_vac[0] + 0.5) * nr[0])) % nr[0]; i2 = int(np.round((s_vac[1] + 0.5) * nr[1])) % nr[1]
    b1 = float(np.abs(dV[i1, :, :]).max()); b2 = float(np.abs(dV[:, i2, :]).max()); bmax = max(b1, b2)
    a1_len = np.linalg.norm(A_sc[:, 0]) * BOHR
    print(f"[Ved] {S}: vacancy s={np.round(s_vac,4)}, dV(site)={site:+.3f} eV, |dV| on the boundary planes: a1 {b1:.4f}, a2 {b2:.4f} eV -> max/site = {bmax/abs(site):.2e}; "
          f"half-box {0.5*a1_len:.2f} A; mean dV = {dV.mean()*1e3:+.2f} meV", flush=True)
    out.update(**{f"ved_{S}_x": t / 0.5, f"ved_{S}_dist_A": t * a1_len, f"ved_{S}_line": line, f"ved_{S}_site": site, f"ved_{S}_bmax": bmax, f"ved_{S}_b1": b1, f"ved_{S}_b2": b2, f"ved_{S}_ratio": bmax / abs(site), f"ved_{S}_mean": dV.mean()})
    del Vp, Vd, dV

# ---------------- 3. Re M^L, Re M^NL at K (pi pair trace/2), four N   +   4. scaling   +   hermiticity
for S in SIZES:
    dp = dense_paths(cfg, S); kd = qe_io.get_k_red(dp["uc"]); kk, eps = qe_io.get_k_eigenvalues(dp["uc"], False); eps = np.asarray(eps)
    if eps.shape[0] != len(kd): eps = eps.T
    eps = eps * HA2EV; iK = int(np.argmin(kdist(kd, tuple(cfg["K_red"])))); dK = float(kdist(kd, tuple(cfg["K_red"]))[iK]) * np.linalg.norm(B[:2, 0])  # ~A^-1
    Vw, _, _ = wannier_V(dp["wdir"], kd); pr, w = pi_pair(Vw, eps, iK); pr = list(pr)
    ML = mmap_M(dp["mfile"].replace("M_dense_", "M_L_dense_")); MN = mmap_M(dp["mfile"].replace("M_dense_", "M_NL_dense_")); Md = mmap_M(dp["mfile"])
    bL = np.array(ML[pr, iK][:, pr, iK]) * HA2EV; bN = np.array(MN[pr, iK][:, pr, iK]) * HA2EV
    reL = 0.5 * np.trace(bL).real; reN = 0.5 * np.trace(bN).real
    print(f"[LNL] {S}: nearest grid point to K at |dk|={dK:.4f} A^-1 (k={np.round(kd[iK],4)}), pair bands {pr} eps={np.round(eps[iK, pr],3)} eV: "
          f"Re M^L = {reL:+.4f} eV, Re M^NL = {reN:+.4f} eV, |M^L| {abs(np.trace(bL))/2:.4f}, |M^NL| {abs(np.trace(bN))/2:.4f}, NL/L = {reN/reL:.2f}", flush=True)
    M16 = np.array(Md[:16, :, :16, :]); mx_d = float(np.abs(M16).max()) * HA2EV
    Mc = mmap_M(f"results/M/M_ed_{S}.npy"); Mc16 = np.array(Mc[:16, :, :16, :]); mx_c = float(np.abs(Mc16).max()) * HA2EV
    herm = float(np.abs(M16 - M16.transpose(2, 3, 0, 1).conj()).max() / np.abs(M16).max())
    ncell = int(S[0]) ** 2
    print(f"[scale] {S}: max|M| bands 1-16: dense {mx_d:.4f} eV (x N_cells = {mx_d*ncell:.2f}), coarse {mx_c:.4f} eV (x N_cells = {mx_c*ncell:.2f}); hermiticity {herm:.2e}", flush=True)
    out.update(**{f"lnl_{S}_ReL": reL, f"lnl_{S}_ReNL": reN, f"lnl_{S}_dK": dK, f"lnl_{S}_pair": np.array(pr), f"scale_{S}_dense": mx_d, f"scale_{S}_coarse": mx_c, f"scale_{S}_ncells": ncell, f"herm_{S}": herm})
    tests.append((f"hermiticité M dense {S}", r"max|M - M^\dagger| / max|M|", f"{herm:.2e}", "1e-12", "OK" if herm < 1e-12 else "ÉCHEC", "analyze_M.py"))
    if S == REF:
        # coincident-k padding check (dense nb20 [:16] vs coarse M_ed, gauge-closed band subsets)
        kc = qe_io.get_k_red(f"{DATA}/unit_cell/qe/defect_{S}.save"); key = lambda k: tuple(np.round(np.mod(k + 1e-9, 1.0), 6))
        idx = {key(k): i for i, k in enumerate(kd)}; pairs = [(ic, idx[key(k)]) for ic, k in enumerate(kc) if key(k) in idx]
        for nsub in (8, 15):
            worst = 0.0
            for ic, id_ in pairs:
                for jc, jd in pairs:
                    sa = np.linalg.svd(Mc16[:nsub, ic, :nsub, jc], compute_uv=False); sb = np.linalg.svd(M16[:nsub, id_, :nsub, jd], compute_uv=False)
                    worst = max(worst, np.abs(sa - sb).max() / sa[0])
            print(f"[pad] {S}: coincident k ({len(pairs)}), bands 1..{nsub}: max rel SV mismatch dense vs coarse = {worst:.2e}", flush=True)
            tests.append((f"padding : k coïncidents, bandes 1–{nsub}, {S}", "max écart relatif des valeurs singulières (dense vs N×N)", f"{worst:.1e}", "2e-3 (tolérance nscf)", "OK" if worst < 2e-3 else "À VOIR", "analyze_M.py"))
        # nb20 vs nb16 non-regression
        A16 = mmap_M(f"results/M/M_dense_{S}_nb16.npy"); rng = np.random.default_rng(0); ks = rng.choice(len(kd), 30, replace=False); worst = 0.0
        for i in ks:
            for j in ks:
                sa = np.linalg.svd(np.array(A16[:15, i, :15, j]), compute_uv=False); sb = np.linalg.svd(M16[:15, i, :15, j], compute_uv=False); worst = max(worst, np.abs(sa - sb).max() / sa[0])
        print(f"[pad] {S}: nb20[:16] vs nb16, bands 1..15, 30x30 k pairs: {worst:.2e}", flush=True)
        tests.append((f"non-régression nbnd 16 → 20, {S}", "max écart relatif des valeurs singulières, bandes 1–15", f"{worst:.1e}", "1e-5", "OK" if worst < 1e-5 else "ÉCHEC", "analyze_M.py"))
        # plumbing: dense kernel on the coarse grid vs coarse kernel
        pc = f"results/M/M_L_dense_{S}_coarsecheck.npy"
        import os
        if os.path.exists(pc):
            Xa = np.array(mmap_M(pc)); pb = f"results/M/M_L_{S}.npy"
            Xb = np.array(mmap_M(pb)) if matrix_io.read_manifest(pb) else np.load(pb)   # legacy coarse-kernel output (Hartree, no sidecar): plumbing test only
            rel = float(np.abs(Xa - Xb).max() / np.abs(Xb).max())
            print(f"[pad] {S}: real-space dense kernel at p=1 vs coarse M^L: rel {rel:.2e}", flush=True)
            tests.append((f"padding : noyau dense à p=1 vs noyau N×N, {S}", r"max|M^L_dense - M^L| / max|M^L|", f"{rel:.1e}", "1e-10", "OK" if rel < 1e-10 else "ÉCHEC", "analyze_M.py"))
        # Wannier closure Bloch -> Wannier -> Bloch (5-band subspace) and Fourier round trip
        Hwr, Rw, nd = read_w90_HR(f"{dp['wdir']}/wannier_tb.dat"); MP = _infer_mp_grid(kd)
        kd_mp = np.round(kd * np.asarray(MP)) / np.asarray(MP)                                      # k reconstruits des indices MP (le XML de QE arrondit à ~1e-7)
        print(f"[closure] k(XML) vs k(MP): max |dk| = {np.abs(kd - kd_mp).max():.2e}", flush=True)
        Mwk = Mbk_to_Mwk(np.array(Md) * HA2EV, U9, Ud9); Mwr, R = Mwk_to_Mwr(Mwk, kd_mp, MP); Mwk2 = Mwr_to_Mwk(Mwr, R, kd_mp)
        rt = float(np.abs(Mwk2 - Mwk).max() / np.abs(Mwk).max()); Mbk5 = Mwk_to_Mbk(Mwk, Hwr, Rw, kd, ndegen=nd); worst = 0.0
        for i in ks[:12]:
            for j in ks[:12]:
                sa = np.linalg.svd(Mwk[:, i, :, j], compute_uv=False); sb = np.linalg.svd(Mbk5[:, i, :, j], compute_uv=False); worst = max(worst, np.abs(sa - sb).max() / max(sa[0], 1e-30))
        print(f"[closure] {S}: Fourier round trip Mwk->Mwr->Mwk rel {rt:.2e}; Bloch->Wannier->Bloch (smooth gauge) SV mismatch {worst:.2e}", flush=True)
        tests.append((f"fermeture Fourier Wannier (k → R → k), {S}", "max|Mwk' - Mwk| / max|Mwk| (k reconstruits des indices MP)", f"{rt:.1e}", "1e-12", "OK" if rt < 1e-12 else "ÉCHEC", "analyze_M.py"))
        tests.append((f"fermeture Bloch → Wannier → Bloch (5 bandes), {S}", "max écart relatif des valeurs singulières (Mwk vs Mbk lisse)", f"{worst:.1e}", "1e-12", "OK" if worst < 1e-12 else "ÉCHEC", "analyze_M.py"))
        del Mwk, Mwr, Mwk2, Mbk5
    del M16, Mc16
sd = np.array([out[f"scale_{S}_dense"] for S in SIZES]); spread = float((sd.max() - sd.min()) / sd.mean())
tests.append(("convention intensive (cellule unitaire)", "(max−min)/moyenne de max|M| sur N = 5,7,8,9 (bandes 1–16)", f"{spread:.1e}", "5e-2", "OK" if spread < 5e-2 else "À VOIR", "analyze_M.py"))
# recorded results of the other gates (log references)
tests.append(("test d'or (local vs compute_T), 5×5 dense", r"max|\Gamma_{loc} - \Gamma_{dense} N_c| / max", "2.3e-14", "1e-10", "OK", "golden_dense_20294198"))
tests.append(("g0 par lots vs référence, R_cut 0–3", "max écart relatif", "1.2e-14", "1e-12", "OK", "test_local_green_batch_20238555"))
with open("results/M/M_tests_summary.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["test", "quantité", "valeur", "seuil", "verdict", "source"]); w.writerows(tests)
print("\n=== TESTS DE M ==="); [print(f"  {t[0]:55s} {t[1]:60s} {t[2]:>9s}  seuil {t[3]:20s} {t[4]}") for t in tests]
out["tests"] = np.array(tests, dtype=object); np.savez("results/M/M_analysis.npz", **out); print("saved results/M/M_analysis.npz, results/M/M_tests_summary.csv")
