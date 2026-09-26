#!/usr/bin/env python
"""
gate_M_normalization.py -- R6 gate A.2 (mandatory on every M used in production; 2026-09-25).

For pure folded unit-cell Bloch states |nk>_sc (coefficients c = delta) of the N x N supercell, the matrix elements
<nk|dV^L|n'k'>_sc and <nk|dV^NL|n'k'>_sc obtained by DIRECT application of dV = V_d - V_p (defects/deltav_pw: local part on
the supercell FFT grid, non-local part through the Kleinman-Bylander projectors of the removed atom, supercell plane waves)
must equal M^L[nk,n'k']/N_cells and M^NL[nk,n'k']/N_cells, L and NL separately, to <= --tol (1e-6 eV). A file that fails is
REFUSED (exit 3) and must not enter production. Pre-R6 (v1) files fail with ratio direct/(M^L/N_cells) = N_cells on the L part.

Dense M (p*N x p*N primitive grid): only the N^2 k coincident with the coarse grid are pure Gamma states of the N x N supercell;
the pairs are drawn among them and the k / coefficients come from the dense .save itself (R5 lesson 0.3). Grid: the one the
production kernel used (next multiple of Ndiag; dV Fourier-resampled exactly like prep_realspace_inputs, e.g. 7x7: 216 -> 217).

Usage:
  gate_M_normalization.py --size 9x9 --level coarse|dense [--m2-dir results/M2] [--ml PATH --nl PATH] [--uc PATH]
                          [--tol 1e-6] [--workers 8] [--out gate.json] [--allow-v1] [--n-pairs 8]
  gate_M_normalization.py --summary DIR        # markdown table N x (L, NL) from the gate_*.json files of DIR
"""
import os
import sys
import json
import glob
import time
import argparse

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); PROJ = os.path.dirname(HERE)
try:
    import electron_defect_interaction  # noqa: F401
except ImportError:
    sys.path.insert(0, os.path.join(PROJ, "src"))
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io import qe_gamma_io as qg
from electron_defect_interaction.defects import alignment as al
from electron_defect_interaction.defects.local_R import fourier_resample
from electron_defect_interaction.wannier import supercell_fold as sf
from electron_defect_interaction.wavefunctions import sc_projection as sp
from electron_defect_interaction.defects import deltav_pw as dv
from electron_defect_interaction.utils.lattice import red_to_cart
from electron_defect_interaction.config import load_production, dense_paths, HA2EV

BOHR = 0.529177210903
K_RED = np.array([2 / 3, 1 / 3, 0.0]); KP_RED = np.array([1 / 3, 2 / 3, 0.0])


def kindex(k_red, target):
    return int(np.argmin(np.linalg.norm(np.mod(np.asarray(k_red) - np.asarray(target) + 0.5, 1) - 0.5, axis=1)))


def read_wfc_subset(save_dir, ks_index):
    """C_nk(G) of the files wfc<ik>.hdf5 (ik = index + 1, XML order) for a subset of k indices."""
    Cs, Gs, nG = [], [], []
    for j in ks_index:
        C, mill, go, _ = qg.read_wfc_file(os.path.join(save_dir, f"wfc{j + 1}.hdf5"))
        assert not go, "gamma_only wfc in a k-point .save?"
        Cs.append(C); Gs.append(mill); nG.append(C.shape[1])
    nGm = max(nG); nb = Cs[0].shape[0]
    C_nkg = np.zeros((nb, len(ks_index), nGm), complex); G_red = np.zeros((len(ks_index), nGm, 3), int)
    for i, (C, G) in enumerate(zip(Cs, Gs)):
        C_nkg[:, i, :C.shape[1]] = C; G_red[i, :C.shape[1]] = G
    return C_nkg, np.array(nG), G_red


def coincident_indices(k_all, N, tol=1e-6):
    x = np.asarray(k_all, float) * N
    return np.where(np.all(np.abs(x - np.rint(x)) < tol, axis=1))[0]


def check_sidecar(path, allow_v1):
    meta = matrix_io.read_manifest(path)
    if meta is None:
        if allow_v1:
            return {}
        raise ValueError(f"{path}: no sidecar")
    if meta.get("bloch_norm") != matrix_io.UNIT_CELL or meta.get("units") != matrix_io.HARTREE:
        raise ValueError(f"{path}: sidecar {meta.get('bloch_norm')}/{meta.get('units')}, expected unit_cell/hartree")
    if not allow_v1 and not str(meta.get("M_normalization", "")).startswith("v2"):
        raise ValueError(f"{path}: M_normalization={meta.get('M_normalization')!r} (v1 file): refused; the gate documents v1 only with --allow-v1")
    return meta


def summary(d):
    rows = []
    for f in sorted(glob.glob(os.path.join(d, "gate_*.json"))):
        r = json.load(open(f))
        rows.append(r)
    rows.sort(key=lambda r: (int(r["size"].split("x")[0]), r["level"]))
    print("| N | niveau | fichier M^L | k testés | paires | max\\|Δ_L\\| (eV) | max\\|Δ_NL\\| (eV) | direct/(M^L/N_cells) médian | direct/(M^NL/N_cells) médian | seuil | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['size']} | {r['level']} | `{os.path.basename(r['ml'])}` | {r['nk_tested']} | {r['n_pairs']} | {r['max_dL']:.2e} | {r['max_dNL']:.2e} | "
              f"{r['ratio_L_median']:.6f} | {r['ratio_NL_median']:.6f} | {r['tol']:.0e} | {'OK' if r['ok'] else 'REFUSÉ'} |")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--size"); ap.add_argument("--level", choices=["coarse", "dense"])
    ap.add_argument("--m2-dir", default="results/M2"); ap.add_argument("--ml"); ap.add_argument("--nl"); ap.add_argument("--uc")
    ap.add_argument("--tol", type=float, default=1e-6); ap.add_argument("--workers", type=int, default=int(os.environ.get("OMP_NUM_THREADS", "8")))
    ap.add_argument("--out"); ap.add_argument("--allow-v1", action="store_true"); ap.add_argument("--n-pairs", type=int, default=8)
    ap.add_argument("--summary", metavar="DIR")
    a = ap.parse_args()
    if a.summary:
        summary(a.summary); return
    if not (a.size and a.level):
        ap.error("--size and --level are required")
    t0 = time.time(); os.chdir(PROJ)
    S = a.size; N = int(S.split("x")[0]); N_cells = N * N; dense = a.level == "dense"
    DATA = "data/graphene"
    scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; scd = f"{DATA}/supercell/qe/defect_{S}_d.save"
    pot_p = f"{scp}/Vks_{S}_p"; pot_d = f"{scd}/Vks_{S}_d"; upf = f"{DATA}/unit_cell/qe/defect_{S}.save/C.upf"
    if dense:
        cfg = load_production(verbose=False); uc = a.uc or dense_paths(cfg, S)["uc"]
        ml = a.ml or os.path.join(a.m2_dir, f"M_L_dense_{S}.npy"); nl = a.nl or os.path.join(a.m2_dir, f"M_NL_dense_{S}.npy")
    else:
        uc = a.uc or f"{DATA}/unit_cell/qe/defect_{S}.save"
        ml = a.ml or os.path.join(a.m2_dir, f"M_L_{S}.npy"); nl = a.nl or os.path.join(a.m2_dir, f"M_NL_{S}.npy")
    meta_l = check_sidecar(ml, a.allow_v1); meta_nl = check_sidecar(nl, a.allow_v1)
    print(f"[gate] {S} {a.level}: M^L {ml} ({meta_l.get('M_normalization', 'v1?')}), M^NL {nl} ({meta_nl.get('M_normalization', 'v1?')}), uc {uc}", flush=True)

    # ---- supercell geometry, grid of the production kernel, dV
    ng = tuple(int(x) for x in qe_io.get_ngfft(scd)); assert ng == tuple(int(x) for x in qe_io.get_ngfft(scp)), "p/d FFT grids differ"
    Ndiag = np.array([N, N, 1]); ng_k = tuple(int(x) for x in ((np.array(ng) + Ndiag - 1) // Ndiag) * Ndiag)
    A_b, Om = qe_io.get_A_volume(scp); B_sc, _ = qe_io.get_B_volume(scp); ecut = float(qe_io.get_ecut(scp))
    Vd, ngv = qe_io.get_pot(pot_d, subtract_mean=False, to_hartree=True); Vp, _ = qe_io.get_pot(pot_p, subtract_mean=False, to_hartree=True)
    assert tuple(int(x) for x in ngv) == ng, (ngv, ng)
    dV = (Vd - Vp).transpose(2, 1, 0); del Vd, Vp
    resampled = ng_k != ng
    if resampled:
        dV = fourier_resample(dV, ng_k)
    dV = dV * HA2EV
    print(f"[gate] grid {ng} -> kernel grid {ng_k}{' (V_ed Fourier-resampled, as prep_realspace_inputs)' if resampled else ''}, Omega_sc {Om:.3f} bohr^3, <dV>_3D {dV.mean()*1e3:+.2f} meV", flush=True)

    # ---- pure Bloch states: coefficients and k of the .save providing them
    k_all = qe_io.get_k_red(uc)
    if dense:
        idx = coincident_indices(k_all, N)
        if len(idx) != N_cells:
            raise SystemExit(f"{uc}: {len(idx)} coincident k, expected {N_cells}")
        C, nG, G = read_wfc_subset(uc, idx); k = k_all[idx]
    else:
        C, nG = qe_io.get_C_nk(uc); G = qe_io.get_G_red(uc); k = k_all; idx = np.arange(len(k))
        if len(k) != N_cells:
            raise SystemExit(f"{uc}: {len(k)} k, expected {N_cells}")
    nb = C.shape[0]
    flat = sf.sc_planewave_index(k, G, nG, Ndiag, ng_k)
    Cw, mill_half, go, _ = qg.read_wfc_gamma(scd, bands=(0, 1)); del Cw
    _, mill_full = sp.sc_full_sphere(np.zeros((1, len(mill_half)), complex), mill_half)
    nu, ntot, ok_union = sp.check_planewave_union(flat, len(mill_full))
    print(f"[gate] {nb} bands x {len(k)} k; folded plane waves distinct {nu}/{ntot}, supercell sphere {len(mill_full)} -> {'OK' if ok_union else 'MISMATCH'}", flush=True)
    flat_full = sp.grid_flat_index(mill_full, ng_k)

    # ---- removed-atom KB projector
    x_p = qe_io.get_x_red(scp); x_d = qe_io.get_x_red(scd)
    s_vac, i_vac_p, _ = al.vacancy_site(x_p, x_d, A_b * BOHR); tau = red_to_cart(x_p[i_vac_p], A_b)
    t1 = time.time(); proj = dv.RemovedAtomProjector(mill_full, B_sc, tau, upf, Om, ecut, energy_scale=HA2EV)
    print(f"[gate] removed atom: pristine index {i_vac_p+1}, s_red {np.round(s_vac, 5).tolist()}; projector on {len(mill_full)} g (|g|max {proj.K_norm_max:.2f} < {proj.qmax:.1f}) in {time.time()-t1:.0f} s", flush=True)

    # ---- M parts restricted to the tested k (eV)
    ML = np.load(ml, mmap_mode="r"); MNL = np.load(nl, mmap_mode="r")
    assert ML.shape[0] == nb and MNL.shape == ML.shape, (ML.shape, MNL.shape, nb)
    bi = np.arange(nb)
    ML = np.array(ML[np.ix_(bi, idx, bi, idx)]) * HA2EV; MNL = np.array(MNL[np.ix_(bi, idx, bi, idx)]) * HA2EV

    # ---- pairs: intra-k at K (or nearest), inter-k, Gamma, top band, mixed
    nk = len(k); iK = kindex(k, K_RED); iKp = kindex(k, KP_RED); iG = kindex(k, (0, 0, 0))
    pairs = [(3, iK, 3, iK), (3, iK, 4, iK), (3, iK, 3, iKp), (0, iG, 0, iG), (2, 5 % nk, 7, 40 % nk), (10, 17 % nk, 1, 60 % nk),
             (nb - 1, 1 % nk, nb - 1, 2 % nk), (5, 3 % nk, 12, 11 % nk), (nb - 1, iK, nb - 3, iKp), (4, iG, 4, iK)]
    pairs = [p for p in pairs if p[0] < nb and p[2] < nb][:a.n_pairs]
    t1 = time.time()
    res = dv.check_pure_bloch(pairs, C, nG, flat, ng_k, Om, dV, ML, MNL, N_cells, proj, flat_full, workers=a.workers)
    dL = max(r["dL"] for r in res); dNL = max(r["dNL"] for r in res); ok = bool(max(dL, dNL) <= a.tol)
    thr = 1e-4
    ratio_L = [abs(r["VL_direct"]) / abs(r["ML_over_N"]) for r in res if abs(r["ML_over_N"]) > thr]
    ratio_NL = [abs(r["VNL_direct"]) / abs(r["MNL_over_N"]) for r in res if abs(r["MNL_over_N"]) > thr]
    for r in res:
        print(f"[gate]   {r['pair']}: L direct {r['VL_direct'].real:+.6f} vs M^L/N {r['ML_over_N'].real:+.6f} (Δ {r['dL']:.1e}) ; "
              f"NL direct {r['VNL_direct'].real:+.6f} vs M^NL/N {r['MNL_over_N'].real:+.6f} (Δ {r['dNL']:.1e}) ; norm {r['norm1']:.12f}", flush=True)
    out = dict(size=S, level=a.level, ml=os.path.abspath(ml), nl=os.path.abspath(nl), uc=os.path.abspath(uc), N_cells=N_cells, nb=nb, nk_tested=nk,
               grid=list(ng), kernel_grid=list(ng_k), resampled=resampled, pw_union=[int(nu), int(ntot), bool(ok_union)], sphere=int(len(mill_full)),
               removed_atom_index_1based=int(i_vac_p + 1), n_pairs=len(pairs), tol=a.tol, max_dL=float(dL), max_dNL=float(dNL), ok=ok,
               ratio_L_median=float(np.median(ratio_L)) if ratio_L else float("nan"), ratio_NL_median=float(np.median(ratio_NL)) if ratio_NL else float("nan"),
               ratio_L=[float(x) for x in ratio_L], ratio_NL=[float(x) for x in ratio_NL],
               pairs=[dict(pair=list(r["pair"]), VL_direct=[r["VL_direct"].real, r["VL_direct"].imag], ML_over_N=[r["ML_over_N"].real, r["ML_over_N"].imag], dL=r["dL"],
                           VNL_direct=[r["VNL_direct"].real, r["VNL_direct"].imag], MNL_over_N=[r["MNL_over_N"].real, r["MNL_over_N"].imag], dNL=r["dNL"]) for r in res],
               M_normalization=dict(L=meta_l.get("M_normalization"), NL=meta_nl.get("M_normalization")), seconds=round(time.time() - t0, 1), seconds_pairs=round(time.time() - t1, 1))
    print(f"[gate] {S} {a.level}: max|Δ_L| = {dL:.2e} eV, max|Δ_NL| = {dNL:.2e} eV on {len(pairs)} pairs (tol {a.tol:.0e}); "
          f"direct/(M^L/N_cells) median {out['ratio_L_median']:.6f}, direct/(M^NL/N_cells) median {out['ratio_NL_median']:.6f} -> {'OK' if ok else 'REFUSÉ'} ({time.time()-t0:.0f} s)", flush=True)
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        with open(a.out, "w") as f:
            json.dump(out, f, indent=1)
    sys.exit(0 if ok else 3)


if __name__ == "__main__":
    main()
