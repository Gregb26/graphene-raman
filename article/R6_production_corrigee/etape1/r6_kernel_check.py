#!/usr/bin/env python
"""
r6_kernel_check.py -- R6 étape 1.1 : test du noyau local corrigé (commit ee051ec, 2026-09-25). Sous-commandes :
  compare  : M^L recalculé (ml/M_L_<S>_v2.npy, noyau corrigé) contre N_cells x M_L_<S>.npy (v1) ; attendu 1e-12 relatif
  vp       : (srun MPI) potentiel périodique V_p sur la super-cellule 9x9, 4 bandes, noyau corrigé : blocs k = k' contre
             N_cells x <u|V_p|u>_maille (reconstruction KS, grille de la maille), blocs k != k' attendus nuls
  serial   : noyau série compute_ML_R (fold_wfk_to_sc) sur la 5x5, 2 bandes, contre le noyau MPI corrigé (écart attendu si bug)
  shared   : noyau partagé compute_ML_R_mpi_shared (--coarse 5x5, 2 bandes) contre le noyau MPI corrigé
  ksrec    : reconstruction KS rejouée (ksrec/results/M/ks_reconstruction.npz) contre results/M/ks_reconstruction.npz (9x9)
Sorties : kernel/*.json, journal r6_log.txt. Aucun fichier de production touché.
"""
import os
import sys
import json
import time
import argparse

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src")); sys.path.insert(0, os.path.join(PROJ, "scripts"))
os.chdir(PROJ)
from electron_defect_interaction.io import qe_io, matrix_io  # noqa: E402
from electron_defect_interaction.config import HA2EV  # noqa: E402

RES = "results/M"; ML = os.path.join(WORK, "ml"); OUT = os.path.join(WORK, "kernel"); DATA = "data/graphene"
os.makedirs(OUT, exist_ok=True)


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"; print(line, flush=True)
    with open(os.path.join(WORK, "r6_log.txt"), "a") as f:
        f.write(line + "\n")


def herm(M):
    nb, nk = M.shape[:2]; O = M.reshape(nb * nk, nb * nk); return float(np.abs(O - O.conj().T).max())


def save(name, d):
    with open(os.path.join(OUT, name), "w") as f:
        json.dump(d, f, indent=1)
    log(f"saved kernel/{name}")


def cmd_compare(a):
    out = {}
    for S in a.sizes.split(","):
        n = int(S.split("x")[0]); N = n * n
        p2 = os.path.join(ML, f"M_L_{S}_v2.npy"); v2 = np.load(p2); meta = matrix_io.read_manifest(p2) or {}
        r = dict(size=S, N_cells=N, v2=p2, shape=list(v2.shape), max_abs_v2=float(np.abs(v2).max()), herm_v2=herm(v2), M_normalization=meta.get("M_normalization"))
        p1 = os.path.join(RES, f"M_L_{S}.npy")
        if os.path.exists(p1):
            v1 = np.load(p1); ref = N * v1
            r.update(v1=p1, max_abs_v1=float(np.abs(v1).max()), rel_max_diff=float(np.abs(v2 - ref).max() / np.abs(ref).max()))
            m = np.abs(v1) > 1e-6 * np.abs(v1).max(); q = np.abs(v2[m]) / np.abs(ref[m])
            r.update(ratio_v2_over_Nv1=dict(median=float(np.median(q)), min=float(q.min()), max=float(q.max()), n=int(m.sum())))
            log(f"[compare] {S}: max|v2 - {N} v1| / max|{N} v1| = {r['rel_max_diff']:.2e} ; v2/({N} v1) médian {r['ratio_v2_over_Nv1']['median']:.6f} "
                f"(min {q.min():.4f}, max {q.max():.4f}) ; max|v2| {r['max_abs_v2']:.4e} Ha ; hermiticité {r['herm_v2']:.1e} -> {'OK' if r['rel_max_diff'] <= 1e-12 else 'ÉCART'}")
        else:
            log(f"[compare] {S}: pas de M_L v1 sur disque ; max|v2| {r['max_abs_v2']:.4e} Ha, max|v2|/N_cells {r['max_abs_v2']/N:.4e}, hermiticité {r['herm_v2']:.1e}")
        pc = os.path.join(RES, f"M_L_dense_{S}_coarsecheck.npy")
        if os.path.exists(pc) and np.load(pc, mmap_mode="r").shape == v2.shape:
            cc = np.load(pc); r["rel_max_diff_vs_N_coarsecheck"] = float(np.abs(v2 - N * cc).max() / np.abs(N * cc).max())
            log(f"[compare] {S}: contre {N} x coarsecheck (noyau partagé, sept.) : {r['rel_max_diff_vs_N_coarsecheck']:.2e}")
        out[S] = r
    save("compare.json", out)


def uc_local_term(uc, scp, pot_p, bands):
    """<u_mk|V_p|u_nk> sur la grille de la maille (norme maille, dV = Omega_uc/N_r), V_p restreint comme test_ks_reconstruction."""
    from test_ks_reconstruction import sc_pot_on_uc_grid
    from electron_defect_interaction.wavefunctions.wfk import compute_psi_nk
    V_uc, spread = sc_pot_on_uc_grid(uc, scp, pot_p)
    C, nG = qe_io.get_C_nk(uc); G = qe_io.get_G_red(uc); k = qe_io.get_k_red(uc); _, Om = qe_io.get_A_volume(uc); ng = qe_io.get_ngfft(uc)
    C = C[bands]; nb = len(bands); nk = len(k)
    psi, _ = compute_psi_nk(C, nG, G, k, Om, ngfft=ng)
    Nr = int(np.prod(ng)); dV = Om / Nr; Vr = np.asarray(V_uc).reshape(Nr); L = np.zeros((nk, nb, nb), complex)
    for ik in range(nk):
        P = psi[:, ik].reshape(nb, Nr); L[ik] = dV * (P.conj() * Vr) @ P.T
    return L, float(spread), int(nk)


def cmd_vp2(a):
    """Même chose que vp mais tout en un (rang 0 fait la référence après le noyau)."""
    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import prep_realspace_inputs, compute_ML_R_mpi
    comm = MPI.COMM_WORLD; rank = comm.Get_rank()
    S = "9x9"; uc = f"{DATA}/unit_cell/qe/defect_{S}.save"; scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; pot_p = f"{scp}/Vks_{S}_p"
    bands = [0, 1, 2, 3]; N = 81
    prep = None
    if rank == 0:
        prep = prep_realspace_inputs(uc, scp, pot_p, pot_p, subtract_mean=False, pristine=True, bands=bands, io=qe_io)
        log(f"[vp] prep : grille {tuple(int(x) for x in prep['ngfft'])}, Ndiag {tuple(int(x) for x in prep['Ndiag'])}, Omega_sc {prep['Omega_sc']:.3f} bohr^3, pristine=True, bandes {bands}")
    prep = comm.bcast(prep, root=0)
    t0 = time.time(); M_p = compute_ML_R_mpi(prep, grid_block=50000)
    if rank != 0:
        return
    log(f"[vp] noyau corrigé sur la super-cellule : M^L(V_p) {M_p.shape} en {time.time()-t0:.0f} s, max|.| {np.abs(M_p).max():.6e} Ha")
    t0 = time.time(); L_uc, spread, nk = uc_local_term(uc, scp, pot_p, bands)
    log(f"[vp] référence maille : <u|V_p|u> ({nk} k, {len(bands)} bandes) en {time.time()-t0:.0f} s ; périodicité de V_p (écart entre images) {spread*HA2EV*1e3:.3f} meV")
    diag = np.array([M_p[:, k, :, k] for k in range(nk)])                    # (nk, nb, nb)
    ref = N * L_uc
    d_diag = float(np.abs(diag - ref).max()); rel_diag = d_diag / float(np.abs(ref).max())
    off = M_p.copy()
    for k in range(nk):
        off[:, k, :, k] = 0.0
    off_max = float(np.abs(off).max()); off_rel = off_max / float(np.abs(M_p).max())
    d_uc = float(np.abs(diag - L_uc).max()) / float(np.abs(L_uc).max())
    r = dict(size=S, bands=bands, N_cells=N, max_abs_Mp=float(np.abs(M_p).max()), max_abs_Luc=float(np.abs(L_uc).max()),
             rel_max_diff_diag_vs_Ncells_Luc=rel_diag, abs_max_diff_diag_vs_Ncells_Luc_Ha=d_diag, rel_max_diff_diag_vs_Luc=d_uc,
             offdiag_k_max_Ha=off_max, offdiag_k_rel=off_rel, spread_Vp_meV=spread * HA2EV * 1e3,
             mean_diag_Mp_eV=float(np.mean([M_p[n, k, n, k].real for n in range(len(bands)) for k in range(nk)]) * HA2EV),
             mean_diag_Luc_eV=float(np.mean([L_uc[k, n, n].real for n in range(len(bands)) for k in range(nk)]) * HA2EV))
    log(f"[vp] blocs k = k' : max|M^L(V_p) - {N} <u|V_p|u>_maille| / max = {rel_diag:.2e} (contre <u|V_p|u> sans facteur : {d_uc:.2e}) ; "
        f"blocs k != k' : max {off_max:.2e} Ha = {off_rel:.2e} x max|M| ; moyenne diagonale M^L(V_p) {r['mean_diag_Mp_eV']:+.4f} eV, maille {r['mean_diag_Luc_eV']:+.4f} eV")
    save("vp.json", r)
    np.save(os.path.join(OUT, "Mp_9x9_4b.npy"), M_p)


def cmd_serial(a):
    from electron_defect_interaction.defects.local_R import compute_ML_R
    S = "5x5"; uc = f"{DATA}/unit_cell/qe/defect_{S}.save"; scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; scd = f"{DATA}/supercell/qe/defect_{S}_d.save"
    t0 = time.time(); Ms = compute_ML_R(uc, scp, f"{scp}/Vks_{S}_p", f"{scd}/Vks_{S}_d", subtract_mean=False, bands=[0, 1], io=qe_io)
    v2 = np.load(os.path.join(ML, f"M_L_{S}_v2.npy"))[:2, :, :2, :]
    m = np.abs(v2) > 1e-6 * np.abs(v2).max(); q = np.abs(Ms[m]) / np.abs(v2[m])
    ng = tuple(int(x) for x in qe_io.get_ngfft(scp)); N_uc = int(np.prod(ng)) // 25
    r = dict(size=S, bands=[0, 1], max_abs_serial=float(np.abs(Ms).max()), max_abs_mpi_v2=float(np.abs(v2).max()), ratio_serial_over_mpi=dict(median=float(np.median(q)), min=float(q.min()), max=float(q.max())),
             N_grid=int(np.prod(ng)), N_uc_grid=N_uc, N_cells=25, rel_max_diff=float(np.abs(Ms - v2).max() / np.abs(v2).max()), seconds=time.time() - t0)
    log(f"[serial] compute_ML_R (fold_wfk_to_sc, ee051ec) 5x5 bandes 0-1 : max|M| {r['max_abs_serial']:.4e} Ha ; noyau MPI corrigé {r['max_abs_mpi_v2']:.4e} ; "
        f"rapport série/MPI médian {r['ratio_serial_over_mpi']['median']:.6g} (min {q.min():.6g}, max {q.max():.6g}) ; N_uc(grille) = {N_uc}, N_cells = 25 ; écart relatif {r['rel_max_diff']:.2e}")
    save("serial.json", r)


def cmd_shared(a):
    p = os.path.join(ML, "M_L_5x5_shared_coarsecheck_v2.npy")
    if not os.path.exists(p):
        log("[shared] fichier absent : compute_ML_R_mpi_shared --coarse 5x5 n'a pas produit de sortie (voir le .err du job)"); save("shared.json", dict(produced=False)); return
    Msh = np.load(p); v2 = np.load(os.path.join(ML, "M_L_5x5_v2.npy"))[:2, :, :2, :]
    rel = float(np.abs(Msh - v2).max() / np.abs(v2).max())
    log(f"[shared] compute_ML_R_mpi_shared --coarse 5x5 bandes 0-1 contre noyau MPI corrigé : écart relatif {rel:.2e}")
    save("shared.json", dict(produced=True, rel_max_diff=rel, max_abs=float(np.abs(Msh).max())))


def cmd_ksrec(a):
    new = os.path.join(WORK, "ksrec", "results", "M", "ks_reconstruction.npz"); old = os.path.join(RES, "ks_reconstruction.npz")
    zn = np.load(new, allow_pickle=True); zo = np.load(old, allow_pickle=True); r = {}
    for key in sorted(k for k in zn.files if k.startswith("9x9_")):
        if key in zo.files:
            d = float(np.abs(np.asarray(zn[key], float) - np.asarray(zo[key], float)).max()); r[key] = dict(new=np.asarray(zn[key]).tolist(), old=np.asarray(zo[key]).tolist(), max_abs_diff=d)
    worst = max(v["max_abs_diff"] for v in r.values())
    log(f"[ksrec] reconstruction KS 9x9 (grossier + dense) rejouée après ee051ec contre results/M/ks_reconstruction.npz : {len(r)} clés, écart max {worst:.2e} ; "
        f"coarse max {float(zn['9x9_coarse_max_meV']):.4f} meV (ancien {float(zo['9x9_coarse_max_meV']):.4f}), dense max {float(zn['9x9_dense_max_meV']):.4f} meV (ancien {float(zo['9x9_dense_max_meV']):.4f})")
    save("ksrec.json", dict(keys=r, max_abs_diff=worst))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["compare", "vp", "serial", "shared", "ksrec"]); ap.add_argument("--sizes", default="9x9,7x7,5x5,6x6,8x8")
    a = ap.parse_args()
    {"compare": cmd_compare, "vp": cmd_vp2, "serial": cmd_serial, "shared": cmd_shared, "ksrec": cmd_ksrec}[a.cmd](a)


if __name__ == "__main__":
    main()
