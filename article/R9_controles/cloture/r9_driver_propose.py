#!/usr/bin/env python
"""
r9_driver.py -- pilote unique de la campagne R9 (contrôles avant le chapitre 4). Lancé par submit_r9.sh (une tâche par job).

Sous-commandes :
  a2c      A.2, grossier : M^L d'un potentiel constant (ΔV = 1 sur toute la super-cellule 9x9, bandes 0-3, noyau compute_ML_R_mpi)
           contre N_cells·𝕀 (seuil 1e-12 relatif). MPI (srun).
  a2d      A.2, dense : M^L[1_boîte] = M^L de ΔV = 1 sur la boîte N×N (zéro ailleurs dans la cellule zéro-paddée p·N), par le noyau
           partagé de production compute_ML_R_mpi_shared tel quel ; ΔV = 1 est injecté par un module io substitut (ConstPotIO), le
           noyau n'est pas modifié (P5 refusée). Sortie cache/ML_box_<S>.npy, sans dimension : M_aligné = M2 − C_N·M^L[1_boîte]. MPI (srun).
  a2dpost  A.2, dense : contrôles de M^L[1_boîte] (hermiticité, diagonale = N_cells, k coïncidents) puis rotation en base de Wannier
           F_W(R, R') = <w_R|1_boîte|w_R'> (V† M V, double TF), étiquettes recentrées sur le R_d de M2 ; F_W contre 𝕀₅ par site
           (boîte, bord, hors boîte, amas R_cut 3 et 4) ; règle de somme Σ_R F_ww(R, R) = N_cells. Série.
  a0       A.0 : porte (niveau 1, 9x9, M2 tel quel, R_cut 3, 240², η 0,02, N_k^int 300 -> 3 132,60 meV, −0,175 eV).
  a1       A.1 : C_N par taille, (i) sphères autour de tous les atomes (P1), plateau à d >= 0,75 r_max ; (ii) diagonales de M_W
           (périodique : M2 dense restreint aux N² k coïncidents ; boîte : M_W de production) ; (iii) <ΔV>_3D ; C_N retenu et étiquette.
  a3       A.3 : C14 de R6 (lecture), convention intensive, tab:rcut_M, niveau 1 et courbe Γ_T (six tailles, R_cut 3 et 4), résonance 9x9,
           variantes brut / aligné (i) / aligné exact (F_W, 9x9 et 5x5).
  a3pole   A.3 : critère de pôle 9x9 (min|det|, min|λ|), blocs π et complet ; 1 fil BLAS + fils Python.
  b        B : position de la résonance contre N_k^int (9x9, 12x12 ; 300, 450, 600, 900), variantes « brut » (M2 tel quel), « aligne »
           et « exact » (C_N de A.1) ; porte B.0 ; tableaux, figure resonance_vs_nkint. btables : tables et figure depuis le json.
  c        C : chaîne repliée (H(R) 27x27 + V_loc 9x9 de production, R_cut 3), porte C.0, C.1 (N = 6…27, complet), C.2 (bloc π, N = 27…81),
           ajustements, C.3 (LDOS de la limite diluée, P3) ; figure folded_vs_R7.
  d        D : grandeur de Kaasbjerg — porte D.0 (107,22), D.2 (bloc 2x2 à (K, K)), D.1 (carte k = K + δx̂, P2) ; figure kaasbjerg_fig3_map.
Clôture (2026-09-28) : option --cn lu | plateau (a3, a3pole, b) — C_N = Lu (1,0 Å) ou moyenne du plateau (i) d'A.1 ; sorties suffixées _lu / _plateau,
  variantes « aligne_<cn> », « exact_<cn> » ; sans --cn, comportement de R9 inchangé.
  r0       porte R.0 : --cn lu redonne R9 (9x9, 12x12, R_cut 3 : médiane de niveau 1 alignée (i), tab:rcut_M).
  synth    synthèse (tel quel, Lu, plateau), figures resonance_vs_nkint (trois variantes ; version R9 gardée sous _R9) et rcut_aligned.
Données de production en lecture seule (results/M2, caches Wannier, .save, Vks, config/production.json) ; sorties dans ce répertoire
(a/, a2/, b/, c/, d/, cache/, fig/) ; journal r9_log.txt. Paramètres de production (config v2) sauf N_k^int en B.
"""
import os
import sys
import json
import time
import argparse
import subprocess

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                   # .../graphene/qe
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src"))
sys.path.insert(0, os.path.join(PROJ, "scripts"))                              # compute_M_dense_stages.paths (mêmes entrées que la production)
os.chdir(PROJ)                                                                 # chemins data/... relatifs au dépôt

from electron_defect_interaction.config import load_production, dense_paths, results_dir, HA2EV  # noqa: E402
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance  # noqa: E402
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order  # noqa: E402
from electron_defect_interaction.defects.many_body import local_tmatrix as lt  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mwr_to_Mwk_pairs  # noqa: E402  (P2, R9)
from electron_defect_interaction.wannier import supercell_fold as sf  # noqa: E402  (repliement ; P4, R9)
from electron_defect_interaction.defects import alignment as al  # noqa: E402  (far_atom_alignment ; P1, R9)

BOHR = 0.529177210903
NW = 5; WF_SIGMA = [0, 1, 2]; WF_PI = [3, 4]
K_RED = np.array([2 / 3, 1 / 3, 0.0]); KP_RED = np.array([1 / 3, 2 / 3, 0.0])
RANK = int(os.environ.get("PMI_RANK", os.environ.get("SLURM_PROCID", "0")))


# ----------------------------------------------------------------------------------------------- utilitaires
def log(msg):
    """Journal horodaté (rang 0 seulement) : stdout + r9_log.txt."""
    if RANK != 0:
        return
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "r9_log.txt"), "a") as f:
        f.write(line + "\n")


def ensure(sub):
    d = os.path.join(WORK, sub)
    os.makedirs(d, exist_ok=True)
    return d


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, complex):
        return [x.real, x.imag]
    return x


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(jsonable(d), f, indent=1, ensure_ascii=False)


def git_head():
    try:
        return subprocess.run(["git", "-C", PROJ, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:
        return "?"


def size_n(S):
    return int(S.split("x")[0])


class ConstPotIO:
    """
    Module io substitut (R9 ; P5 refusée, noyau de production inchangé). Délègue tout à qe_io, sauf get_pot : le vrai fichier pp.x
    est lu (grille et disposition identiques à la production), puis remplacé par des 1 pour le potentiel « défectueux » et des 0
    pour le « parfait ». Le noyau forme ΔV = V_d − V_p = 1 sur toute la grille de la super-cellule N×N, c'est-à-dire 1_boîte dans la
    cellule zéro-paddée du noyau dense (et 1 partout pour le noyau grossier, k commensurables).
    """

    def __init__(self, pot_d):
        self.pot_d = os.path.realpath(pot_d)

    def __getattr__(self, name):
        return getattr(qe_io, name)

    def get_pot(self, filepath, subtract_mean=True, to_hartree=True):
        V, ng = qe_io.get_pot(filepath, subtract_mean=False, to_hartree=to_hartree)
        one = os.path.realpath(filepath) == self.pot_d
        return (np.ones_like(V) if one else np.zeros_like(V)), ng


def coarse_paths(S):
    """Entrées grossières de production (scripts/submit_r6_kernel.sh, compute_M.py --stage ml)."""
    return dict(uc=f"data/graphene/unit_cell/qe/defect_{S}.save", scp=f"data/graphene/supercell/qe/defect_{S}_p.save",
                pot_p=f"data/graphene/supercell/qe/defect_{S}_p.save/Vks_{S}_p",
                pot_d=f"data/graphene/supercell/qe/defect_{S}_d.save/Vks_{S}_d")


def load_wannier(dp, k):
    """U, U_dis réordonnés sur les k du .save dense, H(R) ; porte de jauge (refus si le manifeste ne passe pas)."""
    paths = wannier_provenance.load_wannier_checked(dp["manifest"])
    U, kU = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(kU, k)]
    Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_kpoint_order(kUd, k)]
    Hwr, Rw, nd = read_w90_HR(paths["tb"])
    return U, Ud, Hwr, Rw, nd


def box_geometry(R, MP, n):
    """
    Mailles de la boîte N×N dans les étiquettes R de Mwk_to_Mwr (boîte MP D×D, périodique) : cellule c = R mod D ; dans la boîte si
    0 <= c_a < n ; distance au bord (en mailles) = min_a min(c_a, n − 1 − c_a) (0 = dernière maille), négative hors de la boîte
    (−1 = première maille à l'extérieur).
    """
    D = np.asarray(MP[:2], int)
    c = np.mod(np.asarray(R, int)[:, :2], D)
    inb = np.all(c < n, axis=1)
    edge = np.min(np.minimum(c, n - 1 - c), axis=1)
    out = np.where(inb, edge, 0)
    # hors boîte : distance à la boîte la plus courte, côté périodique compris
    d_lo = np.where(c >= n, np.minimum(c - (n - 1), D - c), 0)
    out = np.where(inb, out, -np.max(d_lo, axis=1))
    return inb, out


# ----------------------------------------------------------------------------------------------- A.2 grossier
def cmd_a2c(a):
    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import prep_realspace_inputs, compute_ML_R_mpi
    comm = MPI.COMM_WORLD; rank = comm.Get_rank()
    S = "9x9"; n = size_n(S); N_cells = n * n; P = coarse_paths(S); bands = [0, 1, 2, 3]
    t0 = time.time(); prep = None
    if rank == 0:
        log(f"[A.2c] {S} grossier, bandes {bands}, ΔV = 1 (ConstPotIO) ; {comm.Get_size()} rangs")
        prep = prep_realspace_inputs(P["uc"], P["scp"], P["pot_p"], P["pot_d"], subtract_mean=False, bands=bands, io=ConstPotIO(P["pot_d"]))
        log(f"[A.2c] ΔV : min {prep['Ved'].min():.3e} max {prep['Ved'].max():.3e} ; grille {prep['ngfft']} ; Ndiag {prep['Ndiag'].tolist()}")
    prep = comm.bcast(prep, root=0)
    M = compute_ML_R_mpi(prep)
    if rank == 0:
        nb, nk = M.shape[:2]; Op = M.reshape(nb * nk, nb * nk)
        dev = float(np.abs(Op - N_cells * np.eye(nb * nk)).max() / N_cells)
        herm = float(np.abs(Op - Op.conj().T).max() / np.abs(Op).max())
        diag = np.diag(Op).real
        d = ensure("a2"); np.save(os.path.join(d, "ML_const_9x9_coarse_4b.npy"), M)
        out = dict(size=S, bands=bands, N_cells=N_cells, shape=list(M.shape), max_rel_dev_from_Ncells_identity=dev, threshold=1e-12,
                   ok=bool(dev < 1e-12), hermiticity_rel=herm, diag_min=float(diag.min()), diag_max=float(diag.max()),
                   max_offdiag_abs=float(np.abs(Op - np.diag(np.diag(Op))).max()), elapsed_s=time.time() - t0, head=git_head())
        save_json(os.path.join(d, "a2c_9x9.json"), out)
        log(f"[A.2c] max|M^L[1] − {N_cells}·𝕀|/{N_cells} = {dev:.2e} (seuil 1e-12 -> {'OK' if dev < 1e-12 else 'ÉCHEC'}) ; diagonale "
            f"[{diag.min():.12f}, {diag.max():.12f}] ; hors diagonale max {out['max_offdiag_abs']:.2e} ; hermiticité {herm:.1e} ; {time.time()-t0:.0f} s")
    comm.Barrier()


# ----------------------------------------------------------------------------------------------- A.2 dense
def cmd_a2d(a):
    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import compute_ML_R_mpi_shared
    import compute_M_dense_stages as cmds
    comm = MPI.COMM_WORLD; rank = comm.Get_rank()
    S = a.size; P = cmds.paths(S); n = size_n(S)
    d = ensure("cache"); out = os.path.join(d, f"ML_box_{S}.npy")
    if os.path.exists(out):
        log(f"[A.2d] {out} existe déjà : sauté"); return
    if rank == 0:
        log(f"[A.2d] {S} dense p={P['p']} D={P['D']} : M^L[1_boîte] par compute_ML_R_mpi_shared (noyau inchangé, HEAD {git_head()}), "
            f"ΔV = 1 (ConstPotIO) ; uc {P['uc_dense']} ; {comm.Get_size()} rangs ; bloc {a.block_size}")
    t0 = time.time()
    M = compute_ML_R_mpi_shared(P["uc_dense"], P["scp"], P["pot_p"], P["pot_d"], subtract_mean=False, bands=None,
                                io=ConstPotIO(P["pot_d"]), grid_block=a.block_size)
    if rank == 0:
        np.save(out, M)
        side = dict(quantity="M^L[1_boite] : noyau compute_ML_R_mpi_shared avec ΔV = 1 sur la grille de la super-cellule NxN "
                             "(zéro hors de la boîte dans la cellule zéro-paddée p·N)",
                    units="sans dimension (ΔV = 1) : M_aligné = M2 − C_N · M^L[1_boîte], C_N dans l'unité de M2",
                    bloch_norm="unit_cell", index="[bra_band, k', ket_band, k]", size=S, p=P["p"], D=P["D"], N_cells=n * n,
                    uc_dense=P["uc_dense"], sc_p=P["scp"], pot_grid_from=P["pot_d"], kernel="local_R.compute_ML_R_mpi_shared (noyau v2, 09280bb)",
                    head=git_head(), date=time.strftime("%Y-%m-%d %H:%M"), elapsed_s=time.time() - t0, shape=list(M.shape), campaign="R9 A.2")
        save_json(out.replace(".npy", ".json"), side)
        log(f"[A.2d] {S} écrit {out} {M.shape} en {(time.time()-t0)/60:.1f} min")
    comm.Barrier()


def cmd_a2dpost(a):
    S = a.size; n = size_n(S); N_cells = n * n; cfg = load_production(verbose=False); dp = dense_paths(cfg, S)
    f = os.path.join(WORK, "cache", f"ML_box_{S}.npy"); t0 = time.time()
    Mb = np.load(f); nb, nk = Mb.shape[:2]
    Op = Mb.reshape(nb * nk, nb * nk)
    out = dict(size=S, N_cells=N_cells, shape=list(Mb.shape))
    out["hermiticity_rel"] = float(np.abs(Op - Op.conj().T).max() / np.abs(Op).max())
    dg = np.diag(Op)
    out["diag_Ncells"] = dict(max_abs_dev=float(np.abs(dg - N_cells).max()), rel=float(np.abs(dg - N_cells).max() / N_cells))
    del Op
    k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k); D = int(MP[0]); p = D // n
    kint = np.mod(np.rint(k * D).astype(int), D)
    co = np.where((kint[:, 0] % p == 0) & (kint[:, 1] % p == 0))[0]
    blk = Mb[:, co][:, :, :, co].reshape(nb * len(co), nb * len(co))
    out["coincident_k"] = dict(n=int(len(co)), max_rel_dev_from_Ncells_identity=float(np.abs(blk - N_cells * np.eye(len(blk))).max() / N_cells))
    log(f"[A.2d post] {S} : hermiticité {out['hermiticity_rel']:.1e} ; diagonale − N_cells max {out['diag_Ncells']['max_abs_dev']:.2e} ; "
        f"{len(co)} k coïncidents : max|M − N_cells 𝕀|/N_cells = {out['coincident_k']['max_rel_dev_from_Ncells_identity']:.2e}")
    # rotation en base de Wannier (même chaîne que la production : k du XML, U de la wannierisation dense, boîte MP simple)
    U, Ud, Hwr, Rw, nd = load_wannier(dp, k)
    Fw, R = Mwk_to_Mwr(Mbk_to_Mwk(Mb, U, Ud), k, MP); del Mb
    # R_d de la production : recentrage de M_W(M2 dense) (le maximum de ||F_W(R,R)|| n'est pas défini pour F ≈ 𝕀)
    M2 = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    Mw2, R2 = Mwk_to_Mwr(Mbk_to_Mwk(M2, U, Ud), k, MP); del M2
    assert np.array_equal(R, R2)
    Rn, R_d = lt.recenter_mwr(Mw2, R2, MP); lt.mwr_locality(Mw2, Rn); del Mw2
    inb, edge = box_geometry(R, MP, n)
    nR = len(R); I5 = np.eye(NW)
    diagblk = np.array([Fw[:, i, :, i] for i in range(nR)])                        # (nR, 5, 5)
    ddiag = np.array([np.abs(diagblk[i] - (I5 if inb[i] else 0 * I5)).max() for i in range(nR)])
    sumrule = [float(diagblk[:, w, w].real.sum()) for w in range(NW)]
    offsite = np.abs(Fw).copy()
    for i in range(nR):
        offsite[:, i, :, i] = 0.0
    off_max = float(offsite.max()); del offsite
    by_edge = {}
    for e in sorted(set(edge.tolist())):
        m = edge == e
        by_edge[int(e)] = dict(n_sites=int(m.sum()), in_box=bool(inb[m][0]), max_abs_dev=float(ddiag[m].max()),
                               diag_min=[float(diagblk[m, w, w].real.min()) for w in range(NW)],
                               diag_max=[float(diagblk[m, w, w].real.max()) for w in range(NW)])
    clusters = {}
    for rc in (3, 4):
        sel = np.where(np.linalg.norm(Rn, axis=1) <= rc + 1e-9)[0]
        Fc = Fw[:, sel][:, :, :, sel]
        Pb = np.zeros((NW, len(sel), NW, len(sel)))
        for j, i in enumerate(sel):
            if inb[i]:
                Pb[:, j, :, j] = I5
        dev = np.abs(Fc - Pb)
        rows = []
        for j, i in enumerate(sel):
            rows.append(dict(Rn=Rn[i].tolist(), in_box=bool(inb[i]), edge=int(edge[i]), F_diag=[float(diagblk[i, w, w].real) for w in range(NW)],
                             max_abs_dev_onsite=float(ddiag[i])))
        clusters[f"R_cut_{rc}"] = dict(n_sites=int(len(sel)), n_out_of_box=int((~inb[sel]).sum()), max_abs_dev_F_minus_Pbox=float(dev.max()),
                                       max_abs_dev_offsite=float(max(dev[:, j, :, jj].max() for j in range(len(sel)) for jj in range(len(sel)) if j != jj)),
                                       sites=rows)
    out.update(R_d=R_d.tolist(), MP=list(MP), sumrule_per_orbital=sumrule, max_abs_offsite=off_max, by_edge_distance=by_edge, clusters=clusters,
               elapsed_s=time.time() - t0)
    d = ensure("a2"); save_json(os.path.join(d, f"a2d_{S}.json"), out)
    np.savez(os.path.join(WORK, "cache", f"Fw_{S}.npz"), Fw=Fw, R=R, Rn=Rn, R_d=R_d, MP=np.array(MP), in_box=inb, edge=edge)
    log(f"[A.2d post] {S} : R_d {R_d.tolist()} ; Σ_R F_ww(R,R) = {np.round(sumrule, 6).tolist()} (attendu {N_cells}) ; max|F(R,R') hors site| {off_max:.2e}")
    for e, v in by_edge.items():
        log(f"[A.2d post] {S} distance au bord {e:+d} ({v['n_sites']} sites, {'boîte' if v['in_box'] else 'hors boîte'}) : max|F(R,R) − {'𝕀' if v['in_box'] else '0'}| "
            f"{v['max_abs_dev']:.2e} ; diag min {np.round(v['diag_min'], 5).tolist()} max {np.round(v['diag_max'], 5).tolist()}")
    for key, v in clusters.items():
        log(f"[A.2d post] {S} {key} : {v['n_sites']} sites dont {v['n_out_of_box']} hors boîte ; max|F − P_boîte| {v['max_abs_dev_F_minus_Pbox']:.2e} "
            f"(hors site {v['max_abs_dev_offsite']:.2e})")


# ----------------------------------------------------------------------------------------------- mise en place commune
SIZES_A = ["5x5", "6x6", "7x7", "8x8", "9x9", "12x12"]
FAMILIES = {"3m": ["6x6", "9x9", "12x12"], "non-3m": ["5x5", "7x7", "8x8"]}
PLATEAU_FRAC = 0.75           # décision 2 (Greg, 2026-09-28) : atomes à distance >= 0,75 r_max de la lacune (vraie image minimale)
PLATEAU_TOL_MEV = 5.0         # « sans plateau » si max|écart| au plateau > 5 meV (critère appliqué ; à confirmer par Greg)
TOL_I_II_MEV = 5.0            # critère (i) contre (ii) du prompt (inchangé)
LU_PUBLISHED_MEV = {"5x5": -122.33, "6x6": -98.62, "7x7": -38.54, "8x8": -22.97, "9x9": -24.65, "12x12": -29.24}   # R4 D5 = R5 C, 1,0 Å
NN_CELLS = [(0, 0, 0), (-1, 0, 0), (0, -1, 0)]          # mailles des trois p_z(B) voisins de la lacune A en Rn = 0 (R4)
WF_PZ_A, WF_PZ_B = 3, 4
CN = None                     # clôture : None = R9 (C_N retenu par A.1) ; "lu" | "plateau" (option --cn)


def cn_suffix():
    """Suffixe des fichiers de sortie : "" (R9), "_lu", "_plateau"."""
    return "" if CN is None else f"_{CN}"


def vkey(var):
    """Clé de variante dans les sorties : « brut » inchangé ; avec --cn : « aligne_lu », « exact_plateau », …"""
    return var if (CN is None or var == "brut") else f"{var}_{CN}"


def hwk_chunked(Hwr, Rw, nd, k, chunk=400_000):
    """H_W(k) par tranches de k (même formule que Hwr_to_Hwk ; une seule tranche, donc identique à la production, pour N_k^int <= 600)."""
    out = np.empty((len(k), Hwr.shape[1], Hwr.shape[1]), dtype=complex)
    for s in range(0, len(k), chunk):
        out[s:s + chunk] = lt.Hwr_to_Hwk(Hwr, Rw, k[s:s + chunk], ndegen=nd)[0]
    return out


def mwr_cached(S, cfg):
    """
    M_W(R, R') du M2 dense de production (V† M V, double TF sur la boîte MP D×D, k du XML : chaîne de compute_spectral_wannier.py et de
    resonance_metrics.py), étiquettes recentrées ; cache/Mwr_<S>.npz, invalidé si la taille ou la date du fichier M2 changent.
    """
    dp = dense_paths(cfg, S); f = os.path.join(ensure("cache"), f"Mwr_{S}.npz"); src = os.stat(dp["mfile"])
    if os.path.exists(f):
        Z = np.load(f)
        if int(Z["src_size"]) == src.st_size and float(Z["src_mtime"]) == src.st_mtime:
            return dict(Mwr=Z["Mwr"], R=Z["R"], Rn=Z["Rn"], R_d=Z["R_d"], MP=tuple(int(x) for x in Z["MP"]))
    t0 = time.time()
    M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k)
    U, Ud, _, _, _ = load_wannier(dp, k)
    Mwr, R = Mwk_to_Mwr(Mbk_to_Mwk(M, U, Ud), k, MP); del M
    Rn, R_d = lt.recenter_mwr(Mwr, R, MP); lt.mwr_locality(Mwr, Rn)
    np.savez(f, Mwr=Mwr, R=R, Rn=Rn, R_d=R_d, MP=np.array(MP), src_size=src.st_size, src_mtime=src.st_mtime)
    log(f"[Mwr] {S} : rotation du M2 dense en {time.time()-t0:.0f} s, R_d {R_d.tolist()} -> {os.path.basename(f)}")
    return dict(Mwr=Mwr, R=R, Rn=Rn, R_d=R_d, MP=tuple(MP))


def wannier_setup(S, cfg, rcuts=(3,)):
    """M_W (cache), U, H(R) (porte de jauge), boîte, V_loc par R_cut (extract_V_loc), E_D (grille 90², comme la production)."""
    dp = dense_paths(cfg, S); n = size_n(S)
    W = mwr_cached(S, cfg)
    k = qe_io.get_k_red(dp["uc"]); U, Ud, Hwr, Rw, nd = load_wannier(dp, k)
    inb, edge = box_geometry(W["R"], W["MP"], n)
    loc = {}
    for rc in rcuts:
        Rloc = W["Rn"][np.linalg.norm(W["Rn"], axis=1) <= rc + 1e-9]
        V, herm = lt.extract_V_loc(W["Mwr"], W["Rn"], Rloc)
        iloc = np.array([int(np.where((W["Rn"] == r).all(axis=1))[0][0]) for r in Rloc])
        loc[rc] = dict(Rloc=Rloc, V=V, herm=herm, iloc=iloc, in_box=inb[iloc], edge=edge[iloc])
    _, E_ref, _ = lt.Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
    gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
    return dict(S=S, n=n, dp=dp, k=k, U=U, Ud=Ud, Hwr=Hwr, Rw=Rw, nd=nd, inb=inb, edge=edge, loc=loc, E_D=E_D, **W)


def resonance_setup(S, cfg, rcuts=(3,)):
    """wannier_setup + états 240², grille en énergie de resonance_metrics.py (« res ») et de rcut_resigma.py / niveau 1 (« sig »)."""
    st = wannier_setup(S, cfg, rcuts)
    N, eta, ew, npe = cfg["grid"], cfg["eta_eV"], cfg["e_window_eV"], cfg["ne_per_eta"]
    k_out = lt.mp_grid(N, N, 1); _, E_out, U_out = lt.Hwr_to_Hwk(st["Hwr"], st["Rw"], k_out, ndegen=st["nd"])
    E_D = st["E_D"]; sel = np.abs(E_out - E_D) <= ew; de = eta / npe
    egrid = np.arange(E_D - ew - eta, E_D + ew + eta + de, de)                                   # resonance_metrics.py l. 62
    E_sel = E_out[sel]; egrid_sig = np.arange(E_sel.min() - eta, E_sel.max() + eta + de, de)    # rcut_resigma.py l. 45, scattering_rate_fast
    iK = int(np.argmin(np.linalg.norm(np.mod(k_out - np.array(cfg["K_red"]) + 0.5, 1) - 0.5, axis=1)))
    A_uc, _ = qe_io.get_A_volume(st["dp"]["uc"]); Bc = 2 * np.pi * np.linalg.inv(A_uc * BOHR).T                 # colonnes b_i (Å⁻¹)
    st.update(ew=ew, k_out=k_out, E_out=E_out, U_out=U_out, sel=sel, de=de, eta=eta, egrid=egrid, egrid_sig=egrid_sig, iK=iK, Bc=Bc, N=N, phi={})
    L3 = st["loc"][min(rcuts)]
    st.update(V=L3["V"], Rloc=L3["Rloc"], loc_in_box=L3["in_box"])                           # compatibilité (R_cut de production = 3)
    for rc in rcuts:
        L = st["loc"][rc]
        log(f"[setup] {S} R_cut {rc} : R_d {st['R_d'].tolist()}, {len(L['Rloc'])} sites (dim {len(L['Rloc'])*NW}, {int((~L['in_box']).sum())} hors boîte), "
            f"herm V_loc {L['herm']:.1e}")
    log(f"[setup] {S} : E_D {E_D:.6f} eV ; {int(sel.sum())} états à ±{ew} eV sur {N}² ; grilles « res » {len(egrid)} et « sig » {len(egrid_sig)} énergies, "
        f"pas {de*1e3:.2f} meV ; K = k_out[{iK}]")
    return st


def phi_of(st, rc):
    """<wR|nk> = U_out[k, w, n] e^{2πi k.R} sur l'amas du R_cut (disposition L*nw + w), mémorisé."""
    if rc not in st["phi"]:
        Rloc = st["loc"][rc]["Rloc"]; nL = len(Rloc)
        st["phi"][rc] = np.einsum("kL,kwn->knLw", lt._phase(st["k_out"], Rloc), st["U_out"], optimize=True).reshape(len(st["k_out"]), NW, nL * NW)
    return st["phi"][rc]


def g0_cached(st, nk_int, rc=3, grid="res"):
    """g₀ sur l'amas du R_cut pour une grille en énergie (local_green_batch, exact) ; cache npy vérifié par ses métadonnées."""
    egrid = st["egrid"] if grid == "res" else st["egrid_sig"]; Rloc = st["loc"][rc]["Rloc"]
    name = f"g0_{st['S']}_nk{nk_int}.npy" if (rc == 3 and grid == "res") else f"g0_{st['S']}_nk{nk_int}_rc{rc}_{grid}.npy"
    f = os.path.join(ensure("cache"), name); fj = f.replace(".npy", ".json")
    meta = dict(S=st["S"], nk_int=nk_int, E_D=st["E_D"], e0=float(egrid[0]), de=st["de"], nE=len(egrid), eta=st["eta"], nL=len(Rloc), Rloc=Rloc.tolist())
    if os.path.exists(f) and os.path.exists(fj):
        old = json.load(open(fj))
        if all(old.get(k) == v if not isinstance(v, float) else abs(old.get(k, 1e9) - v) < 1e-12 for k, v in meta.items()):
            log(f"[g0] {st['S']} {nk_int}² R_cut {rc} « {grid} » : cache {name}"); return np.load(f)
    t0 = time.time()
    k_int = lt.mp_grid(nk_int, nk_int, 1); Hk = hwk_chunked(st["Hwr"], st["Rw"], st["nd"], k_int)
    g0 = lt.local_green_batch(Hk, k_int, Rloc, egrid, st["eta"]); del Hk
    np.save(f, g0); save_json(fj, dict(meta, elapsed_s=time.time() - t0, head=git_head()))
    log(f"[g0] {st['S']} {nk_int}² R_cut {rc} « {grid} » : {len(egrid)} énergies en {(time.time()-t0)/60:.1f} min -> {name}")
    return g0


# ----------------------------------------------------------------------------------------------- variantes de V_loc (alignement)
def C_of(S):
    """C_N retenu par A.1 (eV) et son étiquette (« plateau », ou « sans plateau » / « écart (i)-(ii) » avec la valeur Lu)."""
    f = os.path.join(WORK, "a", "a1_results.json")
    if not os.path.exists(f):
        raise SystemExit(f"C_N de {S} : {f} absent (A.1 non exécutée)")
    r = json.load(open(f))["sizes"][S]
    if CN == "plateau":                                                        # clôture : moyenne du plateau (i), 1,0 Å, d >= 0,75 r_max
        return float(r["C_i_eV"]), "plateau (i)"
    if CN == "lu":                                                             # clôture : Lu (1,0 Å) ; = C_retenu de R9 (six tailles « sans plateau »)
        assert "Lu" in r["etiquette"], f"{S} : C_retenu n'est pas Lu ({r['etiquette']})"
        return float(r["C_retenu_eV"]), "Lu (1,0 Å)"
    return float(r["C_retenu_eV"]), r["etiquette"]


def V_variant(st, rc, var, C_eV):
    """
    brut   : V_loc de production ;
    aligne : approximation (i) de 0.2, V_loc − C_N·P_boîte (𝕀₅ sur les sites de l'amas dont la maille est dans la boîte N×N) ;
    exact  : V_loc − C_N·P F_W P, F_W(R, R') = <w_R|1_boîte|w_R'> de A.2 (cache/Fw_<S>.npz ; 9x9 et 5x5).
    """
    L = st["loc"][rc]; V = L["V"]; var = var.split("_")[0]                # « aligne_plateau » -> « aligne » (clôture)
    if var == "brut":
        return V
    if var == "aligne":
        return V - C_eV * np.diag(np.repeat(L["in_box"].astype(float), NW))
    if var == "exact":
        Z = np.load(os.path.join(WORK, "cache", f"Fw_{st['S']}.npz"))
        assert np.array_equal(Z["R"], st["R"]), "F_W et M_W : étiquettes R différentes"
        idx = L["iloc"]; nL = len(idx)
        F = np.transpose(Z["Fw"][:, idx][:, :, :, idx], (1, 0, 3, 2)).reshape(nL * NW, nL * NW)     # disposition d'extract_V_loc
        return V - C_eV * F
    raise ValueError(var)


def has_exact(S):
    return os.path.exists(os.path.join(WORK, "cache", f"Fw_{S}.npz"))


# ----------------------------------------------------------------------------------------------- observables (formules de production)
def lor(x, eta):
    return (eta / np.pi) / (x * x + eta * eta)


def jn_sig(st):
    """Point de la grille « sig » le plus proche de chaque état (argmin, comme rcut_resigma.py l. 46 et scattering_rate_fast)."""
    if "jn_sig" not in st:
        eg = st["egrid_sig"]; E = st["E_out"]; jn = np.empty(E.shape, int)
        for s in range(0, len(E), 2048):
            jn[s:s + 2048] = np.abs(eg[None, None, :] - E[s:s + 2048, :, None]).argmin(-1)
        st["jn_sig"] = jn
    return st["jn_sig"]


def sigma_states(st, V, g0s, rc):
    """Σ_nk = <nk|t(ε_nk)|nk> (complexe) des états de la fenêtre, t = V[1 − g₀V]⁻¹ au point « sig » le plus proche (niveau 1, rcut_resigma)."""
    phi = phi_of(st, rc); sel = st["sel"]; jn = jn_sig(st); I = np.eye(V.shape[0])
    Sig = np.full(st["E_out"].shape, np.nan + 0j)
    for j in np.unique(jn[sel]):
        t = V @ np.linalg.solve(I - g0s[j] @ V, I)
        m = sel & (jn == j); P = phi[m]
        Sig[m] = np.einsum("mi,mi->m", P.conj() @ t, P)
    return Sig


def level1_stats(st, Sig):
    """Médiane Γ·N_cells et E_res (compute_spectral_wannier.py), Re Σ médian et |Re Σ|/Γ médian (rcut_resigma.py), meV / eV."""
    G = -2 * Sig.imag; Re = Sig.real; m = st["sel"]; E = st["E_out"]; E_D = st["E_D"]
    mm = m & (np.abs(E - E_D) <= 1.5); a_ = int(np.argmax(G[mm])); ik, n = np.argwhere(mm)[a_]
    return dict(median_meV=float(np.median(G[m])) * 1e3, E_res=float(E[mm][a_] - E_D), median_ReSigma_meV=float(np.median(Re[m])) * 1e3,
                median_absRe_over_G=float(np.median(np.abs(Re[m]) / G[m])), min_G=float(G[m].min()), n_states=int(m.sum()),
                E_res_state=(int(ik), int(n)))


def resonance_observables(st, V, g0, rc=3, born=False):
    """Mêmes formules que resonance_metrics.py (l. 67–120) sur la grille « res » : t, Γ par état (rint), courbe lorentzienne, T̄(K), Born."""
    egrid, de, eta, E_out, sel, E_D = st["egrid"], st["de"], st["eta"], st["E_out"], st["sel"], st["E_D"]
    phi = phi_of(st, rc); nE = len(egrid); I = np.eye(V.shape[0])
    t = np.array([V @ np.linalg.solve(I - g0[j] @ V, I) for j in range(nE)])
    j_all = np.rint((E_out - egrid[0]) / de).astype(int)

    def onshell(tc):
        G = np.full(E_out.shape, np.nan)
        for j in np.unique(j_all[sel]):
            m = sel & (j_all == j); P = phi[m]
            G[m] = -2.0 * np.einsum("mi,mi->m", P.conj() @ tc[j], P).imag
        return G
    G = onshell(t); Es = E_out[sel]

    def curve(g, eg):
        w = lor(eg[:, None] - Es[None, :], eta); return (w * g[sel][None, :]).sum(1) / w.sum(1)
    eg5 = egrid[::2]; c5 = curve(G, eg5); c25 = curve(G, egrid)
    PK = phi[st["iK"], 3:5]
    tr = np.array([0.5 * np.trace(PK.conj() @ t[j] @ PK.T) for j in range(nE)])
    mm = np.isfinite(G) & (np.abs(E_out - E_D) <= 1.5)
    ik_all, n_all = np.nonzero(mm); a = int(np.argmax(np.abs(G[mm])))
    atED = sel & (np.abs(E_out - E_D) < 1e-6)
    out = dict(G=G, curve_eg5=c5, curve_eg25=c25, tr=tr,
               peak_GT_prod=float(eg5[int(np.argmax(c5))] - E_D), peak_GT_fine=float(egrid[int(np.argmax(c25))] - E_D), max_GT_fine=float(c25.max()),
               peak_ImTbar=float(egrid[int(np.argmax(-tr.imag))] - E_D), max_mImTbar=float((-tr.imag).max()),
               ReTbar_at_ED=float(np.interp(E_D, egrid, tr.real)), G_ED_meV=float(np.abs(G[atED]).mean() * 1e3), n_ED=int(atED.sum()),
               median_GT_states_meV=float(np.nanmedian(np.abs(G[sel])) * 1e3),
               E_res=float(E_out[ik_all[a], n_all[a]] - E_D), E_res_state=(int(ik_all[a]), int(n_all[a])))
    if born:
        tB = np.array([V + V @ g0[j] @ V for j in range(nE)]); GB = onshell(tB); cB = curve(GB, eg5)
        mw = np.abs(eg5 - E_D) <= st["ew"]                                     # resonance_metrics.py l. 105 : courbes sur ±e_window
        out.update(median_GT_curve_meV=float(np.nanmedian(c5[mw]) * 1e3), median_GBorn_curve_meV=float(np.nanmedian(cB[mw]) * 1e3),
                   Born_over_T_median=float(np.nanmedian(cB[mw] / c5[mw])), Born_over_T_min=float(np.nanmin(cB[mw] / c5[mw])),
                   Born_over_T_max=float(np.nanmax(cB[mw] / c5[mw])), peak_GBorn=float(eg5[int(np.argmax(cB))] - E_D))
    return out


def crown_table(st, band=3, lo=-0.40, hi=-0.02):
    """
    Couronnes de la grille de sortie N² autour de K et K' dans la bande `band` (3 = π) : x = |i b₁ + j b₂|²/|b₁|² (norme de Löschian du
    décalage (i, j)/N à la vallée la plus proche), énergie des états (ε − E_D) par x ; rang = ordre des x distincts > 0 réalisés.
    """
    k, E, N, Bc = st["k_out"], st["E_out"][:, band] - st["E_D"], st["N"], st["Bc"]
    b1 = Bc[:2, 0]; B2 = Bc[:2, :2]
    xs = np.full(len(k), np.inf)
    for Kv in (K_RED, KP_RED):
        dk = np.mod(k - Kv + 0.5, 1.0) - 0.5; ij = np.rint(dk[:, :2] * N)
        q = ij @ B2.T; x = (q * q).sum(1) / (b1 @ b1)
        xs = np.minimum(xs, x)
    xr = np.rint(xs).astype(int); assert np.abs(xs - xr).max() < 1e-6, "x non entier"
    distinct = np.unique(xr[xr > 0]); rank = {int(x): i + 1 for i, x in enumerate(distinct)}
    rows = []
    for x in distinct:
        m = xr == x; e = E[m]
        if e.mean() < lo or e.mean() > hi:
            continue
        rows.append(dict(x=int(x), rank=rank[int(x)], n_states=int(m.sum()), E_mean=float(e.mean()), E_min=float(e.min()), E_max=float(e.max())))
    rows.sort(key=lambda r: r["E_mean"])
    return rows, xr, rank


def crown_of_state(st, ik, n, xr, rank):
    E = st["E_out"][:, n] - st["E_D"]; x = int(xr[ik]); m = xr == x
    return dict(x=x, rank=rank.get(x, 0), crown_E_mean=float(E[m].mean()), crown_E_min=float(E[m].min()), crown_E_max=float(E[m].max()), crown_n=int(m.sum()))


def crowns_containing(rows, v, tol=5e-4):
    """Couronnes dont l'intervalle [min, max] des énergies contient la valeur v (à ±tol : E_res de R6 imprimés au meV)."""
    return [r for r in rows if r["E_min"] - tol <= v <= r["E_max"] + tol]


# ----------------------------------------------------------------------------------------------- B : résonance contre N_k^int
R6_ERES = {"9x9": {"R6 3.3 N_k^int 150": -0.238, "R6 3.3 N_k^int 300 (production)": -0.175, "R6 3.3 N_k^int 450": -0.202,
                   "R6 3.3 N_k^int 600": -0.202, "R6 3.3 R_cut 1": -0.227},
           "12x12": {"R6 3.3 production": -0.175}}


def cmd_b(a):
    cfg = load_production(verbose=True); d = ensure("b")
    sizes = a.sizes.split(","); nks = [int(x) for x in a.nk.split(",")]; variants = a.variants.split(",")
    fres = os.path.join(d, f"b_results{cn_suffix()}.json"); res = json.load(open(fres)) if os.path.exists(fres) else {}
    res.setdefault("meta", {}).update(head=git_head(), config=dict(R_cut=cfg["R_cut"], grid=cfg["grid"], eta=cfg["eta_eV"], ne_per_eta=cfg["ne_per_eta"],
                                                                    e_window=cfg["e_window_eV"]))
    for S in sizes:
        st = resonance_setup(S, cfg)
        rows, xr, rank = crown_table(st)
        res.setdefault(S, {})["crowns_pi_240"] = rows
        res[S]["E_D"] = st["E_D"]; res[S]["R_d"] = st["R_d"].tolist(); res[S]["loc_out_of_box"] = int((~st["loc_in_box"]).sum())
        sp = np.diff([r["E_mean"] for r in rows])
        log(f"[B] {S} couronnes π de la grille {st['N']}² dans [−0,40, −0,02] eV : " + " ; ".join(f"x={r['x']} (rang {r['rank']}) {r['E_mean']:+.4f}" for r in rows)
            + f" ; écarts {np.round(sp*1e3, 1).tolist()} meV")
        r6 = {}
        for lab, v in R6_ERES.get(S, {}).items():
            r6[lab] = dict(E_res=v, crowns=crowns_containing(rows, v))
            log(f"[B] {S} {lab} : E_res {v:+.3f} -> couronnes {[(r['x'], r['rank'], round(r['E_min'], 4), round(r['E_max'], 4)) for r in r6[lab]['crowns']]}")
        res[S]["R6_E_res_crowns"] = r6
        cf = os.path.join(d, f"b_curves_{S}{cn_suffix()}.npz"); curves = dict(np.load(cf)) if os.path.exists(cf) else {}
        run = list(variants)
        if S == "9x9" and "brut" not in run:
            run = ["brut"] + run                                                                  # la porte B.0 précède toujours la partie
        for var in run:
            if var.startswith("exact") and not has_exact(S):
                log(f"[B] {S} : pas de F_W (A.2 exact seulement pour 9x9, 5x5), variante exacte sautée"); continue
            C, tag = (0.0, "—") if var == "brut" else C_of(S)
            V = V_variant(st, 3, var, C)
            for nk in (nks if var in variants else [300]):
                g0 = g0_cached(st, nk); t0 = time.time()
                o = resonance_observables(st, V, g0); del g0
                cr = crown_of_state(st, o["E_res_state"][0], o["E_res_state"][1], xr, rank)
                row = {k: v for k, v in o.items() if k not in ("G", "curve_eg5", "curve_eg25", "tr")}
                row.update(E_res_crown=cr, C_eV=C, C_label=tag, elapsed_s=time.time() - t0)
                if var in variants:
                    res[S].setdefault(vkey(var), {})[str(nk)] = row
                    curves[f"{vkey(var)}_nk{nk}_GT"] = o["curve_eg25"]; curves[f"{vkey(var)}_nk{nk}_tr"] = o["tr"]
                log(f"[B] {S} {var} N_k^int {nk} : pic Γ_T {o['peak_GT_fine']:+.4f} (courbe production pas 5 meV {o['peak_GT_prod']:+.3f}) ; pic −Im T̄(K) "
                    f"{o['peak_ImTbar']:+.4f} (max {o['max_mImTbar']:.2f}) ; Re T̄(E_D) {o['ReTbar_at_ED']:+.3f} ; Γ_T(E_D) {o['G_ED_meV']:.1f} meV ({o['n_ED']} états) ; "
                    f"médiane états {o['median_GT_states_meV']:.2f} meV ; E_res {o['E_res']:+.4f} (couronne x={cr['x']}, rang {cr['rank']}, "
                    f"{cr['crown_E_mean']:+.4f}) ; {time.time()-t0:.0f} s")
                if S == "9x9" and var == "brut" and nk == 300:                                  # porte B.0 contre R6 (resonance_9x9.npz)
                    z = np.load(os.path.join(results_dir(cfg), "resonance_9x9.npz"))
                    dc = float(np.abs(o["curve_eg5"] - z["Gamma_T"]).max() / np.abs(z["Gamma_T"]).max())
                    dt = float(np.abs(o["tr"] - z["Tbar_tr"]).max() / np.abs(z["Tbar_tr"]).max())
                    ok = (abs(o["peak_GT_prod"] - (-0.180)) < 5e-4) and (abs(o["peak_ImTbar"] - (-0.177)) < 1e-3) and dc < 1e-6 and dt < 1e-6
                    res["gate_B0"] = dict(peak_GT_prod=o["peak_GT_prod"], expected_peak_GT=-0.180, peak_ImTbar=o["peak_ImTbar"], expected_peak_ImTbar=-0.177,
                                          curve_rel_dev_vs_R6=dc, Tbar_rel_dev_vs_R6=dt, median_states_meV=o["median_GT_states_meV"],
                                          expected_median_states_meV=float(z["median_GT_states_meV"]), ok=bool(ok), head=git_head())
                    log(f"[B.0] porte : pic Γ_T {o['peak_GT_prod']:+.3f} (R6 −0,180), pic −Im T̄(K) {o['peak_ImTbar']:+.4f} (R6 −0,177) ; courbe vs resonance_9x9.npz "
                        f"{dc:.1e}, T̄ {dt:.1e} ; médiane {o['median_GT_states_meV']:.2f} (R6 {float(z['median_GT_states_meV']):.2f}) -> {'OK' if ok else 'ÉCHEC'}")
                    save_json(fres, res)
                    if not ok:
                        log("[B.0] ÉCHEC de la porte : STOP sur la partie B"); sys.exit(3)
                save_json(fres, res)
        np.savez(cf, egrid=st["egrid"], E_D=st["E_D"], **{k: v for k, v in curves.items() if k not in ("egrid", "E_D")})
    save_json(fres, res)
    b_tables(res); b_figure(res)


def cmd_btables(a):
    res = json.load(open(os.path.join(WORK, "b", f"b_results{cn_suffix()}.json")))
    for S in ("9x9", "12x12"):                                                  # correspondance des E_res de R6 par intervalle de couronne
        if S in res:
            res[S]["R6_E_res_crowns"] = {lab: dict(E_res=v, crowns=crowns_containing(res[S]["crowns_pi_240"], v)) for lab, v in R6_ERES.get(S, {}).items()}
    save_json(os.path.join(WORK, "b", f"b_results{cn_suffix()}.json"), res)
    b_tables(res); b_figure(res)


VAR_LABEL = {"brut": "M2 tel quel", "aligne": "aligné (i)", "exact": "aligné exact (F_W)",
             "aligne_lu": "aligné (i), C_N Lu", "exact_lu": "aligné exact (F_W), C_N Lu",
             "aligne_plateau": "aligné (i), C_N plateau", "exact_plateau": "aligné exact (F_W), C_N plateau"}


def b_tables(res):
    L = ["# R9 — B : position de la résonance contre N_k^int (tables générées par r9_driver.py)", ""]
    if "gate_B0" in res:
        g = res["gate_B0"]
        L += [f"Porte B.0 (9×9, M2 tel quel, 300²) : pic de la courbe Γ_T {g['peak_GT_prod']:+.3f} eV (R6 −0,180), pic de −Im T̄(K) {g['peak_ImTbar']:+.4f} "
              f"(R6 −0,177), écart relatif des courbes à `resonance_9x9.npz` {g['curve_rel_dev_vs_R6']:.1e} (T̄ {g['Tbar_rel_dev_vs_R6']:.1e}) → "
              f"{'OK' if g['ok'] else 'ÉCHEC'}", ""]
    for S in ("9x9", "12x12"):
        if S not in res:
            continue
        for var in [v for v in VAR_LABEL if v in res[S]]:
            r0 = next(iter(res[S][var].values()))
            cl = "" if var == "brut" else f" ; C_N = {r0['C_eV']*1e3:+.2f} meV ({r0['C_label']})"
            L += [f"## {S}, {VAR_LABEL[var]} (énergies en eV relatives à E_D = {res[S]['E_D']:.6f} eV ; courbes au pas de 2,5 meV{cl})", "",
                  "| N_k^int | pic Γ_T | max Γ_T (eV) | pic Γ_T, courbe prod. (5 meV) | pic −Im T̄(K) | max −Im T̄ (eV) | Re T̄(E_D) | Γ_T(E_D) (meV) | médiane états (meV) | E_res | couronne de E_res : x (rang), ε − E_D moyenne [min, max] |",
                  "|---|---|---|---|---|---|---|---|---|---|---|"]
            for nk, r in sorted(res[S][var].items(), key=lambda kv: int(kv[0])):
                c = r["E_res_crown"]
                L.append(f"| {nk} | {r['peak_GT_fine']:+.4f} | {r['max_GT_fine']:.3f} | {r['peak_GT_prod']:+.3f} | {r['peak_ImTbar']:+.4f} | {r['max_mImTbar']:.2f} | "
                         f"{r['ReTbar_at_ED']:+.3f} | {r['G_ED_meV']:.1f} | {r['median_GT_states_meV']:.2f} | {r['E_res']:+.4f} | {c['x']} ({c['rank']}), "
                         f"{c['crown_E_mean']:+.4f} [{c['crown_E_min']:+.4f}, {c['crown_E_max']:+.4f}] |")
            L.append("")
        L += [f"Couronnes π de la grille 240² ({S}), ε − E_D dans [−0,40 ; −0,02] eV :", "", "| x | rang | n états | ε − E_D moyenne | min | max |", "|---|---|---|---|---|---|"]
        for r in res[S]["crowns_pi_240"]:
            L.append(f"| {r['x']} | {r['rank']} | {r['n_states']} | {r['E_mean']:+.4f} | {r['E_min']:+.4f} | {r['E_max']:+.4f} |")
        L += ["", f"E_res de R6 en regard ({S}) — couronnes dont l'intervalle [min, max] contient la valeur à ±0,5 meV :", ""]
        for lab, v in res[S].get("R6_E_res_crowns", {}).items():
            L.append(f"- {lab} : {v['E_res']:+.3f} eV → " + (", ".join(f"x = {c['x']} (rang {c['rank']}, [{c['E_min']:+.4f}, {c['E_max']:+.4f}])" for c in v["crowns"])
                                                          or "aucune couronne à ±0,5 meV"))
        L.append("")
    with open(os.path.join(WORK, "b", f"B_tables{cn_suffix()}.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def fig_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use(os.path.join(PROJ, "figures", "memoire.mplstyle"))
    import _palette as pal                                                     # scripts/_palette.py
    return plt, pal


def savefig(fig, name):
    d = ensure("fig")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(d, f"{name}.{ext}"), dpi=200 if ext == "png" else None)
    log(f"[fig] {name}.pdf/.png")


def b_figure(res):
    plt, pal = fig_style()
    sizes = [S for S in ("9x9", "12x12") if S in res and "brut" in res[S]]
    fig, axes = plt.subplots(1, len(sizes), figsize=(6.5, 3.0), sharey=True, squeeze=False)
    for ax, S, lab in zip(axes[0], sizes, "ab"):
        for var, ls, mk in [(v, l_, m_) for v, l_, m_ in (("brut", "-", "o"), ("aligne", "--", "s"), ("exact", ":", "^"), ("aligne_lu", "--", "s"),
                                                          ("exact_lu", ":", "^"), ("aligne_plateau", "-.", "D"), ("exact_plateau", ":", "v")) if v in res[S]]:
            rr = sorted(res[S][var].items(), key=lambda kv: int(kv[0])); x = [1.0 / int(nk) for nk, _ in rr]
            suf = "" if var == "brut" else f" ({VAR_LABEL[var]})"
            ax.plot(x, [r["peak_GT_fine"] for _, r in rr], ls, marker=mk, ms=3, color=pal.NAVY, label=r"pic de $\Gamma_T$" + suf)
            ax.plot(x, [r["peak_ImTbar"] for _, r in rr], ls, marker=mk, ms=3, color=pal.ORANGE, label=r"pic de $-\mathrm{Im}\,\bar T(K)$" + suf)
            ax.plot(x, [r["E_res"] for _, r in rr], ls, marker=mk, ms=3, color=pal.GREEN, label=r"$E_\mathrm{res}$ (états $240^2$)" + suf)
        ax.set_xlabel(r"$1/N_k^\mathrm{int}$"); ax.set_title(f"({lab}) {S.replace('x', '×')}", fontsize=9)
    axes[0][0].set_ylabel(r"Énergie $\varepsilon - E_D$ (eV)"); axes[0][-1].legend(fontsize=5, loc="lower left")
    fig.tight_layout(); savefig(fig, "resonance_vs_nkint" + cn_suffix()); plt.close(fig)


# ----------------------------------------------------------------------------------------------- A.0 et A.1
def plateau_stats(values, mask):
    v = np.asarray(values)[mask]; mu = float(v.mean())
    return dict(mean=mu, rms=float(np.sqrt(((v - mu) ** 2).mean())), max_abs_dev=float(np.abs(v - mu).max()), n=int(mask.sum()))


def cmd_a0(a):
    """Porte A.0 : niveau 1, 9x9, M2 tel quel, R_cut 3, 240², η 0,02, N_k^int 300 -> médiane 3 132,60 meV, E_res −0,175 eV (R6 3.3)."""
    cfg = load_production(verbose=True); d = ensure("a")
    st = resonance_setup("9x9", cfg, (3,))
    g0s = g0_cached(st, cfg["nk_int"], 3, "sig")
    s = level1_stats(st, sigma_states(st, st["loc"][3]["V"], g0s, 3))
    z = np.load(os.path.join(results_dir(cfg), "specwd_9x9_prod.npz"))["results"]
    row = [r for r in z if r[0] == 3 and r[1] == cfg["grid"] and abs(r[2] - cfg["eta_eV"]) < 1e-12][0]
    ok = abs(s["median_meV"] - 3132.60) < 0.01 and abs(s["E_res"] - (-0.175)) < 5e-4 and abs(s["median_meV"] - row[3]) < 1e-6 * row[3]
    out = dict(level1=s, expected=dict(median_meV=3132.60, E_res=-0.175, specwd_row=list(map(float, row)), ReSigma_R6_meV=-26.65), ok=bool(ok), head=git_head())
    save_json(os.path.join(d, "a0_gate.json"), out)
    log(f"[A.0] porte : médiane {s['median_meV']:.2f} meV (R6 3 132,60 ; specwd {row[3]:.4f}), E_res {s['E_res']:+.4f} (R6 −0,175), Re Σ médian "
        f"{s['median_ReSigma_meV']:.2f} meV (R6 −26,65) -> {'OK' if ok else 'ÉCHEC'}")
    if not ok:
        log("[A.0] ÉCHEC de la porte : STOP sur la partie A"); sys.exit(3)


def load_pot_eV(path):
    V, _ = qe_io.get_pot(path, subtract_mean=False, to_hartree=True)
    return V.transpose(2, 1, 0) * HA2EV


def cmd_a1(a):
    """A.1 : C_N de trois façons par taille — (i) sphères de P1, (ii) diagonales de M_W (périodique et boîte), (iii) <ΔV>_3D."""
    cfg = load_production(verbose=True); d = ensure("a"); out = dict(sizes={}, rule=dict(plateau_frac=PLATEAU_FRAC, plateau_tol_meV=PLATEAU_TOL_MEV,
                                                                                         tol_i_ii_meV=TOL_I_II_MEV), head=git_head())
    prof = {}
    for S in a.sizes.split(","):
        t0 = time.time(); n = size_n(S); N_cells = n * n; P = coarse_paths(S); dp = dense_paths(cfg, S)
        A_b, _ = qe_io.get_A_volume(P["scp"].replace("_p.save", "_d.save")); A_A = A_b * BOHR
        x_p = qe_io.get_x_red(P["scp"]); x_d = qe_io.get_x_red(P["scp"].replace("_p.save", "_d.save"))
        Vp = load_pot_eV(P["pot_p"]); Vd = load_pot_eV(P["pot_d"])
        r10 = al.atom_sphere_shifts(Vd, Vp, x_d, x_p, A_A, 1.0); r05 = al.atom_sphere_shifts(Vd, Vp, x_d, x_p, A_A, 0.5)
        lu = al.far_atom_alignment(Vd, Vp, x_d, x_p, A_A, (0.5, 1.0))
        mean3d = float(Vd.mean() - Vp.mean()); del Vd, Vp
        ifar = r10["i_far_axis"]
        chk_lu = dict(P1_at_far_atom_meV=float(r10["shift"][ifar] * 1e3), far_atom_alignment_meV=float(lu["shifts"][1.0] * 1e3),
                      published_meV=LU_PUBLISHED_MEV.get(S), same_as_function=bool(abs(r10["shift"][ifar] - lu["shifts"][1.0]) < 1e-12),
                      P1_05_at_far_atom_meV=float(r05["shift"][ifar] * 1e3), far_atom_alignment_05_meV=float(lu["shifts"][0.5] * 1e3))
        rmax = float(r10["dist"].max()); mfar = r10["dist"] >= PLATEAU_FRAC * rmax
        pl10 = plateau_stats(r10["shift"], mfar); pl05 = plateau_stats(r05["shift"], mfar)
        # (ii) diagonales de M_W. Périodique : M2 dense restreint aux N² k coïncidents (= M de la super-cellule périodique, k commensurables)
        kd = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(kd); D = int(MP[0]); p = D // n
        kint = np.mod(np.rint(kd * D).astype(int), D); co = np.where((kint[:, 0] % p == 0) & (kint[:, 1] % p == 0))[0]
        matrix_io.check_manifest(dp["mfile"], require_normalization=matrix_io.M_NORM_V2)
        Mm = np.load(dp["mfile"], mmap_mode="r"); nb = Mm.shape[0]
        Mco = np.array(Mm[np.ix_(np.arange(nb), co, np.arange(nb), co)]) * HA2EV; del Mm
        U, Ud, _, _, _ = load_wannier(dp, kd)
        kco = np.rint(kd[co] * n) / n                                                         # k = m/N exacts (le XML arrondit à ~1e-7)
        Mwp, Rp = Mwk_to_Mwr(Mbk_to_Mwk(Mco, U[co], Ud[co]), kco, (n, n, 1)); del Mco
        tau = np.mod(qe_io.get_x_red(dp["uc"]), 1.0)                                          # C1 (1/3,1/3) porte sp² et p_z(A), C2 p_z(B)
        wf_atom = [0, 0, 0, 0, 1]

        def orbital_sites(R, cells_mod):
            """Position réduite (super-cellule) de l'atome de chaque orbitale (R, w) et distance vraie au site de la lacune."""
            c = np.mod(np.asarray(R, int)[:, :2], cells_mod)
            xs = np.array([[(c[i, 0] + tau[wf_atom[w], 0]) / n, (c[i, 1] + tau[wf_atom[w], 1]) / n, x_p[0, 2]] for i in range(len(R)) for w in range(NW)])
            dist = al.true_min_image_dist(xs, r10["s_vac"], A_A).reshape(len(R), NW)
            match = np.array([al.true_min_image_dist(x_p, x, A_A).min() for x in xs]).reshape(len(R), NW)
            return dist, match
        dP, mP = orbital_sites(Rp, n)
        diagP = np.array([[Mwp[w, i, w, i].real for w in range(NW)] for i in range(len(Rp))])
        mpz = np.zeros_like(dP, bool); mpz[:, [WF_PZ_A, WF_PZ_B]] = True                     # masque p_z pour les N² mailles périodiques
        mfarP = (dP >= PLATEAU_FRAC * rmax) & mpz
        plP = plateau_stats(diagP, mfarP); plP_sp2 = plateau_stats(diagP, (dP >= PLATEAU_FRAC * rmax) & ~mpz)
        # boîte : M_W dense de production (cache), mailles de la boîte seulement, distance au bord >= 1
        W = mwr_cached(S, cfg); inb, edge = box_geometry(W["R"], W["MP"], n)
        dB, mB = orbital_sites(W["R"], D)
        diagB = np.array([[W["Mwr"][w, i, w, i].real for w in range(NW)] for i in range(len(W["R"]))])
        mpzB = np.zeros_like(dB, bool); mpzB[:, [WF_PZ_A, WF_PZ_B]] = True                  # masque p_z pour les D² mailles de la boîte MP
        mfarB = (dB >= PLATEAU_FRAC * rmax) & mpzB & (inb & (edge >= 1))[:, None]
        plB = plateau_stats(diagB, mfarB) if mfarB.any() else dict(mean=np.nan, rms=np.nan, max_abs_dev=np.nan, n=0)
        C_i = pl10["mean"]; C_ii = plP["mean"]
        has_plateau = pl10["max_abs_dev"] * 1e3 <= PLATEAU_TOL_MEV
        ok_i_ii = abs(C_i - C_ii) * 1e3 <= TOL_I_II_MEV
        if has_plateau and ok_i_ii:
            C_ret, tag = C_i, "plateau"
        else:
            C_ret = float(lu["shifts"][1.0])
            tag = "sans plateau" if not has_plateau else f"écart (i)-(ii) {abs(C_i - C_ii)*1e3:.1f} meV > {TOL_I_II_MEV:.0f}"
            tag += " : valeur Lu (1,0 Å)"
        out["sizes"][S] = dict(N=n, N_cells=N_cells, r_max_A=rmax, n_atoms=len(x_d), plateau_i_10=pl10, plateau_i_05=pl05, plateau_ii_periodic_pz=plP,
                               plateau_ii_periodic_sp2=plP_sp2, plateau_ii_box_pz=plB, mean3d_meV=mean3d * 1e3, Lu=chk_lu,
                               far_atom=dict(index1=int(ifar) + 1, dist_axis_A=float(r10["dist_axis"][ifar]), dist_true_A=float(r10["dist"][ifar])),
                               orbital_site_match_max_A=float(max(mP.max(), mB[inb].max())), has_plateau=bool(has_plateau), ok_i_ii=bool(ok_i_ii),
                               C_i_eV=C_i, C_ii_eV=C_ii, C_retenu_eV=C_ret, etiquette=tag, elapsed_s=time.time() - t0)
        prof[S] = dict(dist=r10["dist"], shift10=r10["shift"], shift05=r05["shift"], dP=dP, diagP=diagP, dB=dB[inb], diagB=diagB[inb], edgeB=edge[inb], rmax=rmax)
        np.savez(os.path.join(d, f"a1_profiles_{S}.npz"), **{k: v for k, v in prof[S].items()})
        log(f"[A.1] {S} : (i) 1,0 Å plateau {pl10['mean']*1e3:+.2f} meV (rms {pl10['rms']*1e3:.2f}, max|écart| {pl10['max_abs_dev']*1e3:.2f}, {pl10['n']} atomes, "
            f"d >= {PLATEAU_FRAC*rmax:.2f} Å) ; 0,5 Å {pl05['mean']*1e3:+.2f} ; (ii) p_z périodique {plP['mean']*1e3:+.2f} (rms {plP['rms']*1e3:.2f}, max {plP['max_abs_dev']*1e3:.2f}, "
            f"{plP['n']}) ; boîte {plB['mean']*1e3:+.2f} ({plB['n']}) ; (iii) {mean3d*1e3:+.2f} ; Lu {lu['shifts'][1.0]*1e3:+.2f} (P1 {r10['shift'][ifar]*1e3:+.2f}, publié "
            f"{LU_PUBLISHED_MEV.get(S)}) ; site des orbitales {out['sizes'][S]['orbital_site_match_max_A']:.1e} Å ; -> C_N {C_ret*1e3:+.2f} meV ({tag}) ; {time.time()-t0:.0f} s")
        save_json(os.path.join(d, "a1_results.json"), out)
    a1_tables(out); a1_figure(prof, out)


def a1_tables(out):
    L = ["# R9 — A.1 : décalage uniforme C_N par taille (tables générées par r9_driver.py a1)", "",
         f"Plateau (décision 2) : moyenne sur les atomes à distance vraie (image minimale) ≥ {PLATEAU_FRAC} r_max de la lacune ; « sans plateau » si max|écart| > "
         f"{PLATEAU_TOL_MEV:.0f} meV ; critère (i)/(ii) : {TOL_I_II_MEV:.0f} meV. meV sauf mention.", "",
         "| N | r_max (Å) | (i) 1,0 Å : moyenne ; rms ; max\\|écart\\| ; n | (i) 0,5 Å | (ii) p_z, M_W périodique : moyenne ; rms ; max ; n | (ii) sp², périodique | (ii) p_z, boîte (bord ≥ 1) ; n | (iii) ⟨ΔV⟩_3D | Lu 1,0 Å (P1 ; publié) | \\|(i) − (ii)\\| | C_N retenu |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for S, r in out["sizes"].items():
        a, b, c, e, f = r["plateau_i_10"], r["plateau_i_05"], r["plateau_ii_periodic_pz"], r["plateau_ii_periodic_sp2"], r["plateau_ii_box_pz"]
        L.append(f"| {S} | {r['r_max_A']:.2f} | {a['mean']*1e3:+.2f} ; {a['rms']*1e3:.2f} ; {a['max_abs_dev']*1e3:.2f} ; {a['n']} | {b['mean']*1e3:+.2f} | "
                 f"{c['mean']*1e3:+.2f} ; {c['rms']*1e3:.2f} ; {c['max_abs_dev']*1e3:.2f} ; {c['n']} | {e['mean']*1e3:+.2f} | {f['mean']*1e3:+.2f} ; {f['n']} | "
                 f"{r['mean3d_meV']:+.2f} | {r['Lu']['P1_at_far_atom_meV']:+.2f} ; {r['Lu']['published_meV']} | {abs(r['C_i_eV']-r['C_ii_eV'])*1e3:.2f} | "
                 f"{r['C_retenu_eV']*1e3:+.2f} ({r['etiquette']}) |")
    L += ["", "Atome de Lu (far_atom, image minimale axe par axe) : " + " ; ".join(f"{S} atome {r['far_atom']['index1']} à {r['far_atom']['dist_axis_A']:.2f} Å (axe par axe), "
                                                                                    f"{r['far_atom']['dist_true_A']:.2f} Å (vraie)" for S, r in out["sizes"].items())]
    with open(os.path.join(WORK, "a", "A1_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def a1_figure(prof, out):
    plt, pal = fig_style()
    sizes = list(prof)
    fig, axes = plt.subplots(2, 3, figsize=(6.5, 4.6), sharex=False, squeeze=False)
    for ax, S, lab in zip(axes.ravel(), sizes, "abcdef"):
        p = prof[S]; r = out["sizes"][S]
        ax.plot(p["dist"], p["shift10"] * 1e3, "o", ms=2, color=pal.NAVY, label=r"sphères 1,0 Å (i)")
        ax.plot(p["dP"][:, [3, 4]].ravel(), p["diagP"][:, [3, 4]].ravel() * 1e3, "s", ms=1.8, color=pal.ORANGE, mfc="none", label=r"$M_W$ $p_z$, périodique (ii)")
        ax.axhline(r["C_i_eV"] * 1e3, color=pal.NAVY, lw=0.8); ax.axhline(r["Lu"]["P1_at_far_atom_meV"], color=pal.REF, lw=0.8, ls="--", label="Lu (1,0 Å)")
        ax.axvspan(PLATEAU_FRAC * p["rmax"], p["rmax"] * 1.02, color=pal.LIGHT, alpha=0.4, lw=0)
        lo = min(np.percentile(p["shift10"] * 1e3, 5), r["C_i_eV"] * 1e3) - 40; hi = max(np.percentile(p["shift10"] * 1e3, 95), r["C_i_eV"] * 1e3) + 40
        ax.set_ylim(lo, hi); ax.set_title(f"({lab}) {S.replace('x', '×')}", fontsize=8)
        ax.set_xlabel(r"Distance à la lacune (Å)", fontsize=7); ax.set_ylabel(r"$\Delta V$ (meV)", fontsize=7)
    axes[0][0].legend(fontsize=5)
    fig.tight_layout(); savefig(fig, "offset_profiles"); plt.close(fig)


# ----------------------------------------------------------------------------------------------- A.3
def conv_intensive(cfg):
    """Convention intensive : (max − min)/moyenne de max|M| (bandes 1–16, eV), grossiers 5, 7, 8, 9 et denses six tailles, brut et aligné."""
    res = results_dir(cfg); rows = {}
    for S in ["5x5", "7x7", "8x8", "9x9"]:
        C, tag = C_of(S); n2 = size_n(S) ** 2
        M = matrix_io.load_M_checked(f"{res}/M_ed_{S}.npy", require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
        nb, nk = M.shape[:2]; Ma = M.copy()
        for kk in range(nk):
            Ma[np.arange(nb), kk, np.arange(nb), kk] -= n2 * C                        # M_aligné = M2 − N_cells C_N 𝕀 (exact, k commensurables)
        i0 = np.unravel_index(np.abs(M).argmax(), M.shape); i1 = np.unravel_index(np.abs(Ma).argmax(), Ma.shape)
        rows[f"grossier {S}"] = dict(set="grossiers 5, 7, 8, 9", S=S, brut=float(np.abs(M).max()), **{vkey("aligne"): float(np.abs(Ma).max())}, C_eV=C, C_label=tag,
                                    argmax_brut_kk=bool(i0[1] == i0[3]), argmax_aligne_kk=bool(i1[1] == i1[3]), mode_aligne="exact (N_cells C_N 𝕀)")
        del M, Ma
    for S in SIZES_A:
        C, tag = C_of(S); n2 = size_n(S) ** 2; dp = dense_paths(cfg, S)
        matrix_io.check_manifest(dp["mfile"], require_normalization=matrix_io.M_NORM_V2)
        Mm = np.load(dp["mfile"], mmap_mode="r"); M = np.array(Mm[:16, :, :16, :]) * HA2EV; del Mm
        nk = M.shape[1]
        if has_exact(S):
            Mb = np.load(os.path.join(WORK, "cache", f"ML_box_{S}.npy"), mmap_mode="r"); Ma = M - C * np.array(Mb[:16, :, :16, :]); del Mb
            mode = "exact (M2 − C_N M^L[1_boîte])"
        else:
            Ma = M.copy()
            for kk in range(nk):
                Ma[np.arange(16), kk, np.arange(16), kk] -= n2 * C                   # blocs k = k' seulement : M^L[1_boîte]_{n'k,nk} = N_cells δ_nn' exactement
            mode = "blocs k = k' seulement"
        i0 = np.unravel_index(np.abs(M).argmax(), M.shape); i1 = np.unravel_index(np.abs(Ma).argmax(), Ma.shape)
        rows[f"dense {S}"] = dict(set="denses six tailles", S=S, brut=float(np.abs(M).max()), **{vkey("aligne"): float(np.abs(Ma).max())}, C_eV=C, C_label=tag,
                                 argmax_brut_kk=bool(i0[1] == i0[3]), argmax_aligne_kk=bool(i1[1] == i1[3]), mode_aligne=mode)
        del M, Ma
    spreads = {}
    for setname, keys in (("grossiers 5, 7, 8, 9", [f"grossier {S}" for S in ["5x5", "7x7", "8x8", "9x9"]]), ("denses six tailles", [f"dense {S}" for S in SIZES_A]),
                          ("denses 5, 7, 8, 9", [f"dense {S}" for S in ["5x5", "7x7", "8x8", "9x9"]])):
        for var in ("brut", vkey("aligne")):
            v = np.array([rows[k][var] for k in keys]); spreads[f"{setname} | {var}"] = float((v.max() - v.min()) / v.mean())
    return dict(rows=rows, spreads=spreads)


def offdiag_intensive(cfg):
    """
    R.2 (clôture) : max|M| hors k = k' (bandes 1–16, eV) — grossiers 5, 7, 8, 9 (paires k ≠ k') et denses six tailles (hors blocs k = k'),
    tel quel, aligné Lu et aligné plateau (les deux C_N lus dans a1_results.json, quelle que soit l'option --cn). L'alignement ne touche les
    éléments k ≠ k' que pour les denses 5x5 et 9x9 (M2 − C_N M^L[1_boîte], exact) ; ailleurs (grossiers : N_cells C_N 𝕀 ; denses 6, 7, 8, 12 :
    blocs k = k' seulement) les trois valeurs sont égales par construction (colonne « mode »).
    """
    res = results_dir(cfg); a1 = json.load(open(os.path.join(WORK, "a", "a1_results.json")))["sizes"]; rows = {}

    def off_max(M):
        A = np.abs(M)
        for kk in range(M.shape[1]):
            A[:, kk, :, kk] = 0.0
        return float(A.max())
    for S in ["5x5", "7x7", "8x8", "9x9"]:
        M = matrix_io.load_M_checked(f"{res}/M_ed_{S}.npy", require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
        m = off_max(M); del M
        rows[f"grossier {S}"] = dict(set="grossiers 5, 7, 8, 9", S=S, brut=m, lu=m, plateau=m, mode="k ≠ k' inchangés (N_cells C_N 𝕀)")
    for S in SIZES_A:
        dp = dense_paths(cfg, S); matrix_io.check_manifest(dp["mfile"], require_normalization=matrix_io.M_NORM_V2)
        Mm = np.load(dp["mfile"], mmap_mode="r"); M = np.array(Mm[:16, :, :16, :]) * HA2EV; del Mm
        m = off_max(M)
        if has_exact(S):
            Mb = np.array(np.load(os.path.join(WORK, "cache", f"ML_box_{S}.npy"), mmap_mode="r")[:16, :, :16, :])
            lu = off_max(M - a1[S]["C_retenu_eV"] * Mb); pl = off_max(M - a1[S]["C_i_eV"] * Mb); mode = "exact (M2 − C_N M^L[1_boîte])"; del Mb
        else:
            lu = pl = m; mode = "k ≠ k' inchangés (blocs k = k' seulement)"
        rows[f"dense {S}"] = dict(set="denses six tailles", S=S, brut=m, lu=lu, plateau=pl, mode=mode); del M
    spreads = {}
    for setname, keys in (("grossiers 5, 7, 8, 9", [f"grossier {S}" for S in ["5x5", "7x7", "8x8", "9x9"]]), ("denses six tailles", [f"dense {S}" for S in SIZES_A])):
        for var in ("brut", "lu", "plateau"):
            v = np.array([rows[k][var] for k in keys]); spreads[f"{setname} | {var}"] = float((v.max() - v.min()) / v.mean())
    return dict(rows=rows, spreads=spreads)


def fine_pi(Mwr, Rn, nw, Upi, Pm, Pp, mask):
    """π/π* de M sur la grille fine avec les R, R' de `mask` (m_rcut_convergence.py l. 36–44, recopiée telle quelle)."""
    nR = len(Rn); nk = Pm.shape[0]
    Wm = Mwr.copy(); keep = np.asarray(mask, bool); Wm[:, ~keep, :, :] = 0; Wm[:, :, :, ~keep] = 0
    A = Wm.transpose(1, 0, 3, 2).reshape(nR * nw, nR * nw)
    X = A.reshape(nR * nw, nR, nw); Y = np.einsum("aRw,kR->akw", X, Pp, optimize=True).reshape(nR, nw, nk, nw)
    Y = np.einsum("Rwkv,kvn->Rwkn", Y, Upi, optimize=True)
    Z = np.einsum("Rwkn,KR->Kwkn", Y, Pm, optimize=True)
    return np.einsum("Kwkn,Kwm->mKnk", Z, Upi.conj(), optimize=True)


def rcut_table(cfg, S, rcuts=range(7), nf=60):
    """tab:rcut_M (m_rcut_convergence.py : k du XML ramenés sur la grille MP) pour M2 tel quel, aligné (i) et, si A.2 l'a, aligné exact."""
    dp = dense_paths(cfg, S); n = size_n(S); C, tag = C_of(S)
    M = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV, require_normalization=matrix_io.M_NORM_V2)
    k = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k); k = np.round(k * np.asarray(MP)) / np.asarray(MP)
    U, Ud, Hwr, Rw, nd = load_wannier(dp, k)
    Mwr, R = Mwk_to_Mwr(Mbk_to_Mwk(M, U, Ud), k, MP); del M
    Rn, Rd = lt.recenter_mwr(Mwr, R, MP); lt.mwr_locality(Mwr, Rn)
    inb, _ = box_geometry(R, MP, n); nw, nR = Mwr.shape[0], Mwr.shape[1]
    variants = {"brut": Mwr}
    Ma = Mwr.copy()
    for i in np.where(inb)[0]:
        Ma[np.arange(nw), i, np.arange(nw), i] -= C
    variants[vkey("aligne")] = Ma
    if has_exact(S):
        Z = np.load(os.path.join(WORK, "cache", f"Fw_{S}.npz")); assert np.array_equal(Z["R"], R)
        variants[vkey("exact")] = Mwr - C * Z["Fw"]
    kf = lt.mp_grid(nf, nf, 1); nk = len(kf)
    _, Ef, Uf = lt.Hwr_to_Hwk(Hwr, Rw, kf, ndegen=nd); Upi = np.ascontiguousarray(Uf[:, :, 3:5])
    Pm = np.exp(-2j * np.pi * (kf @ Rn.T)); Pp = np.exp(+2j * np.pi * (kf @ Rn.T)); dist = np.linalg.norm(Rn, axis=1)
    out = dict(S=S, C_eV=C, C_label=tag, R_d=Rd.tolist(), rows={})
    for var, Mw in variants.items():
        ref = fine_pi(Mw, Rn, nw, Upi, Pm, Pp, np.ones(nR, bool)); mref = float(np.abs(ref).max())
        rr = []
        for rc in rcuts:
            mask = dist <= rc + 1e-9; Mt = fine_pi(Mw, Rn, nw, Upi, Pm, Pp, mask)
            rr.append(dict(R_cut=int(rc), n_sites=int(mask.sum()), max_dM_over_maxM=float(np.abs(Mt - ref).max() / mref), maxM_pi_eV=mref))
        out["rows"][var] = rr
        log(f"[A.3 rcut] {S} {var} : " + " ; ".join(f"R_cut {r['R_cut']} {r['max_dM_over_maxM']:.3e}" for r in rr) + f" ; max|M_π| {mref:.3f} eV")
    return out


def cmd_a3(a):
    cfg = load_production(verbose=True); d = ensure("a"); f = os.path.join(d, f"a3_results{cn_suffix()}.json")
    out = json.load(open(f)) if os.path.exists(f) else {}
    out["head"] = git_head()
    # C14 de R6 (lecture seule) : énergie de l'écart relatif maximal de Γ_T
    z = np.load(os.path.join(results_dir(cfg), "resonance_9x9_shiftL.npz")); E_D = float(z["E_D"]); eg = z["eg"]; m = np.abs(eg - E_D) <= cfg["e_window_eV"]
    c14 = {}
    for key, C in (("Gamma_T_shiftp25", +25), ("Gamma_T_shiftm25", -25)):
        rel = np.abs(z[key][m] - z["Gamma_T"][m]) / z["Gamma_T"][m]; i = int(np.argmax(rel))
        c14[f"{C:+d} meV"] = dict(max_rel=float(rel[i]), at_eV=float(eg[m][i] - E_D), median_rel=float(np.median(rel)))
    out["C14_R6"] = c14; log(f"[A.3 C14] {c14}")
    save_json(f, out)
    if "conv_intensive" not in out or a.redo:
        out["conv_intensive"] = conv_intensive(cfg); save_json(f, out)
        log(f"[A.3 conv. intensive] {out['conv_intensive']['spreads']}")
    if CN is not None and ("offdiag_intensive" not in out or a.redo):             # R.2 (clôture) : test intensif hors k = k'
        out["offdiag_intensive"] = offdiag_intensive(cfg); save_json(f, out)
        log(f"[R.2 hors k = k'] {out['offdiag_intensive']['spreads']}")
    for S in ("9x9", "12x12"):
        if f"rcut_{S}" not in out or a.redo:
            out[f"rcut_{S}"] = rcut_table(cfg, S); save_json(f, out)
    # niveau 1 et courbe Γ_T, six tailles, R_cut 3 et 4
    out.setdefault("level1", {})
    for S in a.sizes.split(","):
        if S in out["level1"] and not a.redo:
            continue
        st = resonance_setup(S, cfg, (3, 4)); C, tag = C_of(S)
        z = np.load(os.path.join(results_dir(cfg), f"specwd_{S}_prod.npz"))["results"]
        rows_c, xr_c, rank_c = crown_table(st, lo=-1.5, hi=0.0)                  # couronnes π de la grille 240² (E_res avec sa couronne)
        res = dict(C_eV=C, C_label=tag, E_D=st["E_D"], R_d=st["R_d"].tolist(), rc={}, crowns_pi_240=rows_c)
        for rc in (3, 4):
            g0s = g0_cached(st, cfg["nk_int"], rc, "sig"); g0r = g0_cached(st, cfg["nk_int"], rc, "res")
            res["rc"][str(rc)] = dict(n_out_of_box=int((~st["loc"][rc]["in_box"]).sum()), variants={})
            for var in (["brut"] if CN is None else []) + ["aligne"] + (["exact"] if has_exact(S) else []):   # --cn : tel quel déjà dans R9
                V = V_variant(st, rc, var, 0.0 if var == "brut" else C); t0 = time.time()
                s = level1_stats(st, sigma_states(st, V, g0s, rc))                  # grille « sig » : niveau 1 et Re Σ (R6 3.3)
                s["E_res_crown"] = crown_of_state(st, s["E_res_state"][0], s["E_res_state"][1], xr_c, rank_c)
                o = resonance_observables(st, V, g0r, rc, born=(S == "9x9" and rc == 3))   # grille « res » : courbe Γ_T, T̄(K) (R6 3.4)
                for kk_, vv_ in o.items():
                    if kk_ not in ("G", "curve_eg5", "curve_eg25", "tr", "E_res_state"):
                        s["E_res_resgrid" if kk_ == "E_res" else kk_] = vv_
                s["elapsed_s"] = time.time() - t0
                if var == "brut":
                    row = [r for r in z if r[0] == rc and r[1] == cfg["grid"] and abs(r[2] - cfg["eta_eV"]) < 1e-12]
                    if row:
                        s["regression_specwd"] = dict(median_meV=float(row[0][3]), E_res=float(row[0][4]),
                                                      ok=bool(abs(s["median_meV"] - row[0][3]) < 1e-6 * row[0][3] and abs(s["E_res"] - row[0][4]) < 1e-9))
                res["rc"][str(rc)]["variants"][vkey(var)] = s
                log(f"[A.3 niveau 1] {S} R_cut {rc} {var} : médiane {s['median_meV']:.2f} meV, E_res {s['E_res']:+.3f}, Re Σ médian {s['median_ReSigma_meV']:.2f} meV, "
                    f"pic Γ_T {s['peak_GT_prod']:+.3f} (fin {s['peak_GT_fine']:+.4f}), pic −Im T̄(K) {s['peak_ImTbar']:+.4f}, Re T̄(E_D) {s['ReTbar_at_ED']:+.3f}"
                    + (f", Born/T médian {s['Born_over_T_median']:.2f}" if "Born_over_T_median" in s else "")
                    + (f" ; régression specwd {'OK' if s['regression_specwd']['ok'] else 'ÉCART'}" if "regression_specwd" in s else "") + f" ; {s['elapsed_s']:.0f} s")
            del g0s, g0r
        out["level1"][S] = res; save_json(f, out)
    # familles
    fam = {}
    for rc in ("3", "4"):
        for var in sorted({v for S in out["level1"] for v in out["level1"][S]["rc"][rc]["variants"]}):
            for fname, members in FAMILIES.items():
                v = [out["level1"][S]["rc"][rc]["variants"][var]["median_meV"] for S in members if S in out["level1"] and var in out["level1"][S]["rc"][rc]["variants"]]
                if len(v) == len(members):
                    v = np.array(v); fam[f"R_cut {rc} | {var} | {fname}"] = dict(mean=float(v.mean()), spread=float((v.max() - v.min()) / v.mean()))
    out["families"] = fam; save_json(f, out)
    a3_tables(out)


def cmd_a3pole(a):
    """Critère de pôle (min|det|, min|λ|), 9x9, 300², blocs π et complet, trois variantes ; 1 fil BLAS + fils Python (R4_EIG_WORKERS)."""
    from electron_defect_interaction.defects.many_body import pole_criterion as pc
    cfg = load_production(verbose=True); d = ensure("a"); f = os.path.join(d, f"a3_results{cn_suffix()}.json")
    out = json.load(open(f)) if os.path.exists(f) else {}
    st = resonance_setup("9x9", cfg, (3,)); C, tag = C_of("9x9"); g0 = g0_cached(st, cfg["nk_int"], 3, "res")
    eg = st["egrid"] - st["E_D"]; w = np.abs(eg) <= cfg["e_window_eV"]; nL = len(st["loc"][3]["Rloc"])
    blocks = {"pi": pc.block_indices(nL, NW, WF_PI), "complet": None}
    res = {}
    for var in (["brut"] if CN is None else []) + ["aligne"] + (["exact"] if has_exact("9x9") else []):
        V = V_variant(st, 3, var, 0.0 if var == "brut" else C); var = vkey(var); res[var] = {}
        for bname, idx in blocks.items():
            t0 = time.time(); r = pc.det_eig_criterion(V, g0[w], idx=idx, vectors=False)
            rel = np.exp(r["logdet_rel"]); i = int(np.argmin(rel)); j = int(np.argmin(r["minlam"]))
            res[var][bname] = dict(min_absdet_rel=float(rel[i]), at_eV=float(eg[w][i]), min_abs_lambda=float(r["minlam"][j]), lambda_at_eV=float(eg[w][j]),
                                   lambda_min=[float(r["lam_min"][j].real), float(r["lam_min"][j].imag)], elapsed_s=time.time() - t0)
            log(f"[A.3 pôle] 9x9 {var} bloc {bname} : min|det|/max {rel[i]:.3e} à {eg[w][i]:+.3f} ; min|λ| {r['minlam'][j]:.4f} à {eg[w][j]:+.3f} "
                f"(λ {r['lam_min'][j].real:+.4f}{r['lam_min'][j].imag:+.4f}i) ; {time.time()-t0:.0f} s")
    out["pole_9x9"] = dict(variants=res, C_eV=C, C_label=tag, R6=dict(pi=dict(min_absdet_rel=2.33e-2, at_eV=-0.172, min_abs_lambda=0.3967, lambda_at_eV=-0.170),
                                                                          complet=dict(min_absdet_rel=1.32e-4, at_eV=-0.812, min_abs_lambda=0.0019, lambda_at_eV=-0.812)))
    save_json(f, out); a3_tables(out)


def a3_tables(out):
    L = ["# R9 — A.3 : variantes alignées (tables générées par r9_driver.py a3)", ""]
    if "C14_R6" in out:
        L += ["C14 de R6 (`resonance_9x9_shiftL.npz`, lecture) : " + " ; ".join(f"C = {k} : écart relatif maximal de Γ_T {v['max_rel']:.3e} à {v['at_eV']:+.3f} eV (médian {v['median_rel']:.2e})"
                                                                           for k, v in out["C14_R6"].items()), ""]
    if "conv_intensive" in out:
        ci = out["conv_intensive"]
        L += ["## Convention intensive (max|M|, bandes 1–16, eV)", "", "| ensemble | taille | brut | aligné | mode aligné | C_N (meV) | max en k = k' (brut ; aligné) |", "|---|---|---|---|---|---|---|"]
        for k, r in ci["rows"].items():
            ak = next(x for x in r if x.startswith("aligne"))
            L.append(f"| {r['set']} | {r['S']} | {r['brut']:.3f} | {r[ak]:.3f} | {r['mode_aligne']} | {r['C_eV']*1e3:+.2f} ({r['C_label']}) | {r['argmax_brut_kk']} ; {r['argmax_aligne_kk']} |")
        L += ["", "| ensemble | variante | (max − min)/moyenne |", "|---|---|---|"] + [f"| {k.split(' | ')[0]} | {k.split(' | ')[1]} | {v:.3e} |" for k, v in ci["spreads"].items()] + [""]
    for S in ("9x9", "12x12"):
        if f"rcut_{S}" in out:
            r = out[f"rcut_{S}"]; vars_ = list(r["rows"])
            L += [f"## tab:rcut_M, {S} (max|ΔM|/max|M| de la paire π, grille fine 60² ; C_N {r['C_eV']*1e3:+.2f} meV, {r['C_label']})", "",
                  "| R_cut | sites | " + " | ".join(VAR_LABEL[v] for v in vars_) + " |", "|---|---|" + "---|" * len(vars_)]
            for i in range(len(r["rows"]["brut"])):
                L.append(f"| {r['rows']['brut'][i]['R_cut']} | {r['rows']['brut'][i]['n_sites']} | " + " | ".join(f"{r['rows'][v][i]['max_dM_over_maxM']:.3e}" for v in vars_) + " |")
            L.append("")
    if "level1" in out:
        L += ["## Niveau 1 (240², η 0,02, N_k^int 300) et courbe Γ_T (grille de resonance_metrics)", "",
              "| taille | R_cut | variante | médiane Γ·N_cells (meV) | E_res | Re Σ médian (meV) | pic Γ_T (5 meV ; 2,5 meV) | pic −Im T̄(K) | Re T̄(E_D) | sites hors boîte | C_N (meV) | régression R6 |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for S, r in out["level1"].items():
            for rc, rr in r["rc"].items():
                for var, s in rr["variants"].items():
                    reg = "—" if "regression_specwd" not in s else ("OK" if s["regression_specwd"]["ok"] else "ÉCART")
                    L.append(f"| {S} | {rc} | {VAR_LABEL[var]} | {s['median_meV']:.2f} | {s['E_res']:+.3f} | {s['median_ReSigma_meV']:.2f} | {s['peak_GT_prod']:+.3f} ; {s['peak_GT_fine']:+.4f} | "
                             f"{s['peak_ImTbar']:+.4f} | {s['ReTbar_at_ED']:+.3f} | {rr['n_out_of_box']} | {0 if var == 'brut' else r['C_eV']*1e3:+.2f} ({'—' if var == 'brut' else r['C_label']}) | {reg} |")
        L.append("")
    if "families" in out:
        L += ["Familles (médianes Γ·N_cells) : " + " ; ".join(f"{k} : moyenne {v['mean']:.1f} meV, (max − min)/moyenne {v['spread']*100:.2f} %" for k, v in out["families"].items()), ""]
    if "pole_9x9" in out:
        p = out["pole_9x9"]
        L += ["## Critère de pôle, 9×9, 300², fenêtre ±3 eV", "", "| variante | bloc | min\\|det\\|/max (ε − E_D) | min\\|λ\\| (ε − E_D) ; λ |", "|---|---|---|---|"]
        for var, bl in p["variants"].items():
            for b, r in bl.items():
                L.append(f"| {VAR_LABEL[var]} | {b} | {r['min_absdet_rel']:.3e} ({r['at_eV']:+.3f}) | {r['min_abs_lambda']:.4f} ({r['lambda_at_eV']:+.3f}) ; {r['lambda_min'][0]:+.4f}{r['lambda_min'][1]:+.4f}i |")
        L += [f"R6 (M2 tel quel) : π 2,33e-2 (−0,172), |λ| 0,397 (−0,170) ; complet 1,32e-4 (−0,812), |λ| 0,0019 (−0,812).", ""]
    with open(os.path.join(WORK, "a", "A3_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------------- C : chaîne repliée
C1_SIZES = [6, 9, 12, 15, 18, 21, 24, 27]
C2_SIZES = [27, 36, 45, 54, 63, 72, 81]
C2_MAX_S = 3600.0             # une diagonalisation prévue au-delà d'une heure (extrapolation n³) n'est pas lancée
WIN_C = (-3.0, 1.0)


def fold_block(Hwr, Rw, nd, V, Rloc, N, wfs):
    """H_SC = fold(H(R)) + fold(V_loc) sur les orbitales wfs (tout le jeu : 5 ; bloc π : [3, 4]) ; disposition c*nw' + w'."""
    wfs = list(wfs); nwb = len(wfs); nL = len(Rloc)
    H0 = sf.fold_hwr_to_supercell(Hwr[:, wfs][:, :, wfs], Rw, nd, N)
    if V is None:
        return H0
    idx = pc_block(nL, wfs)
    Vb = V[np.ix_(idx, idx)].reshape(nL, nwb, nL, nwb).transpose(1, 0, 3, 2)
    return H0 + sf.fold_mwr_to_supercell(Vb, Rloc, N)


def pc_block(nL, wfs):
    return np.array([L * NW + w for L in range(nL) for w in wfs], dtype=int)


def nn_index(N, wfs):
    """Indices, dans le modèle replié sur `wfs`, des p_z(B) des trois voisins de la lacune (mailles NN_CELLS, recentrées en 0)."""
    wfs = list(wfs); return [int(sf.cell_of(np.array(R), N)) * len(wfs) + wfs.index(WF_PZ_B) for R in NN_CELLS]


def site_index(N, wfs):
    """Indices du « poids site + voisins » de R6 : p_z(A) et sp² de la maille 0, p_z(B) des trois voisins."""
    wfs = list(wfs); c0 = int(sf.cell_of(np.array([0, 0, 0]), N)) * len(wfs)
    return [c0 + wfs.index(w) for w in wfs if w in (0, 1, 2, WF_PZ_A)] + nn_index(N, wfs)


def c_states(e, v, idx_nn, idx_site, E_D, thr):
    """États de la fenêtre : ε − E_D, poids sur les p_z voisins, poids site + voisins ; localisés (poids voisins > seuil) ; quasi-lié = poids max."""
    x = e - E_D; sel = np.where((x >= WIN_C[0]) & (x <= WIN_C[1]))[0]
    wnn = (np.abs(v[idx_nn][:, sel]) ** 2).sum(0); wsite = (np.abs(v[idx_site][:, sel]) ** 2).sum(0)
    if len(sel) == 0:
        return dict(n_window=0)
    j = int(np.argmax(wnn))
    return dict(n_window=int(len(sel)), quasi_bound=dict(x=float(x[sel][j]), w_nn=float(wnn[j]), w_site=float(wsite[j])),
                localized=[dict(x=float(x[sel][i]), w_nn=float(wnn[i]), w_site=float(wsite[i])) for i in np.where(wnn > thr)[0]] if thr is not None else [],
                n_localized=int((wnn > thr).sum()) if thr is not None else None)


def ldos_max(e, v, idx, E_D, eta=0.02, lo=-1.0, hi=0.0):
    eg = E_D + np.arange(lo - 0.2, hi + 0.2 + 1e-9, 0.0025); rho = sf.ldos_from_eigenpairs(e, v, idx, eg, eta).sum(1)
    m = (eg - E_D >= lo) & (eg - E_D <= hi); i = int(np.argmax(rho[m]))
    return dict(at_eV=float(eg[m][i] - E_D), value=float(rho[m][i])), eg - E_D, rho


def cmd_c(a):
    cfg = load_production(verbose=True); d = ensure("c"); f = os.path.join(d, "c_results.json")
    out = json.load(open(f)) if os.path.exists(f) else {}
    out["head"] = git_head()
    st = resonance_setup("9x9", cfg, (3,)); L = st["loc"][3]; Rloc = L["Rloc"]; C9, tag9 = C_of("9x9")
    Vs = {"brut": L["V"], "aligne": V_variant(st, 3, "aligne", C9)}
    Hwr, Rw, nd = st["Hwr"], st["Rw"], st["nd"]
    _, EK, _ = lt.Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); E_DW = float(0.5 * (EK[0, 3] + EK[0, 4]))
    out.update(C9_eV=C9, C9_label=tag9, E_D_wannier_K=E_DW)
    ALL = list(range(NW)); PI = list(WF_PI); parW = lambda N: np.array([1.0 if w in WF_SIGMA else -1.0 for c in range(N * N) for w in range(NW)])
    # ---- C.0 porte : N = 9, V_loc non aligné, (c-3) de R6 : π localisé −0,677 (w₂ 0,227)
    if "C0" not in out:
        out["C0"] = c0_gate(st, Vs["brut"], E_DW); save_json(f, out)
        if not out["C0"]["ok"]:
            log("[C.0] ÉCHEC de la porte : STOP sur la partie C"); sys.exit(3)
    # ---- C.1 : N = 6 … 27, matrice complète
    out.setdefault("C1", {})
    for N in C1_SIZES:
        if str(N) in out["C1"] and not a.redo:
            continue
        t0 = time.time(); r = dict(N=N)
        H0 = fold_block(Hwr, Rw, nd, None, Rloc, N, ALL); e0, v0 = np.linalg.eigh(0.5 * (H0 + H0.conj().T))
        q = np.argsort(np.abs(e0 - E_DW))[:4]; E_Dq = float(e0[q].mean())
        dev1, _, ok1 = sf.wannier_gate(Hwr, Rw, nd, N)
        idx_nn = nn_index(N, ALL); idx_site = site_index(N, ALL)
        # seuil : 3 × poids moyen des états impairs de la fenêtre de la parfaite repliée sur les p_z voisins
        par = parW(N); io = np.where(par < 0)[0]
        x0 = e0 - E_Dq; w0 = (np.abs(v0[idx_nn]) ** 2).sum(0); odd0 = (np.abs(v0[io]) ** 2).sum(0) > 0.5
        msk0 = (x0 >= WIN_C[0]) & (x0 <= WIN_C[1]) & odd0
        thr = 3.0 * float(w0[msk0].mean())
        r.update(E_D_quadruplet=E_Dq, quadruplet_spread=float(np.ptp(e0[q])), E_D_minus_wannier=E_Dq - E_DW, gate1_dev=float(dev1), gate1_ok=bool(ok1), thr_nn=thr,
                 first_crown=first_crown(Hwr, Rw, nd, N, E_DW))
        qe = qe_r7(N)
        r["QE"] = qe
        for var, V in Vs.items():
            H = fold_block(Hwr, Rw, nd, V, Rloc, N, ALL); H = 0.5 * (H + H.conj().T)
            ie = np.where(par > 0)[0]
            cpl = float(np.abs(H[np.ix_(ie, io)]).max())
            eo, vo = np.linalg.eigh(H[np.ix_(io, io)]); ee, ve = np.linalg.eigh(H[np.ix_(ie, ie)])
            pos_o = {g: i for i, g in enumerate(io)}; pos_e = {g: i for i, g in enumerate(ie)}
            so = c_states(eo, vo, [pos_o[i] for i in idx_nn], [pos_o[i] for i in idx_site if i in pos_o], E_Dq, thr)
            se_site = [pos_e[i] for i in idx_site if i in pos_e]
            xe = ee - E_Dq; sel_e = np.where((xe >= WIN_C[0]) & (xe <= WIN_C[1]))[0]
            even = [dict(x=float(xe[i]), w_site=float((np.abs(ve[se_site, i]) ** 2).sum())) for i in sel_e]
            lm, egc, rho = ldos_max(eo, vo, [pos_o[i] for i in idx_nn], E_Dq)
            eqe = qe_r7_spectrum(N)
            s, resid, nst = al.rigid_shift_fit(eqe, np.sort(np.concatenate([eo, ee])) - E_Dq, -4.0) if eqe is not None else (np.nan, np.nan, 0)
            r[var] = dict(coupling_even_odd=cpl, odd=so, even_window=even, ldos_max=lm, rigid_shift_meV=float(s) * 1e3, rigid_resid_eV=float(resid), rigid_n=int(nst))
            np.savez(os.path.join(d, f"c1_ldos_N{N}_{var}.npz"), eg=egc, rho=rho)
            qb = so.get("quasi_bound", {})
            log(f"[C.1] N={N} {var} : π quasi-lié {qb.get('x', np.nan):+.4f} (poids voisins {qb.get('w_nn', np.nan):.3f}), {so.get('n_localized')} localisés (seuil {thr:.4f}) ; "
                f"pairs fenêtre {[round(e_['x'], 3) for e_ in even]} ; max LDOS voisins {lm['at_eV']:+.4f} ; décalage rigide {s*1e3:+.1f} meV ; couplage σ–π {cpl:.1e}")
        r["elapsed_s"] = time.time() - t0; out["C1"][str(N)] = r; save_json(f, out)
    # ---- C.2 : bloc π, N = 27 (contrôle contre le complet) … 81
    out.setdefault("C2", {}); t_prev = None; n_prev = None
    for N in C2_SIZES:
        if str(N) in out["C2"] and not a.redo:
            continue
        n_dim = 2 * N * N
        if t_prev is not None:
            pred = t_prev * (n_dim / n_prev) ** 3
            if pred > C2_MAX_S:
                out["C2"][str(N)] = dict(skipped=True, predicted_s=pred); save_json(f, out)
                log(f"[C.2] N={N} : diagonalisation prévue {pred/60:.0f} min > {C2_MAX_S/60:.0f} min, non lancée"); continue
        r = dict(N=N, dim=n_dim, first_crown=first_crown(Hwr, Rw, nd, N, E_DW))
        idx_nn = nn_index(N, PI)
        for var, V in Vs.items():
            t0 = time.time(); H = fold_block(Hwr, Rw, nd, V, Rloc, N, PI); tb = time.time() - t0
            t1 = time.time(); e, v = np.linalg.eigh(0.5 * (H + H.conj().T)); te = time.time() - t1; del H
            so = c_states(e, v, idx_nn, idx_nn, E_DW, None)
            lm, egc, rho = ldos_max(e, v, idx_nn, E_DW)
            r[var] = dict(odd=so, ldos_max=lm, t_build_s=tb, t_eigh_s=te)
            np.savez(os.path.join(d, f"c2_ldos_N{N}_{var}.npz"), eg=egc, rho=rho)
            t_prev, n_prev = te, n_dim
            log(f"[C.2] N={N} {var} (bloc π, dim {n_dim}) : π quasi-lié {so['quasi_bound']['x']:+.4f} (poids voisins {so['quasi_bound']['w_nn']:.4f}) ; max LDOS "
                f"{lm['at_eV']:+.4f} ; première couronne π {r['first_crown']['pi']:+.4f} ; eigh {te:.0f} s")
            del e, v
        out["C2"][str(N)] = r; save_json(f, out)
    # contrôle N = 27 : bloc π contre matrice complète
    if "27" in out["C1"] and "27" in out["C2"] and not out["C2"]["27"].get("skipped"):
        out["C2_vs_C1_27"] = {var: dict(dx_quasi_bound=out["C2"]["27"][var]["odd"]["quasi_bound"]["x"] - out["C1"]["27"][var]["odd"]["quasi_bound"]["x"],
                                        dldos=out["C2"]["27"][var]["ldos_max"]["at_eV"] - out["C1"]["27"][var]["ldos_max"]["at_eV"]) for var in Vs}
    # ajustements 1/N et 1/N² sur 6 … 27 (comme R7), extrapolation aux N de C.2
    fits = {}
    Ns = np.array([N for N in C1_SIZES if str(N) in out["C1"]], float)
    for var in list(Vs) + ["QE"]:
        for q_ in ("quasi_bound", "ldos"):
            if var == "QE":
                if q_ == "ldos":
                    continue
                y = np.array([out["C1"][str(int(N))]["QE"]["pi_quasi_bound"]["x"] for N in Ns])
            else:
                y = np.array([out["C1"][str(int(N))][var]["odd"]["quasi_bound"]["x"] if q_ == "quasi_bound" else out["C1"][str(int(N))][var]["ldos_max"]["at_eV"] for N in Ns])
            for p in (1, 2):
                A_ = np.vstack([np.ones_like(Ns), Ns ** -p]).T; coef, *_ = np.linalg.lstsq(A_, y, rcond=None); resid = y - A_ @ coef
                ext = {str(N): float(coef[0] + coef[1] * N ** -p) for N in C2_SIZES}
                fits[f"{var} | {q_} | 1/N^{p}"] = dict(eps_inf=float(coef[0]), a=float(coef[1]), rms=float(np.sqrt((resid ** 2).mean())), max_resid=float(np.abs(resid).max()), extrapolated=ext)
    out["fits"] = fits; save_json(f, out)
    # ---- C.3 : limite diluée (matrice T), même V_loc, g₀ 600² (et 900²)
    out.setdefault("C3", {})
    idx_cl = [int(np.where((Rloc == np.array(R)).all(axis=1))[0][0]) * NW + WF_PZ_B for R in NN_CELLS]
    for nk in (600, 900):
        g0 = g0_cached(st, nk); eg = st["egrid"] - st["E_D"]; m = (eg >= -1.0) & (eg <= 0.0); r = {}
        for var, V in Vs.items():
            rho, rho0 = lt.cluster_ldos(g0, V, idx_cl); tot = rho.sum(1); i = int(np.argmax(tot[m]))
            r[var] = dict(at_eV=float(eg[m][i]), value=float(tot[m][i]), per_orbital_at_max=[float(x) for x in rho[m][i]])
            np.savez(os.path.join(d, f"c3_ldos_nk{nk}_{var}.npz"), eg=eg, rho=rho, rho0=rho0)
            log(f"[C.3] {nk}² {var} : max de la LDOS des trois voisins (matrice T) à {r[var]['at_eV']:+.4f} eV ({r[var]['value']:.3f} états/eV)")
        out["C3"][str(nk)] = r; del g0
    save_json(f, out); c_tables(out); c_figure(out)


def cmd_cfig(a):
    """Tables et figure de C depuis c/c_results.json (sans recalcul) + ajustements supplémentaires sur 9 … 27 (N = 6 : maximum de LDOS au bord de la fenêtre)."""
    f = os.path.join(WORK, "c", "c_results.json"); out = json.load(open(f)); sup = {}
    Ns = np.array([N for N in C1_SIZES if N >= 9], float)
    for var in ("brut", "aligne", "QE"):
        for q_ in ("quasi_bound", "ldos"):
            if var == "QE" and q_ == "ldos":
                continue
            if var == "QE":
                y = np.array([out["C1"][str(int(N))]["QE"]["pi_quasi_bound"]["x"] for N in Ns])
            else:
                y = np.array([out["C1"][str(int(N))][var]["odd"]["quasi_bound"]["x"] if q_ == "quasi_bound" else out["C1"][str(int(N))][var]["ldos_max"]["at_eV"] for N in Ns])
            for p in (1, 2):
                A_ = np.vstack([np.ones_like(Ns), Ns ** -p]).T; coef, *_ = np.linalg.lstsq(A_, y, rcond=None); resid = y - A_ @ coef
                sup[f"{var} | {q_} | 1/N^{p}"] = dict(eps_inf=float(coef[0]), a=float(coef[1]), rms=float(np.sqrt((resid ** 2).mean())), max_resid=float(np.abs(resid).max()),
                                                    extrapolated={str(N): float(coef[0] + coef[1] * N ** -p) for N in C2_SIZES})
    out["fits_9_27"] = sup; save_json(f, out); c_tables(out); c_figure(out)


def first_crown(Hwr, Rw, nd, N, E_D):
    """Énergies (ε − E_D) de la première couronne de la grille N×N autour de K (K + b₁/N), bandes π et π*."""
    _, E, _ = lt.Hwr_to_Hwk(Hwr, Rw, (K_RED + np.array([1.0 / N, 0, 0]))[None], ndegen=nd)
    return dict(pi=float(E[0, 3] - E_D), pistar=float(E[0, 4] - E_D))


def qe_r7(N):
    """QE (R7 D1, fenêtre [−3, +1], alignement Lu 1,0 Å + quadruplet de la parfaite) : π quasi-lié (impair de w₂ max), états pairs, localisés."""
    fz = os.path.join(GQ, "defects", "R7_tailles_3m", "d1_8pts", f"window_{N}x{N}.npz")
    if not os.path.exists(fz):
        return None
    z = np.load(fz); x = z["D_x"]; par = z["D_parity"]; w2 = z["D_w2"]; thr = float(z["thr_w2"])
    odd = np.where(par < 0)[0]; ev = np.where(par > 0)[0]; j = odd[int(np.argmax(w2[odd]))]
    return dict(pi_quasi_bound=dict(x=float(x[j]), w2=float(w2[j])), even=[dict(x=float(x[i]), w2=float(w2[i])) for i in ev], thr_w2=thr,
                n_localized_odd=int((w2[odd] > thr).sum()), n_localized_even=int((w2[ev] > thr).sum()), E_D=float(z["E_D"]), Lu_1A=float(z["shift_Lu_1A"]))


def qe_r7_spectrum(N):
    fz = os.path.join(GQ, "defects", "R7_tailles_3m", "d1_8pts", f"window_{N}x{N}.npz")
    if not os.path.exists(fz):
        return None
    z = np.load(fz); return np.asarray(z["e_d_all"]) - float(z["shift_Lu_1A"]) - float(z["E_D"])


def c0_gate(st, V, E_DW):
    """Porte C.0 : (c-3) de R6 2.1 — H(R) replié + V_loc(R_cut 3) replié, N = 9 ; w₂ (disque de 2 Å) par la base dense (route de r6_d4.py)."""
    R4DIR = os.path.join(GQ, "defects", "R4_quasi_lie"); R5DIR = os.path.join(GQ, "defects", "R5_base_vs_M")
    sys.path.insert(0, R4DIR); sys.path.insert(0, R5DIR)
    import r5_driver as r5                                                    # noqa: E402  (aides de R5 et R4 ; journaux redirigés ici)
    r5.log = log; r5.r4.log = log; r4 = r5.r4
    os.chdir(PROJ)
    from electron_defect_interaction.io import qe_gamma_io as qg
    N = 9; nw = NW; Rloc = st["loc"][3]["Rloc"]
    H = fold_block(st["Hwr"], st["Rw"], st["nd"], V, Rloc, N, range(NW))
    parW = np.array([1.0 if w in WF_SIGMA else -1.0 for c in range(N * N) for w in range(NW)])
    bl, cpl, herm = r4.spectrum_blocks(H, parW, E_DW)
    B = r5.load_bases(with_dense=True); b20 = B["bases"][20]; idx81 = B["idx81"]; k81 = B["k81"]; ng = B["ng"]
    z1 = np.load(os.path.join(R4DIR, "d1", "d1_states.npz")); s_vac = z1["s_vac"]; thr = float(z1["thr"])
    mask2 = qg.inplane_disc_mask(ng, B["A_A"], s_vac, r4.R_W2)
    V81 = np.einsum("kbw,kwv->kbv", st["Ud"][idx81], st["U"][idx81]); phase_Rd = np.exp(-2j * np.pi * (k81 @ np.asarray(st["R_d"], float)))
    locs = {}
    for lab in ("even", "odd"):
        x = bl[lab]["x"]; sel = np.where((x >= WIN_C[0]) & (x <= WIN_C[1]))[0]; loc = []
        for j in sel:
            cr = np.zeros(N * N * nw, complex); cr[bl[lab]["idx"]] = bl[lab]["v"][:, j]
            ck = sf.rvec_to_kvec(cr, k81, N, nw).reshape(len(k81), nw) * phase_Rd[:, None]
            dnk = np.einsum("kbw,kw->bk", V81, ck)
            w2 = float(sf.folded_density_2d(dnk, b20["C"], b20["nG"], b20["flat"], ng, workers=int(os.environ.get("OMP_NUM_THREADS", "8")))[mask2].sum())
            if w2 > thr:
                loc.append((float(x[j]), w2))
        locs[lab] = loc
    ref = json.load(open(os.path.join(GQ, "defects", "R6_production_corrigee", "d4", "d4_results.json")))["variants"]["c_3"]["blocks"]
    pi = min(locs["odd"], key=lambda t: abs(t[0] + 0.677)) if locs["odd"] else (np.nan, np.nan)
    ok = abs(pi[0] - (-0.677)) <= 1e-3 and abs(pi[1] - 0.227) <= 1e-3
    out = dict(pi=dict(x=pi[0], w2=pi[1]), expected=dict(x=-0.677, w2=0.227), localized=locs, R6_c3=dict(odd=ref["odd"]["localized"], even=ref["even"]["localized"]),
               coupling_even_odd=cpl, ok=bool(ok))
    log(f"[C.0] porte : π localisé {pi[0]:+.4f} (w₂ {pi[1]:.4f}) contre (c-3) de R6 −0,677 (0,227) -> {'OK' if ok else 'ÉCHEC'} ; localisés {locs} ; R6 {out['R6_c3']}")
    return out


def c_tables(out):
    L = ["# R9 — C : chaîne repliée contre R7 (tables générées par r9_driver.py c)", "",
         f"V_loc : 9×9 de production (R_cut 3, 29 sites) ; aligné = V_loc − C_9·P_boîte, C_9 = {out['C9_eV']*1e3:+.2f} meV ({out['C9_label']}). "
         "Références d'énergie : QE R7 = ε − E_D(quadruplet de la parfaite QE) après alignement Lu 1,0 Å ; modèle = ε − E_D(quadruplet de la parfaite repliée).", ""]
    if "C0" in out:
        c = out["C0"]; L += [f"Porte C.0 (N = 9, non aligné) : π localisé {c['pi']['x']:+.4f} eV, w₂ {c['pi']['w2']:.4f} (R6 (c-3) −0,677, 0,227) → {'OK' if c['ok'] else 'ÉCHEC'}", ""]
    L += ["## C.1 (matrice complète)", "",
          "| N | 1ʳᵉ couronne π (eV) | QE : π quasi-lié ε ; w₂ | QE : pairs (ε ; w₂) | modèle non aligné : π quasi-lié ε ; poids voisins ; localisés | non aligné : pairs de la fenêtre (ε ; poids site) | non aligné : max LDOS voisins | aligné : π quasi-lié ε ; poids ; localisés | aligné : pairs | aligné : max LDOS | décalage rigide (meV) non aligné ; aligné | QE : LDOS |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for N, r in sorted(out.get("C1", {}).items(), key=lambda kv: int(kv[0])):
        q = r["QE"]; b, al_ = r["brut"], r["aligne"]
        qs = f"{q['pi_quasi_bound']['x']:+.4f} ; {q['pi_quasi_bound']['w2']:.3f}" if q else "—"
        qe_even = "; ".join(f"{e['x']:+.3f} ({e['w2']:.3f})" for e in q["even"]) if q else "—"
        fmt = lambda s: f"{s['odd']['quasi_bound']['x']:+.4f} ; {s['odd']['quasi_bound']['w_nn']:.3f} ; {s['odd']['n_localized']}"
        ev = lambda s: "; ".join(f"{e['x']:+.3f} ({e['w_site']:.2f})" for e in s["even_window"]) or "—"
        L.append(f"| {N} | {r['first_crown']['pi']:+.3f} | {qs} | {qe_even} | {fmt(b)} | {ev(b)} | {b['ldos_max']['at_eV']:+.4f} | {fmt(al_)} | {ev(al_)} | {al_['ldos_max']['at_eV']:+.4f} | "
                 f"{b['rigid_shift_meV']:+.1f} ; {al_['rigid_shift_meV']:+.1f} | — |")
    L += ["", "## C.2 (bloc π)", "", "| N | dim | 1ʳᵉ couronne π (eV) | non aligné : π quasi-lié ; poids voisins | non aligné : max LDOS | aligné : π quasi-lié ; poids | aligné : max LDOS | eigh (s) |",
          "|---|---|---|---|---|---|---|---|"]
    for N, r in sorted(out.get("C2", {}).items(), key=lambda kv: int(kv[0])):
        if r.get("skipped"):
            L.append(f"| {N} | {2*int(N)**2} | — | non lancé (prévu {r['predicted_s']/60:.0f} min) | | | | |"); continue
        b, al_ = r["brut"], r["aligne"]
        L.append(f"| {N} | {r['dim']} | {r['first_crown']['pi']:+.3f} | {b['odd']['quasi_bound']['x']:+.4f} ; {b['odd']['quasi_bound']['w_nn']:.4f} | {b['ldos_max']['at_eV']:+.4f} | "
                 f"{al_['odd']['quasi_bound']['x']:+.4f} ; {al_['odd']['quasi_bound']['w_nn']:.4f} | {al_['ldos_max']['at_eV']:+.4f} | {b['t_eigh_s']:.0f} ; {al_['t_eigh_s']:.0f} |")
    if "C2_vs_C1_27" in out:
        L += ["", "Contrôle N = 27, bloc π contre matrice complète : " + " ; ".join(f"{k} : Δε {v['dx_quasi_bound']:+.1e} eV, Δ max LDOS {v['dldos']:+.1e} eV" for k, v in out["C2_vs_C1_27"].items())]
    if "fits" in out:
        L += ["", "## Ajustements ε(N) = ε∞ + a/N^p sur 6 … 27 et valeurs extrapolées", "", "| série | grandeur | loi | ε∞ | a | rms | " + " | ".join(f"N = {N}" for N in C2_SIZES) + " |",
              "|---|---|---|---|---|---|" + "---|" * len(C2_SIZES)]
        for k, v in out["fits"].items():
            s_, q_, p_ = k.split(" | ")
            L.append(f"| {s_} | {q_} | {p_} | {v['eps_inf']:+.4f} | {v['a']:+.3f} | {v['rms']:.4f} | " + " | ".join(f"{v['extrapolated'][str(N)]:+.4f}" for N in C2_SIZES) + " |")
    if "fits_9_27" in out:
        L += ["", "Supplément (hors prompt) : mêmes ajustements sur 9 … 27 (à N = 6, le π quasi-lié est hors de [−1, 0] eV et le maximum de LDOS est au bord de la fenêtre).", "",
              "| série | grandeur | loi | ε∞ | a | rms | " + " | ".join(f"N = {N}" for N in C2_SIZES) + " |", "|---|---|---|---|---|---|" + "---|" * len(C2_SIZES)]
        for k, v in out["fits_9_27"].items():
            s_, q_, p_ = k.split(" | ")
            L.append(f"| {s_} | {q_} | {p_} | {v['eps_inf']:+.4f} | {v['a']:+.3f} | {v['rms']:.4f} | " + " | ".join(f"{v['extrapolated'][str(N)]:+.4f}" for N in C2_SIZES) + " |")
    if "C3" in out:
        L += ["", "## C.3 (limite diluée, matrice T, même V_loc)", "", "| N_k^int | non aligné : max LDOS voisins (eV) ; valeur | aligné |", "|---|---|---|"]
        for nk, r in out["C3"].items():
            L.append(f"| {nk} | {r['brut']['at_eV']:+.4f} ; {r['brut']['value']:.3f} | {r['aligne']['at_eV']:+.4f} ; {r['aligne']['value']:.3f} |")
    with open(os.path.join(WORK, "c", "C_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def c_figure(out):
    plt, pal = fig_style()
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.0))
    N1 = sorted(int(N) for N in out.get("C1", {})); N2 = sorted(int(N) for N in out.get("C2", {}) if not out["C2"][N].get("skipped"))
    qe = [(N, out["C1"][str(N)]["QE"]["pi_quasi_bound"]["x"]) for N in N1 if out["C1"][str(N)]["QE"]]
    ax[0].plot([1 / N for N, _ in qe], [x for _, x in qe], "o", ms=3, color=pal.REF, label="QE (R7)")
    for var, c, mk in (("brut", pal.NAVY, "s"), ("aligne", pal.ORANGE, "^")):
        lab = "modèle replié" + (" (aligné)" if var == "aligne" else "")
        ax[0].plot([1 / N for N in N1] + [1 / N for N in N2 if N > 27], [out["C1"][str(N)][var]["odd"]["quasi_bound"]["x"] for N in N1]
                   + [out["C2"][str(N)][var]["odd"]["quasi_bound"]["x"] for N in N2 if N > 27], mk, ms=3, mfc="none", color=c, label=lab)
        ax[1].plot([1 / N for N in N1] + [1 / N for N in N2 if N > 27], [out["C1"][str(N)][var]["ldos_max"]["at_eV"] for N in N1]
                   + [out["C2"][str(N)][var]["ldos_max"]["at_eV"] for N in N2 if N > 27], mk, ms=3, mfc="none", color=c, label=lab)
        if "C3" in out and "600" in out["C3"]:
            ax[1].axhline(out["C3"]["600"][var]["at_eV"], color=c, lw=0.8, ls="--", label=r"matrice $T$" + (" (aligné)" if var == "aligne" else ""))
    Ns = N1 + [N for N in N2 if N > 27]
    for a_ in ax:
        a_.set_ylim(-1.12, 0.05)                                                # N = 6 : maximum de LDOS au bord de [−1, 0] (état π à −1,04 eV)
    ax[1].plot([1 / N for N in Ns], [(out["C1"].get(str(N)) or out["C2"].get(str(N)))["first_crown"]["pi"] for N in Ns], "-", color=pal.LIGHT, lw=0.8, label="première couronne")
    ax[0].set_xlabel(r"$1/N$"); ax[0].set_ylabel(r"État $\pi$ quasi-lié $\varepsilon - E_D$ (eV)"); ax[0].set_title("(a)", fontsize=9); ax[0].legend(fontsize=5)
    ax[1].set_xlabel(r"$1/N$"); ax[1].set_ylabel(r"Maximum de la LDOS des voisins (eV)"); ax[1].set_title("(b)", fontsize=9); ax[1].legend(fontsize=5)
    fig.tight_layout(); savefig(fig, "folded_vs_R7"); plt.close(fig)


# ----------------------------------------------------------------------------------------------- D : grandeur de Kaasbjerg
def pi_pair_U(Uk):
    """Bandes π, π* (ordre des énergies) : les deux états propres de plus grand poids p_z |U_3n|² + |U_4n|² (analyze_M.pi_pair)."""
    w = np.abs(Uk[..., WF_PZ_A, :]) ** 2 + np.abs(Uk[..., WF_PZ_B, :]) ** 2
    top = np.sort(np.argsort(-w, axis=-1)[..., :2], axis=-1)                  # bandes triées en énergie (eigh : croissantes)
    return top[..., 0], top[..., 1]


def cmd_d(a):
    cfg = load_production(verbose=True); d = ensure("d"); out = dict(head=git_head())
    st = wannier_setup("9x9", cfg, (3,)); dp = st["dp"]; C9, tag9 = C_of("9x9"); n2 = 81
    A_uc, _ = qe_io.get_A_volume(dp["uc"]); A_A = A_uc * BOHR; A_cell = float(np.linalg.norm(np.cross(A_A[:, 0], A_A[:, 1])))
    Bc = 2 * np.pi * np.linalg.inv(A_A).T; b1 = float(np.linalg.norm(Bc[:, 0]))
    # D.0 et D.2 : M2 dense sur la grille (définition d'analyze_M.py §1 et §3)
    kd = st["k"]; kk, eps = qe_io.get_k_eigenvalues(dp["uc"], False); eps = np.asarray(eps)
    if eps.shape[0] != len(kd):
        eps = eps.T
    eps = eps * HA2EV; V9 = np.einsum("kbw,kwv->kbv", st["Ud"], st["U"])
    kdist = lambda k, k0: np.linalg.norm(np.mod(k - np.asarray(k0) + 0.5, 1.0) - 0.5, axis=1)
    iK = int(np.argmin(kdist(kd, cfg["K_red"])))

    def pi_pair(ik):
        w = np.abs(V9[ik][:, 3]) ** 2 + np.abs(V9[ik][:, 4]) ** 2; top = np.argsort(-w)[:2]
        return tuple(sorted(top, key=lambda n_: eps[ik, n_]))
    pairK = pi_pair(iK)
    Mm = np.load(dp["mfile"], mmap_mode="r"); matrix_io.check_manifest(dp["mfile"], require_normalization=matrix_io.M_NORM_V2)
    MK = np.array(Mm[:, :, :, iK]) * HA2EV; nk = MK.shape[1]
    Vpi = np.array([A_cell * np.sqrt(sum(abs(MK[pi_pair(ik)[0], ik, n_]) ** 2 for n_ in pairK)) for ik in range(nk)])
    ok0 = abs(Vpi[iK] - 107.22) < 0.01 and abs(Vpi.max() - 116.46) < 0.01
    out["D0"] = dict(Vpi_at_K=float(Vpi[iK]), Vpi_max=float(Vpi.max()), expected=dict(at_K=107.22, max=116.46), A_cell=A_cell, pairK=list(map(int, pairK)), ok=bool(ok0))
    log(f"[D.0] porte : Ṽ_π(K) (définition d'analyze_M §1) {Vpi[iK]:.2f} eV Å² (R6 107,22), max {Vpi.max():.2f} (R6 116,46) -> {'OK' if ok0 else 'ÉCHEC'}")
    save_json(os.path.join(d, "d_results.json"), out)
    if not ok0:
        log("[D.0] ÉCHEC de la porte : STOP sur la partie D"); sys.exit(3)
    blocks = {}
    parts = {"tot": dp["mfile"], "L": dp["mfile"].replace("M_dense_", "M_L_dense_"), "NL": dp["mfile"].replace("M_dense_", "M_NL_dense_")}
    Mbox = np.load(os.path.join(WORK, "cache", "ML_box_9x9.npy"), mmap_mode="r") if has_exact("9x9") else None
    pk = list(pairK)
    for part, path in parts.items():
        matrix_io.check_manifest(path, require_normalization=matrix_io.M_NORM_V2)
        Bk = np.array(np.load(path, mmap_mode="r")[np.ix_(pk, [iK], pk, [iK])]).reshape(2, 2) * HA2EV
        vars_ = {"brut": Bk}
        if Mbox is not None and part in ("tot", "L"):
            vars_["exact"] = Bk - C9 * np.array(Mbox[np.ix_(pk, [iK], pk, [iK])]).reshape(2, 2)
        for var, Bm in vars_.items():
            ev = np.linalg.eigvalsh(0.5 * (Bm + Bm.conj().T))
            blocks[f"{part} | {var}"] = dict(eigenvalues_eVA2=[float(x) * A_cell for x in ev], half_trace_eVA2=float(0.5 * np.trace(Bm).real) * A_cell,
                                             row_norms_eVA2=[float(np.linalg.norm(Bm[i])) * A_cell for i in range(2)], frobenius_eVA2=float(np.linalg.norm(Bm)) * A_cell,
                                             hermiticity=float(np.abs(Bm - Bm.conj().T).max()))
    out["D2"] = dict(blocks=blocks, note="bloc 2×2 de la paire π/π* dégénérée à (K, K), M2 dense 9×9, × A_cell ; (a) R6 3.2 = norme de la ligne a = pairK[0] ; (b) = ½ Tr")
    log(f"[D.2] {json.dumps(jsonable(blocks), ensure_ascii=False)[:600]}")
    del MK, Mm
    # D.1 : carte de Kaasbjerg, k = K + δ x̂ hors grille, interpolation de Wannier (M_W, R de Mwk_to_Mwr)
    delta = 0.01 * b1; xhat = np.array([1.0, 0.0, 0.0])
    k = K_RED + np.linalg.solve(Bc, delta * xhat); Kc = Bc @ K_RED
    ang = float(np.degrees(np.arccos(np.clip(Kc @ xhat / np.linalg.norm(Kc), -1, 1))))
    kp = lt.mp_grid(240, 240, 1)
    exK = kdist(kp, K_RED) < 1e-9; exKp = kdist(kp, KP_RED) < 1e-9
    Ek, Uk = lt.Hwr_to_Hwk(st["Hwr"], st["Rw"], k[None], ndegen=st["nd"])[1:]
    Ep, Up = lt.Hwr_to_Hwk(st["Hwr"], st["Rw"], kp, ndegen=st["nd"])[1:]
    pk_v, pk_c = pi_pair_U(Uk[0][None])                                          # bandes π, π* à k
    pv, pc_ = pi_pair_U(Up)
    Mwr_var = {"brut": st["Mwr"]}
    if has_exact("9x9"):
        Z = np.load(os.path.join(WORK, "cache", "Fw_9x9.npz")); assert np.array_equal(Z["R"], st["R"]); Mwr_var["exact"] = st["Mwr"] - C9 * Z["Fw"]
    # disques autour de K et K' (image réciproque la plus proche), rayon 0,05 |b₁|
    kc = (Bc[:2, :2] @ kp[:, :2].T).T
    imgs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)]) @ Bc[:2, :2].T

    def disk(K0):
        K0c = Bc[:2, :2] @ K0[:2]; dd = np.min(np.linalg.norm(kc[:, None, :] - K0c[None, None, :] + imgs[None], axis=2), axis=1)
        return dd <= 0.05 * b1
    dK, dKp = disk(K_RED) & ~exK, disk(KP_RED) & ~exKp
    maps = {}
    for var, Mw in Mwr_var.items():
        Mk = Mwr_to_Mwk_pairs(Mw, st["R"], kp, k[None])[:, :, :, 0]                # (5, n', 5)
        Mb = np.einsum("kwb,wkW,WB->bkB", Up.conj(), Mk, Uk[0], optimize=True)      # (bandes k', k', bandes k)
        Vv = A_cell * np.abs(Mb[pv, np.arange(len(kp)), pk_v[0]]); Vc = A_cell * np.abs(Mb[pc_, np.arange(len(kp)), pk_c[0]])
        Vv[exK | exKp] = np.nan; Vc[exK | exKp] = np.nan
        maps[var] = dict(valence=Vv, conduction=Vc)
        stats = {}
        for lab, m in (("K", dK), ("K'", dKp)):
            for band, Vm in (("valence", Vv), ("conduction", Vc)):
                stats[f"{band} | {lab}"] = dict(mean=float(np.nanmean(Vm[m])), min=float(np.nanmin(Vm[m])), max=float(np.nanmax(Vm[m])), n=int(np.isfinite(Vm[m]).sum()))
        out.setdefault("D1", {})[var] = stats
        log(f"[D.1] {var} : " + " ; ".join(f"{kk_} moyenne {v['mean']:.2f} [{v['min']:.2f}, {v['max']:.2f}] eV Å² ({v['n']})" for kk_, v in stats.items()))
    out["D1_meta"] = dict(k_red=k.tolist(), delta_Ainv=delta, xhat_angle_to_GammaK_deg=ang, E_D=st["E_D"],
                          eps_k_minus_ED=[float(Ek[0, pk_v[0]] - st["E_D"]), float(Ek[0, pk_c[0]] - st["E_D"])],
                          grid=240, excluded_K_Kp=int(exK.sum() + exKp.sum()), disk_radius_Ainv=0.05 * b1, C9_eV=C9, C9_label=tag9,
                          Kaasbjerg=dict(V_tilde_eVA2=70, V0_eV=27, supercell="11x11"))
    np.savez(os.path.join(d, "d1_maps.npz"), kc=kc, Bc=Bc, **{f"{v}_{b}": maps[v][b] for v in maps for b in maps[v]})
    save_json(os.path.join(d, "d_results.json"), out); d_tables(out); d_figure(kc, Bc, maps, Kc)


def cmd_dfig(a):
    """Tables et figure de D depuis d/d_results.json et d/d1_maps.npz (sans recalcul) ; corrige l'étiquette des énergies à k du premier passage."""
    out = json.load(open(os.path.join(WORK, "d", "d_results.json"))); m = out["D1_meta"]
    if "E_D" not in m:                                                          # passage du 2026-09-28 09:30 : énergies absolues sous le nom eps_k_minus_ED
        E_D = json.load(open(os.path.join(WORK, "b", "b_results.json")))["9x9"]["E_D"]   # même H(R) 27x27, même définition (grille 90²)
        m["eps_k_abs"] = m["eps_k_minus_ED"]; m["eps_k_minus_ED"] = [e - E_D for e in m["eps_k_abs"]]; m["E_D"] = E_D
        save_json(os.path.join(WORK, "d", "d_results.json"), out)
    Z = np.load(os.path.join(WORK, "d", "d1_maps.npz")); Bc = Z["Bc"]
    maps = {v: dict(valence=Z[f"{v}_valence"], conduction=Z[f"{v}_conduction"]) for v in ("brut", "aligne", "exact") if f"{v}_valence" in Z}
    d_tables(out); d_figure(Z["kc"], Bc, maps, Bc @ K_RED)


def d_tables(out):
    L = ["# R9 — D : grandeur de Kaasbjerg (tables générées par r9_driver.py d)", "",
         f"Porte D.0 : Ṽ_π(K) = {out['D0']['Vpi_at_K']:.2f} eV Å² (R6 107,22), max {out['D0']['Vpi_max']:.2f} (R6 116,46) → {'OK' if out['D0']['ok'] else 'ÉCHEC'} ; A_cell = {out['D0']['A_cell']:.4f} Å².", ""]
    if "D2" in out:
        L += ["## D.2 : bloc 2 × 2 de la paire π/π* à (K, K) (eV Å²)", "", "| partie | variante | valeurs propres | ½ Tr | normes des lignes | ‖·‖_F |", "|---|---|---|---|---|---|"]
        for k, v in out["D2"]["blocks"].items():
            p_, var = k.split(" | ")
            L.append(f"| {p_} | {VAR_LABEL.get(var, var)} | {v['eigenvalues_eVA2'][0]:.2f} ; {v['eigenvalues_eVA2'][1]:.2f} | {v['half_trace_eVA2']:.2f} | "
                     f"{v['row_norms_eVA2'][0]:.2f} ; {v['row_norms_eVA2'][1]:.2f} | {v['frobenius_eVA2']:.2f} |")
    if "D1" in out:
        m = out["D1_meta"]
        L += ["", f"## D.1 : Ṽ = A_cell·|M^nn_k'k|, k = K + δx̂ (δ = {m['delta_Ainv']:.4f} Å⁻¹ ; x̂ à {m['xhat_angle_to_GammaK_deg']:.1f}° de ΓK ; ε(k) − E_D : π {m['eps_k_minus_ED'][0]:+.4f}, π* {m['eps_k_minus_ED'][1]:+.4f} eV), carte 240², "
              f"disques de rayon {m['disk_radius_Ainv']:.4f} Å⁻¹, K et K' exacts exclus ({m['excluded_K_Kp']} points)", "",
              "| variante | bande | disque | moyenne | min | max | n |", "|---|---|---|---|---|---|---|"]
        for var, s in out["D1"].items():
            for k, v in s.items():
                b, lab = k.split(" | ")
                L.append(f"| {VAR_LABEL.get(var, var)} | {b} | {lab} | {v['mean']:.2f} | {v['min']:.2f} | {v['max']:.2f} | {v['n']} |")
        L += ["", "Kaasbjerg (PRB 101, 045433, Fig. 3) : ~70 eV Å² près de K et K′ ; V₀ ≈ 27 eV ; super-cellule 11×11."]
    with open(os.path.join(WORK, "d", "D_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")


def d_figure(kc, Bc, maps, Kc):
    plt, pal = fig_style()
    imgs = np.array([[i, j] for i in (-1, 0, 1) for j in (-1, 0, 1)]) @ Bc[:2, :2].T
    kf = np.array([kc[i] + imgs[np.argmin(np.linalg.norm(kc[i] + imgs, axis=1))] for i in range(len(kc))])
    # coins de la zone = points K et K' de la cellule QE (a₁, a₂ à 60° : K = (2/3, 1/3)) ; la liste d'analyze_M.py (1/3, 1/3), … ne l'est pas
    corners = np.array([Bc[:2, :2] @ np.array(c) for c in [(2/3, 1/3), (1/3, 2/3), (-1/3, 1/3), (-2/3, -1/3), (-1/3, -2/3), (1/3, -1/3)]])
    corners = corners[np.argsort(np.arctan2(corners[:, 1], corners[:, 0]))]
    vars_ = list(maps); fig, axes = plt.subplots(len(vars_), 2, figsize=(6.5, 3.0 * len(vars_)), squeeze=False)
    vmax = max(np.nanmax(maps[v][b]) for v in vars_ for b in ("valence", "conduction"))
    vmin = min(np.nanmin(maps[v][b]) for v in vars_ for b in ("valence", "conduction"))
    tlab = {"brut": "M2 tel quel", "aligne": "aligné (i)", "exact": r"aligné exact ($F_W$)"}
    for i, var in enumerate(vars_):
        for j, (band, lab) in enumerate((("valence", r"valence ($\pi$)"), ("conduction", r"conduction ($\pi^*$)"))):
            ax = axes[i][j]; z = maps[var][band]; m = np.isfinite(z)
            sc = ax.scatter(kf[m, 0], kf[m, 1], c=z[m], s=0.6, cmap=pal.CMAP_SEQ, vmin=vmin, vmax=vmax, rasterized=True)
            ax.plot(np.r_[corners[:, 0], corners[0, 0]], np.r_[corners[:, 1], corners[0, 1]], color=pal.REF, lw=0.6)
            ax.plot(Kc[0], Kc[1], "x", color=pal.ORANGE, ms=4)
            ax.set_aspect("equal"); ax.set_xlabel(r"$k'_x$ (Å$^{-1}$)", fontsize=7); ax.set_ylabel(r"$k'_y$ (Å$^{-1}$)", fontsize=7)
            ax.set_title(f"({'abcd'[2*i+j]}) {lab}, {tlab[var]}", fontsize=8)
            fig.colorbar(sc, ax=ax, label=r"$A_\mathrm{cell}|M_{k'k}|$ (eV Å$^2$)")
    fig.tight_layout(); savefig(fig, "kaasbjerg_fig3_map"); plt.close(fig)



# ----------------------------------------------------------------------------------------------- clôture : porte R.0 et synthèse
def cmd_r0(a):
    """Porte R.0 (clôture) : la variante Lu rejouée par --cn lu redonne R9 — 9x9 et 12x12, R_cut 3 : médiane de niveau 1 alignée (i) et tab:rcut_M."""
    global CN
    CN = "lu"; cfg = load_production(verbose=True); d = ensure("a")
    R9 = json.load(open(os.path.join(d, "a3_results.json"))); out = dict(head=git_head(), sizes={}); ok_all = True
    for S in ("9x9", "12x12"):
        st = resonance_setup(S, cfg, (3,)); C, tag = C_of(S)
        s = level1_stats(st, sigma_states(st, V_variant(st, 3, "aligne", C), g0_cached(st, cfg["nk_int"], 3, "sig"), 3))
        ref_m = R9["level1"][S]["rc"]["3"]["variants"]["aligne"]["median_meV"]
        rt = rcut_table(cfg, S, rcuts=[3])
        got_r = rt["rows"][vkey("aligne")][0]["max_dM_over_maxM"]
        ref_r = [r for r in R9[f"rcut_{S}"]["rows"]["aligne"] if r["R_cut"] == 3][0]["max_dM_over_maxM"]
        ok = abs(s["median_meV"] - ref_m) <= 1e-9 * ref_m and abs(got_r - ref_r) <= 1e-9 * ref_r and C == R9["level1"][S]["C_eV"]
        ok_all &= ok
        out["sizes"][S] = dict(C_eV=C, C_R9_eV=R9["level1"][S]["C_eV"], median_meV=s["median_meV"], median_R9_meV=ref_m, rcut3=got_r, rcut3_R9=ref_r, ok=bool(ok))
        log(f"[R.0] {S} : C_N {C*1e3:+.4f} meV (R9 {R9['level1'][S]['C_eV']*1e3:+.4f}) ; médiane alignée (i) {s['median_meV']:.4f} meV (R9 {ref_m:.4f}) ; "
            f"tab:rcut_M R_cut 3 {got_r:.6e} (R9 {ref_r:.6e}) -> {'OK' if ok else 'ÉCHEC'}")
    out["ok"] = bool(ok_all); save_json(os.path.join(d, "r0_gate.json"), out)
    if not ok_all:
        log("[R.0] ÉCHEC de la porte : STOP sur le rejeu"); sys.exit(3)


def cmd_synth(a):
    """
    R.4 (clôture) : trois colonnes (tel quel, Lu, plateau) pour R.1–R.3, depuis a3_results.json et b_results.json (R9 : brut, aligne = Lu, exact = Lu)
    et a3_results_plateau.json, b_results_plateau.json (clôture) -> cloture/synthese.md ; figures resonance_vs_nkint (trois variantes ; la version
    de R9 est gardée sous resonance_vs_nkint_R9) et rcut_aligned. Aucun calcul.
    """
    import shutil
    A9 = json.load(open(os.path.join(WORK, "a", "a3_results.json"))); AP = json.load(open(os.path.join(WORK, "a", "a3_results_plateau.json")))
    B9 = json.load(open(os.path.join(WORK, "b", "b_results.json"))); BP = json.load(open(os.path.join(WORK, "b", "b_results_plateau.json")))
    a1 = json.load(open(os.path.join(WORK, "a", "a1_results.json")))["sizes"]
    fg = os.path.join(WORK, "a", "r0_gate.json"); G0 = json.load(open(fg)) if os.path.exists(fg) else None
    f3 = lambda x, fmt: "—" if x is None else fmt.format(x)
    L = ["# R9 — Clôture : synthèse tel quel / C_N Lu / C_N plateau (générée par r9_driver.py synth)", "",
         "C_N (meV), Lu ; plateau : " + " ; ".join(f"{S} {a1[S]['C_retenu_eV']*1e3:+.2f} ; {a1[S]['C_i_eV']*1e3:+.2f}" for S in SIZES_A), ""]
    if G0:
        L += ["Porte R.0 : " + " ; ".join(f"{S} médiane {v['median_meV']:.4f} (R9 {v['median_R9_meV']:.4f}), tab:rcut_M {v['rcut3']:.4e} (R9 {v['rcut3_R9']:.4e})"
                                    for S, v in G0["sizes"].items()) + f" → {'OK' if G0['ok'] else 'ÉCHEC'}", ""]

    def crowns(S, e):
        c = crowns_containing(AP["level1"][S]["crowns_pi_240"], e, tol=1e-6)
        return "/".join(f"x={r['x']}" for r in c) or "—"
    L += ["## R.1 niveau 1 (240², η 0,02, N_k^int 300) et courbe Γ_T", "",
          "| taille | R_cut | grandeur | tel quel | aligné (i) Lu | aligné (i) plateau | exact Lu | exact plateau |", "|---|---|---|---|---|---|---|---|"]
    Q = (("médiane Γ·N_cells (meV)", "median_meV", "{:.2f}"), ("E_res (eV)", "E_res", "{:+.4f}"), ("pic Γ_T, 5 meV", "peak_GT_prod", "{:+.3f}"),
         ("pic Γ_T, 2,5 meV", "peak_GT_fine", "{:+.4f}"), ("pic −Im T̄(K)", "peak_ImTbar", "{:+.4f}"), ("Re Σ médian (meV)", "median_ReSigma_meV", "{:.2f}"),
         ("Re T̄(E_D) (eV)", "ReTbar_at_ED", "{:+.3f}"))
    for S in SIZES_A:
        for rc in ("3", "4"):
            v9 = A9["level1"][S]["rc"][rc]["variants"]; vp = AP["level1"][S]["rc"][rc]["variants"]
            cols = [v9.get("brut"), v9.get("aligne"), vp.get("aligne_plateau"), v9.get("exact"), vp.get("exact_plateau")]
            for lab, key, fmt in Q:
                vals = [f3(c[key] if c else None, fmt) for c in cols]
                if key == "E_res":
                    vals = [v if c is None else f"{v} ({crowns(S, c[key])})" for v, c in zip(vals, cols)]
                L.append(f"| {S} | {rc} | {lab} | " + " | ".join(vals) + " |")
    L += ["", "Familles (médianes Γ·N_cells ; moyenne meV, (max − min)/moyenne) :", "", "| R_cut | famille | tel quel | Lu | plateau |", "|---|---|---|---|---|"]
    for rc in ("3", "4"):
        for fname, members in FAMILIES.items():
            cells = []
            for src, var in ((A9, "brut"), (A9, "aligne"), (AP, "aligne_plateau")):
                v = np.array([src["level1"][S]["rc"][rc]["variants"][var]["median_meV"] for S in members])
                cells.append(f"{v.mean():.1f} ; {100*(v.max()-v.min())/v.mean():.2f} %")
            L.append(f"| {rc} | {fname} | " + " | ".join(cells) + " |")
    for S in ("9x9", "12x12"):
        r9 = A9[f"rcut_{S}"]["rows"]; rp = AP[f"rcut_{S}"]["rows"]
        cols = [("tel quel", r9["brut"]), ("Lu", r9["aligne"]), ("plateau", rp["aligne_plateau"])] + \
               ([("exact Lu", r9["exact"]), ("exact plateau", rp["exact_plateau"])] if "exact" in r9 and "exact_plateau" in rp else [])
        L += ["", f"## R.1 tab:rcut_M, {S}", "", "| R_cut | sites | " + " | ".join(c for c, _ in cols) + " |", "|---|---|" + "---|" * len(cols)]
        for i in range(len(r9["brut"])):
            L.append(f"| {r9['brut'][i]['R_cut']} | {r9['brut'][i]['n_sites']} | " + " | ".join(f"{rows[i]['max_dM_over_maxM']:.3e}" for _, rows in cols) + " |")
    c9, cp = A9["conv_intensive"], AP["conv_intensive"]
    L += ["", "## R.1 convention intensive (max|M|, bandes 1–16, eV)", "", "| ensemble | taille | tel quel | Lu | plateau | mode aligné |", "|---|---|---|---|---|---|"]
    for k, r in c9["rows"].items():
        L.append(f"| {r['set']} | {r['S']} | {r['brut']:.3f} | {r['aligne']:.3f} | {cp['rows'][k]['aligne_plateau']:.3f} | {r['mode_aligne']} |")
    L += ["", "| ensemble | tel quel | Lu | plateau |", "|---|---|---|---|"]
    for setname in ("grossiers 5, 7, 8, 9", "denses six tailles", "denses 5, 7, 8, 9"):
        L.append(f"| {setname} | {c9['spreads'][setname + ' | brut']:.3e} | {c9['spreads'][setname + ' | aligne']:.3e} | {cp['spreads'][setname + ' | aligne_plateau']:.3e} |")
    oi = AP["offdiag_intensive"]
    L += ["", "## R.2 test intensif hors k = k' (max|M| sur k ≠ k', bandes 1–16, eV)", "", "| ensemble | taille | tel quel | Lu | plateau | mode |", "|---|---|---|---|---|---|"]
    for k, r in oi["rows"].items():
        L.append(f"| {r['set']} | {r['S']} | {r['brut']:.3f} | {r['lu']:.3f} | {r['plateau']:.3f} | {r['mode']} |")
    L += ["", "| ensemble | tel quel | Lu | plateau |", "|---|---|---|---|"]
    for setname in ("grossiers 5, 7, 8, 9", "denses six tailles"):
        L.append(f"| {setname} | " + " | ".join(f"{oi['spreads'][f'{setname} | {v}']:.3e}" for v in ("brut", "lu", "plateau")) + " |")
    L += ["", "## R.1 résonance 9×9 (R_cut 3, 300²) et critère de pôle", "", "| grandeur | tel quel | Lu | plateau | exact Lu | exact plateau |", "|---|---|---|---|---|---|"]
    v9 = A9["level1"]["9x9"]["rc"]["3"]["variants"]; vp = AP["level1"]["9x9"]["rc"]["3"]["variants"]
    cols = [v9.get("brut"), v9.get("aligne"), vp.get("aligne_plateau"), v9.get("exact"), vp.get("exact_plateau")]
    for lab, key, fmt in (("Born/T médian", "Born_over_T_median", "{:.2f}"), ("max courbe Γ_T (eV)", "max_GT_fine", "{:.3f}"), ("max −Im T̄(K) (eV)", "max_mImTbar", "{:.2f}"),
                          ("|Re Σ|/Γ médian", "median_absRe_over_G", "{:.3f}")):
        L.append(f"| {lab} | " + " | ".join(f3(c.get(key) if c else None, fmt) for c in cols) + " |")
    p9 = A9["pole_9x9"]["variants"]; pp = AP["pole_9x9"]["variants"]
    pc_ = [p9.get("brut"), p9.get("aligne"), pp.get("aligne_plateau"), p9.get("exact"), pp.get("exact_plateau")]
    for b in ("pi", "complet"):
        L.append(f"| bloc {b} : min\\|det\\|/max (ε − E_D) | " + " | ".join("—" if c is None else f"{c[b]['min_absdet_rel']:.3e} ({c[b]['at_eV']:+.3f})" for c in pc_) + " |")
        L.append(f"| bloc {b} : min\\|λ\\| (ε − E_D) | " + " | ".join("—" if c is None else f"{c[b]['min_abs_lambda']:.4f} ({c[b]['lambda_at_eV']:+.3f})" for c in pc_) + " |")
    L += ["", "## R.3 B : résonance contre N_k^int", "", "| taille | N_k^int | grandeur | tel quel | Lu | plateau | exact Lu | exact plateau |", "|---|---|---|---|---|---|---|---|"]
    for S in ("9x9", "12x12"):
        for nk in ("300", "450", "600", "900"):
            cols = [B9[S]["brut"].get(nk), B9[S]["aligne"].get(nk), BP[S]["aligne_plateau"].get(nk),
                    B9[S].get("exact", {}).get(nk), BP[S].get("exact_plateau", {}).get(nk)]
            for lab, key, fmt in (("pic Γ_T", "peak_GT_fine", "{:+.4f}"), ("pic −Im T̄(K)", "peak_ImTbar", "{:+.4f}"), ("Γ_T(E_D) (meV)", "G_ED_meV", "{:.1f}"),
                                  ("E_res", "E_res", "{:+.4f}")):
                vals = [f3(c[key] if c else None, fmt) for c in cols]
                if key == "E_res":
                    vals = [v if c is None else f"{v} (x={c['E_res_crown']['x']})" for v, c in zip(vals, cols)]
                L.append(f"| {S} | {nk} | {lab} | " + " | ".join(vals) + " |")
    d = ensure("cloture")
    with open(os.path.join(d, "synthese.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    log(f"[synth] cloture/synthese.md ({len(L)} lignes)")
    # figures
    for ext in ("pdf", "png"):
        src = os.path.join(WORK, "fig", f"resonance_vs_nkint.{ext}"); dst = os.path.join(WORK, "fig", f"resonance_vs_nkint_R9.{ext}")
        if os.path.exists(src) and not os.path.exists(dst):
            shutil.copy2(src, dst)
    merged = {S: {"brut": B9[S]["brut"], "aligne_lu": B9[S]["aligne"], "aligne_plateau": BP[S]["aligne_plateau"]} for S in ("9x9", "12x12")}
    b_figure(merged)
    plt, pal = fig_style()
    fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.0), sharey=True)
    for ax, S, lab in zip(axes, ("9x9", "12x12"), "ab"):
        r9 = A9[f"rcut_{S}"]["rows"]; rp = AP[f"rcut_{S}"]["rows"]
        series = [(r9["brut"], pal.NAVY, "-", "o", "M2 tel quel"), (r9["aligne"], pal.ORANGE, "--", "s", r"aligné, $C_N$ Lu"),
                  (rp["aligne_plateau"], pal.GREEN, "-.", "D", r"aligné, $C_N$ plateau")]
        if "exact" in r9 and "exact_plateau" in rp:
            series += [(r9["exact"], pal.ORANGE, ":", "^", r"exact, $C_N$ Lu"), (rp["exact_plateau"], pal.GREEN, ":", "v", r"exact, $C_N$ plateau")]
        for rows, c, ls, mk, lb in series:
            ax.semilogy([r["R_cut"] for r in rows], [r["max_dM_over_maxM"] for r in rows], ls, marker=mk, ms=3, color=c, label=lb, mfc="none")
        ax.set_xlabel(r"Rayon de coupure $R_\mathrm{cut}$ (mailles)"); ax.set_title(f"({lab}) {S.replace('x', '×')}", fontsize=9)
    axes[0].set_ylabel(r"$\max|\Delta M|/\max|M|$ (paire $\pi$)"); axes[0].legend(fontsize=5)
    fig.tight_layout(); savefig(fig, "rcut_aligned"); plt.close(fig)


# ----------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="R9 : pilote unique (voir la docstring)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("a2c")
    p = sub.add_parser("a2d"); p.add_argument("--size", required=True); p.add_argument("--block-size", type=int, default=2000)
    p = sub.add_parser("a2dpost"); p.add_argument("--size", required=True)
    sub.add_parser("a0")
    p = sub.add_parser("a1"); p.add_argument("--sizes", default=",".join(SIZES_A))
    p = sub.add_parser("a3"); p.add_argument("--sizes", default=",".join(SIZES_A)); p.add_argument("--redo", action="store_true")
    p.add_argument("--cn", choices=["lu", "plateau"], default=None)
    p = sub.add_parser("a3pole"); p.add_argument("--cn", choices=["lu", "plateau"], default=None)
    p = sub.add_parser("b"); p.add_argument("--sizes", default="9x9,12x12"); p.add_argument("--nk", default="300,450,600,900")
    p.add_argument("--variants", default="brut"); p.add_argument("--cn", choices=["lu", "plateau"], default=None)
    p = sub.add_parser("btables"); p.add_argument("--cn", choices=["lu", "plateau"], default=None)
    sub.add_parser("r0")
    sub.add_parser("synth")
    p = sub.add_parser("c"); p.add_argument("--redo", action="store_true")
    sub.add_parser("d")
    sub.add_parser("dfig")
    sub.add_parser("cfig")
    a = ap.parse_args()
    global CN
    CN = getattr(a, "cn", None)                                                  # clôture : C_N Lu ou plateau (a3, a3pole, b, btables)
    {"r0": cmd_r0, "synth": cmd_synth, "a2c": cmd_a2c, "a2d": cmd_a2d, "a2dpost": cmd_a2dpost, "a0": cmd_a0, "a1": cmd_a1, "a3": cmd_a3, "a3pole": cmd_a3pole, "b": cmd_b,
     "btables": cmd_btables, "c": cmd_c, "d": cmd_d, "dfig": cmd_dfig, "cfig": cmd_cfig}[a.cmd](a)


if __name__ == "__main__":
    main()
