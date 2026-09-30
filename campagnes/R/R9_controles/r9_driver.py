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
Audit (2026-09-28) :
  audit    image minimale dans la cellule à 60° : réduction axe par axe (recopiée) contre vraie image (9 images), sur les géométries réelles
           5…27 non relaxées et 9x9 relaxée R1, les étiquettes R de Wannier, les grilles k ; géométrie seule, nœud de connexion ;
           sortie audit/audit_results.json (rapport audit_image_minimale.md).
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
sys.path.insert(0, os.path.join(PROJ, "scripts")); sys.path.insert(0, os.path.join(PROJ, "scripts", "m"))   # _bands, _palette ; compute_M_dense_stages (scripts/m/ depuis le 2026-09-30)                              # compute_M_dense_stages.paths (mêmes entrées que la production)
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


FIG_LABEL = {"aligne": r"aligné", "exact": r"aligné exact", "aligne_lu": r"aligné, $C_N$ Lu", "exact_lu": r"exact, $C_N$ Lu",
             "aligne_plateau": r"aligné, $C_N$ plateau", "exact_plateau": r"exact, $C_N$ plateau"}


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
    sizes = [S for S in ("9x9", "12x12") if S in res and any(v in res[S] for v in VAR_LABEL)]     # clôture : fichiers --cn sans « brut »
    fig, axes = plt.subplots(1, len(sizes), figsize=(6.5, 3.5), sharey=True, squeeze=False)
    for ax, S, lab in zip(axes[0], sizes, "ab"):
        for var, ls, mk in [(v, l_, m_) for v, l_, m_ in (("brut", "-", "o"), ("aligne", "--", "s"), ("exact", ":", "^"), ("aligne_lu", "--", "s"),
                                                          ("exact_lu", ":", "^"), ("aligne_plateau", "-.", "D"), ("exact_plateau", ":", "v")) if v in res[S]]:
            rr = sorted(res[S][var].items(), key=lambda kv: int(kv[0])); x = [1.0 / int(nk) for nk, _ in rr]
            suf = "" if var == "brut" else " (" + FIG_LABEL.get(var, var) + ")"
            ax.plot(x, [r["peak_GT_fine"] for _, r in rr], ls, marker=mk, ms=3, color=pal.NAVY, label=r"pic de $\Gamma_T$" + suf)
            ax.plot(x, [r["peak_ImTbar"] for _, r in rr], ls, marker=mk, ms=3, color=pal.ORANGE, label=r"pic de $-\mathrm{Im}\,\bar T(K)$" + suf)
            ax.plot(x, [r["E_res"] for _, r in rr], ls, marker=mk, ms=3, color=pal.GREEN, label=r"$E_\mathrm{res}$ (états $240^2$)" + suf)
        ax.set_xlabel(r"$1/N_k^\mathrm{int}$"); ax.set_title(f"({lab}) {S.replace('x', '×')}", fontsize=9)
    axes[0][0].set_ylabel(r"Énergie $\varepsilon - E_D$ (eV)")
    h, l_ = axes[0][0].get_legend_handles_labels()                               # légende commune sous les panneaux (ne couvre aucune courbe)
    fig.legend(h, l_, loc="lower center", ncol=3, fontsize=5, frameon=False)
    fig.tight_layout(rect=(0, 0.035 * ((len(l_) + 2) // 3) + 0.03, 1, 1)); savefig(fig, "resonance_vs_nkint" + cn_suffix()); plt.close(fig)


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
        rows[f"grossier {S}"] = dict(set="grossiers 5, 7, 8, 9", S=S, brut=m, lu=m, plateau=m,
                                    mode="identique par construction (a) : identité exacte, M^L[C] ∝ δ_kk'")
    for S in SIZES_A:
        dp = dense_paths(cfg, S); matrix_io.check_manifest(dp["mfile"], require_normalization=matrix_io.M_NORM_V2)
        Mm = np.load(dp["mfile"], mmap_mode="r"); M = np.array(Mm[:16, :, :16, :]) * HA2EV; del Mm
        m = off_max(M)
        if has_exact(S):
            Mb = np.array(np.load(os.path.join(WORK, "cache", f"ML_box_{S}.npy"), mmap_mode="r")[:16, :, :16, :])
            lu = off_max(M - a1[S]["C_retenu_eV"] * Mb); pl = off_max(M - a1[S]["C_i_eV"] * Mb); del Mb
            mode = "exact (M2 − C_N M^L[1_boîte]) : effet de l'alignement hors diagonale mesuré"
        else:
            lu = pl = m; mode = "identique par construction (b) : alignement approché sur les blocs k = k' seulement"
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


def cmd_a3tables(a):
    """Tables d'A.3 depuis a3_results<suffixe>.json (sans calcul)."""
    a3_tables(json.load(open(os.path.join(WORK, "a", f"a3_results{cn_suffix()}.json"))))


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
    with open(os.path.join(WORK, "a", f"A3_tables{cn_suffix()}.md"), "w") as f:            # clôture : A3_tables_plateau.md (R9 : A3_tables.md)
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
    L += ["", "## R.2 test intensif hors k = k' (max|M| sur k ≠ k', bandes 1–16, eV)", "",
          "Lignes « identique par construction » : (a) grossiers, identité exacte (M^L[C] ∝ δ_kk') ; (b) denses 6, 7, 8, 12, alignement approché sur les seuls blocs "
          "k = k'. Seuls les denses 5×5 et 9×9 (alignement exact) mesurent un effet de l'alignement hors diagonale.", "",
          "| ensemble | taille | tel quel | Lu | plateau | mode |", "|---|---|---|---|---|---|"]
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


# ----------------------------------------------------------------------------------------------- audit : image minimale (2026-09-28)
# Géométrie seule (nœud de connexion) : positions et grilles FFT lues dans les XML, étiquettes R, grilles k. Lectures en plus : a/a1_profiles_*.npz
# et cache/{Mwr,Fw}_*.npz (R9), results/M2/{ved_analysis,mwr_locality}.npz, relax.out R2 12x12, json de R4/R5/R7. Aucune routine de production
# appelée hors lecture : la réduction axe par axe est recopiée ici à l'identique ; la vraie image minimale est une recherche sur 9 images.
AUD_SIZES = [5, 6, 7, 8, 9, 10, 11, 12, 15, 18, 21, 24, 27]
QE_TMP = "/home/gregb26/links/scratch/qe_tmp"
IMG = np.array([(i, j, 0) for i in (-1, 0, 1) for j in (-1, 0, 1)], float)                 # images dans le plan (a3 ⟂ plan)
A_UNIT = np.array([[4.0354919061, -2.3298923383], [4.0354919061, 2.3298923383]]) / 4.6597846766   # lignes a1, a2 en unités de a (mwr_locality)
TIE = 1e-6


def red_axis(d):
    """Réduction axe par axe (min_image_dist, sphere_average, inplane_disc_mask, analyze_Ved, analyze_M) : d -> (d + 1/2) mod 1 − 1/2."""
    return (np.asarray(d, float) + 0.5) % 1.0 - 0.5


def cart_axis(d, A):
    """Vecteur cartésien de la réduction axe par axe (A : colonnes a_i)."""
    return red_axis(d) @ np.asarray(A).T


def cart_true(d, A):
    """
    Vraie image minimale : réduction axe par axe, puis la plus courte des 9 images d + s, s ∈ {−1, 0, 1}² dans le plan (a3 ⟂ plan).
    Suffisant pour a1, a2 à 60° : la cellule de Wigner-Seitz a ses coordonnées réduites dans [−2/3, 2/3], donc s = v* − d ∈ (−7/6, 7/6].
    Retourne (vecteur (n, 3), distance (n,), nombre d'images à égalité à TIE près (n,)).
    """
    dd = red_axis(np.atleast_2d(d))
    r = (dd[:, None, :] + IMG[None]) @ np.asarray(A).T
    L = np.linalg.norm(r, axis=2); m = L.min(1)
    return r[np.arange(len(dd)), L.argmin(1)], m, (L <= m[:, None] + TIE).sum(1)


def aud_geom(N):
    t = f"{N}x{N}"; sd = f"{QE_TMP}/defect_{t}_d/defect_{t}_d.save"; sp_ = f"{QE_TMP}/defect_{t}_p/defect_{t}_p.save"
    A_b, _ = qe_io.get_A_volume(sd); A_bp, _ = qe_io.get_A_volume(sp_)
    assert np.allclose(A_b, A_bp) and abs(A_b[0, 2]) + abs(A_b[1, 2]) + abs(A_b[2, 0]) + abs(A_b[2, 1]) < 1e-10, t
    return dict(N=N, A_b=A_b, A=A_b * BOHR, x_d=qe_io.get_x_red(sd), x_p=qe_io.get_x_red(sp_), ng=tuple(int(v) for v in qe_io.get_ngfft(sd)))


def aud_vacancy(x_p, x_d, A):
    """vacancy_site (axe par axe, cartésien) contre vraie image ; variante d'analyze_M / analyze_Ved / tag_vacancy_sublattice (norme réduite)."""
    dax = np.array([np.linalg.norm(cart_axis(x_d - p, A), axis=1).min() for p in x_p])
    dtr = np.array([cart_true(x_d - p, A)[1].min() for p in x_p])
    dred = np.array([np.linalg.norm(red_axis(np.mod(x_d, 1.0) - np.mod(p, 1.0)), axis=1).min() for p in x_p])
    ia, it, ir = int(np.argmax(dax)), int(np.argmax(dtr)), int(np.argmax(dred))
    out = dict(i_axis_1based=ia + 1, i_true_1based=it + 1, i_rednorm_1based=ir + 1, dmin_axis_A=float(dax[ia]), dmin_true_A=float(dtr[it]),
               next_axis_A=float(np.sort(dax)[-2]), next_true_A=float(np.sort(dtr)[-2]), same=bool(ia == it == ir and abs(dax[ia] - dtr[it]) < 1e-9))
    return out, np.mod(x_p[ia], 1.0)


def aud_far(x_d, x_p, s_vac, A):
    """far_atom (axe par axe) contre l'ensemble des atomes vraiment les plus loin ; partenaire de la parfaite (far_atom_alignment)."""
    dax = np.linalg.norm(cart_axis(x_d - s_vac, A), axis=1); _, dtr, _ = cart_true(x_d - s_vac, A)
    ia = int(np.argmax(dax)); rmax = float(dtr.max()); far = np.where(dtr >= rmax - TIE)[0]
    diff = np.abs(dax - dtr) > TIE

    def partner(i):
        c = np.mod(x_d[i], 1.0)
        pa = int(np.argmin(np.linalg.norm(cart_axis(x_p - c, A), axis=1))); pt = int(np.argmin(cart_true(x_p - c, A)[1]))
        return dict(axis_1based=pa + 1, true_1based=pt + 1, displacement_A=float(cart_true(x_p[pt] - c, A)[1][0]))
    out = dict(n_atoms=len(x_d), i_axis_1based=ia + 1, d_axis_A=float(dax[ia]), d_true_of_axis_atom_A=float(dtr[ia]),
               rank_true_of_axis_atom=int((dtr > dtr[ia] + TIE).sum()) + 1, frac_rmax_of_axis_atom=float(dtr[ia] / rmax),
               r_max_true_A=rmax, far_true_1based=[int(i) + 1 for i in far], a_sc_over_sqrt3_A=float(np.linalg.norm(A[:, 0]) / np.sqrt(3)),
               axis_atom_in_far_true=bool(ia in far), n_atoms_dist_axis_ne_true=int(diff.sum()),
               min_true_dist_where_differs_A=float(dtr[diff].min()) if diff.any() else None,
               min_axis_dist_where_differs_A=float(dax[diff].min()) if diff.any() else None,
               partner_axis_atom=partner(ia), partner_far_true=[partner(i) for i in far])
    return out, dax, dtr, far


def aud_shells(dax, dtr, shells):
    out = {}
    for lab, lo, hi in shells:
        a = np.where((dax > lo) & (dax < hi))[0]; t = np.where((dtr > lo) & (dtr < hi))[0]
        out[lab] = dict(n_axis=int(len(a)), n_true=int(len(t)), same=bool(np.array_equal(a, t)))
    return out


def aud_sphere(c, A, ng, rho):
    """Points de grille dans la sphère (3D) de rayon rho autour de c : critère axe par axe (sphere_average) contre vraie image."""
    n = np.asarray(ng); reach = rho * np.linalg.norm(np.linalg.inv(A), axis=1); ax_ = []
    for i in range(3):
        h = int(np.ceil(reach[i] * n[i])) + 2; ic = int(np.rint(c[i] * n[i]))
        ax_.append(np.arange(n[i]) if 2 * h + 1 >= n[i] else np.mod(np.arange(ic - h, ic + h + 1), n[i]))
    I = np.meshgrid(*ax_, indexing="ij"); d = np.stack([I[k].ravel() / n[k] for k in range(3)], -1) - c
    rax = np.linalg.norm(cart_axis(d, A), axis=1); rtr = cart_true(d, A)[1]
    ma, mt = rax < rho, rtr < rho
    return dict(n_axis=int(ma.sum()), n_true=int(mt.sum()), same=bool(np.array_equal(ma, mt)), n_within_1e9_of_rho=int((np.abs(rax - rho) < 1e-9).sum()))


def aud_disc(c, A, ng, rho):
    """Disque dans le plan (inplane_disc_mask) de rayon rho autour de c : axe par axe contre vraie image, grille n1 × n2."""
    n1, n2 = int(ng[0]), int(ng[1])
    D1, D2 = np.meshgrid(np.arange(n1) / n1, np.arange(n2) / n2, indexing="ij")
    d = np.stack([D1.ravel() - c[0], D2.ravel() - c[1], np.zeros(n1 * n2)], -1)
    dd = red_axis(d); x = dd[:, 0] * A[0, 0] + dd[:, 1] * A[0, 1]; y = dd[:, 0] * A[1, 0] + dd[:, 1] * A[1, 1]
    ma = (x * x + y * y) < rho * rho; mt = cart_true(d, A)[1] < rho
    return dict(n_axis=int(ma.sum()), n_true=int(mt.sum()), same=bool(np.array_equal(ma, mt)))


def aud_radial(S, V, x_d, A_b):
    """
    Profil radial d'analyze_Ved.py (plan z = z_C, moyenne par anneaux de 0,05 Å jusqu'à |a1|/2 ; masqué : r >= 0,5 Å de chaque atome de la
    cellule avec lacune), refait sur le plan stocké dans ved_analysis.npz (carte recentrée, roulée en sens inverse) : axe par axe (porte :
    redonne {S}_rad et {S}_rad_masked) contre vraie image.
    """
    s_vac = np.asarray(V[f"{S}_s_vac"], float); M = np.asarray(V[f"{S}_map"]); n0, n1 = M.shape
    i0 = int(np.round(s_vac[0] * n0)); j0 = int(np.round(s_vac[1] * n1)); plane = np.roll(M, (-(n0 // 2 - i0), -(n1 // 2 - j0)), axis=(0, 1))
    i = np.arange(n0); j = np.arange(n1)
    S1, S2 = np.meshgrid((i / n0 - s_vac[0] + 0.5) % 1 - 0.5, (j / n1 - s_vac[1] + 0.5) % 1 - 0.5, indexing="ij")
    X = (S1 * A_b[0, 0] + S2 * A_b[0, 1]) * BOHR; Y = (S1 * A_b[1, 0] + S2 * A_b[1, 1]) * BOHR; R_ax = np.sqrt(X ** 2 + Y ** 2)
    d = np.stack([S1.ravel(), S2.ravel(), np.zeros(S1.size)], -1); R_tr = cart_true(d, A_b * BOHR)[1].reshape(S1.shape)
    xd = np.mod(x_d, 1.0); mask_ax = np.ones(plane.shape, bool); mask_tr = np.ones(plane.shape, bool)
    for s_at in xd:
        d1 = (i / n0 - s_at[0] + 0.5) % 1 - 0.5; d2 = (j / n1 - s_at[1] + 0.5) % 1 - 0.5; D1, D2 = np.meshgrid(d1, d2, indexing="ij")
        mask_ax &= np.sqrt(((D1 * A_b[0, 0] + D2 * A_b[0, 1]) * BOHR) ** 2 + ((D1 * A_b[1, 0] + D2 * A_b[1, 1]) * BOHR) ** 2) >= 0.5
        e = np.stack([D1.ravel(), D2.ravel(), np.zeros(D1.size)], -1); mask_tr &= (cart_true(e, A_b * BOHR)[1] >= 0.5).reshape(D1.shape)
    rmax = 0.5 * np.linalg.norm(A_b[:, 0]) * BOHR; edges = np.arange(0, rmax + 0.05, 0.05); rc = 0.5 * (edges[1:] + edges[:-1])

    def prof(R, m):
        cnt, _ = np.histogram(R[m], edges); sm, _ = np.histogram(R[m], edges, weights=plane[m]); return np.where(cnt > 0, sm / np.maximum(cnt, 1), np.nan), cnt
    out = dict(n_grid=[n0, n1], half_box_A=rmax, inscribed_parallelogram_A=float(rmax * np.sqrt(3) / 2), core_mask_same=bool(np.array_equal(mask_ax, mask_tr)))
    for lab, key, ma, mt in (("sans masque", "rad", R_ax < rmax, R_tr < rmax), ("masqué", "rad_masked", (R_ax < rmax) & mask_ax, (R_tr < rmax) & mask_tr)):
        pa, ca = prof(R_ax, ma); pt, ct = prof(R_tr, mt); ref = np.asarray(V[f"{S}_{key}"], float)
        gate = float(np.nanmax(np.abs(pa - ref))) if np.array_equal(np.isnan(pa), np.isnan(ref)) else float("inf")
        ch = np.where(ca != ct)[0]; dv = np.abs(pt - pa)
        out[lab] = dict(gate_vs_ved_analysis_max_abs_eV=gate, n_bins=int(len(rc)), n_bins_count_differs=int(len(ch)),
                        r_first_bin_differs_A=float(rc[ch[0]]) if len(ch) else None, n_points_axis=int(ma.sum()), n_points_true=int(mt.sum()),
                        max_abs_diff_meV=float(np.nanmax(dv) * 1e3) if len(ch) else 0.0, r_at_max_diff_A=float(rc[int(np.nanargmax(dv))]) if len(ch) else None,
                        bins=[dict(r_A=float(rc[b]), n_axis=int(ca[b]), n_true=int(ct[b]), axis_meV=float(pa[b] * 1e3), true_meV=float(pt[b] * 1e3)) for b in ch])
    return out


def aud_rgrid(D):
    """Étiquettes R de Mwk_to_Mwr / recenter_mwr (boîte MP D×D, réduction axe par axe) : distance cartésienne (unités de a) axe par axe contre vraie."""
    c = np.arange(D) - D // 2; R = np.array([(i, j) for i in c for j in c], float)
    imgs = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)], float) * D
    dax = np.linalg.norm(R @ A_UNIT, axis=1); L = np.linalg.norm((R[:, None, :] + imgs[None]) @ A_UNIT, axis=2); dtr = L.min(1)
    ties = (L <= dtr[:, None] + TIE).sum(1); rtr = np.linalg.norm(R[:, None, :] + imgs[None], axis=2).min(1); rax = np.linalg.norm(R, axis=1)
    diff = np.abs(dax - dtr) > TIE
    cl = {}
    for rc in range(5):
        a = rax <= rc + 1e-9; t = rtr <= rc + 1e-9; cart = dtr <= rc + 1e-9
        cl[str(rc)] = dict(n_axis_rednorm=int(a.sum()), n_true_rednorm=int(t.sum()), same=bool(np.array_equal(a, t)), n_cartesian_true=int(cart.sum()),
                          only_rednorm=[R[i].astype(int).tolist() for i in np.where(a & ~cart)[0]], only_cartesian=[R[i].astype(int).tolist() for i in np.where(cart & ~a)[0]])
    return dict(D=D, n_R=int(len(R)), n_dist_differs=int(diff.sum()), n_ties_true=int((ties > 1).sum()),
                min_axis_dist_where_differs_a=float(dax[diff].min()) if diff.any() else None, min_true_dist_where_differs_a=float(dtr[diff].min()) if diff.any() else None,
                max_axis_dist_a=float(dax.max()), max_true_dist_a=float(dtr.max()), max_shift_a=float((dax - dtr).max()), clusters=cl), R, dax, dtr


def aud_locality(S, L, D):
    """fig_locality_final / mwr_locality.npz : abscisse |R| (unités de a) axe par axe contre vraie ; porte : redonne dist et w de mwr_locality.npz."""
    out = {}
    if os.path.exists(os.path.join(WORK, "cache", f"Mwr_{S}.npz")):
        Z = np.load(os.path.join(WORK, "cache", f"Mwr_{S}.npz")); Rn = Z["Rn"][:, :2].astype(float); Mw = Z["Mwr"]
        i0 = int(np.where((Rn == 0).all(1))[0][0]); w = np.array([np.linalg.norm(Mw[:, i, :, i0]) for i in range(len(Rn))]); del Mw, Z
        da = np.linalg.norm(Rn @ A_UNIT, axis=1); o = np.argsort(da)
        imgs = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)], float) * D
        dt = np.linalg.norm((Rn[:, None, :] + imgs[None]) @ A_UNIT, axis=2).min(1)
        ref_d, ref_w = L[f"{S}_dense_dist"], L[f"{S}_dense_w"]
        gate = dict(max_abs_dist=float(np.abs(da[o] - ref_d).max()), max_rel_w=float(np.abs(w[o] - ref_w).max() / ref_w.max()))
        ch = np.abs(da - dt) > TIE; shown = (w >= 2e-5) & (w <= 60)
        out["dense"] = dict(D=D, gate_vs_mwr_locality=gate, n_points=int(len(Rn)), n_x_changes=int(ch.sum()), n_x_changes_shown=int((ch & shown).sum()),
                            min_x_axis_changed_a=float(da[ch].min()) if ch.any() else None, max_w_changed_eV=float(w[ch].max()) if ch.any() else None,
                            max_x_shift_a=float((da - dt)[ch].max()) if ch.any() else 0.0, max_x_axis_a=float(da.max()), max_x_true_a=float(dt.max()))
    return out


def aud_locality_coarse(S, L, N):
    """Grille grossière N×N de mwr_locality.npz : abscisse axe par axe -> vraie, par groupe d'égale distance axe par axe (w lus dans le npz)."""
    g, R, dax, dtr = aud_rgrid(N)
    ref_d, ref_w = L[f"{S}_coarse_dist"], L[f"{S}_coarse_w"]
    key = np.round(dax, 9); groups = {}
    for kx, t in zip(key, dtr):
        groups.setdefault(kx, set()).add(round(float(t), 9))
    ok_map = all(len(v) == 1 for v in groups.values())
    new = np.array([next(iter(groups[round(float(x), 9)])) for x in ref_d]) if ok_map else None
    gate = float(np.abs(np.sort(dax) - ref_d).max())
    ch = np.abs(new - ref_d) > TIE if new is not None else None
    return dict(N=N, gate_vs_mwr_locality_dist=gate, map_well_defined=ok_map, n_points=int(len(ref_d)),
                n_x_changes=int(ch.sum()) if ch is not None else None, n_x_changes_shown=int((ch & (ref_w >= 2e-5)).sum()) if ch is not None else None,
                changes=[dict(x_axis_a=float(ref_d[i]), x_true_a=float(new[i]), w_eV=float(ref_w[i])) for i in np.where(ch)[0]] if ch is not None else None)


def aud_interp_bound(C9):
    """
    D.1 (P2, M_W de production 9x9, étiquettes R brutes de Mwk_to_Mwr) : paires (R, R') dont l'étiquette diffère de la vraie image minimale
    autour du défaut R_d (égalités comprises) ; borne |ΔM(k', k)|_F <= 2 Σ ||M(R, R')||_F sur ces paires (hors grille k, la phase change).
    """
    Z = np.load(os.path.join(WORK, "cache", "Mwr_9x9.npz")); R = Z["R"][:, :2].astype(int); Rd = Z["R_d"][:2].astype(int); D = int(Z["MP"][0])
    W = np.sqrt((np.abs(Z["Mwr"]) ** 2).sum(axis=(0, 2))); del Z
    rel_ax = ((R - Rd + D // 2) % D) - D // 2
    imgs = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]) * D
    Lr = np.linalg.norm((rel_ax[:, None, :] + imgs[None]).astype(float) @ A_UNIT, axis=2); m = Lr.min(1); ties = (Lr <= m[:, None] + TIE).sum(1)
    R_true = Rd + rel_ax + imgs[Lr.argmin(1)]
    out = dict(D=D, R_d=Rd.tolist(), n_R=int(len(R)), onsite_norm_eV=float(W.max()))

    def bound(c, Wm):
        return float(2 * (Wm[c].sum() + Wm[~c][:, c].sum()))
    for lab, c in (("brut_vs_vraie", (R != R_true).any(1) | (ties > 1)), ("axe_autour_Rd_vs_vraie", (Rd + rel_ax != R_true).any(1) | (ties > 1))):
        out[lab] = dict(n_R_changed=int(c.sum()), n_ties=int(((ties > 1) & c).sum()), bound_eV=bound(c, W), max_block_eV=float(max(W[c].max(), W[:, c].max())) if c.any() else 0.0,
                        min_true_dist_changed_a=float(m[c].min()) if c.any() else None)
    f = os.path.join(WORK, "cache", "Fw_9x9.npz")
    if os.path.exists(f):
        F = np.load(f); Wf = np.sqrt((np.abs(F["Fw"]) ** 2).sum(axis=(0, 2))); assert np.array_equal(F["R"][:, :2].astype(int), R); del F
        for lab in ("brut_vs_vraie", "axe_autour_Rd_vs_vraie"):
            c = (R != R_true).any(1) | (ties > 1) if lab == "brut_vs_vraie" else (Rd + rel_ax != R_true).any(1) | (ties > 1)
            out[lab]["bound_C_Fw_eV"] = abs(C9) * bound(c, Wf)
    return out


def aud_interp_D1(C9, cfg):
    """
    Carte D.1 refaite (sommes de Fourier recopiées, aucune routine de production) avec trois conventions d'étiquettes R, bra et ket :
      brut    étiquettes de Mwk_to_Mwr (boîte centrée sur l'origine) = convention du pilote D.1 (porte : redonne les moyennes de d_results.json) ;
      axe_Rd  réduction axe par axe autour du défaut R_d (étiquettes de recenter_mwr + R_d) ;
      vraie   vraie image minimale autour de R_d, poids 1/n sur les n images à égalité (comme les dégénérescences de Wigner-Seitz).
    M_W(k', k) = Σ_{R,R'} e^{−2πi k'·R} M_W(R, R') e^{2πi k·R'} ; bandes π par le poids p_z des vecteurs propres de H(k) = Σ e^{2πi k·R} H(R)/ndeg.
    """
    Z = np.load(os.path.join(WORK, "cache", "Mwr_9x9.npz")); R = Z["R"].astype(float); Rd = Z["R_d"].astype(float); D = int(Z["MP"][0]); Mw = Z["Mwr"]; del Z
    F = np.load(os.path.join(WORK, "cache", "Fw_9x9.npz"))["Fw"]
    dd = json.load(open(os.path.join(WORK, "d", "d_results.json"))); meta = dd["D1_meta"]; k = np.array(meta["k_red"], float)
    dp = dense_paths(cfg, "9x9"); A_uc, _ = qe_io.get_A_volume(dp["uc"]); A_A = A_uc * BOHR; A_cell = float(np.linalg.norm(np.cross(A_A[:, 0], A_A[:, 1])))
    Bc = 2 * np.pi * np.linalg.inv(A_A).T; b1 = float(np.linalg.norm(Bc[:, 0]))
    paths = wannier_provenance.load_wannier_checked(dp["manifest"]); Hwr, Rw, nd = read_w90_HR(paths["tb"]); Rw = np.asarray(Rw, float)

    def eig(kk):
        Hk = np.einsum("kr,rwv->kwv", np.exp(2j * np.pi * (kk @ Rw.T)), Hwr / np.asarray(nd)[:, None, None], optimize=True)
        return np.linalg.eigh(Hk)
    n = 240; i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); kp = np.stack([i.ravel() / n, j.ravel() / n, np.zeros(n * n)], -1)
    Ep, Up = eig(kp); Ek, Uk = eig(k[None])
    wz = lambda U: np.abs(U[..., WF_PZ_A, :]) ** 2 + np.abs(U[..., WF_PZ_B, :]) ** 2
    top = np.sort(np.argsort(-wz(Up), axis=-1)[..., :2], axis=-1); pv, pc = top[:, 0], top[:, 1]
    tk = np.sort(np.argsort(-wz(Uk), axis=-1)[..., :2], axis=-1)[0]
    kd = lambda kk, k0: np.linalg.norm(red_axis(kk - k0), axis=1)
    exK, exKp = kd(kp, K_RED) < 1e-9, kd(kp, KP_RED) < 1e-9
    kc = (Bc[:2, :2] @ kp[:, :2].T).T; imgs_k = np.array([[a, b] for a in (-1, 0, 1) for b in (-1, 0, 1)]) @ Bc[:2, :2].T
    disk = lambda K0: np.min(np.linalg.norm(kc[:, None, :] - (Bc[:2, :2] @ K0[:2])[None, None, :] + imgs_k[None], axis=2), axis=1) <= 0.05 * b1
    dK, dKp = disk(K_RED) & ~exK, disk(KP_RED) & ~exKp
    # étiquettes : liste d'images (et poids) par R
    rel_ax = ((R[:, :2] - Rd[:2] + D // 2) % D) - D // 2
    im = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], float) * D
    Lr = np.linalg.norm((rel_ax[:, None, :] + im[None]) @ A_UNIT, axis=2); mn = Lr.min(1)
    conv = {"brut": [(R[:, :2], np.ones(len(R)))], "axe_Rd": [(Rd[:2] + rel_ax, np.ones(len(R)))]}
    tied = Lr <= mn[:, None] + TIE; nt = tied.sum(1)
    conv["vraie"] = [(Rd[:2] + rel_ax + im[s], tied[:, s] / nt) for s in range(len(im))]

    def phases(kk, labels, sign):
        return sum(w[None, :] * np.exp(sign * 2j * np.pi * (kk[:, :2] @ Rl.T)) for Rl, w in labels)
    res = {}; maps = {}
    for cname, labels in conv.items():
        Pk = phases(k[None], labels, +1)[0]                                                     # (nR,)
        for var, M in (("brut", Mw), ("exact", None)):
            X = np.einsum("wrWR,R->wrW", Mw, Pk, optimize=True)
            if var == "exact":
                X = X - C9 * np.einsum("wrWR,R->wrW", F, Pk, optimize=True)
            Mk = np.zeros((len(kp), NW, NW), complex)
            for s in range(0, len(kp), 9600):
                Pb = phases(kp[s:s + 9600], labels, -1); Mk[s:s + 9600] = np.einsum("kr,wrW->kwW", Pb, X, optimize=True)
            Mb = np.einsum("kwb,kwW,WB->kbB", Up.conj(), Mk, Uk[0], optimize=True)
            Vv = A_cell * np.abs(Mb[np.arange(len(kp)), pv, tk[0]]); Vc = A_cell * np.abs(Mb[np.arange(len(kp)), pc, tk[1]])
            Vv[exK | exKp] = np.nan; Vc[exK | exKp] = np.nan
            maps[(cname, var)] = (Vv, Vc)
            res[f"{cname} | {var}"] = {f"{b} | {lab}": dict(mean=float(np.nanmean(V[m])), min=float(np.nanmin(V[m])), max=float(np.nanmax(V[m])))
                                       for lab, m in (("K", dK), ("K'", dKp)) for b, V in (("valence", Vv), ("conduction", Vc))}
    gate = max(abs(res[f"brut | {var}"][key]["mean"] - dd["D1"][var][key]["mean"]) for var in ("brut", "exact") for key in dd["D1"][var])
    diffs = {}
    for a_, b_ in (("brut", "vraie"), ("axe_Rd", "vraie"), ("brut", "axe_Rd")):
        for var in ("brut", "exact"):
            dv = np.nanmax(np.abs(np.r_[maps[(a_, var)][0] - maps[(b_, var)][0], maps[(a_, var)][1] - maps[(b_, var)][1]]))
            dm = max(abs(res[f"{a_} | {var}"][key]["mean"] - res[f"{b_} | {var}"][key]["mean"]) for key in res[f"{a_} | {var}"])
            dd_ = max(np.nanmax(np.abs((maps[(a_, var)][t] - maps[(b_, var)][t])[m])) for t in (0, 1) for m in (dK, dKp))
            diffs[f"{a_} -> {b_} | {var}"] = dict(max_abs_map_eVA2=float(dv), max_abs_in_disks_eVA2=float(dd_), max_abs_disk_mean_eVA2=float(dm))
    return dict(gate_vs_d_results_max_abs_mean_eVA2=float(gate), n_ties_R=int((nt > 1).sum()), A_cell_A2=A_cell, disk_means=res, differences=diffs)


def aud_kgrids(B2):
    """Point de grille le plus proche de K, K', Γ : argmin de la norme réduite axe par axe (kindex, kdist) contre vraie distance cartésienne."""
    grids = dict(coarse=[5, 6, 7, 8, 9, 10, 11, 12], dense=[24, 25, 27, 28, 32], sortie=[240], interne=[300, 450, 600, 900], epw_q=[24])
    out = {}
    for fam, ns in grids.items():
        for n in ns:
            i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); k = np.stack([i.ravel() / n, j.ravel() / n], -1)
            for lab, tg in (("K", K_RED[:2]), ("K'", KP_RED[:2]), ("Γ", np.zeros(2))):
                d = red_axis(k - tg); coded = np.linalg.norm(d, axis=1); ic = int(np.argmin(coded))
                imgs = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], float)
                Lc = np.linalg.norm((d[:, None, :] + imgs[None]) @ B2.T, axis=2).min(1); tset = np.where(Lc <= Lc.min() + 1e-12)[0]
                out[f"{fam} {n} {lab}"] = dict(n=n, coded_ij=[int(round(v * n)) for v in k[ic]], n_coded_ties=int((coded <= coded.min() + 1e-12).sum()),
                                               on_grid=bool(Lc.min() < 1e-9), dist_coded_Ainv=float(Lc[ic]), dist_true_min_Ainv=float(Lc.min()),
                                               n_true_nearest=int(len(tset)), coded_in_true_set=bool(ic in tset),
                                               true_set_ij=[[int(round(v * n)) for v in k[t]] for t in tset])
    return out


def aud_crowns(B2, n=240, xmax=400):
    """Couronnes de crown_table (R9) : décalage axe par axe à K, K', puis norme de Löschian, contre la vraie image (9 images)."""
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij"); k = np.stack([i.ravel() / n, j.ravel() / n], -1); b1 = B2[:, 0]
    xa = np.full(len(k), np.inf); xt = np.full(len(k), np.inf); imgs = np.array([(a, b) for a in (-1, 0, 1) for b in (-1, 0, 1)], float)
    for Kv in (K_RED[:2], KP_RED[:2]):
        dk = red_axis(k - Kv); ij = np.rint(dk * n); q = ij @ B2.T; xa = np.minimum(xa, (q * q).sum(1) / (b1 @ b1))
        qt = (np.rint(dk * n)[:, None, :] + n * imgs[None]) @ B2.T; xt = np.minimum(xt, ((qt * qt).sum(2) / (b1 @ b1)).min(1))
    m = xt <= xmax
    return dict(n=n, x_max=xmax, n_k=int(m.sum()), n_differs=int((np.abs(xa[m] - xt[m]) > 1e-6).sum()))


def aud_transplant():
    """R7c (make_inputs_r7_relax.py) : sites 12x12 transplantés, |r| < R_CUT = 14,5 Å avec r = d − rint(d) (axe par axe), contre la vraie image."""
    import re
    scf = os.path.join(GQ, "defects", "super_cell", "12x12", "defective", "scf.in")
    rel = os.path.join(GQ, "defects", "super_cell_relaxed", "series", "12x12", "nspin1", "relax.out")
    txt = open(scf).read(); _, rest = txt.split("CELL_PARAMETERS bohr"); cell_txt, pos_txt = rest.split("ATOMIC_POSITIONS crystal"); pos_txt = pos_txt.split("K_POINTS")[0]
    A_rows = np.array([[float(v) for v in l.split()] for l in cell_txt.strip().splitlines()]) * BOHR
    X0 = np.array([[float(v) for v in l.split()[1:4]] for l in pos_txt.strip().splitlines()])
    blk = re.search(r"Begin final coordinates(.*?)End final coordinates", open(rel, errors="ignore").read(), re.S).group(1).split("ATOMIC_POSITIONS (crystal)")[1]
    Xf = np.array([[float(v) for v in l.split()[1:4]] for l in blk.strip().splitlines() if re.match(r"\s*C\b", l)])
    s12 = np.array([(5 + 2 / 3) / 12, (5 + 2 / 3) / 12, 0.0]); R_CUT = 14.5
    d = X0 - s12; r_ax = np.linalg.norm((d - np.rint(d)) @ A_rows, axis=1); r_tr = cart_true(d, A_rows.T)[1]
    u = Xf - X0; un = np.linalg.norm((u - np.rint(u)) @ A_rows, axis=1)
    sa, st = r_ax < R_CUT, r_tr < R_CUT; miss = st & ~sa
    L12 = float(np.linalg.norm(A_rows[0]))
    out = dict(n_atoms=len(X0), L12_A=L12, R_CUT_A=R_CUT, inscribed_parallelogram_A=L12 * np.sqrt(3) / 4, n_sel_axis=int(sa.sum()), n_sel_true=int(st.sum()),
               n_missed=int(miss.sum()), n_extra=int((sa & ~st).sum()), missed=[dict(atom_1based=int(i) + 1, r_true_A=float(r_tr[i]), r_axis_A=float(r_ax[i]), u_A=float(un[i]))
                                                                                  for i in np.where(miss)[0]],
               u_max_missed_A=float(un[miss].max()) if miss.any() else 0.0, u_max_selected_12_14p5_A=float(un[sa & (r_ax >= 12.0)].max()))
    return out


def aud_zone():
    """
    Écart 3 : coins d'analyze_M.py §1 (map_corners de M_analysis.npz, formule de la l. 56) contre les points K de la cellule QE (a1, a2 à 60°,
    b1, b2 à 120°), K = (2/3, 1/3) de production.json et ses images. Un coin de la zone est un sommet de Wigner-Seitz : |p| = |p − G| pour
    exactement deux G ≠ 0, et |p| <= |p − G| pour tous.
    """
    Z = np.load(os.path.join(PROJ, "results", "M2", "M_analysis.npz"), allow_pickle=True)
    B = np.asarray(Z["map_B"])[:2, :2]; C = np.asarray(Z["map_corners"]); Kf = np.asarray(Z["map_K"])
    red_now = [(1 / 3, 1 / 3), (1 / 3, -2 / 3), (-2 / 3, 1 / 3), (-1 / 3, -1 / 3), (-1 / 3, 2 / 3), (2 / 3, -1 / 3)]        # analyze_M.py l. 56
    red_new = [(2 / 3, 1 / 3), (1 / 3, 2 / 3), (-1 / 3, 1 / 3), (-2 / 3, -1 / 3), (-1 / 3, -2 / 3), (1 / 3, -1 / 3)]       # K, K' et images
    G = np.array([B @ np.array([i, j]) for i in range(-2, 3) for j in range(-2, 3) if (i, j) != (0, 0)])

    def info(red):
        rows = []
        for c in red:
            p = B @ np.array(c); dG = np.linalg.norm(p - G, axis=1); d0 = float(np.linalg.norm(p)); neq = int((np.abs(dG - d0) < 1e-9).sum())
            rows.append(dict(red=[float(v) for v in c], cart_Ainv=[float(v) for v in p], norm_Ainv=d0, n_G_equidistant=neq,
                             is_WS_vertex=bool(neq == 2 and d0 <= dG.min() + 1e-9)))
        return rows
    now, new = info(red_now), info(red_new)
    P_now = np.array([r["cart_Ainv"] for r in now]); P_now = P_now[np.argsort(np.arctan2(P_now[:, 1], P_now[:, 0]))]
    K_cfg = B @ np.array(load_production(verbose=False)["K_red"][:2])
    # même vérification pour les images de K de scripts/_bands.py (SYM_POINTS, étiquettes des chemins de bandes) : norme de la vraie image
    import _bands
    im = np.array([(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)], float)
    bands_K = [dict(red=[float(v) for v in s[:2]], norm_true_image_Ainv=float(np.linalg.norm((red_axis(np.array(s[:2], float))[None] + im) @ B.T, axis=1).min()))
               for s in _bands.SYM_POINTS["K"]]
    return dict(stored_equals_formula=bool(np.allclose(P_now, C, atol=1e-12)), K_cfg_cart_Ainv=[float(v) for v in K_cfg], K_norm_Ainv=float(np.linalg.norm(K_cfg)),
                map_K_folded_Ainv=[float(v) for v in Kf], b_norm_Ainv=float(np.linalg.norm(B[:, 0])), b1_dot_b2_over_b2=float(B[:, 0] @ B[:, 1] / (B[:, 0] @ B[:, 0])),
                corners_now=now, corners_proposed=new, bands_SYM_POINTS_K=bands_K)


def aud_label():
    """Écart 1 : ligne « convention intensive » d'analyze_M.py (texte lu), ligne de M_tests_summary.csv, étendue recalculée (M_analysis.npz) sur les
    tailles du libellé (5, 7, 8, 9) et sur les tailles effectivement traitées par analyze_M (SIZES : M_dense et M_ed présents)."""
    import csv
    src = open(os.path.join(PROJ, "scripts", "analyze_M.py")).read().splitlines()
    ln = [i + 1 for i, l in enumerate(src) if "convention intensive (cellule unitaire)" in l]
    Z = np.load(os.path.join(PROJ, "results", "M2", "M_analysis.npz"), allow_pickle=True)
    sc = {S: float(Z[f"scale_{S}_dense"]) for S in ("5x5", "6x6", "7x7", "8x8", "9x9", "12x12") if f"scale_{S}_dense" in Z.files}
    spread = lambda v: float((max(v) - min(v)) / np.mean(v))
    row = [r for r in csv.reader(open(os.path.join(PROJ, "results", "M2", "M_tests_summary.csv"))) if "convention intensive" in r[0]]
    return dict(lines=ln, text=[src[i - 1].strip() for i in ln], scale_dense_eV=sc, spread_label_sizes_5789=spread([sc[S] for S in ("5x5", "7x7", "8x8", "9x9")]),
                spread_all_sizes=spread(list(sc.values())), csv_rows=row)


AUD_PARTS = ["atoms", "r1", "rgrids", "d1", "kgrids", "transplant", "zone", "label"]


def cmd_audit(a):
    """Audit : réduction axe par axe contre vraie image minimale sur les géométries réelles (voir audit_image_minimale.md)."""
    t0 = time.time(); d = ensure("audit"); f_out = os.path.join(d, "audit_results.json")
    parts = AUD_PARTS if a.parts == "all" else a.parts.split(",")
    out = json.load(open(f_out)) if (a.parts != "all" and os.path.exists(f_out)) else dict(sizes={})          # partie rejouée : fusion
    out.setdefault("runs", []).append(dict(head=git_head(), date=time.strftime("%Y-%m-%d %H:%M"), parts=parts))
    cfg = load_production(verbose=False); G9 = aud_geom(9)
    V = np.load(os.path.join(PROJ, "results", "M2", "ved_analysis.npz"), allow_pickle=True)
    r7 = json.load(open(os.path.join(GQ, "defects", "R7_tailles_3m", "d1_8pts", "d1_results.json")))["sizes"]
    r5 = json.load(open(os.path.join(GQ, "defects", "R5_base_vs_M", "c", "c_results.json")))["sizes"]
    r4 = json.load(open(os.path.join(GQ, "defects", "R4_quasi_lie", "d5", "d5_results.json")))
    a1 = json.load(open(os.path.join(WORK, "a", "a1_results.json")))["sizes"]
    for N in (AUD_SIZES if "atoms" in parts else []):
        S = f"{N}x{N}"; G = aud_geom(N); A = G["A"]; t1 = time.time()
        vac, s_vac = aud_vacancy(G["x_p"], G["x_d"], A)
        far, dax, dtr, fset = aud_far(G["x_d"], G["x_p"], s_vac, A)
        dpa = np.linalg.norm(cart_axis(G["x_p"] - s_vac, A), axis=1); dpt = cart_true(G["x_p"] - s_vac, A)[1]
        rec = dict(ngfft=list(G["ng"]), a_sc_A=float(np.linalg.norm(A[:, 0])), vacancy=vac, far=far,
                   shells_d=aud_shells(dax, dtr, (("< 1,8 Å", -1.0, 1.8), ("2,2–2,7 Å", 2.2, 2.7))),
                   shells_p_analyze_Ved=aud_shells(dpa, dpt, (("0,5–1,6 Å", 0.5, 1.6),)),
                   R_CUT_14p5_target=dict(n_axis=int((dax < 14.5).sum()), n_true=int((dtr < 14.5).sum()), same=bool(np.array_equal(dax < 14.5, dtr < 14.5))))
        # sphères de Lu (0,5 et 1,0 Å) autour de l'atome axe par axe et des atomes vraiment les plus loin ; porte : nombre de points de R7 (sphère d)
        sph = {}
        for rho in (0.5, 1.0):
            ia = far["i_axis_1based"] - 1
            sph[str(rho)] = dict(axis_atom=aud_sphere(np.mod(G["x_d"][ia], 1.0), A, G["ng"], rho),
                                 far_true=[aud_sphere(np.mod(G["x_d"][i], 1.0), A, G["ng"], rho) for i in fset])
        if str(N) in r7:
            m7 = r7[str(N)]["alignment"]["means"]
            sph["gate_R7_npts_d"] = {k: dict(R7=int(v[2]), audit_axis=sph[k]["axis_atom"]["n_axis"], ok=bool(int(v[2]) == sph[k]["axis_atom"]["n_axis"])) for k, v in m7.items()}
        # toutes les sphères (P1, A.1) pour les tailles de R9
        if S in SIZES_A:
            allsp = [aud_sphere(np.mod(x, 1.0), A, G["ng"], rho) for rho in (0.5, 1.0) for x in G["x_d"]]
            sph["all_atoms"] = dict(n_spheres=len(allsp), n_differ=int(sum(not s["same"] for s in allsp)), n_boundary_points=int(sum(s["n_within_1e9_of_rho"] for s in allsp)))
        rec["spheres"] = sph
        rec["discs"] = {str(r): aud_disc(s_vac, A, G["ng"], r) for r in (1.0, 2.0)}
        if f"{S}_map" in V.files:
            rec["radial_analyze_Ved"] = aud_radial(S, V, G["x_d"], G["A_b"])
        # valeur de Lu : sphère 1,0 Å autour de l'atome axe par axe (publiée) contre atomes vraiment les plus loin (profils P1 d'A.1)
        pf = os.path.join(WORK, "a", f"a1_profiles_{S}.npz")
        if os.path.exists(pf):
            P = np.load(pf); ia = far["i_axis_1based"] - 1
            lu = dict(gate_dist_true_max_abs_A=float(np.abs(P["dist"] - dtr).max()), axis_atom_meV=float(P["shift10"][ia] * 1e3), published_meV=LU_PUBLISHED_MEV.get(S),
                      far_true_meV=[float(P["shift10"][i] * 1e3) for i in fset], axis_atom_05_meV=float(P["shift05"][ia] * 1e3),
                      far_true_05_meV=[float(P["shift05"][i] * 1e3) for i in fset], plateau_i_meV=float(a1[S]["C_i_eV"] * 1e3))
            lu["delta_mean_meV"] = float(np.mean(lu["far_true_meV"]) - lu["axis_atom_meV"]); lu["delta_mean_05_meV"] = float(np.mean(lu["far_true_05_meV"]) - lu["axis_atom_05_meV"])
            rec["Lu"] = lu
        # atome de Lu enregistré par R4 D5, R5 C, R7 D1 (même règle axe par axe) : porte
        rec["recorded_far_atom"] = dict(R4_D5=r4.get(S, {}).get("far_atom_1based"), R5_C=r5.get(str(N), {}).get("alignment", {}).get("i_far_1based"),
                                        R7_D1=r7.get(str(N), {}).get("alignment", {}).get("i_far_1based"))
        qe = {}
        for src, J in (("R5_C", r5), ("R7_D1", r7)):
            if str(N) in J:
                r = J[str(N)]; qe[src] = dict(shift_10_meV=float(r["alignment"]["shift_10"] * 1e3), pi_x=(r["pi_state"] or {}).get("x"),
                                              sigma_x=[q["x"] for q in (r["sigma_doublet"] or [])])
        rec["qe_levels"] = qe
        rec["elapsed_s"] = time.time() - t1
        out["sizes"][S] = rec
        log(f"[audit] {S} : lacune {'identique' if vac['same'] else 'DIFFÈRE'} ; atome de Lu axe par axe {far['i_axis_1based']} "
            f"({far['d_axis_A']:.2f} Å axe, {far['d_true_of_axis_atom_A']:.2f} Å vrai, rang {far['rank_true_of_axis_atom']}) ; vrais plus loin {far['far_true_1based']} "
            f"à {far['r_max_true_A']:.2f} Å ; {time.time()-t1:.0f} s")
    # R1 9x9 relaxée (R4 D1 : far_atom_alignment sur les positions relaxées)
    if "r1" in parts:
        out["R1"] = {}
    R1_SAVES = [("nspin1", f"{QE_TMP}/vacancy_relaxed/nspin1/vac_9x9_relax_nspin1.save"), ("nspin2", f"{QE_TMP}/vacancy_relaxed/nspin2/vac_9x9_relax_nspin2.save")]
    for tag, sv in (R1_SAVES if "r1" in parts else []):
        A_r, _ = qe_io.get_A_volume(sv); assert np.allclose(A_r, G9["A_b"]); xr = qe_io.get_x_red(sv)
        vac, s_vac = aud_vacancy(G9["x_p"], xr, G9["A"]); far, dax, dtr, fset = aud_far(xr, G9["x_p"], s_vac, G9["A"])
        out["R1"][tag] = dict(vacancy=vac, far=far, shells=aud_shells(dax, dtr, (("< 1,8 Å", -1.0, 1.8), ("2,2–2,7 Å", 2.2, 2.7))),
                              spheres={str(r): dict(axis_atom=aud_sphere(np.mod(xr[far["i_axis_1based"] - 1], 1.0), G9["A"], G9["ng"], r)) for r in (0.5, 1.0)})
        log(f"[audit] R1 {tag} : atome axe par axe {far['i_axis_1based']} ({far['d_axis_A']:.2f} / {far['d_true_of_axis_atom_A']:.2f} Å), vrais {far['far_true_1based']}")
    # étiquettes R (Wannier) : grilles denses de production et grossières N×N ; fig_locality_final
    L = np.load(os.path.join(PROJ, "results", "M2", "mwr_locality.npz"))
    if "rgrids" in parts:
        out["R_grids"] = {}
    for S in (SIZES_A if "rgrids" in parts else []):
        D = int(cfg["dense"][S]["D"]); n = size_n(S)
        out["R_grids"][S] = dict(dense=aud_rgrid(D)[0], coarse=aud_rgrid(n)[0], locality=aud_locality(S, L, D))
        if f"{S}_coarse_dist" in L.files:
            out["R_grids"][S]["locality"]["coarse"] = aud_locality_coarse(S, L, n)
        log(f"[audit] R {S} : dense D={D} {out['R_grids'][S]['dense']['n_dist_differs']} étiquettes à distance changée ; "
            f"amas R_cut 3 identiques {out['R_grids'][S]['dense']['clusters']['3']['same']}")
    if "d1" in parts:
        out["interp_D1"] = aud_interp_bound(float(a1["9x9"]["C_retenu_eV"]))
        out["interp_D1"]["carte"] = aud_interp_D1(float(json.load(open(os.path.join(WORK, "d", "d_results.json")))["D1_meta"]["C9_eV"]), cfg)
        log(f"[audit] D.1 : porte {out['interp_D1']['carte']['gate_vs_d_results_max_abs_mean_eVA2']:.2e} eV Å² ; "
            + " ; ".join(f"{k_} {v['max_abs_in_disks_eVA2']:.3f}" for k_, v in out["interp_D1"]["carte"]["differences"].items()))
    # grilles k : point le plus proche de K, K', Γ ; couronnes de crown_table
    A9 = G9["A"] / 9.0; B2 = (2 * np.pi * np.linalg.inv(A9).T)[:2, :2]
    if "kgrids" in parts:
        out["k_grids"] = aud_kgrids(B2); out["crowns_240"] = aud_crowns(B2)
    if "transplant" in parts:
        out["transplant_R7c"] = aud_transplant()
    if "zone" in parts:
        out["zone_analyze_M"] = aud_zone()
    if "label" in parts:
        out["label_analyze_M"] = aud_label()
    out["runs"][-1]["elapsed_s"] = time.time() - t0
    save_json(f_out, out)
    log(f"[audit] écrit audit/audit_results.json en {time.time()-t0:.0f} s")


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
    p = sub.add_parser("a3tables"); p.add_argument("--cn", choices=["lu", "plateau"], default=None)
    sub.add_parser("r0")
    sub.add_parser("synth")
    p = sub.add_parser("c"); p.add_argument("--redo", action="store_true")
    sub.add_parser("d")
    sub.add_parser("dfig")
    sub.add_parser("cfig")
    p = sub.add_parser("audit"); p.add_argument("--parts", default="all", help="all ou liste parmi " + ",".join(AUD_PARTS) + " (fusion dans le json existant)")
    a = ap.parse_args()
    global CN
    CN = getattr(a, "cn", None)                                                  # clôture : C_N Lu ou plateau (a3, a3pole, b, btables)
    {"r0": cmd_r0, "synth": cmd_synth, "a3tables": cmd_a3tables, "a2c": cmd_a2c, "a2d": cmd_a2d, "a2dpost": cmd_a2dpost, "a0": cmd_a0, "a1": cmd_a1, "a3": cmd_a3, "a3pole": cmd_a3pole, "b": cmd_b,
     "btables": cmd_btables, "c": cmd_c, "d": cmd_d, "dfig": cmd_dfig, "cfig": cmd_cfig, "audit": cmd_audit}[a.cmd](a)


if __name__ == "__main__":
    main()
