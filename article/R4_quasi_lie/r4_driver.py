#!/usr/bin/env python
"""
r4_driver.py -- pilote unique de la campagne R4 (diagnostic de l'état quasi-lié de la lacune 9x9).

Sous-commandes (une par job) : prep (J1), d6 (J2), d1 / d5 / r (J3), d4 (J4), d3 (J5).
Données de production en lecture seule ; sorties dans ce répertoire : prep/, d6/, d1/, d4/, d3/, d5/, r/,
cache/ (M_W recentré, blocs V_loc, g0 par grille), fig/. Journal : r4_log.txt.
Routines de production réutilisées telles quelles (config, porte de jauge, matrix_io, rotation, recentrage,
extract_V_loc, local_green_batch, Hwr_to_Hwk, compute_M_NL) ; fonctions nouvelles dans src/electron_defect_interaction/
(io/qe_gamma_io, io/projwfc_io, defects/alignment, defects/many_body/pole_criterion, defects/many_body/tb_models,
wannier/supercell_fold). Aucun fichier de production modifié.
"""
import os
import sys
import json
import time
import argparse

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                   # .../graphene/qe
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
sys.path.insert(0, os.path.join(PROJ, "src"))
os.chdir(PROJ)

from electron_defect_interaction.config import load_production, dense_paths, HA2EV
from electron_defect_interaction.io import qe_io, matrix_io, wannier_provenance
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _infer_mp_grid, _match_kpoint_order
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.defects.many_body import local_tmatrix as lt
from electron_defect_interaction.defects.many_body import pole_criterion as pc
from electron_defect_interaction.defects.many_body import tb_models as tb
from electron_defect_interaction.wannier import supercell_fold as sf
from electron_defect_interaction.io import qe_gamma_io as qg
from electron_defect_interaction.io import projwfc_io as pj
from electron_defect_interaction.defects import alignment as al

BOHR = 0.529177210903
SCRATCH = "/home/gregb26/links/scratch/qe_tmp"
RES = os.path.join(PROJ, "results", "M")
SC_D = f"{SCRATCH}/defect_9x9_d/defect_9x9_d.save"
SC_P = f"{SCRATCH}/defect_9x9_p/defect_9x9_p.save"
UC9 = f"{SCRATCH}/defect_unit_cell_9x9/defect_unit_cell_9x9.save"
R1 = f"{GQ}/defects/super_cell_relaxed/9x9"
R1_N1 = f"{SCRATCH}/vacancy_relaxed/nspin1/vac_9x9_relax_nspin1.save"
R1_N2 = f"{SCRATCH}/vacancy_relaxed/nspin2/vac_9x9_relax_nspin2.save"
VKS = dict(d=f"{GQ}/defects/super_cell/9x9/defective/Vks_9x9_d", p=f"{GQ}/defects/super_cell/9x9/pristine/Vks_9x9_p",
           n1=f"{R1}/nspin1/Vks_R1_nspin1", n2u=f"{R1}/nspin2/Vks_R1_nspin2_up", n2d=f"{R1}/nspin2/Vks_R1_nspin2_dw")
UPF = f"{UC9}/C.upf"
N_SC = 9; N_CELLS = 81
U_LIST = [3.0, 5.0, 6.617, 8.0, 12.0, 20.0, 27.0, 40.0, 80.0, 1e3]
# valeurs attendues à 1200² (prompt R4, D6.1) : U -> (racine de Re lambda près de E_D, argmin |det|, |det| min)
EXPECTED_1200 = {3.0: (None, -1.28, 0.71), 5.0: (None, -0.90, 0.62), 6.617: (-1.88, -0.71, 0.58), 8.0: (-1.31, -0.59, 0.55),
                 12.0: (-0.66, -0.39, 0.50), 20.0: (-0.31, -0.22, 0.46), 27.0: (-0.21, -0.16, 0.46), 40.0: (-0.13, -0.10, 0.49),
                 80.0: (-0.055, -0.044, 0.64), 1e3: (-0.004, -0.004, 6.1)}
FAIL_U = [8.0, 12.0, 20.0, 27.0, 40.0]
WF_SIGMA = [0, 1, 2]; WF_PI = [3, 4]; WF_PZ_A = 3; WF_PZ_B = 4; NW = 5
NN_CELLS = [(0, 0, 0), (-1, 0, 0), (0, -1, 0)]          # mailles des trois voisins B de A (H(R)[3,4] != 0)
K_RED = np.array([2 / 3, 1 / 3, 0.0]); KP_RED = np.array([1 / 3, 2 / 3, 0.0])
RADII_A = (0.5, 1.0)                                    # rayons de moyenne (Å) pour l'alignement Lu 2019
R_W2, R_W1 = 2.0, 1.0                                   # rayons (Å) du poids en plan près de la lacune


# ----------------------------------------------------------------------------------------------- utilitaires
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "r4_log.txt"), "a") as f:
        f.write(line + "\n")


def ensure(sub):
    d = os.path.join(WORK, sub); os.makedirs(d, exist_ok=True); return d


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        return float(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.complexfloating, complex)):
        return [float(np.real(x)), float(np.imag(x))]
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    return x


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(jsonable(d), f, indent=1, ensure_ascii=False)
    log(f"saved {os.path.relpath(path, WORK)}")


def egrid(E_D, lo, hi, eta, npe):
    de = eta / npe
    return np.arange(E_D + lo, E_D + hi + de / 2, de)


def setup(verbose=True):
    """Chaîne de production : config figée, porte de jauge, H(R), grille dense, U / U_dis, E_D Wannier."""
    cfg = load_production(verbose=verbose); dp = dense_paths(cfg, "9x9")
    paths = wannier_provenance.load_wannier_checked(dp["manifest"]); log(f"[gauge] provenance OK: {dp['manifest']}")
    Hwr, Rw, nd = read_w90_HR(paths["tb"])
    k27 = qe_io.get_k_red(dp["uc"]); MP = _infer_mp_grid(k27)
    U, kU = read_w90_mat(paths["u"]); U = U[_match_kpoint_order(kU, k27)]
    Ud, kUd = read_w90_mat(paths["u_dis"]); Ud = Ud[_match_kpoint_order(kUd, k27)]
    _, E_ref, _ = Hwr_to_Hwk(Hwr, Rw, lt.mp_grid(90, 90, 1), ndegen=nd)
    gap = E_ref[:, 4] - E_ref[:, 3]; iD = int(np.argmin(gap)); E_D = float(0.5 * (E_ref[iD, 3] + E_ref[iD, 4]))
    log(f"[setup] E_D (Wannier, 90x90) = {E_D:.6f} eV, gap {gap[iD]*1e3:.3e} meV; H(R): {Hwr.shape}, MP {MP}")
    return dict(cfg=cfg, dp=dp, paths=paths, Hwr=Hwr, Rw=Rw, nd=nd, k27=k27, MP=MP, U=U, Ud=Ud, E_D=E_D)


def cell_index(R_list, R):
    return tb.cell_index(R_list, R)


def load_cache_mwr():
    z = np.load(os.path.join(WORK, "cache", "Mwr_9x9.npz"))
    return {k: z[k] for k in z.files}


def load_cache_vloc():
    z = np.load(os.path.join(WORK, "cache", "Vloc_9x9.npz"))
    return {k: z[k] for k in z.files}


def g0_cached(name, Hwk, k_int, Rloc, eg, eta):
    """g0 sur une grille d'énergie, mis en cache dans cache/<name>.npz (exact, batché : local_green_batch)."""
    p = os.path.join(ensure("cache"), name + ".npz")
    if os.path.exists(p):
        z = np.load(p)
        if z["eg"].shape == eg.shape and np.allclose(z["eg"], eg) and float(z["eta"]) == eta and np.array_equal(z["Rloc"], Rloc):
            log(f"[g0] cache {name}: {z['g0'].shape}"); return z["g0"]
    t0 = time.time(); g0 = lt.local_green_batch(Hwk, k_int, Rloc, eg, eta)
    log(f"[g0] {name}: {g0.shape} computed in {time.time()-t0:.0f} s (nk_int {len(k_int)}, nE {len(eg)})")
    np.savez(p, g0=g0, eg=eg, eta=eta, Rloc=Rloc, nk_int=len(k_int)); return g0


def index_groups(Rloc, nw=NW):
    """Groupes d'indices L*nw+w : p_z lacune, p_z des trois voisins B, sp² de A (R=0), autres p_z, autres sp²."""
    i0 = cell_index(Rloc, (0, 0, 0)); nL = len(Rloc)
    pz_vac = [i0 * nw + WF_PZ_A]
    pz_nn = [cell_index(Rloc, R) * nw + WF_PZ_B for R in NN_CELLS]
    sp2_0 = [i0 * nw + w for w in WF_SIGMA]
    all_pi = set(pc.block_indices(nL, nw, WF_PI).tolist()); all_sig = set(pc.block_indices(nL, nw, WF_SIGMA).tolist())
    return dict(pz_vac=pz_vac, pz_nn=pz_nn, sp2_A0=sp2_0, pz_other=sorted(all_pi - set(pz_vac) - set(pz_nn)),
                sp2_other=sorted(all_sig - set(sp2_0)))


def fig_style():
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        plt.style.use(os.path.join(PROJ, "figures", "memoire.mplstyle"))
    except Exception as e:
        log(f"[fig] style non chargé : {e}")
    plt.rcParams["text.usetex"] = False
    sys.path.insert(0, os.path.join(PROJ, "scripts"))
    import _palette as pal
    return plt, pal


def savefig(fig, name):
    d = ensure("fig")
    fig.savefig(os.path.join(d, name + ".pdf")); fig.savefig(os.path.join(d, name + ".png"), dpi=200)
    log(f"[fig] {name}.pdf/.png")


# ----------------------------------------------------------------------------------------------- J1 : prep
def cmd_prep(a):
    S = setup(); out = {}; d = ensure("prep"); ensure("cache")
    Hwr, Rw, nd, k27, MP, U, Ud, dp = S["Hwr"], S["Rw"], S["nd"], S["k27"], S["MP"], S["U"], S["Ud"], S["dp"]

    # --- 1. M^L grossier de juin vs M^L recalculé (coarsecheck), et M^NL grossier avec la chaîne actuelle
    ML_june = np.load(os.path.join(RES, "M_L_9x9.npy"))                                    # sans sidecar : np.load direct
    ML_cc = matrix_io.load_M_checked(os.path.join(RES, "M_L_dense_9x9_coarsecheck.npy"), matrix_io.UNIT_CELL, matrix_io.HARTREE)
    M_ed = matrix_io.load_M_checked(os.path.join(RES, "M_ed_9x9.npy"), matrix_io.UNIT_CELL, matrix_io.HARTREE)
    nb, nk = ML_june.shape[:2]
    diag = lambda M: np.array([M[n, k, n, k] for n in range(nb) for k in range(nk)])
    ratio_L = float(np.abs(ML_june - ML_cc).max() / np.abs(ML_june).max())
    out["ML_check"] = dict(max_abs_diff_Ha=float(np.abs(ML_june - ML_cc).max()), max_abs_ML_june_Ha=float(np.abs(ML_june).max()),
                           ratio=ratio_L, rule="<= 2e-3 -> M_ed_9x9 (juin) pour D4 (a1)",
                           diag_mean_june_meV=float(diag(ML_june).real.mean() * HA2EV * 1e3),
                           diag_mean_cc_meV=float(diag(ML_cc).real.mean() * HA2EV * 1e3),
                           offdiag_ratio=float(np.abs((ML_june - ML_cc) - np.diag(np.diag((ML_june - ML_cc).reshape(nb * nk, nb * nk))).reshape(nb, nk, nb, nk)).max() / np.abs(ML_june).max()))
    log(f"[J1] M_L juin vs coarsecheck : ratio = {ratio_L:.3e} (règle <= 2e-3) ; diag mean {out['ML_check']['diag_mean_june_meV']:.2f} vs {out['ML_check']['diag_mean_cc_meV']:.2f} meV")
    from electron_defect_interaction.io.pseudo_io import read_upf
    from electron_defect_interaction.defects.non_local import compute_M_NL
    t0 = time.time(); M_NL_new = compute_M_NL(UC9, SC_P, SC_D, UPF, io=qe_io, pseudo_reader=read_upf)
    log(f"[J1] M^NL grossier recalculé (chaîne actuelle) : {M_NL_new.shape}, {time.time()-t0:.0f} s")
    matrix_io.save_M(os.path.join(d, "M_NL_coarse_new_9x9.npy"), M_NL_new, matrix_io.UNIT_CELL, part="M_NL_coarse_R4", note="compute_M_NL sur defect_unit_cell_9x9.save, 2026-09-25")
    M_NL_june = M_ed - ML_june
    out["MNL_check"] = dict(max_abs_diff_Ha=float(np.abs(M_NL_june - M_NL_new).max()), max_abs_MNL_new_Ha=float(np.abs(M_NL_new).max()),
                            ratio=float(np.abs(M_NL_june - M_NL_new).max() / np.abs(M_NL_new).max()))
    log(f"[J1] M^NL juin (M_ed - M_L) vs recalculé : ratio = {out['MNL_check']['ratio']:.3e}")
    if ratio_L <= 2e-3:
        out["M_coarse_choice"] = "M_ed_9x9.npy (juin) ; L = M_L_9x9.npy, NL = M_ed - M_L"
        np.save(os.path.join(d, "M_coarse_tot.npy"), M_ed); np.save(os.path.join(d, "M_coarse_L.npy"), ML_june); np.save(os.path.join(d, "M_coarse_NL.npy"), M_NL_june)
    else:
        out["M_coarse_choice"] = "RECALCULÉ : M_coarse = M_L_coarsecheck + M_NL_new (règle 2e-3 non satisfaite)"
        np.save(os.path.join(d, "M_coarse_tot.npy"), ML_cc + M_NL_new); np.save(os.path.join(d, "M_coarse_L.npy"), ML_cc); np.save(os.path.join(d, "M_coarse_NL.npy"), M_NL_new)
    log(f"[J1] choix M grossier : {out['M_coarse_choice']}")
    del ML_june, ML_cc, M_ed, M_NL_new, M_NL_june

    # --- 2. rotations V^dag M V -> M_W(R,R') (production), recentrage, linéarité
    files = dict(tot=dp["mfile"], L=dp["mfile"].replace("M_dense_", "M_L_dense_"), NL=dp["mfile"].replace("M_dense_", "M_NL_dense_"))
    Mwr = {}; Rn = None; Rd = None; argmax_cells = {}
    for tag in ("tot", "L", "NL"):
        t0 = time.time(); M = matrix_io.load_M_checked(files[tag], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)
        Mwk = Mbk_to_Mwk(M, U, Ud); del M
        Mw, R = Mwk_to_Mwr(Mwk, k27, MP); del Mwk
        w = np.array([np.linalg.norm(Mw[:, i, :, i]) for i in range(len(R))]); argmax_cells[tag] = R[int(np.argmax(w))].tolist()
        if tag == "tot":
            Rn, Rd = lt.recenter_mwr(Mw, R, MP); dist, wt = lt.mwr_locality(Mw, Rn)
            out["recenter"] = dict(R_d=Rd.tolist(), locality_first=[(float(x), float(y)) for x, y in zip(dist[:8], wt[:8])])
        Mwr[tag] = Mw; log(f"[J1] rotation {tag}: {time.time()-t0:.0f} s ; argmax on-site cell {argmax_cells[tag]}")
    out["argmax_onsite_cell"] = argmax_cells
    lin = float(np.linalg.norm(Mwr["tot"] - Mwr["L"] - Mwr["NL"]) / np.linalg.norm(Mwr["tot"]))
    out["linearity"] = dict(ratio=lin, expected="<= 1e-12"); log(f"[J1] linéarité ||tot-L-NL||/||tot|| = {lin:.3e}")
    np.savez(os.path.join(WORK, "cache", "Mwr_9x9.npz"), Mwr_tot=Mwr["tot"], Mwr_L=Mwr["L"], Mwr_NL=Mwr["NL"], R=Rn, R_d=Rd, MP=np.array(MP))
    log("[J1] cache/Mwr_9x9.npz écrit")

    # --- 3. V_loc (R_cut production) et blocs miroir
    rc = S["cfg"]["R_cut"]; Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]; nL = len(Rloc)
    V = {}; herm = {}
    for tag in ("tot", "L", "NL"):
        V[tag], herm[tag] = lt.extract_V_loc(Mwr[tag], Rn, Rloc)
    ip = pc.block_indices(nL, NW, WF_PI); isg = pc.block_indices(nL, NW, WF_SIGMA)
    ob = {tag: dict(max_offblock=float(max(np.abs(V[tag][np.ix_(isg, ip)]).max(), np.abs(V[tag][np.ix_(ip, isg)]).max())),
                    max_abs=float(np.abs(V[tag]).max()), herm_residual=float(herm[tag])) for tag in V}
    for tag in ob:
        ob[tag]["ratio"] = ob[tag]["max_offblock"] / ob[tag]["max_abs"]
    out["Vloc"] = dict(R_cut=rc, nL=nL, dim=nL * NW, offblock=ob, D3_approx_flag=bool(ob["tot"]["ratio"] > 1e-3))
    log(f"[J1] V_loc dim {nL*NW}; hors bloc sigma-pi (tot) = {ob['tot']['max_offblock']:.3e} eV / max {ob['tot']['max_abs']:.3f} eV (ratio {ob['tot']['ratio']:.2e}; seuil 1e-3)")
    np.savez(os.path.join(WORK, "cache", "Vloc_9x9.npz"), V_tot=V["tot"], V_L=V["L"], V_NL=V["NL"], Rloc=Rloc, idx_pi=ip, idx_sigma=isg)

    # --- 4. D2, base de Wannier
    i0 = cell_index(Rn, (0, 0, 0)); iNN = [cell_index(Rn, R) for R in NN_CELLS]
    iR0 = cell_index(Rw, (0, 0, 0))
    t_nn = {str(R): complex(Hwr[cell_index(Rw, R), WF_PZ_A, WF_PZ_B] / nd[cell_index(Rw, R)]) for R in NN_CELLS}
    d2w = {}
    for tag in ("tot", "L", "NL"):
        Mw = Mwr[tag]
        d2w[tag] = dict(pz_vac_onsite=complex(Mw[WF_PZ_A, i0, WF_PZ_A, i0]), pz_B_onsite_R0=complex(Mw[WF_PZ_B, i0, WF_PZ_B, i0]),
                        pz_vac_pz_nn={str(R): complex(Mw[WF_PZ_A, i0, WF_PZ_B, i]) for R, i in zip(NN_CELLS, iNN)},
                        pz_nn_diag={str(R): complex(Mw[WF_PZ_B, i, WF_PZ_B, i]) for R, i in zip(NN_CELLS, iNN)},
                        sp2_diag=[complex(Mw[w, i0, w, i0]) for w in WF_SIGMA],
                        sp2_offdiag_R0=complex(Mw[0, i0, 1, i0]), onsite_norm=float(np.linalg.norm(Mw[:, i0, :, i0])))
    A_uc, _ = qe_io.get_A_volume(dp["uc"]); A_cell = float(np.linalg.norm(np.cross(A_uc[:, 0], A_uc[:, 1])) * BOHR ** 2)
    out["D2_wannier"] = dict(elements_eV=d2w, t_nn_HR_eV=t_nn, HR_onsite_pz_eV=float(Hwr[iR0, WF_PZ_A, WF_PZ_A].real),
                             U_c_eV=5.94, kaasbjerg=dict(V0_eVA2=70.0, A_cell_A2=A_cell, per_cell_eV=70.0 / A_cell, per_atom_eV=70.0 / (A_cell / 2), A_cell_kaasbjerg=5.24))
    log(f"[J1] D2 Wannier : pz-pz lacune tot {d2w['tot']['pz_vac_onsite']:.4f} (L {d2w['L']['pz_vac_onsite']:.4f}, NL {d2w['NL']['pz_vac_onsite']:.4f}) eV ; t_nn {t_nn}")

    # --- 5. D2, base de Bloch à K et K' (27x27 dense), convention intensive eV Å²
    eps27 = qe_io.get_eigenvalues(dp["uc"]) * HA2EV                                # (20, 729)
    V27 = np.einsum("kbw,kwv->kbv", Ud, U); wpz = np.abs(V27[:, :, WF_PZ_A]) ** 2 + np.abs(V27[:, :, WF_PZ_B]) ** 2
    def kidx(kref):
        dd = np.linalg.norm(np.mod(k27 - kref + 0.5, 1.0) - 0.5, axis=1); i = int(np.argmin(dd)); return i, float(dd[i])
    iK, dK = kidx(K_RED); iKp, dKp = kidx(KP_RED)
    def pair(ik):
        top = np.argsort(-wpz[ik])[:2]; return sorted(top.tolist(), key=lambda n: eps27[n, ik]), wpz[ik][top].tolist()
    pK, wK = pair(iK); pKp, wKp = pair(iKp)
    A_sc = N_CELLS * A_cell
    d2b = dict(iK=iK, iKp=iKp, dist_K=dK, dist_Kp=dKp, pair_K=pK, pzweight_K=wK, pair_Kp=pKp, pzweight_Kp=wKp,
               eps_pair_K_eV=[float(eps27[n, iK]) for n in pK], eps_pair_Kp_eV=[float(eps27[n, iKp]) for n in pKp],
               A_cell_A2=A_cell, A_sc_A2=A_sc, blocks={})
    for tag in ("tot", "L", "NL"):
        Mm = np.load(files[tag], mmap_mode="r")
        colK = np.array(Mm[:, :, :, iK]) * HA2EV                                     # (nb, k', nb) M[m, k', n, K]
        bKK = colK[np.ix_(pK, [iK], pK)][:, 0, :]; bKpK = colK[np.ix_(pKp, [iKp], pK)][:, 0, :]
        d2b["blocks"][tag] = {}
        for lab, B in (("intra_KK", bKK), ("inter_KpK", bKpK)):
            d2b["blocks"][tag][lab] = dict(half_trace_eV=complex(0.5 * np.trace(B)), half_trace_eVA2=complex(0.5 * np.trace(B) * A_cell),
                                           frobenius_eV=float(np.linalg.norm(B)), frobenius_eVA2=float(np.linalg.norm(B) * A_cell),
                                           block_eV=B.tolist())
        del Mm, colK
    Rz = np.load(os.path.join(RES, "resonance_9x9.npz")); mlm = float(Rz["ML_diag_mean"])
    d2b["ML_diag_mean_production"] = dict(meV=mlm * 1e3, times_A_sc_eVA2=mlm * A_sc, source="results/M/resonance_9x9.npz ML_diag_mean (C14)")
    out["D2_bloch"] = d2b
    log(f"[J1] D2 Bloch : K idx {iK} paire {pK} (poids pz {np.round(wK,3)}), K' idx {iKp} paire {pKp}; A_cell {A_cell:.3f} Å², A_sc {A_sc:.1f} Å²")
    save_json(os.path.join(d, "prep_results.json"), out)


# ----------------------------------------------------------------------------------------------- J2 : D6
def single_site_report(V, supp, g0, eg, E_D, nontrivial_only=True, tag=""):
    """Pour un V à support restreint : |det| non normalisé (slogdet), valeurs propres non triviales, racines de Re lambda."""
    I = np.eye(V.shape[0]); A = I[None] - V[None] @ g0
    sign, logabs = np.linalg.slogdet(A); det = np.exp(logabs)
    lam = pc.nontrivial_eigenvalues(V, g0, supp)                                     # (nE, |S|)
    x = eg - E_D
    j = int(np.argmin(det))
    ilam = np.abs(lam).argmin(1); lam_min = lam[np.arange(len(eg)), ilam]
    # racines : changements de signe de Re de la valeur propre de plus petit module (suivie en identité par tri du module)
    roots = pc.sign_changes(lam_min.real, x)[1]
    mins = pc.local_minima(np.abs(lam_min), x)
    rep = dict(argmin_det=float(x[j]), det_min=float(det[j]), lam_at_argmin_det=[complex(v) for v in lam[j]],
               roots_Re_lam_min=[float(r) for r in roots], nroots=int(len(roots)),
               root_near_ED=(float(roots[np.argmin(np.abs(roots))]) if len(roots) else None),
               minima_abs_lam=[(float(x[i]), float(np.abs(lam_min[i]))) for i in mins],
               global_min_abs_lam=(float(x[int(np.argmin(np.abs(lam_min)))]), float(np.abs(lam_min).min())),
               det_at_m2p5=float(np.interp(-2.5, x, det)), lam_min_at_ED=complex(np.interp(0.0, x, lam_min.real) + 1j * np.interp(0.0, x, lam_min.imag)))
    return rep, det, lam


def cmd_d6(a):
    S = setup(); out = dict(D61={}, D62={}); d = ensure("d6"); cfg = S["cfg"]
    eta, npe = cfg["eta_eV"], cfg["ne_per_eta"]
    C = load_cache_vloc(); Rloc = C["Rloc"]; nL = len(Rloc); ip = C["idx_pi"]
    curves = {}
    # ---------------- D6.1 : modèle synthétique
    t = 2.7; Hs, Rs, nds = tb.graphene_pz_tb(t)
    chk = tb.check_pz_model(Hs, Rs, nds, t, eta=eta, nk=300)
    out["D61"]["model_check"] = chk
    log(f"[D6.1] modèle : écart bandes {chk['path_max_dev']:.2e} eV, E(K) = {chk['E_K_pi']:.2e}/{chk['E_K_pistar']:.2e}, E(M) = {chk['E_M_pi']:.4f}/{chk['E_M_pistar']:.4f}, E_D grille {chk['E_D']:.2e}, g0 batch {chk['g0_batch']:.6f} vs direct {chk['g0_direct']:.6f}")
    if chk["path_max_dev"] > 1e-8 or abs(chk["E_K_pi"]) > 1e-8 or abs(chk["g0_batch"] - chk["g0_direct"]) > 1e-10:
        raise RuntimeError("D6.1 : le modèle synthétique n'a pas passé sa vérification (erreur d'exécution, à corriger)")
    E_D0 = 0.0; eg0 = egrid(E_D0, -3.0, 1.0, eta, npe)
    i0 = cell_index(Rloc, (0, 0, 0)); iv = i0 * NW + WF_PZ_A
    fail = []
    for nk in (300, 1200):
        k_int = lt.mp_grid(nk, nk, 1); Hwk, _, _ = Hwr_to_Hwk(Hs, Rs, k_int, ndegen=nds)
        g0 = g0_cached(f"g0_toy_nk{nk}", Hwk, k_int, Rloc, eg0, eta)
        res = dict(g0_pz_onsite_at_m1=complex(g0[int(np.argmin(np.abs(eg0 + 1.0))), iv, iv]),
                   Re_g0_max=(float(eg0[int(np.argmax(g0[:, iv, iv].real))]), float(g0[:, iv, iv].real.max())), U={})
        for Uv in U_LIST:
            V, supp = tb.onsite_vloc(Rloc, NW, WF_PZ_A, Uv)
            rep, det, lam = single_site_report(V, supp, g0, eg0, E_D0); res["U"][str(Uv)] = rep
            curves[f"toy_nk{nk}_U{Uv}"] = np.stack([eg0, det, lam[:, 0].real, lam[:, 0].imag])
            exp = EXPECTED_1200.get(Uv); rep["expected_1200"] = exp
            if nk == 1200 and exp is not None:
                rep["dev_root"] = (None if exp[0] is None or rep["root_near_ED"] is None else float(rep["root_near_ED"] - exp[0]))
                rep["dev_argmin_det"] = float(rep["argmin_det"] - exp[1]); rep["dev_det_min"] = float(rep["det_min"] - exp[2])
                if Uv in FAIL_U and (rep["root_near_ED"] is None or abs(rep["dev_root"]) > 0.02):
                    fail.append((Uv, rep["root_near_ED"], exp[0]))
            log(f"[D6.1] nk {nk} U {Uv}: racines Re lam {np.round(rep['roots_Re_lam_min'],3)}, argmin|det| {rep['argmin_det']:+.3f} (|det| {rep['det_min']:.3f}), min|lam| {rep['global_min_abs_lam']}")
        V, supp = tb.removed_site_vloc(Hs, Rs, nds, Rloc, NW, WF_PZ_A, 1e3)
        rep, det, lam = single_site_report(V, supp, g0, eg0, E_D0); rep["support"] = [int(s) for s in supp]
        rep["min_abs_each_lambda_over_energies"] = [float(np.abs(lam[:, i]).min()) for i in range(lam.shape[1])]
        rep["min_over_energies_of_min_abs_lam"] = float(np.abs(lam).min(1).min())
        res["removed_site"] = rep
        curves[f"toy_nk{nk}_removed"] = np.concatenate([eg0[None], det[None], lam.real.T, lam.imag.T])
        log(f"[D6.1] nk {nk} site retiré : argmin|det| {rep['argmin_det']:+.4f} (|det| {rep['det_min']:.3f}), min|lam| {rep['min_over_energies_of_min_abs_lam']:.3f}")
        out["D61"][f"nk{nk}"] = res
    out["D61"]["fail"] = fail; out["D61"]["verdict"] = "FAIL" if fail else "PASS"
    log(f"[D6.1] verdict {out['D61']['verdict']} {fail}")
    if fail:
        save_json(os.path.join(d, "d6_results.json"), out); np.savez(os.path.join(d, "d6_curves.npz"), **curves)
        open(os.path.join(d, "D6_FAIL"), "w").write(str(fail)); sys.exit(1)

    # ---------------- D6.2 : vrai H(R), U sur la p_z du site lacunaire ; bloc pi de M_loc complet
    Hwr, Rw, nd, E_D = S["Hwr"], S["Rw"], S["nd"], S["E_D"]
    for nk in (300, 1200):
        k_int = lt.mp_grid(nk, nk, 1); Hwk, _, _ = Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=nd)
        egF = egrid(E_D, -3.0, 3.0, eta, npe) if nk == 300 else egrid(E_D, -3.0, 1.0, eta, npe)
        g0F = g0_cached(f"g0_real_nk{nk}", Hwk, k_int, Rloc, egF, eta)
        win = egF <= E_D + 1.0 + 1e-9; eg = egF[win]; g0 = g0F[win]
        j1 = int(np.argmin(np.abs(eg - (E_D - 1.0))))
        res = dict(g0_pz_onsite_at_ED_m1=complex(g0[j1, iv, iv]), toy_g0_pz_onsite_at_m1=out["D61"][f"nk{nk}"]["g0_pz_onsite_at_m1"],
                   Re_g0_max=(float(eg[int(np.argmax(g0[:, iv, iv].real))] - E_D), float(g0[:, iv, iv].real.max())), U={})
        for Uv in U_LIST:
            V, supp = tb.onsite_vloc(Rloc, NW, WF_PZ_A, Uv)
            rep, det, lam = single_site_report(V, supp, g0, eg, E_D); res["U"][str(Uv)] = rep
            curves[f"real_nk{nk}_U{Uv}"] = np.stack([eg - E_D, det, lam[:, 0].real, lam[:, 0].imag])
            log(f"[D6.2] nk {nk} U {Uv}: racines {np.round(rep['roots_Re_lam_min'],3)}, argmin|det| {rep['argmin_det']:+.3f} (|det| {rep['det_min']:.3f}), min|lam| {rep['global_min_abs_lam']}")
        # bloc pi du M_loc de production
        r = pc.det_eig_criterion(C["V_tot"], g0, idx=ip, vectors=True); x = eg - E_D
        ov = pc.branch_overlaps(r["vec_min"]); ii, xc = pc.sign_changes(r["lam_min"].real, x)
        jd = int(np.argmin(r["logdet_rel"])); jl = int(np.argmin(r["minlam"]))
        res["Mloc_pi_block"] = dict(dim=int(len(ip)), argmin_det=float(x[jd]), det_rel_min=float(np.exp(r["logdet_rel"][jd])), det_abs_min=float(np.exp(r["logabs"][jd])),
                                    minima_det=[(float(x[i]), float(np.exp(r["logdet_rel"][i]))) for i in pc.local_minima(r["logdet_rel"], x)],
                                    global_min_abs_lam=(float(x[jl]), float(r["minlam"][jl]), complex(r["lam_min"][jl])),
                                    minima_abs_lam=[(float(x[i]), float(r["minlam"][i]), complex(r["lam_min"][i])) for i in pc.local_minima(r["minlam"], x)],
                                    roots_Re_lam_min=[(float(e), float(ov[i]), bool(ov[i] < 0.9)) for i, e in zip(ii, xc)])
        curves[f"real_nk{nk}_Mloc_pi"] = np.stack([x, np.exp(r["logdet_rel"]), r["minlam"], r["lam_min"].real, r["lam_min"].imag])
        log(f"[D6.2] nk {nk} bloc pi M_loc : min|det|/max {res['Mloc_pi_block']['det_rel_min']:.3e} à {res['Mloc_pi_block']['argmin_det']:+.3f}, min|lam| {res['Mloc_pi_block']['global_min_abs_lam']}")
        out["D62"][f"nk{nk}"] = res
    save_json(os.path.join(d, "d6_results.json"), out); np.savez(os.path.join(d, "d6_curves.npz"), **curves)
    # figures
    plt, pal = fig_style()
    for kind, ED in (("toy", 0.0), ("real", E_D)):
        fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.4))
        for Uv, c in zip(U_LIST, pal.CYCLE * 2):
            z = curves[f"{kind}_nk1200_U{Uv}"]
            ax[0].plot(z[0], z[1], color=c, lw=1.0, label=f"U = {Uv:g}"); ax[1].plot(z[0], z[2], color=c, lw=1.0)
        ax[0].set_yscale("log"); ax[0].set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax[0].set_ylabel(r"$|\det[1 - V g^{(0)}]|$"); ax[0].set_title("(a) 1200², " + ("modèle" if kind == "toy" else "H(R) de production"))
        ax[1].axhline(0, color=pal.MUTED, lw=0.6); ax[1].set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax[1].set_ylabel(r"Re $\lambda$ (valeur propre non triviale)"); ax[1].set_ylim(-3, 3); ax[1].set_title("(b)")
        ax[0].legend(fontsize=6, ncol=2); fig.tight_layout(); savefig(fig, f"d6_{kind}_1200")
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    for nk, ls in ((300, "--"), (1200, "-")):
        for Uv, c in zip([6.617, 12.0, 27.0], [pal.NAVY, pal.ORANGE, pal.GREEN]):
            z = curves[f"toy_nk{nk}_U{Uv}"]; ax.plot(z[0], z[2], color=c, ls=ls, lw=1.0, label=f"U = {Uv:g}, {nk}²")
    ax.axhline(0, color=pal.MUTED, lw=0.6); ax.set_ylim(-2, 2); ax.set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax.set_ylabel(r"Re $\lambda = 1 - U g_0$"); ax.legend(fontsize=6, ncol=2)
    fig.tight_layout(); savefig(fig, "d6_toy_grid")


# ----------------------------------------------------------------------------------------------- J3 : D1
GEOMS = [dict(tag="P", name="9x9 parfaite", save=SC_P, spin=None), dict(tag="D", name="9x9 lacune non relaxée", save=SC_D, spin=None),
         dict(tag="N1", name="R1 nspin1", save=R1_N1, spin=None), dict(tag="N2u", name="R1 nspin2 ↑", save=R1_N2, spin="up"),
         dict(tag="N2d", name="R1 nspin2 ↓", save=R1_N2, spin="dw")]


def load_pot_eV(path):
    V, ng = qe_io.get_pot(path, subtract_mean=False, to_hartree=True)
    return V.transpose(2, 1, 0) * HA2EV


def second_neighbours(x_d, s_vac, A_A):
    """Indices (0-based) des atomes à 2.2-2.7 Å de la lacune (6 seconds voisins, sous-réseau A) et des 3 premiers (< 1.8 Å)."""
    dd = al.min_image_dist(x_d, s_vac, A_A)
    return np.where(dd < 1.8)[0].tolist(), np.where((dd > 2.2) & (dd < 2.7))[0].tolist(), dd


def cmd_d1(a):
    d = ensure("d1"); out = {}; arrays = {}
    A_b, _ = qe_io.get_A_volume(SC_P); A_A = A_b * BOHR
    x_p = qe_io.get_x_red(SC_P); x_d = qe_io.get_x_red(SC_D); ng = qe_io.get_ngfft(SC_D)
    s_vac, i_vac_p, _ = al.vacancy_site(x_p, x_d, A_A)
    nn1, nn2, dd = second_neighbours(x_d, s_vac, A_A)
    out["geometry"] = dict(s_vac=s_vac.tolist(), missing_atom_p_index0=i_vac_p, ngfft=list(ng), A_A=A_A.tolist(),
                           first_neighbours_index1=[i + 1 for i in nn1], second_neighbours_index1=[i + 1 for i in nn2],
                           d_first=[float(dd[i]) for i in nn1], d_second=[float(dd[i]) for i in nn2])
    log(f"[D1] s_vac {np.round(s_vac,5)}, premiers voisins (1-based) {[i+1 for i in nn1]}, seconds {[i+1 for i in nn2]}")
    # E_D de la super-cellule parfaite : quadruplet dégénéré à Gamma
    eP, efP, _ = qg.get_eigenvalues_spin(SC_P); eP = eP[0] * HA2EV; efP *= HA2EV
    q = np.argsort(np.abs(eP - efP))[:4]; E_D = float(eP[q].mean()); spread = float(eP[q].max() - eP[q].min())
    Sw = setup(); E_D_W = Sw["E_D"]
    out["E_D_SC"] = dict(E_D=E_D, quadruplet_bands_1based=sorted((q + 1).tolist()), quadruplet=sorted(eP[q].tolist()), spread_eV=spread, E_F_P=efP,
                         E_D_minus_EF=E_D - efP, E_D_wannier=E_D_W, E_D_SC_minus_wannier=E_D - E_D_W)
    log(f"[D1] E_D(SC) = {E_D:.5f} eV (quadruplet {np.round(eP[q],5)}, étalement {spread*1e3:.2f} meV), E_F(P) {efP:.5f}, E_D Wannier {E_D_W:.5f} (écart {E_D-E_D_W:+.4f})")
    # comparaison parfaite 9x9 vs maille unitaire repliée (16 bandes x 81 k)
    eU = np.sort(qe_io.get_eigenvalues(UC9).reshape(-1)) * HA2EV
    nlo = min(len(eU), len(eP)); out["P_vs_uc_folded"] = dict(n=int(nlo), max_abs_dev_lowest_eV=float(np.abs(np.sort(eP)[:nlo] - eU[:nlo]).max()),
                                                            max_abs_dev_window_eV=float(np.abs(np.sort(eP)[:nlo] - eU[:nlo])[(np.sort(eP)[:nlo] - E_D > -3) & (np.sort(eP)[:nlo] - E_D < 1)].max()))
    # potentiels et alignement (Lu 2019)
    V_P = load_pot_eV(VKS["p"]); pots = dict(D=load_pot_eV(VKS["d"]))
    have_r1 = all(os.path.exists(VKS[k]) for k in ("n1", "n2u", "n2d"))
    if have_r1:
        pots["N1"] = load_pot_eV(VKS["n1"]); Vu = load_pot_eV(VKS["n2u"]); Vd = load_pot_eV(VKS["n2d"]); pots["N2"] = 0.5 * (Vu + Vd); pots["N2u"] = Vu; pots["N2d"] = Vd
    else:
        log("[D1] potentiels pp.x de R1 absents : alignement des géométries relaxées non calculé (J0 non exécuté)")
    xs = dict(D=x_d, N1=qe_io.get_x_red(R1_N1), N2=qe_io.get_x_red(R1_N2))
    align = {}
    for key in ("D", "N1", "N2", "N2u", "N2d"):
        if key not in pots: continue
        xg = xs["N2"] if key.startswith("N2") else xs[key]
        r = al.far_atom_alignment(pots[key], V_P, xg, x_p, A_A, RADII_A)
        align[key] = dict(i_far_1based=r["i_far_d"] + 1, dist_far_A=r["dist_far"], far_displacement_A=r["far_displacement"],
                          shifts_eV={str(k): float(v) for k, v in r["shifts"].items()}, means={str(k): [float(v[0]), float(v[1]), int(v[2]), int(v[3])] for k, v in r["means"].items()},
                          mean3d_diff_eV=float(pots[key].mean() - V_P.mean()))
        log(f"[D1] alignement {key}: atome loin {r['i_far_d']+1} à {r['dist_far']:.3f} Å ; décalage {', '.join(f'{k} Å: {v*1e3:+.2f} meV' for k, v in r['shifts'].items())} ; <V_d>-<V_p> 3D {align[key]['mean3d_diff_eV']*1e3:+.2f} meV")
    out["alignment"] = align; out["alignment_radius_used_A"] = 1.0
    del pots, V_P
    # états par géométrie
    mask2 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W2); mask1 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W1)
    out["disc_area_fraction"] = dict(r2=float(mask2.mean()), r1=float(mask1.mean()))
    states = {}
    for g in GEOMS:
        tag = g["tag"]; save = g["save"]
        eigs, ef, lsda = qg.get_eigenvalues_spin(save); ef *= HA2EV
        e = eigs[{None: 0, "up": 0, "dw": 1}[g["spin"]]] * HA2EV
        akey = {"P": None, "D": "D", "N1": "N1", "N2u": "N2", "N2d": "N2"}[tag]
        shift = 0.0 if akey is None or akey not in align else align[akey]["shifts_eV"]["1.0"]
        e_al = e - shift; x = e_al - E_D
        sel = np.where((x >= -3.0) & (x <= 1.0))[0]
        if len(sel) == 0:
            log(f"[D1] {tag}: aucun état dans la fenêtre"); continue
        b0, b1 = int(sel[0]), int(sel[-1]) + 1
        xr = qe_io.get_x_red(save); z0 = float(np.mean(xr[:, 2])); zdev = float(np.abs(xr[:, 2] - z0).max() * A_A[2, 2])
        t0 = time.time(); Cw, mill, go, at = qg.read_wfc_gamma(save, spin=g["spin"], bands=(b0, b1))
        par = qg.mirror_parity_z(Cw, mill, z0_red=z0, gamma_only=go)
        rho = qg.density_2d(Cw, mill, ng, gamma_only=go, workers=int(os.environ.get("OMP_NUM_THREADS", "8")))
        w2 = np.array([r[mask2].sum() for r in rho]); w1 = np.array([r[mask1].sum() for r in rho]); del Cw, rho
        rec = dict(name=g["name"], nbnd=int(len(e)), E_F=float(ef), shift_eV=float(shift), bands_1based=[b0 + 1, b1], n_window=int(b1 - b0),
                   z0_red=z0, max_z_dev_A=zdev, gamma_only=bool(go), max_dev_parity=float(np.abs(1 - np.abs(par)).max()),
                   n_even=int((par.real > 0).sum()), n_odd=int((par.real < 0).sum()), read_s=time.time() - t0,
                   nelec=None)
        arrays[f"{tag}_e"] = e[b0:b1]; arrays[f"{tag}_e_al"] = e_al[b0:b1]; arrays[f"{tag}_x"] = x[b0:b1]; arrays[f"{tag}_relEF"] = e[b0:b1] - ef
        arrays[f"{tag}_parity"] = par.real; arrays[f"{tag}_w2"] = w2; arrays[f"{tag}_w1"] = w1; arrays[f"{tag}_band1"] = np.arange(b0 + 1, b1 + 1)
        states[tag] = rec
        log(f"[D1] {tag}: bandes {b0+1}-{b1} ({b1-b0} états), parité max|1-|<s>|| = {rec['max_dev_parity']:.2e}, pairs {rec['n_even']} / impairs {rec['n_odd']}, w2 max {w2.max():.3f}, lu en {rec['read_s']:.0f} s")
    out["states"] = states
    # seuil de localisation : 3 x moyenne de w2 dans la parfaite
    thr = 3.0 * float(arrays["P_w2"].mean()); out["localization_threshold"] = dict(thr_w2=thr, mean_w2_P=float(arrays["P_w2"].mean()), max_w2_P=float(arrays["P_w2"].max()), rule="w2 > 3 x moyenne(parfaite)")
    # projwfc (R1 nspin2 à Gamma)
    pw = os.path.join(R1, "projwfc", "gamma_nspin2", "projwfc.out"); proj = {}
    if os.path.exists(pw):
        st, bands = pj.read_projwfc_states(pw)
        nn1_1 = [i + 1 for i in nn1]; nn2_1 = [i + 1 for i in nn2]
        iso = [81]; pairat = [i for i in nn1_1 if i != 81]
        groups = dict(nn_pz=(set(nn1_1), {"pz"}), iso_s_pxy=(set(iso), {"s", "px", "py"}), pair_s_pxy=(set(pairat), {"s", "px", "py"}),
                      nnn_pz=(set(nn2_1), {"pz"}), nn_all=(set(nn1_1), {"s", "px", "py", "pz"}))
        for b in bands:
            key = f"{'N2u' if b['kblock'] == 0 else 'N2d'}"
            proj.setdefault(key, {})[b["band"]] = dict(e_eV=b["e_eV"], psi2=b["psi2"], **pj.weights_by_group(b, st, groups))
        out["projwfc"] = dict(file=pw, nstates=len(st), nkblocks=max(b["kblock"] for b in bands) + 1, groups={k: (sorted(v[0]), sorted(v[1])) for k, v in groups.items()})
        log(f"[D1] projwfc.out lu : {len(st)} états atomiques, {len(bands)} bandes, {out['projwfc']['nkblocks']} blocs k")
    # tableau des états localisés
    loc = {}
    for tag in states:
        rows = []
        for j in range(len(arrays[f"{tag}_x"])):
            if arrays[f"{tag}_w2"][j] > thr:
                b = int(arrays[f"{tag}_band1"][j])
                row = dict(band=b, x=float(arrays[f"{tag}_x"][j]), relEF=float(arrays[f"{tag}_relEF"][j]), parity=float(arrays[f"{tag}_parity"][j]),
                           w2=float(arrays[f"{tag}_w2"][j]), w1=float(arrays[f"{tag}_w1"][j]))
                if tag in proj and b in proj[tag]:
                    row["projwfc"] = proj[tag][b]
                rows.append(row)
        loc[tag] = rows
    out["localized"] = loc
    for tag in loc:
        log(f"[D1] {tag}: {len(loc[tag])} états avec w2 > {thr:.3f} : " + "; ".join(f"b{r['band']} x={r['x']:+.3f} p={r['parity']:+.2f} w2={r['w2']:.3f}" for r in loc[tag]))
    # entrées projwfc.x préparées (non lancées)
    for tag, save, prefix, outdir in (("defect_9x9_d", SC_D, "defect_9x9_d", f"{SCRATCH}/defect_9x9_d"), ("defect_9x9_p", SC_P, "defect_9x9_p", f"{SCRATCH}/defect_9x9_p"),
                                      ("vac_9x9_relax_nspin1", R1_N1, "vac_9x9_relax_nspin1", f"{SCRATCH}/vacancy_relaxed/nspin1")):
        pd = os.path.join(d, "projwfc_inputs", tag); os.makedirs(pd, exist_ok=True)
        with open(os.path.join(pd, "projwfc.in"), "w") as f:
            f.write(f"&PROJWFC\n  prefix   = '{prefix}'\n  outdir   = '{outdir}'\n  filproj  = 'proj_{tag}'\n  filpdos  = 'pdos_{tag}'\n"
                    "  lwrite_overlaps = .false.\n  Emin = -30.0, Emax = 10.0, DeltaE = 1.0\n  degauss = 0.01, ngauss = -1\n/\n")
    np.savez(os.path.join(d, "d1_states.npz"), E_D=E_D, thr=thr, s_vac=s_vac, **arrays)
    save_json(os.path.join(d, "d1_results.json"), out)
    # figure
    plt, pal = fig_style()
    tags = [t for t in ("P", "D", "N1", "N2u", "N2d") if t in states]
    fig, axs = plt.subplots(1, len(tags), figsize=(6.5, 3.6), sharey=True)
    for ax, tag in zip(np.atleast_1d(axs), tags):
        x = arrays[f"{tag}_x"]; w = arrays[f"{tag}_w2"]; p = arrays[f"{tag}_parity"]
        ax.scatter(w[p > 0], x[p > 0], s=9, color=pal.NAVY, label="paire (σ)"); ax.scatter(w[p < 0], x[p < 0], s=9, color=pal.ORANGE, label="impaire (π)")
        ax.axvline(thr, color=pal.MUTED, lw=0.6, ls="--"); ax.set_title(states[tag]["name"], fontsize=8)
        ax.set_xlim(0, max(0.05, float(w.max()) * 1.1))
    np.atleast_1d(axs)[0].set_ylabel(r"$\varepsilon - E_D$ (eV)"); np.atleast_1d(axs)[0].legend(fontsize=6)
    fig.supxlabel(r"poids $w_2$ à moins de 2 Å de la lacune (seuil en tirets)", fontsize=9); fig.tight_layout(); savefig(fig, "d1_states")


# ----------------------------------------------------------------------------------------------- J3 : D5 et R
def ved_metrics(dV, A_A, s_vac, shift=0.0):
    """Moyenne 3D, valeur au site, ligne le long de a1 (lacune -> demi-boîte), max |dV - c| sur les plans frontière."""
    from scipy.ndimage import map_coordinates
    nr = np.array(dV.shape); mean3d = float(dV.mean())
    t = np.linspace(0, 0.5, 401); s = s_vac[None, :] + t[:, None] * np.array([1.0, 0, 0])[None, :]
    line = map_coordinates(dV, (s * nr[None, :]).T, order=1, mode="wrap")
    i1 = int(np.round((s_vac[0] + 0.5) * nr[0])) % nr[0]; i2 = int(np.round((s_vac[1] + 0.5) * nr[1])) % nr[1]
    iz = int(np.round(s_vac[2] * nr[2])) % nr[2]
    bnd = lambda c: float(max(np.abs(dV[i1] - c).max(), np.abs(dV[:, i2] - c).max()))
    bnd_pl = lambda c: float(max(np.abs(dV[i1, :, iz] - c).max(), np.abs(dV[:, i2, iz] - c).max()))
    site = float(map_coordinates(dV, (s_vac * nr)[:, None], order=1, mode="wrap")[0])
    return dict(mean3d=mean3d, site=site, end_raw=float(line[-1]), end_sub=float(line[-1] - mean3d), end_aligned=float(line[-1] - shift),
                boundary_raw=bnd(0.0), boundary_sub=bnd(mean3d), boundary_aligned=bnd(shift),
                boundary_plane_raw=bnd_pl(0.0), boundary_plane_sub=bnd_pl(mean3d), boundary_plane_aligned=bnd_pl(shift), line=line, t=t)


def cmd_d5(a):
    d = ensure("d5"); out = {}; lines = {}
    sizes = ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]
    d1 = json.load(open(os.path.join(WORK, "d1", "d1_results.json"))) if os.path.exists(os.path.join(WORK, "d1", "d1_results.json")) else None
    for Sz in sizes:
        scp = f"{SCRATCH}/defect_{Sz}_p/defect_{Sz}_p.save"; scd = f"{SCRATCH}/defect_{Sz}_d/defect_{Sz}_d.save"
        vp = f"{GQ}/defects/super_cell/{Sz}/pristine/Vks_{Sz}_p"; vd = f"{GQ}/defects/super_cell/{Sz}/defective/Vks_{Sz}_d"
        if not (os.path.exists(vp) and os.path.exists(vd) and os.path.isdir(scp) and os.path.isdir(scd)):
            log(f"[D5] {Sz}: fichiers absents, sauté"); continue
        t0 = time.time(); A_b, _ = qe_io.get_A_volume(scd); A_A = A_b * BOHR; x_p = qe_io.get_x_red(scp); x_d = qe_io.get_x_red(scd)
        Vp = load_pot_eV(vp); Vd = load_pot_eV(vd)
        r = al.far_atom_alignment(Vd, Vp, x_d, x_p, A_A, RADII_A); s_vac = r["s_vac"]
        m = ved_metrics(Vd - Vp, A_A, s_vac, shift=float(r["shifts"][1.0])); lines[Sz] = np.stack([m.pop("t"), m.pop("line")])
        N = int(Sz.split("x")[0]); A_cell = float(np.linalg.norm(np.cross(A_A[:, 0], A_A[:, 1])) / N ** 2)
        out[Sz] = dict(N=N, A_sc_A2=A_cell * N * N, s_vac=s_vac.tolist(), far_atom_1based=r["i_far_d"] + 1, dist_far_A=r["dist_far"],
                       shift_Lu_eV={str(k): float(v) for k, v in r["shifts"].items()}, G0_component_eV=m["mean3d"], **m)
        log(f"[D5] {Sz}: <dV>_3D {m['mean3d']*1e3:+.2f} meV, site {m['site']:+.3f} eV, bord (ligne a1) brut {m['end_raw']*1e3:+.1f} / -moyenne {m['end_sub']*1e3:+.1f} / aligné {m['end_aligned']*1e3:+.1f} meV ; "
            f"plans frontière max|dV| brut {m['boundary_raw']*1e3:.1f} / -moy {m['boundary_sub']*1e3:.1f} / aligné {m['boundary_aligned']*1e3:.1f} meV ; Lu {r['shifts'][0.5]*1e3:+.2f} / {r['shifts'][1.0]*1e3:+.2f} meV ({time.time()-t0:.0f} s)")
        del Vp, Vd
    Rz = np.load(os.path.join(RES, "resonance_9x9.npz"))
    out["production_subtraction"] = dict(M_dense_chain="subtract_mean=False (compute_M_dense_stages.py, stage ml)", C14_shift_meV=float(Rz["ML_diag_mean"]) * 1e3,
                                         note="aucune moyenne soustraite dans la construction de M dense ; 67.0 meV soustraits seulement dans le contrôle C14")
    np.savez(os.path.join(d, "d5_lines.npz"), **{k: v for k, v in lines.items()})
    save_json(os.path.join(d, "d5_results.json"), out)
    plt, pal = fig_style()
    fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.2))
    Ns = [out[s]["N"] for s in sizes if s in out]
    for key, lab, c in (("end_raw", "brut", pal.REF), ("end_sub", "moyenne soustraite", pal.SKY), ("end_aligned", "aligné (Lu, 1 Å)", pal.NAVY)):
        ax[0].plot(Ns, [out[s][key] * 1e3 for s in sizes if s in out], "o-", color=c, label=lab, ms=3)
    ax[0].set_xlabel("N"); ax[0].set_ylabel(r"$\Delta V^L$ au bord, ligne $a_1$ (meV)"); ax[0].legend(fontsize=6)
    for Sz, c in zip(sizes, pal.CYCLE + [pal.REF, pal.MUTED]):
        if Sz in lines: ax[1].plot(lines[Sz][0] * out[Sz]["N"] * 2.4659, lines[Sz][1], color=c, lw=0.9, label=Sz)
    ax[1].set_xlabel(r"distance à la lacune le long de $a_1$ (Å)"); ax[1].set_ylabel(r"$\Delta V^L$ (eV)"); ax[1].set_ylim(-3, 3); ax[1].legend(fontsize=6, ncol=2)
    fig.tight_layout(); savefig(fig, "d5_boundary")


def cmd_r(a):
    d = ensure("r"); out = {}
    out["code_reading"] = dict(
        M_NL="defects/non_local.py compute_M_NL (l. 109-201) : positions de TOUS les atomes lues séparément dans sc_p (tau_s_p) et sc_d (tau_s_d), phases par atome, M_NL = M_d - M_p -> une géométrie relaxée est traitée par construction ; aucun indice de spin ; les projecteurs KB ne dépendent pas du spin -> ΔV^NL identique pour ↑ et ↓ à géométrie donnée",
        M_L="defects/local_R.py compute_ML_R (l. 15-97) et local_G.py : ΔV = V_d - V_p lus par qe_io.get_pot (pp.x plot_num 1) ; positions non utilisées (le potentiel porte la géométrie) ; aucun indice de spin : nspin = 2 exige deux potentiels pp.x (spin_component 1 et 2) et deux passages de la partie LOCALE seulement",
        pot="io/qe_io.py get_pot (l. 393-453) : un fichier filplot, aucune notion de spin",
        wfc="io/qe_io.py _read_all_wfc (l. 342-391) : npol = 1 exigé ; glob wfc*.hdf5 trié par l'attribut ik -> les fichiers wfcup1/wfcdw1 d'un run lsda (même ik) ne sont pas distingués ; gamma_only non géré (demi-sphère) -- ne concerne que les super-cellules, jamais lues par la chaîne M",
        missing="manquant pour R1 : (i) deux potentiels ΔV_↑, ΔV_↓ (pp.x de J0 les fournit) ; (ii) un pilote qui appelle compute_ML_R deux fois (up, down) et compute_M_NL une fois avec sc_d = .save relaxé, puis somme par spin ; (iii) le sidecar matrix_io doit porter le spin ; (iv) la géométrie relaxée entre dans M_NL par les positions du .save relaxé (déjà supporté)")
    if not all(os.path.exists(VKS[k]) for k in ("n1", "n2u", "n2d")):
        out["dV_R1"] = "potentiels pp.x de R1 absents (J0 non exécuté)"; save_json(os.path.join(d, "r_results.json"), out); return
    A_b, _ = qe_io.get_A_volume(SC_P); A_A = A_b * BOHR; x_p = qe_io.get_x_red(SC_P); x_d = qe_io.get_x_red(SC_D); x_n1 = qe_io.get_x_red(R1_N1); x_n2 = qe_io.get_x_red(R1_N2)
    d1 = json.load(open(os.path.join(WORK, "d1", "d1_results.json"))); al1 = d1["alignment"]
    Vp = load_pot_eV(VKS["p"]); s_vac, _, _ = al.vacancy_site(x_p, x_d, A_A)
    res = {}
    for key, path, sh in (("D_unrelaxed", VKS["d"], al1["D"]["shifts_eV"]["1.0"]), ("N1", VKS["n1"], al1["N1"]["shifts_eV"]["1.0"]),
                          ("N2u", VKS["n2u"], al1["N2u"]["shifts_eV"]["1.0"]), ("N2d", VKS["n2d"], al1["N2d"]["shifts_eV"]["1.0"])):
        V = load_pot_eV(path); m = ved_metrics(V - Vp, A_A, s_vac, shift=float(sh)); m.pop("t"); line = m.pop("line"); res[key] = dict(shift_used_eV=float(sh), **m); res[key]["line_end_points"] = line[[0, 100, 200, 300, 400]].tolist()
        log(f"[R] dV {key}: <dV> {m['mean3d']*1e3:+.2f} meV, site {m['site']:+.3f} eV, bord ligne brut {m['end_raw']*1e3:+.1f} / aligné {m['end_aligned']*1e3:+.1f} meV, plans frontière brut {m['boundary_raw']*1e3:.1f} / aligné {m['boundary_aligned']*1e3:.1f} meV")
        if key == "N2u": Vu = V
        elif key == "N2d":
            dif = Vu - V; res["N2_up_minus_dw"] = dict(max_abs_eV=float(np.abs(dif).max()), mean3d_eV=float(dif.mean()), boundary_max_abs_eV=ved_metrics(dif, A_A, s_vac)["boundary_raw"], site_eV=ved_metrics(dif, A_A, s_vac)["site"])
            log(f"[R] V_up - V_dw : max {res['N2_up_minus_dw']['max_abs_eV']:.3f} eV, moyenne {res['N2_up_minus_dw']['mean3d_eV']*1e3:+.2f} meV, frontière {res['N2_up_minus_dw']['boundary_max_abs_eV']*1e3:.2f} meV")
            del Vu
        del V
    out["dV_R1"] = res
    save_json(os.path.join(d, "r_results.json"), out)


# ----------------------------------------------------------------------------------------------- J4 : D4
def read_wfc_subset(save_dir, ks_index):
    """C_nk(G) des fichiers wfc<ik>.hdf5 (ik = index + 1, ordre du XML) pour un sous-ensemble d'indices k."""
    Cs = []; Gs = []; nG = []
    for j in ks_index:
        C, mill, go, at = qg.read_wfc_file(os.path.join(save_dir, f"wfc{j+1}.hdf5"))
        assert not go; Cs.append(C); Gs.append(mill); nG.append(C.shape[1])
    nGm = max(nG); nb = Cs[0].shape[0]
    C_nkg = np.zeros((nb, len(ks_index), nGm), complex); G_red = np.zeros((len(ks_index), nGm, 3), int)
    for i, (C, G) in enumerate(zip(Cs, Gs)):
        C_nkg[:, i, :C.shape[1]] = C; G_red[i, :C.shape[1]] = G
    return C_nkg, np.array(nG), G_red


def uc_parity(C_nkg, nG, G_red):
    p = np.zeros(C_nkg.shape[:2])
    for k in range(C_nkg.shape[1]):
        p[:, k] = qg.mirror_parity_z(C_nkg[:, k, :nG[k]], G_red[k, :nG[k]], 0.0, gamma_only=False).real
    return p


def spectrum_blocks(H, parity_flat, E_Dm, groups=None):
    """Diagonalisation par bloc de parité (sign de parity_flat) ; renvoie par bloc (e - E_Dm, vecteurs) et le couplage résiduel."""
    res = {}
    Hh = 0.5 * (H + H.conj().T); herm = float(np.abs(H - H.conj().T).max())
    ie = np.where(parity_flat > 0)[0]; io = np.where(parity_flat < 0)[0]
    coupling = float(np.abs(Hh[np.ix_(ie, io)]).max())
    for lab, idx in (("even", ie), ("odd", io)):
        e, v = np.linalg.eigh(Hh[np.ix_(idx, idx)])
        res[lab] = dict(x=e - E_Dm, v=v, idx=idx)
    return res, coupling, herm


def cmd_d4(a):
    S = setup(); d = ensure("d4"); out = dict(gates={}, variants={}); Hwr, Rw, nd, k27, U, Ud, dp = S["Hwr"], S["Rw"], S["nd"], S["k27"], S["U"], S["Ud"], S["dp"]
    z1 = np.load(os.path.join(WORK, "d1", "d1_states.npz")); d1 = json.load(open(os.path.join(WORK, "d1", "d1_results.json")))
    thr = float(z1["thr"]); s_vac = z1["s_vac"]; E_D_SC = float(z1["E_D"])
    qe = dict(x=z1["D_x"], parity=z1["D_parity"], w2=z1["D_w2"], e_abs=z1["D_e"])
    ng = qe_io.get_ngfft(SC_D); A_b, _ = qe_io.get_A_volume(SC_D); A_A = A_b * BOHR
    mask2 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W2)
    # bases
    k81 = qe_io.get_k_red(UC9); eps16 = qe_io.get_eigenvalues(UC9) * HA2EV
    C16, nG16 = qe_io.get_C_nk(UC9); G16 = qe_io.get_G_red(UC9); par16 = uc_parity(C16, nG16, G16)
    idx81 = np.array([int(np.where(np.all(np.abs(np.mod(k27 - k + 0.5, 1) - 0.5) < 1e-6, axis=1))[0][0]) for k in k81])
    eps20 = (qe_io.get_eigenvalues(dp["uc"]) * HA2EV)[:, idx81]
    C20, nG20, G20 = read_wfc_subset(dp["uc"], idx81); par20 = uc_parity(C20, nG20, G20)
    out["bases"] = dict(k81_subset_first=idx81[:5].tolist(), parity16_dev=float(np.abs(1 - np.abs(par16)).max()), parity20_dev=float(np.abs(1 - np.abs(par20)).max()),
                        n_even16=int((par16 > 0).sum()), n_odd16=int((par16 < 0).sum()), n_even20=int((par20 > 0).sum()), n_odd20=int((par20 < 0).sum()))
    iK = int(np.argmin(np.linalg.norm(np.mod(k81 - K_RED + 0.5, 1) - 0.5, axis=1)))
    E_D16 = float(0.5 * (eps16[3, iK] + eps16[4, iK])); E_D20 = float(0.5 * (eps20[3, iK] + eps20[4, iK]))
    _, EK, _ = Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); E_DW = float(0.5 * (EK[0, 3] + EK[0, 4]))
    out["E_D"] = dict(uc16_at_K=E_D16, uc20dense_at_K=E_D20, wannier_at_K=E_DW, SC_P=E_D_SC, uc16_bands34_K=[float(eps16[3, iK]), float(eps16[4, iK])], uc20_bands34_K=[float(eps20[3, iK]), float(eps20[4, iK])])
    log(f"[D4] E_D : uc16 {E_D16:.5f}, dense20 {E_D20:.5f}, Wannier {E_DW:.5f}, SC parfaite {E_D_SC:.5f} eV")
    Nd = (N_SC, N_SC, 1); flat16 = sf.sc_planewave_index(k81, G16, nG16, Nd, ng); flat20 = sf.sc_planewave_index(k81, G20, nG20, Nd, ng)
    V81 = np.einsum("kbw,kwv->kbv", Ud[idx81], U[idx81])                                      # (81, 20, 5)
    nw = NW; nk = len(k81)

    def weights_bloch(blocks, C, nG, flat, nb, lo=-3.0, hi=1.0):
        """poids w2 (disque 2 Å) des états de la fenêtre, par bloc ; d_nk reconstruit sur la base (n, k)."""
        W = {}
        for lab in blocks:
            x = blocks[lab]["x"]; sel = np.where((x >= lo) & (x <= hi))[0]; w = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                dfull = np.zeros(nb * nk, complex); dfull[blocks[lab]["idx"]] = blocks[lab]["v"][:, j]
                w[jj] = sf.folded_density_2d(dfull.reshape(nb, nk), C, nG, flat, ng, workers=int(os.environ.get("OMP_NUM_THREADS", "8")))[mask2].sum()
            W[lab] = dict(sel=sel, w2=w)
        return W

    def record(tag, blocks, W, E_Dm, extra=None):
        rec = dict(E_D=E_Dm, blocks={}, **(extra or {}))
        for lab in ("even", "odd"):
            x = blocks[lab]["x"]; sel = W[lab]["sel"]; w = W[lab]["w2"]
            locs = [(float(x[j]), float(w[i])) for i, j in enumerate(sel) if w[i] > thr]
            rec["blocks"][lab] = dict(n_window=int(len(sel)), x=[float(v) for v in x[sel]], w2=[float(v) for v in w], localized=locs)
        out["variants"][tag] = rec
        log(f"[D4] {tag}: fenêtre pairs {rec['blocks']['even']['n_window']} (localisés {rec['blocks']['even']['localized']}), impairs {rec['blocks']['odd']['n_window']} (localisés {rec['blocks']['odd']['localized']})")
        return rec

    # ---- (a1) M grossier 16 bandes
    for var in ("tot", "L", "NL"):
        M = np.load(os.path.join(WORK, "prep", f"M_coarse_{var}.npy")) * HA2EV
        H = sf.bloch_folded_hamiltonian(eps16, M, N_CELLS); blocks, cpl, herm = spectrum_blocks(H, par16.reshape(-1), E_D16)
        W = weights_bloch(blocks, C16, nG16, flat16, 16); record(f"a1_{var}", blocks, W, E_D16, dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
        del M, H
    # ---- (a2), (a3) M dense restreint aux 81 k
    Mm = np.load(dp["mfile"], mmap_mode="r"); Msub = np.array(Mm[np.ix_(np.arange(20), idx81, np.arange(20), idx81)]) * HA2EV; del Mm
    Mc = np.load(os.path.join(WORK, "prep", "M_coarse_tot.npy")) * HA2EV; Md16 = Msub[:16, :, :16, :]
    sv_c = np.linalg.svd(Mc.reshape(16 * nk, 16 * nk), compute_uv=False); sv_d = np.linalg.svd(Md16.reshape(16 * nk, 16 * nk), compute_uv=False)
    diag_c = np.array([[Mc[n, k, n, k].real for k in range(nk)] for n in range(16)]); diag_d = np.array([[Md16[n, k, n, k].real for k in range(nk)] for n in range(16)])
    out["coarse_vs_dense16"] = dict(max_abs_M_diff_eV=float(np.abs(np.abs(Mc) - np.abs(Md16)).max()), max_abs_M_eV=float(np.abs(Mc).max()),
                                    max_rel_SV_mismatch=float(np.abs(sv_c - sv_d).max() / sv_c.max()), diag_diff_max_eV=float(np.abs(diag_c - diag_d).max()),
                                    diag_diff_max_band_k=[int(x) for x in np.unravel_index(np.abs(diag_c - diag_d).argmax(), diag_c.shape)],
                                    diag_diff_max_bands_1to8_eV=float(np.abs(diag_c - diag_d)[:8].max()), eps16_vs_eps20_max_eV=float(np.abs(eps16 - eps20[:16]).max()),
                                    eps20_gap_16_17_min_eV=float((eps20[16] - eps20[15]).min()), n_C16_parity_dev_gt_1em4=int((np.abs(1 - np.abs(par16)) > 1e-4).sum()))
    log(f"[D4] M grossier (16 bandes) vs M dense restreint : max||M_c|-|M_d|| {out['coarse_vs_dense16']['max_abs_M_diff_eV']:.3f} eV, VS rel {out['coarse_vs_dense16']['max_rel_SV_mismatch']:.2e}, diag bandes 1-8 {out['coarse_vs_dense16']['diag_diff_max_bands_1to8_eV']:.2e} eV, gap 16-17 dense min {out['coarse_vs_dense16']['eps20_gap_16_17_min_eV']:.1e} eV")
    del Mc, Md16
    H = sf.bloch_folded_hamiltonian(eps20[:16], Msub[:16, :, :16, :], N_CELLS); blocks, cpl, herm = spectrum_blocks(H, par20[:16].reshape(-1), E_D20)
    W = weights_bloch(blocks, C20[:16], nG20, flat20, 16); record("a2_16", blocks, W, E_D20, dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
    H = sf.bloch_folded_hamiltonian(eps20, Msub, N_CELLS); blocks, cpl, herm = spectrum_blocks(H, par20.reshape(-1), E_D20)
    W = weights_bloch(blocks, C20, nG20, flat20, 20); record("a3_20", blocks, W, E_D20, dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
    # ---- (b) projection sur 5 WF, partie diagonale V^dag eps V
    Hk = np.einsum("kbw,kb,kbv->kwv", np.conj(V81), eps20.T, V81)                                 # (81, 5, 5)
    Mwk = Mbk_to_Mwk(Msub, U[idx81], Ud[idx81]); Hb = sf.kbasis_matrix(Hk, Mwk, k81, N_CELLS)
    parW = np.array([1.0 if w in WF_SIGMA else -1.0 for k in range(nk) for w in range(nw)])
    eK, vK = np.linalg.eigh(Hk[iK]); wpzK = np.abs(vK[WF_PZ_A]) ** 2 + np.abs(vK[WF_PZ_B]) ** 2      # poids p_z de chaque état propre à K
    pairK = np.argsort(-wpzK)[:2]; ED_b = float(eK[pairK].mean())
    blocks_b, cpl, herm = spectrum_blocks(Hb, parW, ED_b)
    def w_from_k(blocks):
        Wb = {}
        for lab in blocks:
            x = blocks[lab]["x"]; sel = np.where((x >= -3.0) & (x <= 1.0))[0]; w = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                c = np.zeros(nk * nw, complex); c[blocks[lab]["idx"]] = blocks[lab]["v"][:, j]
                dnk = np.einsum("kbw,kw->bk", V81, c.reshape(nk, nw))
                w[jj] = sf.folded_density_2d(dnk, C20, nG20, flat20, ng)[mask2].sum()
            Wb[lab] = dict(sel=sel, w2=w)
        return Wb
    Wb = w_from_k(blocks_b); record("b_5wf", blocks_b, Wb, ED_b, dict(coupling_even_odd=cpl, herm=herm, dim=int(Hb.shape[0]), E_D_def="paire centrale de V^dag eps V à K"))
    # ---- (c) repliement de M_W dense et de H(R)
    Z = load_cache_mwr(); Mwr = Z["Mwr_tot"]; Rn = Z["R"]; R_d = Z["R_d"]; C = load_cache_vloc(); Rloc = C["Rloc"]
    dev1, herm1, ok1 = sf.wannier_gate(Hwr, Rw, nd, N_SC); out["gates"]["gate1_fold_HR"] = dict(max_dev_eV=dev1, herm=herm1, ok=bool(ok1), threshold=1e-8)
    log(f"[D4] porte 1 (H(R) replié vs Wannier 81 k) : {dev1:.2e} eV -> {'OK' if ok1 else 'ÉCHEC'}")
    HS = sf.fold_hwr_to_supercell(Hwr, Rw, nd, N_SC); Hk1r = sf.kblocks_to_rbasis(Hk, k81, N_SC)
    Mall = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC); M3 = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC, R_local=Rloc)
    e_b = np.linalg.eigvalsh(0.5 * (Hb + Hb.conj().T)); e_c1 = np.linalg.eigvalsh(0.5 * ((Hk1r + Mall) + (Hk1r + Mall).conj().T))
    dev2 = float(np.abs(e_b - e_c1).max())
    # Porte 2 bis : les k du .save dense portent un bruit ~1e-7 (XML à 8 chiffres) qui entache les phases de la double TF de
    # production (R jusqu'à 13) ; M_W est recalculé ici avec les k ramenés exactement sur m/27 pour tester le repliement seul.
    knoise = dict(k27_max_dev=float(np.abs(27 * k27 - np.rint(27 * k27)).max()), k81_max_dev=float(np.abs(9 * k81 - np.rint(9 * k81)).max()))
    Mfull = matrix_io.load_M_checked(dp["mfile"], require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV)
    Mwk_full = Mbk_to_Mwk(Mfull, U, Ud); del Mfull
    Mwr_snap, R_snap = Mwk_to_Mwr(Mwk_full, np.rint(27 * k27) / 27.0, S["MP"]); del Mwk_full
    Rn_snap, Rd_snap = lt.recenter_mwr(Mwr_snap, R_snap, S["MP"])
    Mall_snap = sf.fold_mwr_to_supercell(Mwr_snap, Rn_snap, N_SC)
    e_c1s = np.linalg.eigvalsh(0.5 * ((Hk1r + Mall_snap) + (Hk1r + Mall_snap).conj().T)); dev2s = float(np.abs(e_b - e_c1s).max())
    Vsnap, _ = lt.extract_V_loc(Mwr_snap, Rn_snap, Rloc)
    out["gates"]["gate2_fold_Mwr"] = dict(max_dev_eV_production_k=dev2, max_dev_eV_snapped_k=dev2s, ok=bool(dev2s < 1e-9), literal_ok=bool(dev2 < 1e-9), threshold=1e-9,
                                          definition="eig[F^+ (V^+ eps V) F + fold(M_W)] vs eig[(b)] (même partie à un corps) ; M_W de production (k du XML dense) et M_W recalculé avec k = m/27 exacts",
                                          k_noise=knoise, R_d_snap=Rd_snap.tolist(),
                                          Vloc_prod_vs_snap=dict(max_abs_diff_eV=float(np.abs(C["V_tot"] - Vsnap).max()), max_abs_Vloc_eV=float(np.abs(C["V_tot"]).max()),
                                                                 onsite_pz_prod=float(C["V_tot"][cell_index(Rloc, (0, 0, 0)) * nw + WF_PZ_A, cell_index(Rloc, (0, 0, 0)) * nw + WF_PZ_A].real),
                                                                 onsite_pz_snap=float(Vsnap[cell_index(Rloc, (0, 0, 0)) * nw + WF_PZ_A, cell_index(Rloc, (0, 0, 0)) * nw + WF_PZ_A].real)),
                                          Mall_prod_vs_snap_max_abs_diff_eV=float(np.abs(Mall - Mall_snap).max()))
    del Mwr_snap, Mall_snap, Vsnap
    Hwk81, Ew81, _ = Hwr_to_Hwk(Hwr, Rw, k81, ndegen=nd)
    out["gates"]["one_body_residual"] = dict(max_abs_Hwk_minus_VepsV_eV=float(np.abs(Hwk81 - Hk).max()), max_eig_dev_eV=float(np.abs(np.sort(np.linalg.eigvalsh(HS)) - np.sort(np.linalg.eigvalsh(Hk1r))).max()))
    log(f"[D4] porte 2 (repliement de M_W) : {dev2:.2e} eV avec les k du XML (bruit {knoise['k27_max_dev']:.1e}), {dev2s:.2e} eV avec k = m/27 exacts -> {'OK' if dev2s < 1e-9 else 'ÉCHEC'} ; "
        f"V_loc production vs k exacts : max|diff| {out['gates']['gate2_fold_Mwr']['Vloc_prod_vs_snap']['max_abs_diff_eV']:.2e} eV ; résidu H(R) vs V^+epsV aux 81 k : {out['gates']['one_body_residual']['max_abs_Hwk_minus_VepsV_eV']:.2e} eV")
    if not (ok1 and dev2s < 1e-9):
        save_json(os.path.join(d, "d4_results.json"), out); open(os.path.join(d, "D4_FAIL"), "w").write("gates"); sys.exit(1)
    e_call_literal = np.linalg.eigvalsh(0.5 * ((HS + Mall) + (HS + Mall).conj().T))
    out["gates"]["literal_b_vs_call"] = dict(max_dev_eV=float(np.abs(e_b - e_call_literal).max()), note="(b) avec V^+epsV contre (c-all) avec H(R) : inclut le résidu d'interpolation")
    def w_from_r(blocks):
        Wr = {}
        for lab in blocks:
            x = blocks[lab]["x"]; sel = np.where((x >= -3.0) & (x <= 1.0))[0]; w = np.zeros(len(sel)); wsite = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                cr = np.zeros(N_CELLS * nw, complex); cr[blocks[lab]["idx"]] = blocks[lab]["v"][:, j]
                ck = sf.rvec_to_kvec(cr, k81, N_SC, nw).reshape(nk, nw) * np.exp(-2j * np.pi * (k81 @ np.asarray(R_d, float)))[:, None]   # étiquettes recentrées r = r_réel - R_d
                dnk = np.einsum("kbw,kw->bk", V81, ck)
                w[jj] = sf.folded_density_2d(dnk, C20, nG20, flat20, ng)[mask2].sum()
                p = np.abs(cr.reshape(N_CELLS, nw)) ** 2; c0 = sf.cell_of(np.array([0, 0, 0]), N_SC)
                wsite[jj] = p[c0, WF_PZ_A] + sum(p[sf.cell_of(np.array(R), N_SC), WF_PZ_B] for R in NN_CELLS) + p[c0, :3].sum()
            Wr[lab] = dict(sel=sel, w2=w, wsite=wsite)
        return Wr
    for tag, Msc in (("c_all", Mall), ("c_3", M3)):
        H = HS + Msc; blocks, cpl, herm = spectrum_blocks(H, parW, E_DW); Wr = w_from_r(blocks)
        rec = record(tag, blocks, Wr, E_DW, dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
        for lab in ("even", "odd"): rec["blocks"][lab]["w_site"] = [float(v) for v in Wr[lab]["wsite"]]
    # ---- comparaisons : QE (D1 aligné) et échelle (a1) -> (a2) -> (a3) -> (b) ; (c-all) -> (c-3)
    qe_rec = dict(blocks={})
    for lab, sgn in (("even", 1), ("odd", -1)):
        m = (qe["parity"] * sgn > 0); qe_rec["blocks"][lab] = dict(n_window=int(m.sum()), x=[float(v) for v in qe["x"][m]], w2=[float(v) for v in qe["w2"][m]],
                                                                  localized=[(float(x), float(w)) for x, w in zip(qe["x"][m], qe["w2"][m]) if w > thr])
    out["variants"]["QE_D1"] = qe_rec
    # décalage rigide ajusté hors fenêtre (états < E_D - 4 eV), sur l'échelle epsilon - E_D
    fits = {}
    eqe = qe_io.get_eigenvalues(SC_D)[:, 0] * HA2EV - d1["alignment"]["D"]["shifts_eV"]["1.0"] - E_D_SC     # tous les états QE, alignés Lu, rel. E_D
    # spectres complets par variante (recalcul bon marché) pour l'ajustement du décalage rigide hors fenêtre
    spectra = {}
    M = np.load(os.path.join(WORK, "prep", "M_coarse_tot.npy")) * HA2EV; H = sf.bloch_folded_hamiltonian(eps16, M, N_CELLS); spectra["a1_tot"] = np.linalg.eigvalsh(0.5 * (H + H.conj().T)) - E_D16; del M, H
    H = sf.bloch_folded_hamiltonian(eps20[:16], Msub[:16, :, :16, :], N_CELLS); spectra["a2_16"] = np.linalg.eigvalsh(0.5 * (H + H.conj().T)) - E_D20
    H = sf.bloch_folded_hamiltonian(eps20, Msub, N_CELLS); spectra["a3_20"] = np.linalg.eigvalsh(0.5 * (H + H.conj().T)) - E_D20
    spectra["b_5wf"] = e_b - ED_b; spectra["c_all"] = np.linalg.eigvalsh(0.5 * ((HS + Mall) + (HS + Mall).conj().T)) - E_DW; spectra["c_3"] = np.linalg.eigvalsh(0.5 * ((HS + M3) + (HS + M3).conj().T)) - E_DW
    for tag, sp in spectra.items():
        s, resid, n = al.rigid_shift_fit(eqe, sp, -4.0); fits[tag] = dict(shift_eV=s, max_residual_eV=resid, n_states=n)
        log(f"[D4] décalage rigide résiduel {tag} (états QE < E_D - 4 eV, n = {n}) : {s*1e3:+.1f} meV (résidu max {resid*1e3:.1f} meV)")
    out["rigid_shift_fit"] = dict(definition="médiane(x_QE - x_modèle) sur les n états les plus bas, x = epsilon - E_D (QE aligné Lu 1 Å + E_D parfaite ; modèle : son propre E_D à K) ; = écart entre l'alignement de D1 et l'alignement sur E_D", fits=fits, shift_Lu_D_eV=d1["alignment"]["D"]["shifts_eV"]["1.0"])
    # marches
    ladder = ["QE_D1", "a1_tot", "a2_16", "a3_20", "b_5wf"]; steps = {}
    def match(loc_a, loc_b):
        rows = []
        for xa, wa in loc_a:
            if loc_b:
                j = int(np.argmin([abs(xb - xa) for xb, _ in loc_b])); xb, wb = loc_b[j]; rows.append(dict(x_from=xa, w_from=wa, x_to=xb, w_to=wb, dx=xb - xa, moved_gt_0p1=bool(abs(xb - xa) > 0.1)))
            else:
                rows.append(dict(x_from=xa, w_from=wa, x_to=None, w_to=None, dx=None, moved_gt_0p1=None))
        return rows
    for i in range(len(ladder) - 1):
        A_, B_ = ladder[i], ladder[i + 1]; steps[f"{A_}->{B_}"] = {}
        for lab in ("even", "odd"):
            la = out["variants"][A_]["blocks"][lab]["localized"]; lb = out["variants"][B_]["blocks"][lab]["localized"]
            xa = np.array(out["variants"][A_]["blocks"][lab]["x"]); xb = np.array(out["variants"][B_]["blocks"][lab]["x"])
            steps[f"{A_}->{B_}"][lab] = dict(n_from=len(xa), n_to=len(xb), localized_match=match(la, lb),
                                             median_abs_shift_sorted_eV=(float(np.median(np.abs(np.sort(xa) - np.sort(xb)))) if len(xa) == len(xb) else None))
    for lab in ("even", "odd"):
        la = out["variants"]["c_all"]["blocks"][lab]["localized"]; lb = out["variants"]["c_3"]["blocks"][lab]["localized"]
        steps.setdefault("c_all->c_3", {})[lab] = dict(localized_match=match(la, lb), n_from=out["variants"]["c_all"]["blocks"][lab]["n_window"], n_to=out["variants"]["c_3"]["blocks"][lab]["n_window"])
        la = out["variants"]["b_5wf"]["blocks"][lab]["localized"]; lb = out["variants"]["c_all"]["blocks"][lab]["localized"]
        steps.setdefault("b_5wf->c_all", {})[lab] = dict(localized_match=match(la, lb))
    out["steps"] = steps
    np.savez(os.path.join(d, "d4_spectra.npz"), **{k: v for k, v in spectra.items()}, qe_x_all=eqe)
    save_json(os.path.join(d, "d4_results.json"), out)
    # figure : échelle des niveaux dans la fenêtre, par parité, poids w2 en taille
    plt, pal = fig_style()
    order = ["QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"]
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 4.2), sharey=True)
    for ax, lab, ttl in zip(axs, ("odd", "even"), ("(a) impairs (π)", "(b) pairs (σ)")):
        for i, tag in enumerate(order):
            if tag not in out["variants"]: continue
            b = out["variants"][tag]["blocks"][lab]; x = np.array(b["x"]); w = np.array(b["w2"])
            ax.scatter(np.full(len(x), i), x, s=2 + 60 * w, color=pal.NAVY if tag != "QE_D1" else pal.ORANGE, alpha=0.7, lw=0)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, fontsize=6); ax.set_title(ttl); ax.set_ylim(-3, 1); ax.axhline(0, color=pal.MUTED, lw=0.5)
    axs[0].set_ylabel(r"$\varepsilon - E_D$ (eV)"); fig.tight_layout(); savefig(fig, "d4_ladder")


# ----------------------------------------------------------------------------------------------- J5 : D3
def cmd_d3(a):
    S = setup(); d = ensure("d3"); out = dict(gate={}, scan={}); cfg = S["cfg"]; eta, npe = cfg["eta_eV"], cfg["ne_per_eta"]
    Hwr, Rw, nd, E_D = S["Hwr"], S["Rw"], S["nd"], S["E_D"]
    C = load_cache_vloc(); V = dict(tot=C["V_tot"], L=C["V_L"], NL=C["V_NL"]); Rloc = C["Rloc"]; nL = len(Rloc)
    ip, isg = C["idx_pi"], C["idx_sigma"]; grp = index_groups(Rloc)
    blocks = dict(pi=ip, sigma=isg, full=None)
    # phi_K (paire pi, bandes 3-4) comme resonance_metrics.py l. 50, 86-87
    _, EK, UK = Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); ph = lt._phase(K_RED[None], Rloc)
    phi_K = np.einsum("kL,kwn->knLw", ph, UK, optimize=True).reshape(1, NW, nL * NW)[0][3:5]
    out["K"] = dict(E_pi_pistar=[float(EK[0, 3]), float(EK[0, 4])], E_D=E_D)
    # ---- porte de régression : alpha = 1, matrice complète, 300², fenêtre +-3 eV
    k300 = lt.mp_grid(300, 300, 1); Hwk300, _, _ = Hwr_to_Hwk(Hwr, Rw, k300, ndegen=nd)
    egF = egrid(E_D, -3.0, 3.0, eta, npe); g0F = g0_cached("g0_real_nk300", Hwk300, k300, Rloc, egF, eta)
    r = pc.det_eig_criterion(V["tot"], g0F, None, vectors=True); xF = egF - E_D
    jd = int(np.argmin(r["logdet_rel"])); jl = int(np.argmin(r["minlam"]))
    gate = dict(min_det_rel=float(np.exp(r["logdet_rel"][jd])), at=float(xF[jd]), min_abs_lam=float(r["minlam"][jl]), lam=complex(r["lam_min"][jl]), lam_at=float(xF[jl]),
                expected=dict(min_det_rel=2.091e-4, at=-2.530, lam=[0.0003, 0.0108]),
                minima_det=[(float(xF[i]), float(np.exp(r["logdet_rel"][i])), float(r["minlam"][i])) for i in pc.local_minima(r["logdet_rel"], xF)])
    okg = (abs(gate["min_det_rel"] / 2.091e-4 - 1) < 5e-3) and (abs(gate["at"] + 2.530) < 1e-3) and (abs(gate["lam"] - (0.0003 + 0.0108j)) < 1e-4)
    gate["ok"] = bool(okg); out["gate"] = gate
    log(f"[D3] porte : min|det|/max {gate['min_det_rel']:.4e} à {gate['at']:+.3f} (attendu 2.091e-4 à -2.530) ; lambda {gate['lam']:.5f} à {gate['lam_at']:+.3f} -> {'OK' if okg else 'ÉCHEC'}")
    if not okg:
        save_json(os.path.join(d, "d3_results.json"), out); open(os.path.join(d, "D3_FAIL"), "w").write("gate"); sys.exit(1)
    # caractérisation du minimum à -2.530 (alpha = 1, complète, 300²)
    vm = r["vec_min"][jd]; wpi = float((np.abs(vm[ip]) ** 2).sum()); wsg = float((np.abs(vm[isg]) ** 2).sum())
    out["gate"]["minimum_m2p530"] = dict(weights=pc.eigvec_weights(vm, grp), weight_pi_block=wpi, weight_sigma_block=wsg, lam=complex(r["lam_min"][jd]))
    log(f"[D3] minimum -2.530 : poids bloc pi {wpi:.4f}, sigma {wsg:.4f}, groupes {out['gate']['minimum_m2p530']['weights']}")
    # ---- balayage alpha
    alphas = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0]; curves = {}
    for nk in (300, 600):
        if nk == 300:
            win = egF <= E_D + 1.0 + 1e-9; eg = egF[win]; g0 = g0F[win]
        else:
            k600 = lt.mp_grid(600, 600, 1); Hwk600, _, _ = Hwr_to_Hwk(Hwr, Rw, k600, ndegen=nd); eg = egrid(E_D, -3.0, 1.0, eta, npe)
            g0 = g0_cached("g0_real_nk600", Hwk600, k600, Rloc, eg, eta)
        x = eg - E_D; out["scan"][f"nk{nk}"] = {}
        for var in ("tot", "L", "NL"):
            for al_ in alphas:
                Va = {"tot": al_ * V["tot"], "L": al_ * V["L"] + V["NL"], "NL": V["L"] + al_ * V["NL"]}[var]
                for bl, idx in blocks.items():
                    t0 = time.time(); rr = pc.det_eig_criterion(Va, g0, idx, vectors=True)
                    jd = int(np.argmin(rr["logdet_rel"])); jl = int(np.argmin(rr["minlam"]))
                    vfull = np.zeros(nL * NW, complex); vfull[idx if idx is not None else np.arange(nL * NW)] = rr["vec_min"][jl]
                    ov = pc.branch_overlaps(rr["vec_min"]); ii, xc = pc.sign_changes(rr["lam_min"].real, x)
                    rec = dict(min_det_rel=float(np.exp(rr["logdet_rel"][jd])), det_at=float(x[jd]),
                               minima_det=[(float(x[i]), float(np.exp(rr["logdet_rel"][i]))) for i in pc.local_minima(rr["logdet_rel"], x)],
                               min_abs_lam=float(rr["minlam"][jl]), lam=complex(rr["lam_min"][jl]), lam_at=float(x[jl]), weights=pc.eigvec_weights(vfull, grp),
                               weight_pi=float((np.abs(vfull[ip]) ** 2).sum()), weight_sigma=float((np.abs(vfull[isg]) ** 2).sum()),
                               roots_Re_lam_min=[(float(e), float(ov[i]), bool(ov[i] < 0.9)) for i, e in zip(ii, xc)], n_branch_jumps=int((ov < 0.9).sum()))
                    if bl in ("pi", "full"):
                        Vb = Va if idx is None else Va[np.ix_(idx, idx)]; gb = g0 if idx is None else g0[:, idx][:, :, idx]
                        tc = pc.local_t_cache(Vb, gb); PK = phi_K if idx is None else phi_K[:, idx]
                        _, tr = pc.tbar_pair(tc, PK); jt = int(np.argmax(-tr.imag))
                        rec["peak_minus_Im_Tbar"] = float(x[jt]); rec["Tbar_at_ED"] = complex(np.interp(0.0, x, tr.real) + 1j * np.interp(0.0, x, tr.imag)); rec["max_minus_Im_Tbar"] = float(-tr.imag[jt])
                        curves[f"nk{nk}_{var}_a{al_}_{bl}_ImTbar"] = -tr.imag
                        del tc
                    curves[f"nk{nk}_{var}_a{al_}_{bl}_detrel"] = np.exp(rr["logdet_rel"]); curves[f"nk{nk}_{var}_a{al_}_{bl}_minlam"] = rr["minlam"]
                    out["scan"][f"nk{nk}"][f"{var}_a{al_}_{bl}"] = rec
                    log(f"[D3] nk{nk} {var} α={al_:g} {bl}: min|det|/max {rec['min_det_rel']:.2e} à {rec['det_at']:+.3f}; min|λ| {rec['min_abs_lam']:.4f} ({rec['lam']:.4f}) à {rec['lam_at']:+.3f}; racines {[round(e,3) for e,_,_ in rec['roots_Re_lam_min']]}"
                        + (f"; pic -Im T̄ {rec['peak_minus_Im_Tbar']:+.3f}" if 'peak_minus_Im_Tbar' in rec else "") + f" ({time.time()-t0:.0f} s)")
        curves[f"nk{nk}_x"] = x
    np.savez(os.path.join(d, "d3_curves.npz"), **curves); save_json(os.path.join(d, "d3_results.json"), out)
    # figures : |det|/max et min|lambda| vs energie pour les alpha (complète et pi), 600²
    plt, pal = fig_style()
    for bl in ("full", "pi"):
        fig, ax = plt.subplots(1, 3, figsize=(6.5, 3.0))
        x = curves["nk600_x"]
        for al_, c in zip(alphas, pal.CYCLE + [pal.REF, pal.MUTED]):
            ax[0].plot(x, curves[f"nk600_tot_a{al_}_{bl}_detrel"], color=c, lw=0.8, label=f"α = {al_:g}"); ax[1].plot(x, curves[f"nk600_tot_a{al_}_{bl}_minlam"], color=c, lw=0.8)
            ax[2].plot(x, curves[f"nk600_tot_a{al_}_{bl}_ImTbar"], color=c, lw=0.8)
        ax[0].set_yscale("log"); ax[0].set_ylabel(r"$|\det|/\max$"); ax[1].set_yscale("log"); ax[1].set_ylabel(r"$\min_i |\lambda_i|$"); ax[2].set_ylabel(r"$-\mathrm{Im}\,\bar T(K)$ (eV)")
        for a_ in ax: a_.set_xlabel(r"$\varepsilon - E_D$ (eV)")
        ax[0].legend(fontsize=5, ncol=2); ax[0].set_title({"full": "matrice complète", "pi": "bloc π"}[bl] + ", 600²", fontsize=8)
        fig.tight_layout(); savefig(fig, f"d3_alpha_{bl}")


# ----------------------------------------------------------------------------------------------- J5a : g0 600²
def cmd_g0(a):
    """Calcule et met en cache g0 du vrai H(R) à 600² sur [E_D - 3, E_D + 1] (BLAS multi-fils) ; d3 le relit avec 1 fil BLAS."""
    S = setup(); cfg = S["cfg"]; eta, npe = cfg["eta_eV"], cfg["ne_per_eta"]; C = load_cache_vloc(); Rloc = C["Rloc"]
    k600 = lt.mp_grid(600, 600, 1); Hwk600, _, _ = Hwr_to_Hwk(S["Hwr"], S["Rw"], k600, ndegen=S["nd"])
    g0_cached("g0_real_nk600", Hwk600, k600, Rloc, egrid(S["E_D"], -3.0, 1.0, eta, npe), eta)




# ----------------------------------------------------------------------------------------------- figures D6 (relecture des courbes)
def cmd_figs(a):
    """Regénère les figures D6 depuis d6/d6_curves.npz avec dix teintes distinctes (CMAP_SEQ) ; aucune donnée recalculée."""
    curves = dict(np.load(os.path.join(WORK, "d6", "d6_curves.npz"))); plt, pal = fig_style()
    cols = [pal.CMAP_SEQ(t) for t in np.linspace(0.3, 1.0, len(U_LIST))]
    for kind in ("toy", "real"):
        fig, ax = plt.subplots(1, 2, figsize=(6.5, 3.4))
        for Uv, c in zip(U_LIST, cols):
            z = curves[f"{kind}_nk1200_U{Uv}"]
            ax[0].plot(z[0], z[1], color=c, lw=1.0, label=f"U = {Uv:g}"); ax[1].plot(z[0], z[2], color=c, lw=1.0)
        ax[0].set_yscale("log"); ax[0].set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax[0].set_ylabel(r"$|\det[1 - V g^{(0)}]|$"); ax[0].set_title("(a) 1200², " + ("modèle" if kind == "toy" else "H(R) de production"), fontsize=9)
        ax[1].axhline(0, color=pal.MUTED, lw=0.6); ax[1].set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax[1].set_ylabel(r"Re $\lambda$ (valeur propre non triviale)"); ax[1].set_ylim(-3, 3); ax[1].set_title("(b)", fontsize=9)
        ax[0].legend(fontsize=6, ncol=2); fig.tight_layout(); savefig(fig, f"d6_{kind}_1200")
    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    for nk, ls in ((300, "--"), (1200, "-")):
        for Uv, c in zip([6.617, 12.0, 27.0], [pal.NAVY, pal.ORANGE, pal.GREEN]):
            z = curves[f"toy_nk{nk}_U{Uv}"]; ax.plot(z[0], z[2], color=c, ls=ls, lw=1.0, label=f"U = {Uv:g}, {nk}²")
    ax.axhline(0, color=pal.MUTED, lw=0.6); ax.set_ylim(-2, 2); ax.set_xlabel(r"$\varepsilon - E_D$ (eV)"); ax.set_ylabel(r"Re $\lambda = 1 - U g_0$"); ax.legend(fontsize=6, ncol=2)
    fig.tight_layout(); savefig(fig, "d6_toy_grid")

# ----------------------------------------------------------------------------------------------- tables du rapport
def _f(x, nd=3):
    return "—" if x is None else (f"{x:+.{nd}f}" if isinstance(x, (int, float)) else str(x))


def _c(z, nd=4):
    if z is None: return "—"
    if isinstance(z, (list, tuple)): z = complex(z[0], z[1])
    return f"{z.real:+.{nd}f} {z.imag:+.{nd}f} i"


def cmd_tables(a):
    """Tableaux Markdown des sections D6, D1, D5, R, D4, D3 à partir des json (sections/*.tables.md)."""
    d = ensure("sections"); L = []
    # ---- D6
    z = json.load(open(os.path.join(WORK, "d6", "d6_results.json"))); L.append("### D6.1 — jouet, 1200² (valeurs attendues du prompt entre parenthèses ; écarts = calculé − attendu)\n")
    L.append("| U (eV) | racines de Re λ (eV) | racine près de E_D (attendu) | écart | argmin \\|det\\| (attendu) | écart | \\|det\\| min (attendu) | écart | \\|det\\| à −2,5 eV | λ à E_D |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    for U in U_LIST:
        v = z["D61"]["nk1200"]["U"][str(U)]; e = v["expected_1200"]
        L.append(f"| {U:g} | {', '.join(f'{r:+.3f}' for r in v['roots_Re_lam_min']) or 'aucune'} | {_f(v['root_near_ED'])} ({'aucune' if e[0] is None else f'{e[0]:+.3f}'}) | {_f(v.get('dev_root'), 4)} | {v['argmin_det']:+.4f} ({e[1]:+.3f}) | {_f(v.get('dev_argmin_det'), 4)} | {v['det_min']:.3f} ({e[2]}) | {_f(v.get('dev_det_min'), 4)} | {v['det_at_m2p5']:.3f} | {_c(v['lam_min_at_ED'])} |")
    rs = z["D61"]["nk1200"]["removed_site"]
    L.append(f"\nSite retiré (U = 1e3 + −H sur les trois liaisons ; support {rs['support']}) à 1200² : argmin \\|det\\| = {rs['argmin_det']:+.1e} eV, \\|det\\| = {rs['det_min']:.3f} ; min sur les énergies de \\|λ_i\\| pour les quatre valeurs propres non triviales : {', '.join(f'{x:.4f}' for x in rs['min_abs_each_lambda_over_energies'])} ; aucun changement de signe de Re λ_min. À 300² : \\|det\\| min {z['D61']['nk300']['removed_site']['det_min']:.3f} à {z['D61']['nk300']['removed_site']['argmin_det']:+.1e} eV, min \\|λ\\| {z['D61']['nk300']['removed_site']['min_over_energies_of_min_abs_lam']:.4f}.\n")
    L.append("### D6.1 — jouet, 300² (bruit de grille : toutes les racines rapportées)\n")
    L.append("| U (eV) | racines de Re λ (eV) | argmin \\|det\\| | \\|det\\| min | min \\|λ\\| (position) |"); L.append("|---|---|---|---|---|")
    for U in U_LIST:
        v = z["D61"]["nk300"]["U"][str(U)]; L.append(f"| {U:g} | {', '.join(f'{r:+.3f}' for r in v['roots_Re_lam_min']) or 'aucune'} | {v['argmin_det']:+.4f} | {v['det_min']:.3f} | {v['global_min_abs_lam'][1]:.3f} ({v['global_min_abs_lam'][0]:+.3f}) |")
    for nk in ("300", "1200"):
        r = z["D62"][f"nk{nk}"]; L.append(f"\n### D6.2 — vrai H(R), {nk}² (U sur la p_z du site lacunaire ; énergies ε − E_D, E_D Wannier)\n")
        L.append(f"g⁽⁰⁾ p_z–p_z sur site à E_D − 1 eV : vrai H(R) {_c(r['g0_pz_onsite_at_ED_m1'], 5)} eV⁻¹ ; jouet à −1 eV, même grille : {_c(r['toy_g0_pz_onsite_at_m1'], 5)} eV⁻¹. Maximum de Re g⁽⁰⁾ (vrai) : {r['Re_g0_max'][1]:.4f} eV⁻¹ à {r['Re_g0_max'][0]:+.4f} eV (U_c = {1/r['Re_g0_max'][1]:.2f} eV) ; jouet : {z['D61'][f'nk{nk}']['Re_g0_max'][1]:.4f} à {z['D61'][f'nk{nk}']['Re_g0_max'][0]:+.4f} eV (U_c = {1/z['D61'][f'nk{nk}']['Re_g0_max'][1]:.2f} eV).\n")
        L.append("| U (eV) | racines de Re λ (eV) | argmin \\|det\\| | \\|det\\| min | \\|det\\| à −2,5 eV | λ à E_D | jouet même grille : racine près de E_D / argmin \\|det\\| |"); L.append("|---|---|---|---|---|---|---|")
        for U in U_LIST:
            v = r["U"][str(U)]; t = z["D61"][f"nk{nk}"]["U"][str(U)]
            L.append(f"| {U:g} | {', '.join(f'{x:+.3f}' for x in v['roots_Re_lam_min']) or 'aucune'} | {v['argmin_det']:+.4f} | {v['det_min']:.3f} | {v['det_at_m2p5']:.3f} | {_c(v['lam_min_at_ED'])} | {_f(t['root_near_ED'])} / {t['argmin_det']:+.4f} |")
        m = r["Mloc_pi_block"]
        roots_txt = ", ".join("%+.3f (recouvrement %.2f%s)" % (e, o, ", saut" if j else "") for e, o, j in m["roots_Re_lam_min"]) or "aucun"
        mind = ", ".join("%+.3f (%.4f)" % (x, y) for x, y in m["minima_det"]); minl = ", ".join("%+.3f (%.4f)" % (x, y) for x, y, _ in m["minima_abs_lam"])
        L.append(f"\nBloc π de M_loc complet (58 × 58), {nk}² : min \\|det\\|/max = {m['det_rel_min']:.4f} à {m['argmin_det']:+.4f} eV (\\|det\\| non normalisé {m['det_abs_min']:.4f}) ; minima locaux de \\|det\\|/max : {mind} ; min \\|λ\\| = {m['global_min_abs_lam'][1]:.4f} à {m['global_min_abs_lam'][0]:+.4f} eV, λ = {_c(m['global_min_abs_lam'][2])} ; minima locaux de \\|λ_min\\| : {minl} ; changements de signe de Re λ_min : {roots_txt}.\n")
    # ---- D1
    z = json.load(open(os.path.join(WORK, "d1", "d1_results.json"))); L.append("\n### D1 — alignement (Lu 2019) et fenêtres\n")
    L.append("| géométrie | atome le plus loin (1-based) / distance (Å) / déplacement vs parfaite (Å) | décalage 0,5 Å (meV) | décalage 1,0 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ 3D (meV) | points dans la sphère (d / p) |"); L.append("|---|---|---|---|---|---|")
    for k, v in z["alignment"].items():
        L.append(f"| {k} | {v['i_far_1based']} / {v['dist_far_A']:.3f} / {v['far_displacement_A']:.2e} | {v['shifts_eV']['0.5']*1e3:+.2f} | {v['shifts_eV']['1.0']*1e3:+.2f} | {v['mean3d_diff_eV']*1e3:+.2f} | {v['means']['1.0'][2]} / {v['means']['1.0'][3]} (1,0 Å) ; {v['means']['0.5'][2]} / {v['means']['0.5'][3]} (0,5 Å) |")
    L.append("\n| géométrie | E_F (eV) | décalage appliqué (1,0 Å, eV) | bandes de la fenêtre | états | pairs (σ) | impairs (π) | max\\|1 − \\|⟨σ_h⟩\\|\\| | max\\|z − z₀\\| (Å) |"); L.append("|---|---|---|---|---|---|---|---|---|")
    for t, s_ in z["states"].items():
        L.append(f"| {s_['name']} | {s_['E_F']:.5f} | {s_['shift_eV']:+.5f} | {s_['bands_1based'][0]}–{s_['bands_1based'][1]} | {s_['n_window']} | {s_['n_even']} | {s_['n_odd']} | {s_['max_dev_parity']:.1e} | {s_['max_z_dev_A']:.1e} |")
    th = z["localization_threshold"]; L.append(f"\nSeuil de localisation : w₂ > {th['thr_w2']:.4f} (3 × moyenne des {z['states']['P']['n_window']} états de la fenêtre de la parfaite, {th['mean_w2_P']:.4f} ; max dans la parfaite {th['max_w2_P']:.4f} ; fraction d'aire du disque {z['disc_area_fraction']['r2']:.4f}).\n")
    L.append("| géométrie | bande | ε − E_D (eV, aligné) | ε − E_F(défaut) (eV) | ⟨σ_h⟩ | w₂ (2 Å) | w₁ (1 Å) | projwfc : p_z 3 voisins | s+p_x,y atome 81 | s+p_x,y atomes 64+80 | p_z 6 seconds voisins | \\|ψ\\|² projeté |"); L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for t, rows in z["localized"].items():
        for r in rows:
            pj = r.get("projwfc"); L.append(f"| {z['states'][t]['name']} | {r['band']} | {r['x']:+.3f} | {r['relEF']:+.3f} | {r['parity']:+.0f} | {r['w2']:.3f} | {r['w1']:.3f} | " + (f"{pj['nn_pz']:.3f} | {pj['iso_s_pxy']:.3f} | {pj['pair_s_pxy']:.3f} | {pj['nnn_pz']:.3f} | {pj['psi2']:.3f} |" if pj else "— | — | — | — | — |"))
    # ---- D5
    z = json.load(open(os.path.join(WORK, "d5", "d5_results.json"))); L.append("\n### D5 — ΔV^L par taille (eV sauf mention ; ligne a₁ de la lacune à la demi-boîte ; plans frontière à la demi-boîte)\n")
    L.append("| N | A_sc (Å²) | atome loin / distance (Å) | ⟨ΔV⟩_3D = composante G = 0 (meV) | ΔV au site (eV) | bord ligne a₁ brut / − moyenne / aligné (meV) | plans frontière max\\|ΔV − c\\|, tous z : brut / − moyenne / aligné (meV) | idem, plan des atomes (meV) | décalage Lu 0,5 / 1,0 Å (meV) |"); L.append("|---|---|---|---|---|---|---|---|---|")
    for Sz in ["5x5", "6x6", "7x7", "8x8", "9x9", "10x10", "11x11", "12x12"]:
        if Sz not in z: continue
        r = z[Sz]; L.append(f"| {r['N']} | {r['A_sc_A2']:.1f} | {r['far_atom_1based']} / {r['dist_far_A']:.2f} | {r['mean3d']*1e3:+.2f} | {r['site']:+.3f} | {r['end_raw']*1e3:+.1f} / {r['end_sub']*1e3:+.1f} / {r['end_aligned']*1e3:+.1f} | {r['boundary_raw']*1e3:.1f} / {r['boundary_sub']*1e3:.1f} / {r['boundary_aligned']*1e3:.1f} | {r['boundary_plane_raw']*1e3:.1f} / {r['boundary_plane_sub']*1e3:.1f} / {r['boundary_plane_aligned']*1e3:.1f} | {r['shift_Lu_eV']['0.5']*1e3:+.2f} / {r['shift_Lu_eV']['1.0']*1e3:+.2f} |")
    # ---- R
    z = json.load(open(os.path.join(WORK, "r", "r_results.json"))); L.append("\n### R — ΔV = V(géométrie) − V_p(9×9 parfaite non relaxée) (eV)\n")
    L.append("| géométrie | décalage Lu utilisé (eV) | ⟨ΔV⟩_3D (meV) | ΔV au site (eV) | ligne a₁ : t = 0 / 0,125 / 0,25 / 0,375 / 0,5 (eV) | bord ligne brut / aligné (meV) | plans frontière max\\|ΔV\\| brut / aligné (meV) | idem plan des atomes (meV) |"); L.append("|---|---|---|---|---|---|---|---|")
    for k in ("D_unrelaxed", "N1", "N2u", "N2d"):
        r = z["dV_R1"][k]; L.append(f"| {k} | {r['shift_used_eV']:+.4f} | {r['mean3d']*1e3:+.2f} | {r['site']:+.3f} | {' / '.join(f'{x:+.3f}' for x in r['line_end_points'])} | {r['end_raw']*1e3:+.1f} / {r['end_aligned']*1e3:+.1f} | {r['boundary_raw']*1e3:.1f} / {r['boundary_aligned']*1e3:.1f} | {r['boundary_plane_raw']*1e3:.1f} / {r['boundary_plane_aligned']*1e3:.1f} |")
    u = z["dV_R1"]["N2_up_minus_dw"]; L.append(f"\nV_↑ − V_↓ (R1 nspin2) : max\\|·\\| = {u['max_abs_eV']:.3f} eV, moyenne 3D {u['mean3d_eV']*1e3:+.2f} meV, au site {u['site_eV']:+.3f} eV, plans frontière max {u['boundary_max_abs_eV']*1e3:.2f} meV.\n")
    # ---- D4
    p4 = os.path.join(WORK, "d4", "d4_results.json")
    if os.path.exists(p4):
        z = json.load(open(p4)); g = z["gates"]; L.append("\n### D4 — portes et références\n")
        L.append(f"Porte 1 (H(R) replié contre valeurs de Wannier aux 81 k) : {g['gate1_fold_HR']['max_dev_eV']:.2e} eV (seuil 1e-8) → {'OK' if g['gate1_fold_HR']['ok'] else 'ÉCHEC'}. "
                 f"Porte 2 (repliement de M_W, même partie à un corps) : {g['gate2_fold_Mwr']['max_dev_eV_production_k']:.2e} eV avec M_W de production (k du XML dense, bruit {g['gate2_fold_Mwr']['k_noise']['k27_max_dev']:.1e} sur 27k), "
                 f"{g['gate2_fold_Mwr']['max_dev_eV_snapped_k']:.2e} eV avec M_W recalculé à k = m/27 exacts (seuil 1e-9) → {'OK' if g['gate2_fold_Mwr']['ok'] else 'ÉCHEC'} ; V_loc production contre k exacts : max\\|diff\\| {g['gate2_fold_Mwr']['Vloc_prod_vs_snap']['max_abs_diff_eV']:.2e} eV (p_z–p_z {g['gate2_fold_Mwr']['Vloc_prod_vs_snap']['onsite_pz_prod']:.9f} contre {g['gate2_fold_Mwr']['Vloc_prod_vs_snap']['onsite_pz_snap']:.9f}). "
                 f"Résidu H(R) contre V†εV aux 81 k : {g['one_body_residual']['max_abs_Hwk_minus_VepsV_eV']:.2e} eV (valeurs propres {g['one_body_residual']['max_eig_dev_eV']:.2e}) ; écart littéral (b) contre (c-all) avec H(R) : {g['literal_b_vs_call']['max_dev_eV']:.2e} eV.\n")
        E = z["E_D"]; L.append(f"E_D par base : maille grossière 16 bandes à K {E['uc16_at_K']:.5f} (bandes 3–4 : {E['uc16_bands34_K'][0]:.5f}, {E['uc16_bands34_K'][1]:.5f}) ; maille dense 20 bandes {E['uc20dense_at_K']:.5f} ; H(R) {E['wannier_at_K']:.5f} ; super-cellule parfaite {E['SC_P']:.5f} eV. Parité des états de maille : 16 bandes max\\|1 − \\|⟨σ⟩\\|\\| = {z['bases']['parity16_dev']:.1e} ({z['bases']['n_even16']} pairs, {z['bases']['n_odd16']} impairs), 20 bandes {z['bases']['parity20_dev']:.1e} ({z['bases']['n_even20']} / {z['bases']['n_odd20']}).\n")
        if "coarse_vs_dense16" in z:
            c = z["coarse_vs_dense16"]; L.append(f"M grossier (16 bandes, juin) contre M dense restreint aux 81 k et tronqué à 16 bandes : max\\|\\|M_c\\| − \\|M_d\\|\\| = {c['max_abs_M_diff_eV']:.3f} eV (max\\|M\\| {c['max_abs_M_eV']:.3f}) ; valeurs singulières : écart relatif max {c['max_rel_SV_mismatch']:.2e} ; diagonale, bandes 1–8 : {c['diag_diff_max_bands_1to8_eV']:.2e} eV, toutes bandes {c['diag_diff_max_eV']:.3f} eV (bande {c['diag_diff_max_band_k'][0]+1}, k {c['diag_diff_max_band_k'][1]}) ; ε 16 bandes contre 20 bandes : {c['eps16_vs_eps20_max_eV']:.2e} eV ; gap minimal bandes 16–17 (dense) {c['eps20_gap_16_17_min_eV']:.1e} eV ; états de la maille 16 bandes à parité non entière (> 1e-4) : {c['n_C16_parity_dev_gt_1em4']}.\n")
        L.append("| variante | dim | couplage résiduel pair–impair (eV) | E_D (eV) | pairs (σ) dans la fenêtre | impairs (π) | états localisés (ε − E_D ; w₂) pairs | impairs | décalage rigide résiduel (meV, états < E_D − 4) / résidu max (eV) |"); L.append("|---|---|---|---|---|---|---|---|---|")
        for t in ("QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"):
            if t not in z["variants"]: continue
            v = z["variants"][t]; f = z["rigid_shift_fit"]["fits"].get(t)
            cpl = "—" if v.get("coupling_even_odd") is None else "%.1e" % v["coupling_even_odd"]
            rig = "—" if f is None else "%+.1f / %.2f" % (f["shift_eV"] * 1e3, f["max_residual_eV"])
            loc_e = "; ".join("%+.3f ; %.3f" % (x, w) for x, w in v["blocks"]["even"]["localized"]) or "aucun"
            loc_o = "; ".join("%+.3f ; %.3f" % (x, w) for x, w in v["blocks"]["odd"]["localized"]) or "aucun"
            L.append(f"| {t} | {v.get('dim', '—')} | {cpl} | {_f(v.get('E_D'), 5)} | {v['blocks']['even']['n_window']} | {v['blocks']['odd']['n_window']} | {loc_e} | {loc_o} | {rig} |")
        L.append("\nMarches (états localisés de l'étape i → état localisé le plus proche de l'étape i + 1, même parité) :\n")
        L.append("| marche | parité | états dans la fenêtre (i → i+1) | médiane \\|Δε\\| des spectres triés (eV) | correspondances (ε_i ; w₂,i → ε_{i+1} ; w₂,i+1 ; Δε) |"); L.append("|---|---|---|---|---|")
        for k, v in z["steps"].items():
            for lab in ("even", "odd"):
                s_ = v[lab]; L.append(f"| {k} | {lab} | {s_.get('n_from', '—')} → {s_.get('n_to', '—')} | {_f(s_.get('median_abs_shift_sorted_eV'), 4)} | " + ("; ".join(f"{m['x_from']:+.3f} ; {m['w_from']:.3f} → " + ("absent" if m['x_to'] is None else f"{m['x_to']:+.3f} ; {m['w_to']:.3f} ; {m['dx']:+.3f}{' (> 0,1)' if m['moved_gt_0p1'] else ''}") for m in s_["localized_match"]) or "aucun état localisé") + " |")
        L.append("\nTous les états de la fenêtre par variante et parité (ε − E_D ; w₂) :\n")
        for t in ("QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"):
            if t not in z["variants"]: continue
            for lab in ("even", "odd"):
                b = z["variants"][t]["blocks"][lab]; L.append(f"- {t}, {lab} : " + ", ".join(f"{x:+.3f} ({w:.3f})" for x, w in sorted(zip(b['x'], b['w2']))) + (f" ; poids Wannier site+voisins : " + ", ".join(f"{w:.3f}" for w in b["w_site"]) if "w_site" in b else ""))
    # ---- D3
    p3 = os.path.join(WORK, "d3", "d3_results.json")
    if os.path.exists(p3):
        z = json.load(open(p3)); g = z["gate"]; L.append("\n### D3 — porte de régression et balayage α\n")
        L.append(f"Porte (α = 1, matrice complète, 300², fenêtre ±3 eV) : min \\|det\\|/max = {g['min_det_rel']:.4e} à {g['at']:+.3f} eV (attendu 2,091e-4 à −2,530) ; min \\|λ\\| = {g['min_abs_lam']:.4f}, λ = {_c(g['lam'])} à {g['lam_at']:+.3f} eV (attendu +0,0003 + 0,0108 i) → {'OK' if g['ok'] else 'ÉCHEC'}. Minima locaux de \\|det\\|/max : {', '.join(f'{x:+.3f} ({y:.3e}, min|λ| {m:.4f})' for x, y, m in g['minima_det'])}.\n")
        mm = g.get("minimum_m2p530", {}); L.append(f"Vecteur propre de λ_min à −2,530 eV : poids bloc π {mm.get('weight_pi_block', float('nan')):.4f}, bloc σ {mm.get('weight_sigma_block', float('nan')):.4f} ; groupes : " + ", ".join(f"{k} {v:.4f}" for k, v in mm.get("weights", {}).items()) + f" ; λ = {_c(mm.get('lam'))}.\n")
        for nk in ("nk300", "nk600"):
            if nk not in z["scan"]: continue
            L.append(f"\n#### Balayage α, {nk[2:]}² (fenêtre [−3, +1] eV ; \\|det\\| normalisé sur la fenêtre)\n")
            L.append("| variante | α | bloc | min \\|det\\|/max (position) | minima locaux de \\|det\\|/max | min \\|λ\\| (position) ; λ | poids de v_min : p_z lacune / p_z 3 voisins / sp² A(0) / autres p_z / autres sp² | bloc π / σ | racines de Re λ_min (recouvrement ; saut) | pic de −Im T̄(K) (max) | T̄(E_D) |"); L.append("|---|---|---|---|---|---|---|---|---|---|---|")
            for key, r in z["scan"][nk].items():
                var, al_, bl = key.split("_"); w = r["weights"]
                roots_txt = ", ".join("%+.3f (%.2f%s)" % (e, o, "; saut" if j else "") for e, o, j in r["roots_Re_lam_min"]) or "aucune"
                mind = ", ".join("%+.3f (%.2e)" % (x, y) for x, y in r["minima_det"])
                tb = ("%+.3f (%.3f eV)" % (r["peak_minus_Im_Tbar"], r["max_minus_Im_Tbar"])) if "peak_minus_Im_Tbar" in r else "—"
                L.append(f"| {var} | {al_[1:]} | {bl} | {r['min_det_rel']:.3e} ({r['det_at']:+.3f}) | {mind} | {r['min_abs_lam']:.4f} ({r['lam_at']:+.3f}) ; {_c(r['lam'])} | {w['pz_vac']:.3f} / {w['pz_nn']:.3f} / {w['sp2_A0']:.3f} / {w['pz_other']:.3f} / {w['sp2_other']:.3f} | {r['weight_pi']:.3f} / {r['weight_sigma']:.3f} | {roots_txt} | {tb} | " + (_c(r['Tbar_at_ED']) if 'Tbar_at_ED' in r else "—") + " |")
    with open(os.path.join(d, "tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log(f"[tables] sections/tables.md écrit ({len(L)} lignes)")

# ----------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("cmd", choices=["prep", "d6", "d1", "d5", "r", "d4", "d3", "g0", "tables", "figs"])
    a = ap.parse_args()
    log(f"=== r4_driver {a.cmd} (PROJ {PROJ}, WORK {WORK}, OMP {os.environ.get('OMP_NUM_THREADS')}) ===")
    {"prep": cmd_prep, "d6": cmd_d6, "d1": cmd_d1, "d5": cmd_d5, "r": cmd_r, "d4": cmd_d4, "d3": cmd_d3, "g0": cmd_g0, "tables": cmd_tables, "figs": cmd_figs}[a.cmd](a)
    log(f"=== {a.cmd} terminé ===")


if __name__ == "__main__":
    main()
