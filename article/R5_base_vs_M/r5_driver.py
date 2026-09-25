#!/usr/bin/env python
"""
r5_driver.py -- pilote unique de la campagne R5 (base ou M ? ; vérité DFT selon la taille ; convergence en bandes).

Sous-commandes (une par job) : a (J1 : A.1-A.5), c (J2 : C), b1 (J5 : vérifications du .save 128 bandes), b2ml (J4, srun MPI :
M^L 128 bandes), b2nl (J4, processus frais : M^NL), b2sum (J4 : somme + vérifications), b3 (J5 : B.3-B.4).
Données de production en lecture seule ; sorties dans ce répertoire : a/, c/, b/, fig/ ; journal r5_log.txt.
Modules de campagne (répertoire de travail, hors dépôt) : r5_sc_projection, r5_deltav_pw, r5_basis_diagnostics, r5_alignment_ext.
Modules R4 réutilisés (src/, non commités) : qe_gamma_io, supercell_fold, alignment ; aides du pilote R4 (r4_driver : setup,
read_wfc_subset, uc_parity, spectrum_blocks, caches M_W / V_loc, load_pot_eV, jsonable). Routines de production réutilisées
telles quelles : qe_io, matrix_io, config, Mbk_to_Mwk, Hwr_to_Hwk, compute_ML_R_mpi_shared, compute_M_NL. Aucun fichier de
production modifié ; le pilote R4 n'est pas modifié (son journal n'est pas écrit : r4.log est redirigé ici).
"""
import os
import sys
import json
import time
import glob
import argparse

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__))
GQ = os.path.dirname(os.path.dirname(WORK))                                    # .../graphene/qe
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
R4DIR = os.path.join(GQ, "defects", "R4_quasi_lie")
sys.path.insert(0, os.path.join(PROJ, "src")); sys.path.insert(0, R4DIR); sys.path.insert(0, WORK)
os.chdir(PROJ)

from electron_defect_interaction.config import load_production, dense_paths, HA2EV
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.io import qe_gamma_io as qg
from electron_defect_interaction.defects import alignment as al
from electron_defect_interaction.wannier import supercell_fold as sf
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.utils.lattice import red_to_cart
import r4_driver as r4
import r5_sc_projection as sp
import r5_basis_diagnostics as bd
import r5_alignment_ext as ax

BOHR = r4.BOHR
SCRATCH = r4.SCRATCH
RES = r4.RES
SC_D, SC_P, UC9, UPF = r4.SC_D, r4.SC_P, r4.UC9, r4.UPF
VKS = dict(d=r4.VKS["d"], p=r4.VKS["p"])
N_SC, N_CELLS = r4.N_SC, r4.N_CELLS
K_RED, KP_RED = r4.K_RED, r4.KP_RED
R_W2, R_W1 = r4.R_W2, r4.R_W1
WF_SIGMA, WF_PI, WF_PZ_A, WF_PZ_B, NW, NN_CELLS = r4.WF_SIGMA, r4.WF_PI, r4.WF_PZ_A, r4.WF_PZ_B, r4.NW, r4.NN_CELLS
ND = (N_SC, N_SC, 1)
NB128_OUT = f"{SCRATCH}/R5_uc9x9_nb128"
NB128 = f"{NB128_OUT}/defect_unit_cell_9x9.save"
NLIST = [16, 24, 32, 48, 64, 96, 128]
WIN = (-3.0, 1.0)
NT = int(os.environ.get("OMP_NUM_THREADS", "8"))
TOL_PURE = 1e-6                                   # eV, A.2 : états purs contre M/N_cells
UC_SAVE = {5: f"{SCRATCH}/defect_unit_cell_5x5/defect_unit_cell_5x5.save", 6: f"{SCRATCH}/defect_unit_cell_6x6/defect_unit_cell_6x6.save",
           7: f"{SCRATCH}/defect_unit_cell_7x7/defect_7x7.save", 8: f"{SCRATCH}/defect_unit_cell_8x8/defect_unit_cell_8x8.save",
           9: UC9, 10: f"{SCRATCH}/defect_unit_cell_10x10/defect_unit_cell_10x10.save",
           11: f"{SCRATCH}/defect_unit_cell_11x11/defect_unit_cell_11x11.save", 12: f"{SCRATCH}/defect_unit_cell_12x12/defect_unit_cell_12x12.save"}
UC_NOTE = {11: "run bands (chemin Γ-K-M-Γ, 181 k) : la grille 121 k a été écrasée ; Δ_fold par l'état le plus bas"}
FAM_3M = [6, 9, 12]; FAM_N3M = [5, 7, 8, 10, 11]


# ----------------------------------------------------------------------------------------------- utilitaires
def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(os.path.join(WORK, "r5_log.txt"), "a") as f:
        f.write(line + "\n")


r4.log = log                                      # les aides de R4 écrivent dans le journal R5, jamais dans r4_log.txt


def ensure(sub):
    d = os.path.join(WORK, sub); os.makedirs(d, exist_ok=True); return d


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(r4.jsonable(d), f, indent=1, ensure_ascii=False)
    log(f"saved {os.path.relpath(path, WORK)}")


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


def kindex(k_red, target):
    return int(np.argmin(np.linalg.norm(np.mod(np.asarray(k_red) - np.asarray(target) + 0.5, 1) - 0.5, axis=1)))


def load_pot_eV(path):
    return r4.load_pot_eV(path)


def _f(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:+.{nd}f}"


def _re(z):
    """Partie réelle d'un complexe (en mémoire) ou d'une paire [re, im] (relu d'un json)."""
    return float(z[0]) if isinstance(z, (list, tuple)) else float(np.real(z))


# ----------------------------------------------------------------------------------------------- bases repliées (9x9)
def load_bases(with_dense=True, nb128=False):
    """Base 16 bandes (maille grossière), 20 bandes (dense restreint aux 81 k) et, si nb128, 128 bandes (.save neuf)."""
    ng = qe_io.get_ngfft(SC_D); A_b, Om = qe_io.get_A_volume(SC_D)
    B = dict(ng=tuple(int(x) for x in ng), A_b=A_b, A_A=A_b * BOHR, Om=float(Om), bases={})

    def one(save, k_sub=None, tag=""):
        k = qe_io.get_k_red(save); eps = qe_io.get_eigenvalues(save) * HA2EV
        if k_sub is None:
            C, nG = qe_io.get_C_nk(save); G = qe_io.get_G_red(save); ks = k
        else:
            k27 = k; ks = k27[k_sub]; eps = eps[:, k_sub]
            C, nG, G = r4.read_wfc_subset(save, k_sub)
        par = r4.uc_parity(C, nG, G)
        flat = sf.sc_planewave_index(ks, G, nG, ND, ng)
        nu, ntot, ok = sp.check_planewave_union(flat, 2 * 385710 - 1)
        iK = kindex(ks, K_RED); iKp = kindex(ks, KP_RED); iG = kindex(ks, (0, 0, 0))
        E_D = float(0.5 * (eps[3, iK] + eps[4, iK]))
        log(f"[bases] {tag}: {C.shape[0]} bandes x {C.shape[1]} k, ondes planes distinctes {nu}/{ntot} (attendu 771419) -> {'OK' if ok else 'ÉCHEC'}, E_D(K) {E_D:.6f} eV")
        return dict(save=save, k=ks, eps=eps, C=C, nG=nG, G=G, par=par, flat=flat, iK=iK, iKp=iKp, iG=iG, E_D=E_D, pw_union=(nu, ntot, ok), nb=int(C.shape[0]))

    B["bases"][16] = one(UC9, tag="16 (grossier)")
    B["k81"] = B["bases"][16]["k"]
    if with_dense:
        cfg = load_production(verbose=False); dp = dense_paths(cfg, "9x9"); B["dp"] = dp
        k27 = qe_io.get_k_red(dp["uc"])
        idx81 = np.array([int(np.where(np.all(np.abs(np.mod(k27 - k + 0.5, 1) - 0.5) < 1e-6, axis=1))[0][0]) for k in B["k81"]])
        B["idx81"] = idx81; B["k27"] = k27
        B["bases"][20] = one(dp["uc"], k_sub=idx81, tag="20 (dense ⊂ 81 k)")
        # bug 0.3 : index construit avec k81 au lieu de k27[idx81]
        b20 = B["bases"][20]
        flat_wrong = sf.sc_planewave_index(B["k81"], b20["G"], b20["nG"], ND, ng)
        nbad = int(sum(not np.array_equal(flat_wrong[i], b20["flat"][i]) for i in range(len(flat_wrong))))
        try:
            ax.check_k_reference(B["k81"], b20["k"]); msg = "aucune exception (inattendu)"
        except ValueError as e:
            msg = str(e)
        B["bug03"] = dict(n_k_index_differs=nbad, offsets=np.unique(np.rint(b20["k"] - B["k81"]).astype(int), axis=0).tolist(), check_k_reference=msg)
        log(f"[bases] bug 0.3 : index d'ondes planes différent pour {nbad}/81 k avec k81 ; check_k_reference : {msg[:90]}")
    if nb128:
        B["bases"][128] = one(NB128, tag="128 (nscf neuf)")
        ax.check_k_reference(B["k81"], B["bases"][128]["k"])
    return B


def dense_restricted(dp, idx81, part):
    """M dense (20, 729, 20, 729) restreint aux 81 k, en eV ; sidecar vérifié (unit_cell, hartree) sans charger les 3,4 G."""
    path = {"tot": dp["mfile"], "L": os.path.join(RES, "M_L_dense_9x9.npy"), "NL": os.path.join(RES, "M_NL_dense_9x9.npy")}[part]
    meta = matrix_io.read_manifest(path)
    if meta is None or meta.get("bloch_norm") != matrix_io.UNIT_CELL or meta.get("units") != matrix_io.HARTREE:
        raise ValueError(f"{path}: sidecar absent ou convention inattendue ({meta})")
    Mm = np.load(path, mmap_mode="r"); nb = Mm.shape[0]
    M = np.array(Mm[np.ix_(np.arange(nb), idx81, np.arange(nb), idx81)]) * HA2EV
    del Mm
    return M


def coarse_M(part):
    return np.load(os.path.join(R4DIR, "prep", f"M_coarse_{part}.npy")) * HA2EV


def qe_states_9x9():
    """États QE de la 9x9 lacune dans la fenêtre (R4 D1) : bandes 297-326, coefficients (demi-sphère), parité, alignement."""
    z1 = np.load(os.path.join(R4DIR, "d1", "d1_states.npz")); d1 = json.load(open(os.path.join(R4DIR, "d1", "d1_results.json")))
    shift = float(d1["alignment"]["D"]["shifts_eV"]["1.0"]); E_D_SC = float(z1["E_D"]); thr = float(z1["thr"]); s_vac = z1["s_vac"]
    b1 = z1["D_band1"]; b0 = int(b1[0]) - 1; bend = int(b1[-1])
    Cw, mill, go, at = qg.read_wfc_gamma(SC_D, bands=(b0, bend))
    e_all = qg.get_eigenvalues_spin(SC_D)[0][0] * HA2EV
    par = qg.mirror_parity_z(Cw, mill, z0_red=0.0, gamma_only=go).real
    nrm = np.abs(Cw[:, 0]) ** 2 + 2 * (np.abs(Cw[:, 1:]) ** 2).sum(1) if (mill[0] == 0).all() else None
    return dict(Cw=Cw, mill=mill, gamma_only=bool(go), bands1=b1, b0=b0, bend=bend, e=e_all[b0:bend], e_all=e_all, shift=shift, E_D_SC=E_D_SC, thr=thr,
                s_vac=s_vac, x=z1["D_x"], parity=par, parity_d1=z1["D_parity"], w2_d1=z1["D_w2"], norm=nrm)


def highlighted(Q):
    """Indices (dans la fenêtre de 30) des états mis en évidence : σ +0,101 x2 (324, 325), π -0,737 (320), -1,759 x2 (318, 319), +0,269 (326)."""
    b1 = list(Q["bands1"])
    return {lab: b1.index(b) for lab, b in (("sigma_324", 324), ("sigma_325", 325), ("pi_320", 320), ("pi_318", 318), ("pi_319", 319), ("pi_326", 326)) if b in b1}


def window_weights(blocks, C, nG, flat, nb, ng, mask2, lo=WIN[0], hi=WIN[1]):
    """w2 (disque 2 Å) des états propres de la fenêtre, par bloc de parité ; vecteurs (n*nk + k) reconstruits sur la base (n, k)."""
    W = {}
    nk = C.shape[1]
    for lab in blocks:
        x = blocks[lab]["x"]; sel = np.where((x >= lo) & (x <= hi))[0]; w = np.zeros(len(sel))
        for jj, j in enumerate(sel):
            dfull = np.zeros(nb * nk, complex); dfull[blocks[lab]["idx"]] = blocks[lab]["v"][:, j]
            w[jj] = sf.folded_density_2d(dfull.reshape(nb, nk), C, nG, flat, ng, workers=NT)[mask2].sum()
        W[lab] = dict(sel=sel, w2=w)
    return W


def record_variant(out, tag, blocks, W, E_Dm, thr, extra=None):
    rec = dict(E_D=E_Dm, blocks={}, **(extra or {}))
    for lab in ("even", "odd"):
        x = blocks[lab]["x"]; sel = W[lab]["sel"]; w = W[lab]["w2"]
        locs = [(float(x[j]), float(w[i])) for i, j in enumerate(sel) if w[i] > thr]
        rec["blocks"][lab] = dict(n_window=int(len(sel)), x=[float(v) for v in x[sel]], w2=[float(v) for v in w], localized=locs)
    out[tag] = rec
    log(f"[variante] {tag}: pairs {rec['blocks']['even']['n_window']} (localisés {[(round(a,3), round(b,3)) for a, b in rec['blocks']['even']['localized']]}), "
        f"impairs {rec['blocks']['odd']['n_window']} (localisés {[(round(a,3), round(b,3)) for a, b in rec['blocks']['odd']['localized']]})")
    return rec


def deltav_grid_and_projector(B, mill_half):
    """ΔV^L = V_d - V_p (eV) sur la grille [ix, iy, iz] du .save, et le projecteur KB de l'atome retiré sur la sphère complète."""
    import r5_deltav_pw as dv
    Vd, ngv = qe_io.get_pot(VKS["d"], subtract_mean=False, to_hartree=True); Vp, _ = qe_io.get_pot(VKS["p"], subtract_mean=False, to_hartree=True)
    if tuple(int(x) for x in ngv) != tuple(B["ng"]):
        raise ValueError(f"grille pp.x {ngv} != grille du .save {B['ng']}")
    dV = (Vd - Vp).transpose(2, 1, 0) * HA2EV; del Vd, Vp
    x_p = qe_io.get_x_red(SC_P); x_d = qe_io.get_x_red(SC_D)
    s_vac, i_vac_p, _ = al.vacancy_site(x_p, x_d, B["A_A"])
    tau_vac = red_to_cart(x_p[i_vac_p], B["A_b"])
    B_sc, _ = qe_io.get_B_volume(SC_P); ecut = float(qe_io.get_ecut(SC_P))
    _, mill_full = sp.sc_full_sphere(np.zeros((1, len(mill_half)), complex), mill_half)
    flat_full = sp.grid_flat_index(mill_full, B["ng"])
    t0 = time.time(); proj = dv.RemovedAtomProjector(mill_full, B_sc, tau_vac, UPF, B["Om"], ecut, energy_scale=HA2EV)
    log(f"[ΔV] grille {B['ng']}, <ΔV>_3D {dV.mean()*1e3:+.2f} meV, ΔV au site {dV[tuple(np.rint(s_vac * np.array(B['ng'])).astype(int) % np.array(B['ng']))]:+.3f} eV ; "
        f"projecteur KB de l'atome retiré (indice parfaite {i_vac_p+1}, {len(mill_full)} g, |g|max {proj.K_norm_max:.2f} < qmax {proj.qmax:.1f}) en {time.time()-t0:.0f} s")
    return dV, proj, flat_full, mill_full, dict(s_vac=s_vac.tolist(), i_vac_p_1based=i_vac_p + 1, tau_vac_bohr=np.asarray(tau_vac).tolist(), mean3d_meV=float(dV.mean() * 1e3))


# ----------------------------------------------------------------------------------------------- J1 : A
def cmd_a(a):
    import r5_deltav_pw as dv
    d = ensure("a"); out = dict(definitions={}); t_start = time.time()
    B = load_bases(with_dense=True); ng = B["ng"]; Om = B["Om"]; A_A = B["A_A"]
    Q = qe_states_9x9(); ns = len(Q["e"]); hl = highlighted(Q)
    mask2 = qg.inplane_disc_mask(ng, A_A, Q["s_vac"], R_W2)
    out["qe_states"] = dict(bands=Q["bands1"].tolist(), n=ns, shift_Lu_eV=Q["shift"], E_D_SC=Q["E_D_SC"], thr_w2=Q["thr"], gamma_only=Q["gamma_only"],
                            max_parity_dev_vs_d1=float(np.abs(Q["parity"] - Q["parity_d1"]).max()), norm_min=float(Q["norm"].min()), norm_max=float(Q["norm"].max()),
                            highlighted={k: int(v) for k, v in hl.items()})
    out["bug03"] = B["bug03"]; out["pw_union"] = {str(nb): B["bases"][nb]["pw_union"] for nb in B["bases"]}
    log(f"[A] {ns} états QE (bandes {Q['bands1'][0]}-{Q['bands1'][-1]}), parité vs D1 max dev {out['qe_states']['max_parity_dev_vs_d1']:.1e}, norme [{Q['norm'].min():.12f}, {Q['norm'].max():.12f}]")
    # ---- M et H par base
    Ms = {16: {p: coarse_M(p) for p in ("tot", "L", "NL")}, 20: {p: dense_restricted(B["dp"], B["idx81"], p) for p in ("tot", "L", "NL")}}
    out["M_linearity"] = {str(nb): float(np.abs(Ms[nb]["tot"] - Ms[nb]["L"] - Ms[nb]["NL"]).max()) for nb in Ms}
    # convention diagnostique (hors prompt) : M^L x N_cells + M^NL, les deux parties alors en norme unit_cell (voir A.2 états purs)
    for nb in (16, 20):
        Ms[nb]["tot_Lx81"] = N_CELLS * Ms[nb]["L"] + Ms[nb]["NL"]
    out["definitions"]["Lx81"] = "variante diagnostique : M^L multiplié par N_cells = 81 avant la somme avec M^NL (A.2 : M^L de production = élément entre états normés sur la super-cellule)"
    CONV = ("prod", "Lx81"); MKEY = {"prod": "tot", "Lx81": "tot_Lx81"}
    H = {c: {} for c in CONV}; blocks = {c: {} for c in CONV}
    for conv in CONV:
        for nb in (16, 20):
            b = B["bases"][nb]; H[conv][nb] = sf.bloch_folded_hamiltonian(b["eps"], Ms[nb][MKEY[conv]], N_CELLS)
            blocks[conv][nb], cpl, herm = r4.spectrum_blocks(H[conv][nb], b["par"].reshape(-1), b["E_D"])
            log(f"[A] H({nb}, {conv}) : dim {H[conv][nb].shape[0]}, couplage pair-impair {cpl:.1e} eV, hermiticité {herm:.1e}, E_D {b['E_D']:.6f}")
    # ---- recouvrements c_nk (0.2)
    c = {16: np.zeros((ns, 16, 81), complex), 20: np.zeros((ns, 20, 81), complex)}
    t0 = time.time()
    for s in range(ns):
        A = sp.sc_state_grid(Q["Cw"][s], Q["mill"], ng, gamma_only=Q["gamma_only"])
        for nb in (16, 20):
            b = B["bases"][nb]; c[nb][s] = sp.bloch_overlaps(A, b["C"], b["nG"], b["flat"])
    log(f"[A.1] recouvrements des {ns} états sur les bases 16 et 20 en {time.time()-t0:.0f} s")
    np.savez(os.path.join(d, "a_overlaps.npz"), c16=c[16], c20=c[20], x_qe=Q["x"], parity=Q["parity"], e_qe=Q["e"], bands=Q["bands1"])
    # ---- A.1 (convention de production, puis diagnostique Lx81)
    for conv in CONV:
        A1 = {}
        for nb in (16, 20):
            b = B["bases"][nb]; Hh = 0.5 * (H[conv][nb] + H[conv][nb].conj().T); rows = []
            for s in range(ns):
                lab = "even" if Q["parity"][s] > 0 else "odd"; other = "odd" if lab == "even" else "even"
                idx = blocks[conv][nb][lab]["idx"]; cf = c[nb][s].reshape(-1); cb = cf[idx]
                e_ref = (Q["e"][s] - Q["shift"]) - Q["E_D_SC"] + b["E_D"]
                r = bd.basis_diagnostics(Hh[np.ix_(idx, idx)], cb, e_ref, evals=blocks[conv][nb][lab]["x"] + b["E_D"], evecs=blocks[conv][nb][lab]["v"])
                r["delta"] = float(1 - np.vdot(cf, cf).real); r["leak_other_parity"] = float(np.sum(np.abs(cf[blocks[conv][nb][other]["idx"]]) ** 2))
                r["rayleigh_minus_eps"] = r["rayleigh"] - e_ref; r["top_x"] = [(e - b["E_D"], w) for e, w in r["top"]]
                r.update(band=int(Q["bands1"][s]), x_qe=float(Q["x"][s]), parity=lab, w2_qe=float(Q["w2_d1"][s]))
                rows.append(r)
            A1[str(nb)] = rows
            for lab in hl:
                r = rows[hl[lab]]; log(f"[A.1 {conv}] base {nb} {lab}: δ {r['delta']:.4f}, R-ε {r['rayleigh_minus_eps']*1e3:+.1f} meV, résidu {r['residual']:.4f} eV, top {[(round(e,3), round(w,3)) for e, w in r['top_x'][:3]]}")
        out["A1" if conv == "prod" else "A1_Lx81"] = A1
    # ---- A.2 : vérification sur états purs, puis les 30 états
    dV, proj, flat_full, mill_full, geo = deltav_grid_and_projector(B, Q["mill"]); out["deltaV"] = geo
    A2 = dict(pure={}, states={}, verified={}, pure_Lx81={}, verified_Lx81={}, states_Lx81={})
    for nb in (16, 20):
        b = B["bases"][nb]; iK, iKp = b["iK"], b["iKp"]
        pairs = [(3, iK, 3, iK), (3, iK, 4, iK), (3, iK, 3, iKp), (0, 0, 0, 0), (2, 5, 7, 40), (10, 17, 1, 60)]
        if nb == 20:
            pairs.append((17, iK, 19, iKp))
        t0 = time.time()
        res = dv.check_pure_bloch(pairs, b["C"], b["nG"], b["flat"], ng, Om, dV, Ms[nb]["L"], Ms[nb]["NL"], N_CELLS, proj, flat_full, workers=NT)
        dmax = max(max(r["dL"], r["dNL"]) for r in res); okv = bool(dmax <= TOL_PURE)
        ratios = [abs(r["VL_direct"]) / abs(r["ML_over_N"]) for r in res if abs(r["ML_over_N"]) > 1e-4]
        A2["pure"][str(nb)] = res; A2["verified"][str(nb)] = dict(max_abs_diff_eV=float(dmax), tol=TOL_PURE, ok=okv, n_pairs=len(pairs), seconds=time.time() - t0,
                                                                     max_dL=float(max(r["dL"] for r in res)), max_dNL=float(max(r["dNL"] for r in res)),
                                                                     ratio_L_direct_over_ML_over_N=[float(x) for x in ratios])
        log(f"[A.2] base {nb} états purs : max|direct - M/81| = {dmax:.2e} eV sur {len(pairs)} paires (L {A2['verified'][str(nb)]['max_dL']:.1e}, NL {A2['verified'][str(nb)]['max_dNL']:.1e}) -> {'OK' if okv else 'ÉCHEC : STOP A.2 (production)'} ; rapport L direct / (M_L/81) : {[round(x, 4) for x in ratios]}")
        for r in res:
            log(f"       {r['pair']}: L direct {r['VL_direct'].real:+.6f} vs M_L/81 {r['ML_over_N'].real:+.6f} (Δ {r['dL']:.1e}) ; NL direct {r['VNL_direct'].real:+.6f} vs M_NL/81 {r['MNL_over_N'].real:+.6f} (Δ {r['dNL']:.1e})")
        # diagnostique : même vérification avec M^L x N_cells (M^L de production = élément entre états normés sur la super-cellule)
        res_c = dv.check_pure_bloch(pairs, b["C"], b["nG"], b["flat"], ng, Om, dV, N_CELLS * Ms[nb]["L"], Ms[nb]["NL"], N_CELLS, proj, flat_full, workers=NT)
        dmax_c = max(max(r["dL"], r["dNL"]) for r in res_c); okc = bool(dmax_c <= TOL_PURE)
        A2["pure_Lx81"][str(nb)] = res_c; A2["verified_Lx81"][str(nb)] = dict(max_abs_diff_eV=float(dmax_c), tol=TOL_PURE, ok=okc, max_dL=float(max(r["dL"] for r in res_c)), max_dNL=float(max(r["dNL"] for r in res_c)))
        log(f"[A.2] base {nb} états purs, M^L x 81 (diagnostique) : max|direct - M/81| = {dmax_c:.2e} eV (L {A2['verified_Lx81'][str(nb)]['max_dL']:.1e}, NL {A2['verified_Lx81'][str(nb)]['max_dNL']:.1e}) -> {'OK' if okc else 'ÉCHEC'}")
        for conv, okk in (("prod", okv), ("Lx81", okc)):
            if not okk:
                continue
            rows = []; t0 = time.time()
            Mf = {"tot": Ms[nb][MKEY[conv]].reshape(nb * 81, nb * 81), "L": (Ms[nb]["L"] if conv == "prod" else N_CELLS * Ms[nb]["L"]).reshape(nb * 81, nb * 81), "NL": Ms[nb]["NL"].reshape(nb * 81, nb * 81)}
            A2_states_loop(conv, nb, b, Mf, rows, c, Q, ns, ng, Om, dV, proj, flat_full, hl, A2, geo)
    out["A2"] = A2
    del dV, proj

    # ---- A.3 : parité de tous les états QE, entrelacement
    t0 = time.time(); par_all = np.zeros(len(Q["e_all"]))
    for b0 in range(0, len(Q["e_all"]), 32):
        Cs, mill, go, _ = qg.read_wfc_gamma(SC_D, bands=(b0, min(b0 + 32, len(Q["e_all"]))))
        par_all[b0:b0 + Cs.shape[0]] = qg.mirror_parity_z(Cs, mill, 0.0, gamma_only=go).real; del Cs
    x_all = Q["e_all"] - Q["shift"] - Q["E_D_SC"]
    A3 = dict(n_qe=int(len(x_all)), parity_dev_max=float(np.abs(1 - np.abs(par_all)).max()), n_qe_even=int((par_all > 0).sum()), n_qe_odd=int((par_all < 0).sum()), models={})
    for tag, nb in (("a1_16", 16), ("a3_20", 20)):
        A3["models"][tag] = {}
        for lab, sgn in (("even", 1), ("odd", -1)):
            xr = x_all[par_all * sgn > 0]; xm = blocks["prod"][nb][lab]["x"]
            A3["models"][tag][lab] = bd.interlacing(xr, xm, WIN[0], WIN[1]); A3["models"][tag][lab]["n_model_below_lo_m4"] = int((xm < -4).sum())
            A3["models"][tag][lab]["ref_window_m4"] = [float(v) for v in np.sort(xr[(xr >= -4) & (xr <= 1)])]; A3["models"][tag][lab]["model_window_m4"] = [float(v) for v in np.sort(xm[(xm >= -4) & (xm <= 1)])]
            r = A3["models"][tag][lab]
            log(f"[A.3] {tag} {lab}: QE sous E_D-3 : {r['n_ref_below_lo']} ; modèle : {r['n_model_below_lo']} ; fenêtre QE {len(r['ref_window'])} / modèle {len(r['model_window'])} ; "
                f"min_k[λ_k(mod)-λ_(k-s)(QE)] s=0,1,2 (fenêtre) : {[round(r['shifts'][s].get('min_window', np.nan), 3) for s in (0, 1, 2)]} ; (tout) : {[round(r['shifts'][s]['min_all'], 3) for s in (0, 1, 2)]}")
    out["A3"] = A3; np.savez(os.path.join(d, "a_qe_all.npz"), x_all=x_all, parity_all=par_all, e_all=Q["e_all"])
    log(f"[A.3] parité de {len(x_all)} états QE en {time.time()-t0:.0f} s (dev max {A3['parity_dev_max']:.1e})")
    # ---- A.4 : w2 corrigés de toutes les variantes de R4
    t0 = time.time(); V = {}
    b16 = B["bases"][16]; b20 = B["bases"][20]
    for part in ("tot", "L", "NL"):
        Hp = sf.bloch_folded_hamiltonian(b16["eps"], Ms[16][part], N_CELLS); bl, cpl, herm = r4.spectrum_blocks(Hp, b16["par"].reshape(-1), b16["E_D"])
        W = window_weights(bl, b16["C"], b16["nG"], b16["flat"], 16, ng, mask2); record_variant(V, f"a1_{part}", bl, W, b16["E_D"], Q["thr"], dict(coupling_even_odd=cpl, dim=int(Hp.shape[0])))
        if part == "tot":
            blocks_a1 = bl
    Hp = sf.bloch_folded_hamiltonian(b20["eps"][:16], Ms[20]["tot"][:16, :, :16, :], N_CELLS); bl, cpl, herm = r4.spectrum_blocks(Hp, b20["par"][:16].reshape(-1), b20["E_D"])
    W = window_weights(bl, b20["C"][:16], b20["nG"], b20["flat"], 16, ng, mask2); record_variant(V, "a2_16", bl, W, b20["E_D"], Q["thr"], dict(coupling_even_odd=cpl, dim=int(Hp.shape[0])))
    W = window_weights(blocks["prod"][20], b20["C"], b20["nG"], b20["flat"], 20, ng, mask2); record_variant(V, "a3_20", blocks["prod"][20], W, b20["E_D"], Q["thr"], dict(dim=int(H["prod"][20].shape[0])))
    # diagnostique : M^L x 81 + M^NL (hors prompt)
    W = window_weights(blocks["Lx81"][16], b16["C"], b16["nG"], b16["flat"], 16, ng, mask2); record_variant(V, "a1_tot_Lx81", blocks["Lx81"][16], W, b16["E_D"], Q["thr"], dict(dim=int(H["Lx81"][16].shape[0]), diagnostic="M^L x N_cells + M^NL"))
    W = window_weights(blocks["Lx81"][20], b20["C"], b20["nG"], b20["flat"], 20, ng, mask2); record_variant(V, "a3_20_Lx81", blocks["Lx81"][20], W, b20["E_D"], Q["thr"], dict(dim=int(H["Lx81"][20].shape[0]), diagnostic="M^L x N_cells + M^NL"))
    # (b), (c-all), (c-3) : machinerie R4 (porte de jauge, H(R), U, caches M_W et V_loc), index d'ondes planes corrigé
    S = r4.setup(verbose=False); Hwr, Rw, nd, U, Ud = S["Hwr"], S["Rw"], S["nd"], S["U"], S["Ud"]; idx81 = B["idx81"]; k81 = B["k81"]; nk = 81; nw = NW
    V81 = np.einsum("kbw,kwv->kbv", Ud[idx81], U[idx81]); Hk = np.einsum("kbw,kb,kbv->kwv", np.conj(V81), b20["eps"].T, V81)
    Mwk = Mbk_to_Mwk(Ms[20]["tot"], U[idx81], Ud[idx81]); Hb = sf.kbasis_matrix(Hk, Mwk, k81, N_CELLS)
    parW = np.array([1.0 if w in WF_SIGMA else -1.0 for k in range(nk) for w in range(nw)])
    eK, vK = np.linalg.eigh(Hk[b20["iK"]]); wpzK = np.abs(vK[WF_PZ_A]) ** 2 + np.abs(vK[WF_PZ_B]) ** 2; pairK = np.argsort(-wpzK)[:2]; ED_b = float(eK[pairK].mean())
    blocks_b, cpl, herm = r4.spectrum_blocks(Hb, parW, ED_b)

    def w_from_k(bl):
        Wb = {}
        for lab in bl:
            x = bl[lab]["x"]; sel = np.where((x >= WIN[0]) & (x <= WIN[1]))[0]; w = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                cvec = np.zeros(nk * nw, complex); cvec[bl[lab]["idx"]] = bl[lab]["v"][:, j]
                dnk = np.einsum("kbw,kw->bk", V81, cvec.reshape(nk, nw))
                w[jj] = sf.folded_density_2d(dnk, b20["C"], b20["nG"], b20["flat"], ng, workers=NT)[mask2].sum()
            Wb[lab] = dict(sel=sel, w2=w)
        return Wb
    record_variant(V, "b_5wf", blocks_b, w_from_k(blocks_b), ED_b, Q["thr"], dict(coupling_even_odd=cpl, dim=int(Hb.shape[0])))
    Z = r4.load_cache_mwr(); Mwr = Z["Mwr_tot"]; Rn = Z["R"]; R_d = Z["R_d"]; Cc = r4.load_cache_vloc(); Rloc = Cc["Rloc"]
    HS = sf.fold_hwr_to_supercell(Hwr, Rw, nd, N_SC); Mall = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC); M3 = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC, R_local=Rloc)
    _, EK, _ = Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); E_DW = float(0.5 * (EK[0, 3] + EK[0, 4]))
    phase_Rd = np.exp(-2j * np.pi * (k81 @ np.asarray(R_d, float)))

    def w_from_r(bl):
        Wr = {}
        for lab in bl:
            x = bl[lab]["x"]; sel = np.where((x >= WIN[0]) & (x <= WIN[1]))[0]; w = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                cr = np.zeros(N_CELLS * nw, complex); cr[bl[lab]["idx"]] = bl[lab]["v"][:, j]
                ck = sf.rvec_to_kvec(cr, k81, N_SC, nw).reshape(nk, nw) * phase_Rd[:, None]
                dnk = np.einsum("kbw,kw->bk", V81, ck)
                w[jj] = sf.folded_density_2d(dnk, b20["C"], b20["nG"], b20["flat"], ng, workers=NT)[mask2].sum()
            Wr[lab] = dict(sel=sel, w2=w)
        return Wr
    for tag, Msc in (("c_all", Mall), ("c_3", M3)):
        Hc = HS + Msc; bl, cpl, herm = r4.spectrum_blocks(Hc, parW, E_DW); record_variant(V, tag, bl, w_from_r(bl), E_DW, Q["thr"], dict(coupling_even_odd=cpl, dim=int(Hc.shape[0])))
    # QE (D1) et R4 pour comparaison
    d4 = json.load(open(os.path.join(R4DIR, "d4", "d4_results.json")))
    V["QE_D1"] = d4["variants"]["QE_D1"]
    out["A4"] = dict(variants=V, r4_w2={t: {lab: dict(x=d4["variants"][t]["blocks"][lab]["x"], w2=d4["variants"][t]["blocks"][lab]["w2"]) for lab in ("even", "odd")} for t in d4["variants"] if t != "QE_D1"})
    # vérification : w2(a2) = w2(a1) à 1e-3 pour les états localisés (appariement au plus proche en énergie)
    ver = []
    for lab in ("even", "odd"):
        xa = np.array(V["a1_tot"]["blocks"][lab]["x"]); wa = np.array(V["a1_tot"]["blocks"][lab]["w2"]); xb = np.array(V["a2_16"]["blocks"][lab]["x"]); wb = np.array(V["a2_16"]["blocks"][lab]["w2"])
        for x, w in zip(xa, wa):
            if w > Q["thr"]:
                j = int(np.argmin(np.abs(xb - x))); ver.append(dict(parity=lab, x_a1=float(x), w2_a1=float(w), x_a2=float(xb[j]), w2_a2=float(wb[j]), dw2=float(abs(wb[j] - w)), dx=float(xb[j] - x)))
    out["A4"]["verification_a2_vs_a1"] = dict(rows=ver, max_dw2=float(max(r["dw2"] for r in ver)) if ver else None, ok=bool(ver and max(r["dw2"] for r in ver) <= 1e-3), tol=1e-3)
    log(f"[A.4] vérification w2(a2) = w2(a1) (états localisés de a1) : max|Δw2| = {out['A4']['verification_a2_vs_a1']['max_dw2']:.4f} -> {'OK' if out['A4']['verification_a2_vs_a1']['ok'] else 'ÉCART'} ({time.time()-t0:.0f} s)")
    # ---- A.5 : moyenne diagonale de M^L soustraite
    ML = Ms[16]["L"]; nbk = 16 * 81; diag_mean = float(np.mean([ML[n, k, n, k].real for n in range(16) for k in range(81)]))
    e_ref = np.linalg.eigvalsh(0.5 * (H["prod"][16] + H["prod"][16].conj().T))
    A5 = dict(diag_mean_ML_eV=diag_mean, r4_value_eV=0.07658, variants={})
    for tag, cst in (("M_L_minus_mean_identity", diag_mean), ("dV_minus_mean_equiv", N_CELLS * diag_mean)):
        Mm = Ms[16]["tot"].copy().reshape(nbk, nbk); Mm[np.diag_indices(nbk)] -= cst; Hm = sf.bloch_folded_hamiltonian(b16["eps"], Mm.reshape(16, 81, 16, 81), N_CELLS)
        e_m = np.linalg.eigvalsh(0.5 * (Hm + Hm.conj().T)); dshift = e_m - e_ref
        blm, _, _ = r4.spectrum_blocks(Hm, b16["par"].reshape(-1), b16["E_D"])
        loc_shift = []
        for lab in ("even", "odd"):
            for (x, w) in V["a1_tot"]["blocks"][lab]["localized"]:
                xm = blm[lab]["x"]; j = int(np.argmin(np.abs(xm - (x - cst / N_CELLS)))); loc_shift.append(dict(parity=lab, x_a1=x, w2=w, x_new=float(xm[j]), dx=float(xm[j] - x)))
        A5["variants"][tag] = dict(constant_subtracted_eV=cst, rigid_shift_median_eV=float(np.median(dshift)), rigid_shift_min_eV=float(dshift.min()), rigid_shift_max_eV=float(dshift.max()),
                                   max_abs_dev_after_shift_eV=float(np.abs(dshift - np.median(dshift)).max()), expected_shift_eV=-cst / N_CELLS, localized=loc_shift,
                                   max_abs_dev_localized_after_shift_eV=float(max(abs(r["dx"] + cst / N_CELLS) for r in loc_shift)) if loc_shift else None)
        log(f"[A.5] {tag}: constante {cst*1e3:.2f} meV soustraite de la diagonale de M (norme unit_cell) -> décalage rigide médian {np.median(dshift)*1e3:+.3f} meV (attendu {-cst/N_CELLS*1e3:+.3f}), max|Δε - médiane| {np.abs(dshift-np.median(dshift)).max():.1e} eV")
    out["A5"] = A5
    out["timing_s"] = time.time() - t_start
    save_json(os.path.join(d, "a_results.json"), out)
    np.savez(os.path.join(d, "a_spectra.npz"), **{f"{tag}_{lab}_x": np.array(V[tag]["blocks"][lab]["x"]) for tag in V for lab in ("even", "odd")},
             **{f"{tag}_{lab}_w2": np.array(V[tag]["blocks"][lab]["w2"]) for tag in V for lab in ("even", "odd")})
    write_tables_a(out, Q); figs_a(out, Q)


def A2_states_loop(conv, nb, b, Mf, rows, c, Q, ns, ng, Om, dV, proj, flat_full, hl, A2, geo):
    """A.2, les 30 états : (i) c†(M/81)c, (ii) application directe de ΔV à Ψ_in, ⟨Ψ_in|ΔV|Ψ_out⟩, identité et constante d'alignement."""
    import r5_deltav_pw as dv
    t0 = time.time()
    for s in range(ns):
        cf = c[nb][s].reshape(-1); cc = float(np.vdot(cf, cf).real)
        A_in = sp.bloch_superposition_grid(c[nb][s], b["C"], b["nG"], b["flat"], ng); psi_in = sp.grid_to_real(A_in, ng, Om, NT)
        A = sp.sc_state_grid(Q["Cw"][s], Q["mill"], ng, gamma_only=Q["gamma_only"]); psi = sp.grid_to_real(A, ng, Om, NT)
        i_ = {p: complex(np.vdot(cf, Mf[p] @ cf)) / N_CELLS for p in ("tot", "L", "NL")}
        VL_in_in = dv.expect_local(psi_in, psi_in, dV, Om); VL_in_psi = dv.expect_local(psi_in, psi, dV, Om)
        B_in = proj.projections(A_in[flat_full][None])[0]; B_psi = proj.projections(A[flat_full][None])[0]
        VNL_in_in = proj.expect(B_in, B_in); VNL_in_psi = proj.expect(B_in, B_psi)
        VL_in_out = VL_in_psi - VL_in_in; VNL_in_out = VNL_in_psi - VNL_in_in
        lhs = float(np.sum(b["eps"].reshape(-1) * np.abs(cf) ** 2) + i_["tot"].real)
        e_abs = float(Q["e"][s] - Q["shift"])
        rhs = e_abs * cc - (VL_in_out + VNL_in_out).real
        rows.append(dict(band=int(Q["bands1"][s]), x_qe=float(Q["x"][s]), parity="even" if Q["parity"][s] > 0 else "odd", norm_in=cc,
                         cMc_tot=i_["tot"].real, cMc_L=i_["L"].real, cMc_NL=i_["NL"].real,
                         VL_in_in=VL_in_in.real, VNL_in_in=VNL_in_in.real, dL_in_in=abs(VL_in_in - i_["L"]), dNL_in_in=abs(VNL_in_in - i_["NL"]),
                         VL_in_out=VL_in_out.real, VNL_in_out=VNL_in_out.real, VL_in_out_im=VL_in_out.imag, VNL_in_out_im=VNL_in_out.imag,
                         lhs_cHc=lhs, rhs=rhs, s_align_eV=(rhs - lhs) / cc, s_align_unaligned_eV=(rhs - lhs) / cc + Q["shift"]))
        del A_in, psi_in, A, psi
    A2["states" if conv == "prod" else "states_Lx81"][str(nb)] = rows
    sal = np.array([r["s_align_eV"] for r in rows])
    log(f"[A.2 {conv}] base {nb} : {ns} états en {time.time()-t0:.0f} s ; max|(i)-(ii)| in-in L {max(r['dL_in_in'] for r in rows):.1e}, NL {max(r['dNL_in_in'] for r in rows):.1e} eV ; "
        f"s_align médiane {np.median(sal)*1e3:+.2f} meV (min {sal.min()*1e3:+.2f}, max {sal.max()*1e3:+.2f}) ; décalage Lu {Q['shift']*1e3:+.2f} meV, <ΔV>_3D {geo['mean3d_meV']:+.2f} meV")
    for lab in hl:
        r = rows[hl[lab]]; log(f"       {lab}: Σ|c|² {r['norm_in']:.4f} ; c†(M/81)c tot {r['cMc_tot']:+.4f} (L {r['cMc_L']:+.4f}, NL {r['cMc_NL']:+.4f}) ; <in|ΔV|out> L {r['VL_in_out']:+.4f}, NL {r['VNL_in_out']:+.4f} ; s {r['s_align_eV']*1e3:+.2f} meV")


def write_tables_a(out, Q):
    L = []
    L.append("#### A.1 — recouvrements des états QE (9×9 lacune, bandes 297–326) sur les bases repliées\n")
    for key, name in (("A1", "convention de production"), ("A1_Lx81", "diagnostique M^L × 81 (hors prompt)")):
        L.append(f"\n*{name}* :\n")
        L.append("| bande | parité | ε − E_D QE (eV) | w₂ QE | base | δ = 1 − Σ\\|c\\|² | fuite autre parité | R − ε_QE (meV) | résidu (eV) | 5 plus grands poids : (ε_j − E_D ; \\|⟨φ_j\\|c⟩\\|²/c†c) | Σ 5 |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for s in range(len(out[key]["16"])):
            for nb in ("16", "20"):
                r = out[key][nb][s]
                L.append(f"| {r['band']} | {'σ' if r['parity']=='even' else 'π'} | {r['x_qe']:+.3f} | {r['w2_qe']:.3f} | {nb} | {r['delta']:.4f} | {r['leak_other_parity']:.1e} | {r['rayleigh_minus_eps']*1e3:+.1f} | {r['residual']:.4f} | "
                         + ", ".join(f"({e:+.3f} ; {w:.3f})" for e, w in r["top_x"]) + f" | {r['sum_top']:.3f} |")
    L.append("\n#### A.2 — vérification sur états de Bloch purs : ⟨nk\\|ΔV\\|n′k′⟩ direct contre M/81 (eV)\n")
    L.append("| base | (n, k, n′, k′) | ΔV^L direct | M^L/81 | Δ | ΔV^NL direct | M^NL/81 | Δ |")
    L.append("|---|---|---|---|---|---|---|---|")
    for nb in ("16", "20"):
        for r in out["A2"]["pure"].get(nb, []):
            L.append(f"| {nb} | {tuple(r['pair'])} | {_re(r['VL_direct']):+.6f} | {_re(r['ML_over_N']):+.6f} | {r['dL']:.1e} | {_re(r['VNL_direct']):+.6f} | {_re(r['MNL_over_N']):+.6f} | {r['dNL']:.1e} |")
        v = out["A2"]["verified"].get(nb)
        if v:
            L.append(f"\nBase {nb} : max\\|Δ\\| = {v['max_abs_diff_eV']:.2e} eV (L {v['max_dL']:.1e}, NL {v['max_dNL']:.1e} ; seuil {v['tol']:.0e}) → {'OK' if v['ok'] else 'ÉCHEC (A.2 arrêté pour cette base, convention de production)'} ; rapport L direct / (M^L/81) : {', '.join(f'{x:.4f}' for x in v['ratio_L_direct_over_ML_over_N'])}.\n")
    L.append("\nDiagnostique (hors prompt) : même vérification avec M^L × N_cells (= 81) + M^NL :\n")
    L.append("| base | (n, k, n′, k′) | ΔV^L direct | 81·M^L/81 | Δ | ΔV^NL direct | M^NL/81 | Δ |")
    L.append("|---|---|---|---|---|---|---|---|")
    for nb in ("16", "20"):
        for r in out["A2"]["pure_Lx81"].get(nb, []):
            L.append(f"| {nb} | {tuple(r['pair'])} | {_re(r['VL_direct']):+.6f} | {_re(r['ML_over_N']):+.6f} | {r['dL']:.1e} | {_re(r['VNL_direct']):+.6f} | {_re(r['MNL_over_N']):+.6f} | {r['dNL']:.1e} |")
        v = out["A2"]["verified_Lx81"].get(nb)
        if v:
            L.append(f"\nBase {nb}, M^L × 81 : max\\|Δ\\| = {v['max_abs_diff_eV']:.2e} eV (L {v['max_dL']:.1e}, NL {v['max_dNL']:.1e}) → {'OK' if v['ok'] else 'ÉCHEC'}.\n")
    L.append("\n#### A.2 — les 30 états : ⟨Ψ_in\\|ΔV\\|Ψ_in⟩ par (i) c†(M/81)c et (ii) application directe ; ⟨Ψ_in\\|ΔV\\|Ψ_out⟩ ; constante d'alignement (eV)\n")
    L.append("| bande | parité | ε − E_D QE | base | Σ\\|c\\|² | (i) tot | (i) L | (i) NL | (ii) L | (ii) NL | \\|(i)−(ii)\\| L / NL | ⟨in\\|ΔV^L\\|out⟩ | ⟨in\\|ΔV^NL\\|out⟩ | c†H_in c | ε_QE c†c − ⟨in\\|ΔV\\|out⟩ | s (meV) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for key, name in (("states", "convention de production"), ("states_Lx81", "diagnostique M^L × 81")):
        L.append(f"\n*{name}* :\n")
        L.append("| bande | parité | ε − E_D QE | base | Σ\\|c\\|² | (i) tot | (i) L | (i) NL | (ii) L | (ii) NL | \\|(i)−(ii)\\| L / NL | ⟨in\\|ΔV^L\\|out⟩ | ⟨in\\|ΔV^NL\\|out⟩ | c†H_in c | ε_QE c†c − ⟨in\\|ΔV\\|out⟩ | s (meV) |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for nb in ("16", "20"):
            for r in out["A2"][key].get(nb, []):
                L.append(f"| {r['band']} | {'σ' if r['parity']=='even' else 'π'} | {r['x_qe']:+.3f} | {nb} | {r['norm_in']:.4f} | {r['cMc_tot']:+.4f} | {r['cMc_L']:+.4f} | {r['cMc_NL']:+.4f} | {r['VL_in_in']:+.4f} | {r['VNL_in_in']:+.4f} | "
                         f"{r['dL_in_in']:.1e} / {r['dNL_in_in']:.1e} | {r['VL_in_out']:+.4f} | {r['VNL_in_out']:+.4f} | {r['lhs_cHc']:+.4f} | {r['rhs']:+.4f} | {r['s_align_eV']*1e3:+.2f} |")
        if not any(out["A2"][key].get(nb) for nb in ("16", "20")):
            L.append("(non calculé : vérification sur états purs en échec)")
    L.append("\n#### A.3 — comptages et entrelacement (énergies ε − E_D ; QE aligné Lu 1,0 Å ; modèles sur leur E_D)\n")
    A3 = out["A3"]
    L.append(f"QE défaut : {A3['n_qe']} états, {A3['n_qe_even']} pairs, {A3['n_qe_odd']} impairs (parité max dev {A3['parity_dev_max']:.1e}).\n")
    L.append("| modèle | parité | QE sous E_D − 3 | modèle sous E_D − 3 | modèle sous E_D − 4 | fenêtre [−3, +1] QE / modèle | min_k[λ_k − λ_{k}(QE)] (fenêtre ; tout) | s = 1 | s = 2 |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for tag in A3["models"]:
        for lab in ("even", "odd"):
            r = A3["models"][tag][lab]; sh = r["shifts"]
            def ms(s):
                q = sh.get(s, sh.get(str(s))); return f"{q.get('min_window', float('nan')):+.3f} ; {q['min_all']:+.3f}"
            L.append(f"| {tag} | {'σ' if lab=='even' else 'π'} | {r['n_ref_below_lo']} | {r['n_model_below_lo']} | {r['n_model_below_lo_m4']} | {len(r['ref_window'])} / {len(r['model_window'])} | {ms(0)} | {ms(1)} | {ms(2)} |")
    L.append("\nSpectres triés sur [−4, +1] côte à côte (QE ; modèle), par parité :\n")
    for tag in A3["models"]:
        for lab in ("even", "odd"):
            r = A3["models"][tag][lab]
            L.append(f"- {tag}, {'σ' if lab=='even' else 'π'} : QE ({len(r['ref_window_m4'])}) " + ", ".join(f"{v:+.3f}" for v in r["ref_window_m4"]))
            L.append(f"- {tag}, {'σ' if lab=='even' else 'π'} : modèle ({len(r['model_window_m4'])}) " + ", ".join(f"{v:+.3f}" for v in r["model_window_m4"]))
    L.append("\n#### A.4 — tableau D4 de R4 refait avec l'index d'ondes planes corrigé (ε − E_D ; w₂ ; seuil %.4f)\n" % Q["thr"])
    L.append("| variante | pairs (σ) localisés : (ε − E_D ; w₂ corrigé ; w₂ R4) | impairs (π) localisés |")
    L.append("|---|---|---|")
    V = out["A4"]["variants"]; r4w = out["A4"]["r4_w2"]
    for tag in ("QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3", "a1_tot_Lx81", "a3_20_Lx81"):
        cells = []
        for lab in ("even", "odd"):
            b = V[tag]["blocks"][lab]; items = []
            for x, w in b["localized"]:
                old = ""
                if tag in r4w:
                    xr = np.array(r4w[tag][lab]["x"]); wr = np.array(r4w[tag][lab]["w2"])
                    if len(xr):
                        j = int(np.argmin(np.abs(xr - x))); old = f" ; {wr[j]:.3f}"
                items.append(f"{x:+.3f} ({w:.3f}{old})")
            cells.append(", ".join(items) if items else "aucun")
        L.append(f"| {tag} | {cells[0]} | {cells[1]} |")
    L.append("\nTous les états de la fenêtre par variante et parité (ε − E_D ; w₂ corrigé) :\n")
    for tag in ("a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3", "a1_tot_Lx81", "a3_20_Lx81"):
        for lab in ("even", "odd"):
            b = V[tag]["blocks"][lab]
            L.append(f"- {tag}, {lab} : " + ", ".join(f"{x:+.3f} ({w:.3f})" for x, w in zip(b["x"], b["w2"])))
    ver = out["A4"]["verification_a2_vs_a1"]
    L.append(f"\nVérification w₂(a2) = w₂(a1) sur les états localisés de (a1) : max\\|Δw₂\\| = {ver['max_dw2']:.4f} (seuil 1e-3) → {'OK' if ver['ok'] else 'ÉCART'} ; " + " ; ".join(f"{r['parity']} {r['x_a1']:+.3f} : {r['w2_a1']:.4f} → {r['w2_a2']:.4f} (Δx {r['dx']*1e3:+.1f} meV)" for r in ver["rows"]) + "\n")
    L.append("\n#### A.5 — moyenne diagonale de M^L soustraite (base 16 bandes)\n")
    A5 = out["A5"]
    L.append(f"Moyenne diagonale de M^L grossier (norme unit_cell) : {A5['diag_mean_ML_eV']*1e3:.2f} meV (R4 J1 : 76,58 meV).\n")
    L.append("| variante | constante soustraite de la diagonale de M (meV) | décalage rigide observé médian / min / max (meV) | attendu = −c/81 (meV) | max\\|Δε − médiane\\| (eV) | états localisés (a1) : Δx après retrait du décalage (meV) |")
    L.append("|---|---|---|---|---|---|")
    for tag, r in A5["variants"].items():
        L.append(f"| {tag} | {r['constant_subtracted_eV']*1e3:.2f} | {r['rigid_shift_median_eV']*1e3:+.3f} / {r['rigid_shift_min_eV']*1e3:+.3f} / {r['rigid_shift_max_eV']*1e3:+.3f} | {r['expected_shift_eV']*1e3:+.3f} | {r['max_abs_dev_after_shift_eV']:.1e} | "
                 + ", ".join(f"{q['parity']} {q['x_a1']:+.3f}: {(q['dx'] + r['constant_subtracted_eV']/N_CELLS)*1e3:+.4f}" for q in r["localized"]) + " |")
    with open(os.path.join(WORK, "a", "A_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log("[A] tables a/A_tables.md")


def figs_a(out, Q):
    plt, pal = fig_style()
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.4))
    for ax, nb in zip(axs, ("16", "20")):
        rows = out["A1"][nb]; x = np.array([r["x_qe"] for r in rows]); dl = np.array([r["delta"] for r in rows]); p = np.array([r["parity"] == "even" for r in rows])
        ax.scatter(x[~p], dl[~p], s=14, color=pal.ORANGE, label="impaire (π)"); ax.scatter(x[p], dl[p], s=14, color=pal.NAVY, label="paire (σ)")
        ax.set_yscale("log"); ax.set_xlabel(r"$\varepsilon - E_D$ (QE, eV)"); ax.set_title(f"base {nb} bandes", fontsize=9)
    axs[0].set_ylabel(r"$\delta = 1 - \sum |c_{nk}|^2$"); axs[0].legend(fontsize=6); fig.tight_layout(); savefig(fig, "a_delta")
    V = out["A4"]["variants"]; order = ["QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3", "a1_tot_Lx81", "a3_20_Lx81"]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.2), sharey=True)
    for ax, lab, ttl in zip(axs, ("odd", "even"), ("(a) impairs (π)", "(b) pairs (σ)")):
        for i, tag in enumerate(order):
            b = V[tag]["blocks"][lab]; x = np.array(b["x"]); w = np.array(b["w2"])
            ax.scatter(np.full(len(x), i), x, s=2 + 60 * w, color=pal.NAVY if tag != "QE_D1" else pal.ORANGE, alpha=0.7, lw=0)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, fontsize=6); ax.set_title(ttl); ax.set_ylim(-3, 1); ax.axhline(0, color=pal.MUTED, lw=0.5)
    axs[0].set_ylabel(r"$\varepsilon - E_D$ (eV)"); fig.tight_layout(); savefig(fig, "a_ladder_corrected")


# ----------------------------------------------------------------------------------------------- J2 : C
def sc_paths(N):
    t = f"{N}x{N}"
    return dict(N=N, tag=t, d=f"{SCRATCH}/defect_{t}_d/defect_{t}_d.save", p=f"{SCRATCH}/defect_{t}_p/defect_{t}_p.save",
                vd=f"{GQ}/defects/super_cell/{t}/defective/Vks_{t}_d", vp=f"{GQ}/defects/super_cell/{t}/pristine/Vks_{t}_p",
                scf_d=f"{GQ}/defects/super_cell/{t}/defective/scf.out", scf_p=f"{GQ}/defects/super_cell/{t}/pristine/scf.out")


def scf_wall(path):
    import re
    t = None
    for line in open(path, errors="ignore"):
        if "PWSCF" in line and "WALL" in line:
            m = re.search(r"CPU\s+(?:(\d+)h)?\s*(?:(\d+)m)?\s*([\d.]+)s\s+WALL", line)
            if m:
                t = 3600 * float(m.group(1) or 0) + 60 * float(m.group(2) or 0) + float(m.group(3))
    return t


def window_states_gamma(save, E_D, shift, ng, mask2, mask1, lo=WIN[0], hi=WIN[1], z0=0.0):
    eigs, ef, _ = qg.get_eigenvalues_spin(save); e = eigs[0] * HA2EV; ef *= HA2EV
    x = e - shift - E_D; sel = np.where((x >= lo) & (x <= hi))[0]
    if len(sel) == 0:
        return dict(e=e, ef=ef, x=x, n=0)
    b0, b1 = int(sel[0]), int(sel[-1]) + 1
    t0 = time.time(); Cw, mill, go, at = qg.read_wfc_gamma(save, bands=(b0, b1))
    par = qg.mirror_parity_z(Cw, mill, z0_red=z0, gamma_only=go).real
    rho = qg.density_2d(Cw, mill, ng, gamma_only=go, workers=NT); del Cw
    w2 = np.array([r[mask2].sum() for r in rho]); w1 = np.array([r[mask1].sum() for r in rho]); del rho
    return dict(e=e, ef=ef, x=x[b0:b1], b0=b0, b1=b1, n=b1 - b0, parity=par, w2=w2, w1=w1, bands1=np.arange(b0 + 1, b1 + 1), read_s=time.time() - t0,
                parity_dev=float(np.abs(1 - np.abs(par)).max()), n_even=int((par > 0).sum()), n_odd=int((par < 0).sum()))


def cmd_c(a):
    d = ensure("c"); out = dict(sizes={}); t_start = time.time()
    # références de la maille 9x9 et v_F du dense 27x27
    e9 = qe_io.get_eigenvalues(UC9) * HA2EV; k9 = qe_io.get_k_red(UC9); iK9 = kindex(k9, K_RED); iG9 = kindex(k9, (0, 0, 0))
    E_K9 = float(0.5 * (e9[3, iK9] + e9[4, iK9])); B_uc, _ = qe_io.get_B_volume(UC9)
    cfg = load_production(verbose=False); dp = dense_paths(cfg, "9x9"); k27 = qe_io.get_k_red(dp["uc"]); e27 = qe_io.get_eigenvalues(dp["uc"]) * HA2EV
    iKd = kindex(k27, K_RED); iKn = kindex(k27, K_RED + np.array([1 / 27, 0, 0])); dk = np.linalg.norm(B_uc[:, 0]) / 27 / BOHR
    vF = float((e27[4, iKn] - e27[3, iKn]) / (2 * dk))
    out["references"] = dict(E_K_uc9=E_K9, iK9=iK9, e_uc9_K_bands34=[float(e9[3, iK9]), float(e9[4, iK9])], vF_eVA=vF, vF_def=f"(ε_5 − ε_4)/(2|δk|) au point K + b1/27 de la maille dense 27×27, |δk| = {dk:.4f} Å⁻¹",
                             uc9_gamma_lowest8=[float(v) for v in np.sort(e9[:, iG9])[:8]])
    log(f"[C] E_K(maille 9x9) {E_K9:.6f} eV ; ħv_F {vF:.3f} eV Å (|δk| {dk:.4f} Å⁻¹) ; K cart {np.linalg.norm(B_uc @ K_RED)/BOHR:.4f} Å⁻¹")
    S = {}
    for N in range(5, 13):
        P = sc_paths(N); t0 = time.time(); r = dict(N=N, c=1.0 / (2 * N * N))
        A_b, Om = qe_io.get_A_volume(P["p"]); A_A = A_b * BOHR; x_p = qe_io.get_x_red(P["p"]); x_d = qe_io.get_x_red(P["d"])
        ng = tuple(int(v) for v in qe_io.get_ngfft(P["d"])); assert ng == tuple(int(v) for v in qe_io.get_ngfft(P["p"]))
        s_vac, i_vac_p, _ = al.vacancy_site(x_p, x_d, A_A); z0 = float(np.mean(x_d[:, 2]))
        eP, efP, _ = qg.get_eigenvalues_spin(P["p"]); eP = eP[0] * HA2EV; efP *= HA2EV
        eD, efD, _ = qg.get_eigenvalues_spin(P["d"]); eD = eD[0] * HA2EV; efD *= HA2EV
        r.update(nat_d=int(len(x_d)), nat_p=int(len(x_p)), nbnd_p=int(len(eP)), nbnd_d=int(len(eD)), E_F_p=float(efP), E_F_d=float(efD), s_vac=s_vac.tolist(), ngfft=list(ng), A_sc_A2=float(abs(np.linalg.det(A_A[:2, :2]))))
        # E_D : quadruplet (3m) et route « maille »
        if N % 3 == 0:
            E_Dq, qidx, spread = ax.dirac_quadruplet(eP, efP); r["E_D_quadruplet"] = dict(E_D=E_Dq, bands_1based=[i + 1 for i in qidx], spread_eV=spread, E_D_minus_EF=E_Dq - efP)
        uc = UC_SAVE[N]; e_uc = qe_io.get_eigenvalues(uc) * HA2EV; k_uc = qe_io.get_k_red(uc); iGu = kindex(k_uc, (0, 0, 0))
        dG, dG_res, nlow = ax.gamma_shift(e_uc[:, iGu], e9[:, iG9], 8)
        if N != 11:
            folded = np.sort(e_uc.reshape(-1)); dF, dF_res, nF = al.rigid_shift_fit(eP, folded, E_K9 + dG - 4.0); dF_note = f"médiane sur {nF} états sous E_D − 4"
        else:
            dF = ax.lowest_state_shift(eP, e_uc[:, iGu]); dF_res = float("nan"); nF = 1; dF_note = UC_NOTE[11]
        E_D_uc = E_K9 + dG + dF
        iKu = kindex(k_uc, K_RED); dKu = float(np.linalg.norm(np.mod(k_uc[iKu] - K_RED + 0.5, 1) - 0.5))
        r["E_D_uc_route"] = dict(E_K_uc9=E_K9, delta_gamma=dG, delta_gamma_resid=dG_res, nlow=nlow, delta_fold=dF, delta_fold_resid=dF_res, n_fold=nF, note=dF_note, E_D=E_D_uc, uc_save=uc, uc_nk=int(len(k_uc)))
        if dKu < 1e-4:
            E_Ku = float(0.5 * (e_uc[3, iKu] + e_uc[4, iKu])); r["E_D_uc_route"]["E_K_uc_direct"] = E_Ku; r["E_D_uc_route"]["E_D_from_uc_K"] = E_Ku + dF; r["E_D_uc_route"]["k_dev_K"] = dKu
        if N % 3 == 0:
            E_D = r["E_D_quadruplet"]["E_D"]; r["E_D_source"] = "quadruplet"
        elif "E_D_from_uc_K" in r["E_D_uc_route"]:
            E_D = r["E_D_uc_route"]["E_D_from_uc_K"]; r["E_D_source"] = "K de la maille N×N + Δ_fold (même run)"
        else:
            E_D = E_D_uc; r["E_D_source"] = "maille (E_K9 + Δ_Γ + Δ_fold)"
        r["E_D"] = E_D
        # gap à Γ de la parfaite ; distance de la grille à K
        # gap à Γ de la parfaite : n_occ = nelec/2 (comptage d'électrons, pas E_F : pour N = 3m le quadruplet est à demi rempli)
        import xml.etree.ElementTree as ET
        nelec = float(ET.parse(os.path.join(P["p"], "data-file-schema.xml")).getroot().find(".//output/band_structure/nelec").text)
        n_occ = int(round(nelec / 2)); eS = np.sort(eP); homo = float(eS[n_occ - 1]); lumo = float(eS[n_occ])
        r["gamma_gap_p"] = dict(homo=homo, lumo=lumo, gap=lumo - homo, n_occ=n_occ, nelec=nelec, homo_minus_ED=homo - E_D, lumo_minus_ED=lumo - E_D,
                                n_below_EF=int((eP <= efP).sum()), homo_EF_def=float(eP[eP <= efP].max()), lumo_EF_def=float(eP[eP > efP].min()))
        kg = np.array([(i / N, j / N, 0.0) for i in range(N) for j in range(N)]); dmin = np.inf
        for G1 in (-1, 0, 1):
            for G2 in (-1, 0, 1):
                dmin = min(dmin, np.linalg.norm((B_uc @ (kg + np.array([G1, G2, 0]) - K_RED).T).T, axis=1).min())
        dmin /= BOHR; r["grid_to_K"] = dict(min_dk_A=float(dmin), hbar_vF_dk_eV=float(vF * dmin))
        # alignement Lu
        V_P = load_pot_eV(P["vp"]); V_D = load_pot_eV(P["vd"]); assert V_P.shape == ng
        ra = al.far_atom_alignment(V_D, V_P, x_d, x_p, A_A, (0.5, 1.0)); shift = float(ra["shifts"][1.0])
        r["alignment"] = dict(i_far_1based=ra["i_far_d"] + 1, dist_far_A=ra["dist_far"], shift_05=float(ra["shifts"][0.5]), shift_10=shift, mean3d_diff_meV=float((V_D.mean() - V_P.mean()) * 1e3),
                              dV_site_eV=float((V_D - V_P)[tuple(np.rint(s_vac * np.array(ng)).astype(int) % np.array(ng))]))
        del V_P, V_D
        # états de la fenêtre
        mask2 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W2); mask1 = qg.inplane_disc_mask(ng, A_A, s_vac, R_W1)
        r["disc_area_fraction"] = dict(r2=float(mask2.mean()), r1=float(mask1.mean()))
        WP = window_states_gamma(P["p"], E_D, 0.0, ng, mask2, mask1, z0=z0); WD = window_states_gamma(P["d"], E_D, shift, ng, mask2, mask1, z0=z0)
        thr = 3.0 * float(WP["w2"].mean()); r["threshold"] = dict(thr=thr, mean_w2_P=float(WP["w2"].mean()), max_w2_P=float(WP["w2"].max()))
        for tag, W in (("P", WP), ("D", WD)):
            r[f"window_{tag}"] = dict(n=int(W["n"]), bands=[int(W["b0"]) + 1, int(W["b1"])], n_even=W["n_even"], n_odd=W["n_odd"], parity_dev=W["parity_dev"], x=[float(v) for v in W["x"]],
                                     parity=[float(v) for v in W["parity"]], w2=[float(v) for v in W["w2"]], w1=[float(v) for v in W["w1"]], read_s=W["read_s"])
        loc = [dict(band=int(WD["bands1"][j]), x=float(WD["x"][j]), parity=float(WD["parity"][j]), w2=float(WD["w2"][j]), w1=float(WD["w1"][j]), relEF=float(WD["e"][WD["b0"] + j] - efD)) for j in range(WD["n"]) if WD["w2"][j] > thr]
        r["localized"] = loc; odd = [q for q in loc if q["parity"] < 0]; even = [q for q in loc if q["parity"] > 0]
        r["n_localized"] = dict(even=len(even), odd=len(odd))
        r["pi_state"] = max(odd, key=lambda q: q["w2"]) if odd else None
        r["sigma_doublet"] = sorted(sorted(even, key=lambda q: -q["w2"])[:2], key=lambda q: q["x"]) if even else None
        if r["pi_state"] is not None:
            e_pi = r["pi_state"]["x"] + E_D; r["pi_in_gamma_gap"] = bool(homo < e_pi < lumo) if N % 3 != 0 else None; r["pi_state"]["e_aligned_abs"] = e_pi
        r["seconds"] = time.time() - t0; S[N] = r
        log(f"[C] N={N}: E_D {E_D:.5f} ({r['E_D_source']}) ; quadruplet {r.get('E_D_quadruplet', {}).get('E_D', float('nan')):.5f} ; route maille {E_D_uc:.5f} (Δ_Γ {dG*1e3:+.2f} meV, Δ_fold {dF*1e3:+.2f} meV, résidu {dF_res:.1e}) ; "
            f"gap Γ parfaite {lumo-homo:.3f} eV ; ħv_F·min|k−K| {vF*dmin:.3f} eV ; Lu {shift*1e3:+.2f} meV ; fenêtre P {WP['n']} ({WP['n_even']}σ/{WP['n_odd']}π), D {WD['n']} ({WD['n_even']}σ/{WD['n_odd']}π) ; seuil {thr:.4f} ; "
            f"localisés σ {len(even)} / π {len(odd)} ; π {r['pi_state'] and (round(r['pi_state']['x'],3), round(r['pi_state']['w2'],3))} ; σ {r['sigma_doublet'] and [(round(q['x'],3), round(q['w2'],3)) for q in r['sigma_doublet']]} ({time.time()-t0:.0f} s)")
    out["sizes"] = {str(N): S[N] for N in S}
    # ajustements par famille
    fits = {}
    for fam, Ns in (("3m", FAM_3M), ("non3m", FAM_N3M)):
        for what in ("pi", "sigma"):
            pts = [(N, S[N]["pi_state"]["x"]) for N in Ns if S[N]["pi_state"] is not None] if what == "pi" else [(N, float(np.mean([q["x"] for q in S[N]["sigma_doublet"]]))) for N in Ns if S[N]["sigma_doublet"]]
            if len(pts) >= 2:
                fits[f"{what}_{fam}"] = dict(N=[p[0] for p in pts], x=[p[1] for p in pts], fits=ax.size_fits([p[0] for p in pts], [p[1] for p in pts]))
                log(f"[C] ajustement {what} {fam} : N {[p[0] for p in pts]} x {[round(p[1],3) for p in pts]} ; 1/N : ε_∞ {fits[f'{what}_{fam}']['fits'][1]['eps_inf']:+.3f} (rms {fits[f'{what}_{fam}']['fits'][1]['rms']:.3f}) ; 1/N² : ε_∞ {fits[f'{what}_{fam}']['fits'][2]['eps_inf']:+.3f} (rms {fits[f'{what}_{fam}']['fits'][2]['rms']:.3f})")
    out["fits"] = fits
    # coût 15x15 / 18x18
    Ns = np.arange(5, 13); td = np.array([scf_wall(sc_paths(N)["scf_d"]) for N in Ns]); tp = np.array([scf_wall(sc_paths(N)["scf_p"]) for N in Ns])
    cost = dict(N=Ns.tolist(), wall_d_s=td.tolist(), wall_p_s=tp.tolist(), mpi_ranks=64)
    for lab, t in (("d", td), ("p", tp)):
        p = np.polyfit(np.log(Ns), np.log(t), 1)
        cost[f"fit_{lab}"] = dict(exponent=float(p[0]), t15_min=float(np.exp(np.polyval(p, np.log(15))) / 60), t18_min=float(np.exp(np.polyval(p, np.log(18))) / 60),
                                  from_9x9_exp_fit_15_min=float(t[4] * (15 / 9) ** p[0] / 60), from_9x9_exp_fit_18_min=float(t[4] * (18 / 9) ** p[0] / 60),
                                  from_9x9_exp3_15_min=float(t[4] * (15 / 9) ** 3 / 60), from_9x9_exp3_18_min=float(t[4] * (18 / 9) ** 3 / 60))
    w12 = os.path.getsize(f"{SCRATCH}/defect_12x12_d/defect_12x12_d.save/wfc1.hdf5") / 1e9
    cost["wfc1_GB"] = dict(N12=w12, N15=w12 * (15 / 12) ** 4, N18=w12 * (18 / 12) ** 4, law="∝ N⁴ (nbnd × ondes planes)")
    cost["atoms"] = dict(N15=2 * 225 - 1, N18=2 * 324 - 1)
    out["cost_15_18"] = cost
    log(f"[C] coût : exposants d {cost['fit_d']['exponent']:.2f}, p {cost['fit_p']['exponent']:.2f} ; 15x15 {cost['fit_d']['t15_min']:.0f} / {cost['fit_p']['t15_min']:.0f} min ; 18x18 {cost['fit_d']['t18_min']:.0f} / {cost['fit_p']['t18_min']:.0f} min (64 rangs)")
    out["timing_s"] = time.time() - t_start
    save_json(os.path.join(d, "c_results.json"), out)
    write_tables_c(out); figs_c(out)


def write_tables_c(out):
    S = out["sizes"]; L = []
    L.append("#### C — E_D, alignement et états localisés par taille (ε − E_D en eV ; alignement Lu 1,0 Å ; seuil = 3 × ⟨w₂⟩ de la parfaite N×N)\n")
    L.append("| N | c = 1/(2N²) | E_D (eV) [source] | E_D quadruplet | E_D maille (Δ_Γ ; Δ_fold, résidu) | E_D par K de la maille N×N | décalage Lu 1,0 / 0,5 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ (meV) | seuil w₂ | état π localisé (ε − E_D ; w₂ ; w₁) | doublet σ (ε − E_D ; w₂) | localisés σ / π | fenêtre P (σ/π) ; D (σ/π) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for N in range(5, 13):
        r = S[str(N)]; u = r["E_D_uc_route"]; q = r.get("E_D_quadruplet")
        pi = r["pi_state"]; sg = r["sigma_doublet"]
        qcell = f"{q['E_D']:.5f} (étalement {q['spread_eV']:.1e})" if q else "—"
        L.append(f"| {N} | {r['c']:.5f} | {r['E_D']:.5f} [{r['E_D_source']}] | {qcell} | "
                 f"{u['E_D']:.5f} ({u['delta_gamma']*1e3:+.2f} meV ; {u['delta_fold']*1e3:+.2f} meV, {u['delta_fold_resid']:.1e}) | {_f(u.get('E_D_from_uc_K'), 5)} | {r['alignment']['shift_10']*1e3:+.2f} / {r['alignment']['shift_05']*1e3:+.2f} | {r['alignment']['mean3d_diff_meV']:+.2f} | {r['threshold']['thr']:.4f} | "
                 + (f"{pi['x']:+.3f} ; {pi['w2']:.3f} ; {pi['w1']:.3f} (bande {pi['band']})" if pi else "aucun") + " | " + (", ".join(f"{s['x']:+.3f} ; {s['w2']:.3f}" for s in sg) if sg else "aucun")
                 + f" | {r['n_localized']['even']} / {r['n_localized']['odd']} | {r['window_P']['n']} ({r['window_P']['n_even']}/{r['window_P']['n_odd']}) ; {r['window_D']['n']} ({r['window_D']['n_even']}/{r['window_D']['n_odd']}) |")
    L.append("\n#### C — gap à Γ de la parfaite N×N et distance de la grille à K\n")
    L.append("| N | HOMO Γ − E_D (eV) | LUMO Γ − E_D (eV) | gap Γ (eV) | n_occ = nelec/2 (états ≤ E_F) | min\\|k − K\\| de la grille (Å⁻¹) | ħv_F·min\\|k − K\\| (eV) | π localisé dans le gap ? | E_F p / d (eV) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for N in range(5, 13):
        r = S[str(N)]; g = r["gamma_gap_p"]; k = r["grid_to_K"]
        L.append(f"| {N} | {g['homo_minus_ED']:+.4f} | {g['lumo_minus_ED']:+.4f} | {g['gap']:.4f} | {g['n_occ']} ({g['n_below_EF']}) | {k['min_dk_A']:.4f} | {k['hbar_vF_dk_eV']:.3f} | {('oui' if r.get('pi_in_gamma_gap') else 'non') if r.get('pi_in_gamma_gap') is not None else '— (3m)'} | {r['E_F_p']:.5f} / {r['E_F_d']:.5f} |")
    L.append(f"\nħv_F = {out['references']['vF_eVA']:.3f} eV Å ({out['references']['vF_def']}).\n")
    L.append("\n#### C — tous les états localisés (w₂ > seuil) par taille\n")
    for N in range(5, 13):
        r = S[str(N)]
        L.append(f"- N = {N} : " + ("; ".join(f"b{q['band']} {'σ' if q['parity']>0 else 'π'} {q['x']:+.3f} (w₂ {q['w2']:.3f}, w₁ {q['w1']:.3f})" for q in r["localized"]) or "aucun"))
    L.append("\n#### C — ajustements ε(N) = ε_∞ + a/N^p par famille (eV)\n")
    L.append("| état | famille | N | ε − E_D | p = 1 : ε_∞ ; a ; rms ; max résidu | p = 2 : ε_∞ ; a ; rms ; max résidu |")
    L.append("|---|---|---|---|---|---|")
    for key, f in out["fits"].items():
        what, fam = key.split("_"); f1 = f["fits"].get(1, f["fits"].get("1")); f2 = f["fits"].get(2, f["fits"].get("2"))
        L.append(f"| {what} | {fam} | {f['N']} | {', '.join(f'{x:+.3f}' for x in f['x'])} | {f1['eps_inf']:+.4f} ; {f1['a']:+.3f} ; {f1['rms']:.4f} ; {f1['max_resid']:.4f} | {f2['eps_inf']:+.4f} ; {f2['a']:+.3f} ; {f2['rms']:.4f} ; {f2['max_resid']:.4f} |")
    c = out["cost_15_18"]
    L.append("\n#### C — coût estimé des super-cellules 15×15 (449 atomes) et 18×18 (647 atomes), rien lancé\n")
    L.append("| série | murs scf 5…12 (s, 64 rangs) | exposant t ∝ N^p | 15×15 (min) | 18×18 (min) | depuis la 9×9 seule, même exposant | depuis la 9×9, exposant 3 |")
    L.append("|---|---|---|---|---|---|---|")
    for lab, name in (("d", "défaut"), ("p", "parfaite")):
        f = c[f"fit_{lab}"]
        L.append(f"| {name} | {', '.join(f'{t:.0f}' for t in c[f'wall_{lab}_s'])} | {f['exponent']:.2f} | {f['t15_min']:.0f} | {f['t18_min']:.0f} | {f['from_9x9_exp_fit_15_min']:.0f} / {f['from_9x9_exp_fit_18_min']:.0f} | {f['from_9x9_exp3_15_min']:.0f} / {f['from_9x9_exp3_18_min']:.0f} |")
    L.append(f"\n`wfc1.hdf5` : 12×12 {c['wfc1_GB']['N12']:.2f} Go → 15×15 {c['wfc1_GB']['N15']:.1f} Go, 18×18 {c['wfc1_GB']['N18']:.1f} Go par super-cellule ({c['wfc1_GB']['law']}).\n")
    with open(os.path.join(WORK, "c", "C_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log("[C] tables c/C_tables.md")


def figs_c(out):
    plt, pal = fig_style(); S = out["sizes"]
    fig, ax = plt.subplots(figsize=(4.2, 3.4))
    for N in range(5, 13):
        r = S[str(N)]
        for q in r["localized"]:
            ax.scatter(1 / N, q["x"], s=6 + 80 * q["w2"], color=pal.NAVY if q["parity"] > 0 else pal.ORANGE, alpha=0.75, lw=0, marker="o" if N % 3 == 0 else "s")
    Nf = np.linspace(4.5, 40, 200)
    for key, ls in (("pi_3m", "-"), ("pi_non3m", "--")):
        if key in out["fits"]:
            ff = out["fits"][key]["fits"]; f1 = ff.get(1, ff.get("1")); ax.plot(1 / Nf, f1["eps_inf"] + f1["a"] / Nf, ls=ls, lw=0.8, color=pal.ORANGE, label=f"π {key.split('_')[1]} : 1/N")
            f2 = ff.get(2, ff.get("2")); ax.plot(1 / Nf, f2["eps_inf"] + f2["a"] / Nf ** 2, ls=ls, lw=0.8, color=pal.GREEN, label=f"π {key.split('_')[1]} : 1/N²")
    ax.axhline(0, color=pal.MUTED, lw=0.5); ax.set_xlabel("1/N"); ax.set_ylabel(r"$\varepsilon - E_D$ (eV)"); ax.set_xlim(0, 0.22); ax.set_ylim(-3, 1)
    ax.legend(fontsize=5); ax.set_title("états localisés : ronds N = 3m, carrés N ≠ 3m ; bleu σ, orange π", fontsize=7); fig.tight_layout(); savefig(fig, "size")


# ----------------------------------------------------------------------------------------------- J4 : B.2 (M à 128 bandes)
def cmd_b2ml(a):
    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import compute_ML_R_mpi_shared
    comm = MPI.COMM_WORLD; rank = comm.Get_rank(); t0 = time.time()
    if rank == 0:
        log(f"[B.2] M^L 128 bandes : compute_ML_R_mpi_shared({NB128}), {comm.Get_size()} rangs")
    M = compute_ML_R_mpi_shared(NB128, SC_P, VKS["p"], VKS["d"], subtract_mean=False, bands=None, io=qe_io, grid_block=int(os.environ.get("R5_BLOCK", "2000")))
    if rank == 0:
        d = ensure("b"); p = os.path.join(d, "M_L_9x9_nb128.npy")
        matrix_io.save_M(p, M, matrix_io.UNIT_CELL, part="M_L_coarse_nb128", nbnd=128, nk=81, uc_save=NB128, sc_p=SC_P, pot_p=VKS["p"], pot_d=VKS["d"], subtract_mean=False,
                         kernel="compute_ML_R_mpi_shared", date=time.strftime("%Y-%m-%d"), note="R5 B.2 : M^L grossier 128 bandes x 81 k, norme unit_cell, Hartree")
        log(f"[B.2] M^L sauvé {p} {M.shape} en {time.time()-t0:.0f} s ; max|M^L| {np.abs(M).max():.4e} Ha")


def cmd_b2nl(a):
    from electron_defect_interaction.defects.non_local import compute_M_NL
    from electron_defect_interaction.io.pseudo_io import read_upf
    t0 = time.time(); log(f"[B.2] M^NL 128 bandes : compute_M_NL({NB128})")
    M = compute_M_NL(NB128, SC_P, SC_D, UPF, io=qe_io, pseudo_reader=read_upf, bands=None)
    d = ensure("b"); p = os.path.join(d, "M_NL_9x9_nb128.npy")
    matrix_io.save_M(p, M, matrix_io.UNIT_CELL, part="M_NL_coarse_nb128", nbnd=128, nk=81, uc_save=NB128, sc_p=SC_P, sc_d=SC_D, upf=UPF, kernel="compute_M_NL",
                     date=time.strftime("%Y-%m-%d"), note="R5 B.2 : M^NL grossier 128 bandes x 81 k, norme unit_cell, Hartree")
    log(f"[B.2] M^NL sauvé {p} {M.shape} en {time.time()-t0:.0f} s ; max|M^NL| {np.abs(M).max():.4e} Ha")


def cmd_b2sum(a):
    d = ensure("b"); out = {}
    ML = matrix_io.load_M_checked(os.path.join(d, "M_L_9x9_nb128.npy"), require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    MNL = matrix_io.load_M_checked(os.path.join(d, "M_NL_9x9_nb128.npy"), require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    M = ML + MNL; p = os.path.join(d, "M_ed_9x9_nb128.npy")
    matrix_io.save_M(p, M, matrix_io.UNIT_CELL, part="M_ed_coarse_nb128", nbnd=128, nk=81, uc_save=NB128, date=time.strftime("%Y-%m-%d"), note="R5 B.2 : M^L + M^NL, 128 bandes x 81 k")
    nbk = 128 * 81
    for tag, X in (("L", ML), ("NL", MNL), ("tot", M)):
        Xf = X.reshape(nbk, nbk); out[f"herm_{tag}_Ha"] = float(np.abs(Xf - Xf.conj().T).max()); out[f"max_{tag}_Ha"] = float(np.abs(X).max())
    # bloc 16 x 16 contre M_ed_9x9 (juin) : valeurs singulières (jauge) et valeurs propres de H(a1) (physique, indépendantes de la jauge)
    eps_new = qe_io.get_eigenvalues(NB128) * HA2EV; eps16 = qe_io.get_eigenvalues(UC9) * HA2EV
    out["eps16_vs_new_max_eV"] = float(np.abs(eps16 - eps_new[:16]).max())
    cmp = {}
    for tag, X, ref in (("L", ML, coarse_M("L") / HA2EV), ("NL", MNL, coarse_M("NL") / HA2EV), ("tot", M, coarse_M("tot") / HA2EV)):
        X16 = X[:16, :, :16, :].reshape(16 * 81, 16 * 81); R16 = ref.reshape(16 * 81, 16 * 81)
        sv_x = np.linalg.svd(X16, compute_uv=False); sv_r = np.linalg.svd(R16, compute_uv=False)
        Hx = sf.bloch_folded_hamiltonian(eps_new[:16], X[:16, :, :16, :] * HA2EV, N_CELLS); Hr = sf.bloch_folded_hamiltonian(eps16, ref * HA2EV, N_CELLS)
        ex = np.linalg.eigvalsh(0.5 * (Hx + Hx.conj().T)); er = np.linalg.eigvalsh(0.5 * (Hr + Hr.conj().T))
        diag_x = np.array([X[n, k, n, k].real for n in range(16) for k in range(81)]); diag_r = np.array([ref[n, k, n, k].real for n in range(16) for k in range(81)])
        cmp[tag] = dict(sv_rel_max=float(np.abs(sv_x - sv_r).max() / sv_r.max()), sv_max=float(sv_r.max()), eig_H_a1_max_dev_eV=float(np.abs(ex - er).max()),
                        diag_max_dev_Ha=float(np.abs(diag_x - diag_r).max()), diag_mean_Ha=float(diag_x.mean()), diag_mean_ref_Ha=float(diag_r.mean()))
        log(f"[B.2] bloc 16x16 {tag} vs M grossier juin : VS rel {cmp[tag]['sv_rel_max']:.2e}, eig H(a1) {cmp[tag]['eig_H_a1_max_dev_eV']:.2e} eV, diag {cmp[tag]['diag_max_dev_Ha']:.2e} Ha")
    out["block16_vs_june"] = cmp; out["files"] = dict(L=os.path.join(d, "M_L_9x9_nb128.npy"), NL=os.path.join(d, "M_NL_9x9_nb128.npy"), tot=p, GB_each=M.nbytes / 1e9)
    save_json(os.path.join(d, "b2_results.json"), out)


# ----------------------------------------------------------------------------------------------- J5 : B.1 et B.3
def cmd_b1(a):
    d = ensure("b"); out = {}
    k81 = qe_io.get_k_red(UC9); kn = qe_io.get_k_red(NB128)
    out["k_list"] = dict(max_abs_diff=float(np.abs(kn - k81).max()), same_order=bool(np.abs(kn - k81).max() < 1e-8), nk=int(len(kn)))
    eps16 = qe_io.get_eigenvalues(UC9) * HA2EV; en = qe_io.get_eigenvalues(NB128) * HA2EV
    out["eps16_vs_new"] = dict(max_abs_eV=float(np.abs(en[:16] - eps16).max()), tol=1e-6, ok=bool(np.abs(en[:16] - eps16).max() <= 1e-6), nbnd_new=int(en.shape[0]))
    iK = kindex(kn, K_RED); iG = kindex(kn, (0, 0, 0)); iM = kindex(kn, (0.5, 0.0, 0.0)); E_D = float(0.5 * (en[3, iK] + en[4, iK]))
    out["E_D"] = E_D; out["band128"] = {lab: dict(e=float(en[127, i]), x=float(en[127, i] - E_D), k=kn[i].tolist()) for lab, i in (("Gamma", iG), ("K", iK), ("M", iM))}
    out["band_range"] = dict(min_x=float(en.min() - E_D), max_x=float(en.max() - E_D), band_by_band_min_x=[float(v) for v in en.min(1) - E_D], band_by_band_max_x=[float(v) for v in en.max(1) - E_D])
    # parité des 128 x 81 états
    t0 = time.time(); C, nG = qe_io.get_C_nk(NB128); G = qe_io.get_G_red(NB128); par = r4.uc_parity(C, nG, G); del C
    out["parity"] = dict(max_dev=float(np.abs(1 - np.abs(par)).max()), n_even=int((par > 0).sum()), n_odd=int((par < 0).sum()), n_undetermined=int((np.abs(par) < 0.5).sum()),
                         even_per_band=[int(v) for v in (par > 0).sum(1)], odd_per_band=[int(v) for v in (par < 0).sum(1)], seconds=time.time() - t0)
    nsc = open(os.path.join(WORK, "b", "nscf_nb128.out"), errors="ignore").read() if os.path.exists(os.path.join(WORK, "b", "nscf_nb128.out")) else ""
    out["nscf_out"] = dict(wall=scf_wall(os.path.join(WORK, "b", "nscf_nb128.out")) if nsc else None, n_not_converged=nsc.count("not converged"), n_k=nsc.count("k =") if nsc else None)
    out["save_size_GB"] = sum(os.path.getsize(p) for p in glob.glob(f"{NB128}/wfc*.hdf5")) / 1e9
    log(f"[B.1] k identiques {out['k_list']['same_order']} (max {out['k_list']['max_abs_diff']:.1e}) ; ε16 vs neuf max {out['eps16_vs_new']['max_abs_eV']:.2e} eV ; E_D {E_D:.6f} ; bande 128 : Γ {out['band128']['Gamma']['x']:+.3f}, K {out['band128']['K']['x']:+.3f}, M {out['band128']['M']['x']:+.3f} eV ; "
        f"parité dev max {out['parity']['max_dev']:.1e}, pairs {out['parity']['n_even']} / impairs {out['parity']['n_odd']} ; .save {out['save_size_GB']:.2f} Go ; nscf {out['nscf_out']}")
    save_json(os.path.join(d, "b1_results.json"), out); np.savez(os.path.join(d, "b1_parity.npz"), parity=par, eps=en, k=kn)


def cmd_b3(a):
    import r5_deltav_pw as dv
    d = ensure("b"); out = dict(n_list=NLIST); t_start = time.time()
    B = load_bases(with_dense=False, nb128=True); ng = B["ng"]; Om = B["Om"]; A_A = B["A_A"]; b = B["bases"][128]
    Q = qe_states_9x9(); ns = len(Q["e"]); hl = highlighted(Q); mask2 = qg.inplane_disc_mask(ng, A_A, Q["s_vac"], R_W2)
    Ms = {p: matrix_io.load_M_checked(os.path.join(d, f"M_{'ed' if p == 'tot' else p}_9x9_nb128.npy"), require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.EV) for p in ("tot", "L", "NL")}
    out["E_D_128"] = b["E_D"]; out["M_linearity_eV"] = float(np.abs(Ms["tot"] - Ms["L"] - Ms["NL"]).max())
    # recouvrements sur les 128 bandes (une fois)
    t0 = time.time(); c128 = np.zeros((ns, 128, 81), complex); Agrids = {}
    for s in range(ns):
        A = sp.sc_state_grid(Q["Cw"][s], Q["mill"], ng, gamma_only=Q["gamma_only"]); c128[s] = sp.bloch_overlaps(A, b["C"], b["nG"], b["flat"])
        if s in hl.values():
            Agrids[s] = A
    np.savez(os.path.join(d, "b3_overlaps.npz"), c128=c128, x_qe=Q["x"], parity=Q["parity"], bands=Q["bands1"])
    log(f"[B.3] recouvrements 128 bandes en {time.time()-t0:.0f} s ; δ(128) des états en évidence : " + ", ".join(f"{lab} {1-np.sum(np.abs(c128[s])**2):.4f}" for lab, s in hl.items()))
    dV, proj, flat_full, mill_full, geo = deltav_grid_and_projector(B, Q["mill"])
    psi_hl = {s: sp.grid_to_real(Agrids[s], ng, Om, NT) for s in Agrids}; B_hl = {s: proj.projections(Agrids[s][flat_full][None])[0] for s in Agrids}
    # référence QE (D1)
    d4 = json.load(open(os.path.join(R4DIR, "d4", "d4_results.json"))); qe_rec = d4["variants"]["QE_D1"]
    x_all = np.load(os.path.join(WORK, "a", "a_qe_all.npz")) if os.path.exists(os.path.join(WORK, "a", "a_qe_all.npz")) else None
    per_n = {}; prev = None
    for n in NLIST:
        t0 = time.time(); Mn = Ms["tot"][:n, :, :n, :]; eps_n = b["eps"][:n]; par_n = b["par"][:n].reshape(-1)
        H = sf.bloch_folded_hamiltonian(eps_n, Mn, N_CELLS); blocks, cpl, herm = r4.spectrum_blocks(H, par_n, b["E_D"]); Hh = 0.5 * (H + H.conj().T)
        W = window_weights(blocks, b["C"][:n], b["nG"], b["flat"], n, ng, mask2); V = {}; rec = record_variant(V, f"n{n}", blocks, W, b["E_D"], Q["thr"], dict(coupling_even_odd=cpl, dim=int(H.shape[0])))
        odd = rec["blocks"]["odd"]["localized"]; even = rec["blocks"]["even"]["localized"]
        pi = max(odd, key=lambda t: t[1]) if odd else None
        sig_near = [t for t in even if abs(t[0] - 0.101) < 0.3]
        r = dict(n=n, dim=int(H.shape[0]), coupling=cpl, blocks=rec["blocks"], pi_state=pi, sigma_localized=even, sigma_near_0p101=sig_near,
                 median_dx_vs_qe={lab: (float(np.median(np.abs(np.sort(np.array(rec["blocks"][lab]["x"])) - np.sort(np.array(qe_rec["blocks"][lab]["x"]))))) if len(rec["blocks"][lab]["x"]) == len(qe_rec["blocks"][lab]["x"]) else None) for lab in ("even", "odd")},
                 n_window={lab: rec["blocks"][lab]["n_window"] for lab in ("even", "odd")})
        if prev is not None:
            r["dpi_vs_prev"] = (pi[0] - prev["pi_state"][0]) if (pi and prev["pi_state"]) else None
            r["median_dx_vs_prev"] = {lab: (float(np.median(np.abs(np.sort(np.array(rec["blocks"][lab]["x"])) - np.sort(np.array(prev["blocks"][lab]["x"]))))) if len(rec["blocks"][lab]["x"]) == len(prev["blocks"][lab]["x"]) else None) for lab in ("even", "odd")}
        # A.1 à ce n
        rows = []
        for s in range(ns):
            lab = "even" if Q["parity"][s] > 0 else "odd"; idx = blocks[lab]["idx"]; cf = c128[s, :n, :].reshape(-1); cb = cf[idx]
            e_ref = (Q["e"][s] - Q["shift"]) - Q["E_D_SC"] + b["E_D"]
            q = bd.basis_diagnostics(Hh[np.ix_(idx, idx)], cb, e_ref, evals=blocks[lab]["x"] + b["E_D"], evecs=blocks[lab]["v"])
            q["delta"] = float(1 - np.vdot(cf, cf).real); q["rayleigh_minus_eps"] = q["rayleigh"] - e_ref; q["top_x"] = [(e - b["E_D"], w) for e, w in q["top"]]
            q.update(band=int(Q["bands1"][s]), x_qe=float(Q["x"][s]), parity=lab); rows.append(q)
        r["A1"] = rows
        # A.3 à ce n
        if x_all is not None:
            r["A3"] = {}
            for lab, sgn in (("even", 1), ("odd", -1)):
                xr = x_all["x_all"][x_all["parity_all"] * sgn > 0]; xm = blocks[lab]["x"]
                r["A3"][lab] = bd.interlacing(xr, xm, WIN[0], WIN[1])
        # B.4 : δ et <Ψ_in|ΔV|Ψ_out> pour les états en évidence
        r["B4"] = {}
        Mf = {p: Ms[p][:n, :, :n, :].reshape(n * 81, n * 81) for p in ("L", "NL")}
        for lab, s in hl.items():
            cf = c128[s, :n, :].reshape(-1); cc = float(np.vdot(cf, cf).real)
            A_in = sp.bloch_superposition_grid(c128[s, :n, :], b["C"][:n], b["nG"], b["flat"], ng); psi_in = sp.grid_to_real(A_in, ng, Om, NT)
            VL_in_in = dv.expect_local(psi_in, psi_in, dV, Om); VL_in_psi = dv.expect_local(psi_in, psi_hl[s], dV, Om)
            B_in = proj.projections(A_in[flat_full][None])[0]; VNL_in_in = proj.expect(B_in, B_in); VNL_in_psi = proj.expect(B_in, B_hl[s])
            cML = complex(np.vdot(cf, Mf["L"] @ cf)) / N_CELLS; cMNL = complex(np.vdot(cf, Mf["NL"] @ cf)) / N_CELLS
            r["B4"][lab] = dict(delta=1 - cc, VL_in_out=(VL_in_psi - VL_in_in).real, VNL_in_out=(VNL_in_psi - VNL_in_in).real, VL_in_in=VL_in_in.real, VNL_in_in=VNL_in_in.real,
                                cML=cML.real, cMNL=cMNL.real, dL=abs(VL_in_in - cML), dNL=abs(VNL_in_in - cMNL))
            del A_in, psi_in
        # diagnostique : M^L x 81 + M^NL à ce n (spectre et w2 seulement)
        Hc = sf.bloch_folded_hamiltonian(eps_n, N_CELLS * Ms["L"][:n, :, :n, :] + Ms["NL"][:n, :, :n, :], N_CELLS); blc, cplc, _ = r4.spectrum_blocks(Hc, par_n, b["E_D"])
        Vc = {}; recc = record_variant(Vc, f"n{n}_Lx81", blc, window_weights(blc, b["C"][:n], b["nG"], b["flat"], n, ng, mask2), b["E_D"], Q["thr"], dict(coupling_even_odd=cplc))
        oddc = recc["blocks"]["odd"]["localized"]; evenc = recc["blocks"]["even"]["localized"]
        r["Lx81"] = dict(blocks=recc["blocks"], pi_state=max(oddc, key=lambda t: t[1]) if oddc else None, sigma_localized=evenc, sigma_near_0p101=[t for t in evenc if abs(t[0] - 0.101) < 0.3],
                         n_window={lab: recc["blocks"][lab]["n_window"] for lab in ("even", "odd")})
        del Hc, blc
        per_n[str(n)] = r; prev = r
        log(f"[B.3] n={n}: dim {H.shape[0]}, fenêtre σ {r['n_window']['even']} / π {r['n_window']['odd']}, π localisé {pi and (round(pi[0],3), round(pi[1],3))}, σ localisés {[(round(x,3), round(w,3)) for x, w in even]}, "
            f"σ près de +0,101 {[(round(x,3), round(w,3)) for x, w in sig_near]} ; médiane|Δε| vs QE σ {r['median_dx_vs_qe']['even']} π {r['median_dx_vs_qe']['odd']} ; "
            f"δ(π_320) {r['A1'][hl['pi_320']]['delta']:.4f}, δ(σ_324) {r['A1'][hl['sigma_324']]['delta']:.4f} ; B4 π_320 <in|ΔV|out> L {r['B4']['pi_320']['VL_in_out']:+.4f} NL {r['B4']['pi_320']['VNL_in_out']:+.4f} ({time.time()-t0:.0f} s)")
        del H, Hh, blocks
    out["per_n"] = per_n
    out["first_n_sigma_near_0p101"] = next((n for n in NLIST if per_n[str(n)]["sigma_near_0p101"]), None)
    out["first_n_sigma_near_0p101_Lx81"] = next((n for n in NLIST if per_n[str(n)]["Lx81"]["sigma_near_0p101"]), None)
    out["pi_vs_n_Lx81"] = [(n, per_n[str(n)]["Lx81"]["pi_state"][0] if per_n[str(n)]["Lx81"]["pi_state"] else None) for n in NLIST]
    out["pi_vs_n"] = [(n, per_n[str(n)]["pi_state"][0] if per_n[str(n)]["pi_state"] else None, per_n[str(n)]["pi_state"][1] if per_n[str(n)]["pi_state"] else None) for n in NLIST]
    out["timing_s"] = time.time() - t_start
    save_json(os.path.join(d, "b3_results.json"), out)
    write_tables_b(out, Q, hl); figs_b(out, Q, hl)


def write_tables_b(out, Q, hl):
    L = []; P = out["per_n"]
    b1 = json.load(open(os.path.join(WORK, "b", "b1_results.json"))) if os.path.exists(os.path.join(WORK, "b", "b1_results.json")) else None
    b2 = json.load(open(os.path.join(WORK, "b", "b2_results.json"))) if os.path.exists(os.path.join(WORK, "b", "b2_results.json")) else None
    if b1:
        L.append("#### B.1 — vérifications du `.save` à 128 bandes\n")
        L.append(f"- Liste des 81 k identique à la maille grossière : {b1['k_list']['same_order']} (max\\|Δk\\| {b1['k_list']['max_abs_diff']:.1e}) ; 16 premières bandes : max\\|Δε\\| = {b1['eps16_vs_new']['max_abs_eV']:.2e} eV (seuil 1e-6) → {'OK' if b1['eps16_vs_new']['ok'] else 'ÉCART'} ; E_D(K) = {b1['E_D']:.6f} eV.")
        L.append(f"- Bande 128 (ε − E_D) : Γ {b1['band128']['Gamma']['x']:+.3f}, K {b1['band128']['K']['x']:+.3f}, M {b1['band128']['M']['x']:+.3f} eV ; étendue de la base : [{b1['band_range']['min_x']:+.2f}, {b1['band_range']['max_x']:+.2f}] eV.")
        L.append(f"- Parité des 128 × 81 états : max\\|1 − \\|⟨σ_h⟩\\|\\| = {b1['parity']['max_dev']:.1e}, {b1['parity']['n_even']} pairs, {b1['parity']['n_odd']} impairs, {b1['parity']['n_undetermined']} indéterminés ; pairs par bande : {b1['parity']['even_per_band']}.")
        L.append(f"- nscf : mur {b1['nscf_out']['wall']} s, « not converged » : {b1['nscf_out']['n_not_converged']} ; `.save` {b1['save_size_GB']:.2f} Go.\n")
    if b2:
        L.append("#### B.2 — M à 128 bandes (chaîne grossière de production, `.save` neuf)\n")
        L.append("| partie | max\\|M\\| (Ha) | hermiticité (Ha) | bloc 16 × 16 contre M grossier de juin : VS écart rel. max | valeurs propres de H(a1) écart max (eV) | diagonale écart max (Ha) | moyenne diag. neuf / juin (Ha) |")
        L.append("|---|---|---|---|---|---|---|")
        for tag in ("L", "NL", "tot"):
            c = b2["block16_vs_june"][tag]
            L.append(f"| {tag} | {b2[f'max_{tag}_Ha']:.4e} | {b2[f'herm_{tag}_Ha']:.1e} | {c['sv_rel_max']:.2e} | {c['eig_H_a1_max_dev_eV']:.2e} | {c['diag_max_dev_Ha']:.2e} | {c['diag_mean_Ha']:.6f} / {c['diag_mean_ref_Ha']:.6f} |")
        L.append(f"\nε(16 bandes, grossier) contre ε(neuf)[:16] : {b2['eps16_vs_new_max_eV']:.2e} eV ; fichiers de {b2['files']['GB_each']:.2f} Go chacun (hors dépôt).\n")
    L.append("#### B.3 — (a1) à n bandes (sous-blocs du M à 128 bandes ; ε − E_D ; w₂ ; seuil %.4f)\n" % Q["thr"])
    L.append("| n | dim | fenêtre σ / π | état π localisé (ε − E_D ; w₂) | \\|Δε_π\\| vs n précédent (meV) | σ localisés (ε − E_D ; w₂) | σ à moins de 0,3 eV de +0,101 | médiane \\|Δε\\| vs QE σ / π (eV) | médiane \\|Δε\\| vs n précédent σ / π (eV) |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for n in out["n_list"]:
        r = P[str(n)]; pi = r["pi_state"]
        L.append(f"| {n} | {r['dim']} | {r['n_window']['even']} / {r['n_window']['odd']} | " + (f"{pi[0]:+.3f} ; {pi[1]:.3f}" if pi else "aucun") + f" | {_f(r.get('dpi_vs_prev') and r['dpi_vs_prev']*1e3, 1)} | "
                 + (", ".join(f"{x:+.3f} ; {w:.3f}" for x, w in r["sigma_localized"]) or "aucun") + " | " + (", ".join(f"{x:+.3f} ; {w:.3f}" for x, w in r["sigma_near_0p101"]) or "aucun")
                 + f" | {_f(r['median_dx_vs_qe']['even'])} / {_f(r['median_dx_vs_qe']['odd'])} | {_f(r.get('median_dx_vs_prev', {}).get('even'))} / {_f(r.get('median_dx_vs_prev', {}).get('odd'))} |")
    L.append(f"\nPremier n avec un état pair localisé à moins de 0,3 eV de +0,101 : {out['first_n_sigma_near_0p101']}.\n")
    L.append("\n*Diagnostique (hors prompt) : M^L × 81 + M^NL, sous-blocs à n bandes*\n")
    L.append("| n | fenêtre σ / π | état π localisé (ε − E_D ; w₂) | σ localisés (ε − E_D ; w₂) | σ à moins de 0,3 eV de +0,101 |")
    L.append("|---|---|---|---|---|")
    for n in out["n_list"]:
        q = P[str(n)]["Lx81"]; pi = q["pi_state"]
        L.append(f"| {n} | {q['n_window']['even']} / {q['n_window']['odd']} | " + (f"{pi[0]:+.3f} ; {pi[1]:.3f}" if pi else "aucun") + " | " + (", ".join(f"{x:+.3f} ; {w:.3f}" for x, w in q["sigma_localized"]) or "aucun") + " | " + (", ".join(f"{x:+.3f} ; {w:.3f}" for x, w in q["sigma_near_0p101"]) or "aucun") + " |")
    L.append(f"\nPremier n (M^L × 81) avec un état pair localisé à moins de 0,3 eV de +0,101 : {out['first_n_sigma_near_0p101_Lx81']}.\n")
    L.append("\n#### B.3 — A.1 à chaque n pour les états en évidence (δ, R − ε_QE, résidu, poids principal)\n")
    L.append("| état | n | δ | R − ε_QE (meV) | résidu (eV) | 3 plus grands poids (ε_j − E_D ; poids) |")
    L.append("|---|---|---|---|---|---|")
    for lab, s in hl.items():
        for n in out["n_list"]:
            q = P[str(n)]["A1"][s]
            L.append(f"| {lab} ({q['x_qe']:+.3f}) | {n} | {q['delta']:.4f} | {q['rayleigh_minus_eps']*1e3:+.1f} | {q['residual']:.4f} | " + ", ".join(f"({e:+.3f} ; {w:.3f})" for e, w in q["top_x"][:3]) + " |")
    L.append("\n#### B.3 — A.3 à chaque n : comptages sous E_D − 3 et min_k[λ_k(n) − λ_{k−s}(QE)] (fenêtre ; tout)\n")
    L.append("| n | parité | QE sous −3 | modèle sous −3 | fenêtre QE / modèle | s = 0 | s = 1 | s = 2 |")
    L.append("|---|---|---|---|---|---|---|---|")
    for n in out["n_list"]:
        r = P[str(n)]
        if "A3" not in r:
            continue
        for lab in ("even", "odd"):
            q = r["A3"][lab]; sh = q["shifts"]
            L.append(f"| {n} | {'σ' if lab=='even' else 'π'} | {q['n_ref_below_lo']} | {q['n_model_below_lo']} | {len(q['ref_window'])} / {len(q['model_window'])} | " + " | ".join(f"{sh.get(s, sh.get(str(s))).get('min_window', float('nan')):+.3f} ; {sh.get(s, sh.get(str(s)))['min_all']:+.3f}" for s in (0, 1, 2)) + " |")
    L.append("\n#### B.4 — états en évidence : δ(n) et ⟨Ψ_in\\|ΔV\\|Ψ_out⟩ (eV)\n")
    L.append("| état | n | δ | ⟨in\\|ΔV^L\\|in⟩ (direct ; c†M^Lc/81 ; Δ) | ⟨in\\|ΔV^NL\\|in⟩ (direct ; c†M^NLc/81 ; Δ) | ⟨in\\|ΔV^L\\|out⟩ | ⟨in\\|ΔV^NL\\|out⟩ |")
    L.append("|---|---|---|---|---|---|---|")
    for lab in hl:
        for n in out["n_list"]:
            q = P[str(n)]["B4"][lab]
            L.append(f"| {lab} | {n} | {q['delta']:.4f} | {q['VL_in_in']:+.4f} ; {q['cML']:+.4f} ; {q['dL']:.1e} | {q['VNL_in_in']:+.4f} ; {q['cMNL']:+.4f} ; {q['dNL']:.1e} | {q['VL_in_out']:+.4f} | {q['VNL_in_out']:+.4f} |")
    L.append("\nTous les états de la fenêtre par n et parité (ε − E_D ; w₂) :\n")
    for n in out["n_list"]:
        for lab in ("even", "odd"):
            bl = P[str(n)]["blocks"][lab]
            L.append(f"- n = {n}, {lab} : " + ", ".join(f"{x:+.3f} ({w:.3f})" for x, w in zip(bl["x"], bl["w2"])))
            bl = P[str(n)]["Lx81"]["blocks"][lab]
            L.append(f"- n = {n}, {lab}, M^L × 81 : " + ", ".join(f"{x:+.3f} ({w:.3f})" for x, w in zip(bl["x"], bl["w2"])))
    with open(os.path.join(WORK, "b", "B_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log("[B] tables b/B_tables.md")


def figs_b(out, Q, hl):
    plt, pal = fig_style(); P = out["per_n"]; ns_ = out["n_list"]
    a = json.load(open(os.path.join(WORK, "a", "a_results.json"))) if os.path.exists(os.path.join(WORK, "a", "a_results.json")) else None
    fig, ax = plt.subplots(figsize=(4.2, 3.2))
    for i, (lab, s) in enumerate(hl.items()):
        dl = [P[str(n)]["A1"][s]["delta"] for n in ns_]; col = pal.CYCLE[i % len(pal.CYCLE)]
        ax.plot(ns_, dl, "o-", ms=3, lw=0.8, color=col, label=f"{lab} ({Q['x'][s]:+.2f} eV)")
        if a:
            ax.plot([16, 20], [a["A1"]["16"][s]["delta"], a["A1"]["20"][s]["delta"]], "x", ms=4, color=col)
    ax.set_yscale("log"); ax.set_xlabel("n bandes (× 81 k)"); ax.set_ylabel(r"$\delta = 1 - \sum|c|^2$"); ax.legend(fontsize=5); ax.set_title("croix : bases 16 (grossier) et 20 (dense) de A", fontsize=7)
    fig.tight_layout(); savefig(fig, "delta_vs_n")
    fig, axs = plt.subplots(1, 2, figsize=(8.5, 4.2), sharey=True); d4 = json.load(open(os.path.join(R4DIR, "d4", "d4_results.json")))
    cols = ["QE"] + [f"n{n}" for n in ns_] + [f"n{n}×" for n in ns_]
    for ax, lab, ttl in zip(axs, ("odd", "even"), ("(a) impairs (π)", "(b) pairs (σ)")):
        b = d4["variants"]["QE_D1"]["blocks"][lab]; ax.scatter(np.zeros(len(b["x"])), b["x"], s=2 + 60 * np.array(b["w2"]), color=pal.ORANGE, alpha=0.7, lw=0)
        for i, n in enumerate(ns_):
            bl = P[str(n)]["blocks"][lab]; x = np.array(bl["x"]); w = np.array(bl["w2"])
            ax.scatter(np.full(len(x), i + 1), x, s=2 + 60 * w, color=pal.NAVY, alpha=0.7, lw=0)
            bl = P[str(n)]["Lx81"]["blocks"][lab]; x = np.array(bl["x"]); w = np.array(bl["w2"])
            ax.scatter(np.full(len(x), i + 1 + len(ns_)), x, s=2 + 60 * w, color=pal.GREEN, alpha=0.7, lw=0)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=60, fontsize=6); ax.set_title(ttl); ax.set_ylim(-3, 1); ax.axhline(0, color=pal.MUTED, lw=0.5)
    axs[0].set_ylabel(r"$\varepsilon - E_D$ (eV)"); axs[0].set_title("(a) impairs (π) ; vert : M^L × 81", fontsize=8); fig.tight_layout(); savefig(fig, "ladder_bands")


# ----------------------------------------------------------------------------------------------- tables et figures depuis les json
def cmd_atables(a):
    out = json.load(open(os.path.join(WORK, "a", "a_results.json"))); Q = dict(thr=out["qe_states"]["thr_w2"])
    write_tables_a(out, Q); figs_a(out, Q)


def cmd_ctables(a):
    out = json.load(open(os.path.join(WORK, "c", "c_results.json"))); write_tables_c(out); figs_c(out)


def cmd_btables(a):
    out = json.load(open(os.path.join(WORK, "b", "b3_results.json"))); z = np.load(os.path.join(WORK, "b", "b3_overlaps.npz"))
    a_out = json.load(open(os.path.join(WORK, "a", "a_results.json")))
    Q = dict(thr=a_out["qe_states"]["thr_w2"], x=z["x_qe"], bands1=z["bands"]); hl = {k: int(v) for k, v in a_out["qe_states"]["highlighted"].items()}
    write_tables_b(out, Q, hl); figs_b(out, Q, hl)


# ----------------------------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="R5 driver"); p.add_argument("cmd", choices=["a", "c", "b1", "b2ml", "b2nl", "b2sum", "b3", "atables", "ctables", "btables"])
    a = p.parse_args()
    {"a": cmd_a, "c": cmd_c, "b1": cmd_b1, "b2ml": cmd_b2ml, "b2nl": cmd_b2nl, "b2sum": cmd_b2sum, "b3": cmd_b3, "atables": cmd_atables, "ctables": cmd_ctables, "btables": cmd_btables}[a.cmd](a)


if __name__ == "__main__":
    main()
