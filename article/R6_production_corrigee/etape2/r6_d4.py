#!/usr/bin/env python
"""
r6_d4.py -- R6 étape 2 (GO 2) : validation de la chaîne Wannier corrigée (M2) contre la super-cellule 9x9 (escalier D4 de R4).
Sous-commandes :
  prep : M2 dense 9x9 (tot, L, NL) -> M_W(R,R') (rotation V^dag M V, double TF, recentrage) -> V_loc (R_cut 3) ; caches cache/ ;
         comparaison à la rotation v1 de R4 (V_loc(M2) = 81 V_L(v1) + V_NL(v1)) ; éléments p_z-p_z du site de la lacune.
  d4   : escalier (a1 tot/L/NL, a2_16, a3_20, b_5wf, c_all, c_3) avec M2 ; portes 1 et 2 de R4 ; états localisés par parité contre QE (D1) ;
         attendu R5 (variante M^L x 81) : (a1) -0,727, (a3) -0,734 ; marches ; décalage rigide ; figure.
  d3   : critère de pôle avec V_loc(M2), alpha = 1, R_cut 3, 300² (et 600² sur [-3, +1]) : min |det| / min |lambda| par bloc (complet, pi, sigma),
         vecteur propre, racines, pic de -Im Tbar(K) ; contre v1 (R4 D3).
  b3   : (a1) à n = 16 ... 128 bandes avec le M2 nb128 (états localisés, delta des états en évidence) ; contre R5 B.3 (variante M^L x 81).
Réutilise r5_driver (bases, états QE, poids w2) et r4_driver (setup, spectrum_blocks, g0 en cache) ; sorties dans d4/, d3/, b3/, fig/ ;
journal r6_log.txt. Données de production en lecture seule ; M2 lus dans results/M2/ (sidecar v2 exigé).
"""
import os
import sys
import json
import time
import argparse

import numpy as np

WORK = os.path.dirname(os.path.abspath(__file__)); GQ = os.path.dirname(os.path.dirname(WORK))
PROJ = os.environ.get("GRAPHENE_RAMAN") or os.path.join(os.path.dirname(os.path.dirname(GQ)), "graphene-raman")
R5DIR = os.path.join(GQ, "defects", "R5_base_vs_M"); R4DIR = os.path.join(GQ, "defects", "R4_quasi_lie")
sys.path.insert(0, os.path.join(PROJ, "src")); sys.path.insert(0, R4DIR); sys.path.insert(0, R5DIR)
import r5_driver as r5  # noqa: E402  (chdir PROJ)
r4 = r5.r4
from electron_defect_interaction.io import qe_io, matrix_io  # noqa: E402
from electron_defect_interaction.io import qe_gamma_io as qg  # noqa: E402
from electron_defect_interaction.config import HA2EV  # noqa: E402
from electron_defect_interaction.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr  # noqa: E402
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk  # noqa: E402
from electron_defect_interaction.defects.many_body import local_tmatrix as lt  # noqa: E402
from electron_defect_interaction.defects.many_body import pole_criterion as pc  # noqa: E402
from electron_defect_interaction.wannier import supercell_fold as sf  # noqa: E402
from electron_defect_interaction.defects import alignment as al  # noqa: E402
import r5_sc_projection as sp  # noqa: E402

M2 = os.path.join(PROJ, "results", "M2")
N_SC, N_CELLS, NW = r4.N_SC, r4.N_CELLS, r4.NW
WF_SIGMA, WF_PI, WF_PZ_A, WF_PZ_B, NN_CELLS, K_RED = r4.WF_SIGMA, r4.WF_PI, r4.WF_PZ_A, r4.WF_PZ_B, r4.NN_CELLS, r4.K_RED
WIN = r5.WIN; NT = r5.NT; NLIST = r5.NLIST
FILES_DENSE = dict(tot=os.path.join(M2, "M_dense_9x9.npy"), L=os.path.join(M2, "M_L_dense_9x9.npy"), NL=os.path.join(M2, "M_NL_dense_9x9.npy"))
FILES_COARSE = dict(tot=os.path.join(M2, "M_ed_9x9.npy"), L=os.path.join(M2, "M_L_9x9.npy"), NL=os.path.join(M2, "M_NL_9x9.npy"))
FILES_128 = dict(tot=os.path.join(M2, "M_ed_9x9_nb128.npy"), L=os.path.join(M2, "M_L_9x9_nb128.npy"), NL=os.path.join(M2, "M_NL_9x9_nb128.npy"))


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"; print(line, flush=True)
    with open(os.path.join(WORK, "r6_log.txt"), "a") as f:
        f.write(line + "\n")


r5.log = log; r4.log = log


def ensure(sub):
    d = os.path.join(WORK, sub); os.makedirs(d, exist_ok=True); return d


def save_json(path, d):
    with open(path, "w") as f:
        json.dump(r4.jsonable(d), f, indent=1, ensure_ascii=False)
    log(f"saved {os.path.relpath(path, WORK)}")


def load_v2(path, units=matrix_io.EV):
    return matrix_io.load_M_checked(path, require_bloch_norm=matrix_io.UNIT_CELL, units=units, require_normalization=matrix_io.M_NORM_V2)


def dense_restricted(part, idx81):
    p = FILES_DENSE[part]; meta = matrix_io.read_manifest(p)
    if meta is None or meta.get("bloch_norm") != matrix_io.UNIT_CELL or meta.get("units") != matrix_io.HARTREE or not str(meta.get("M_normalization", "")).startswith("v2"):
        raise ValueError(f"{p}: sidecar v2 requis ({meta})")
    Mm = np.load(p, mmap_mode="r"); nb = Mm.shape[0]
    M = np.array(Mm[np.ix_(np.arange(nb), idx81, np.arange(nb), idx81)]) * HA2EV; del Mm
    return M


def fig_style():
    return r5.fig_style()


def savefig(fig, name):
    d = ensure("fig"); fig.savefig(os.path.join(d, name + ".pdf")); fig.savefig(os.path.join(d, name + ".png"), dpi=200); log(f"[fig] {name}.pdf/.png")


# ----------------------------------------------------------------------------------------------- prep : M_W(M2), V_loc(M2)
def cmd_prep(a):
    S = r4.setup(verbose=False); out = {}; ensure("cache"); d = ensure("d4")
    Hwr, Rw, nd, k27, MP, U, Ud, dp = S["Hwr"], S["Rw"], S["nd"], S["k27"], S["MP"], S["U"], S["Ud"], S["dp"]
    Mwr = {}; Rn = Rd = None
    for tag in ("tot", "L", "NL"):
        t0 = time.time(); M = load_v2(FILES_DENSE[tag]); Mwk = Mbk_to_Mwk(M, U, Ud); del M
        Mw, R = Mwk_to_Mwr(Mwk, k27, MP); del Mwk
        if tag == "tot":
            Rn, Rd = lt.recenter_mwr(Mw, R, MP); dist, wt = lt.mwr_locality(Mw, Rn)
            out["recenter"] = dict(R_d=Rd.tolist(), locality_first=[(float(x), float(y)) for x, y in zip(dist[:10], wt[:10])])
        Mwr[tag] = Mw; log(f"[prep] rotation M2 {tag} : {time.time()-t0:.0f} s ; max|M_W| {np.abs(Mw).max():.4f} eV")
    lin = float(np.linalg.norm(Mwr["tot"] - Mwr["L"] - Mwr["NL"]) / np.linalg.norm(Mwr["tot"])); out["linearity"] = lin
    np.savez(os.path.join(WORK, "cache", "Mwr_M2_9x9.npz"), Mwr_tot=Mwr["tot"], Mwr_L=Mwr["L"], Mwr_NL=Mwr["NL"], R=Rn, R_d=Rd, MP=np.array(MP))
    rc = S["cfg"]["R_cut"]; Rloc = Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9]; nL = len(Rloc)
    V = {}; herm = {}
    for tag in ("tot", "L", "NL"):
        V[tag], herm[tag] = lt.extract_V_loc(Mwr[tag], Rn, Rloc)
    ip = pc.block_indices(nL, NW, WF_PI); isg = pc.block_indices(nL, NW, WF_SIGMA)
    np.savez(os.path.join(WORK, "cache", "Vloc_M2_9x9.npz"), V_tot=V["tot"], V_L=V["L"], V_NL=V["NL"], Rloc=Rloc, idx_pi=ip, idx_sigma=isg)
    # comparaison à R4 (v1) : même R_loc, V_loc(M2) = 81 V_L(v1) + V_NL(v1) ; M_W idem
    C1 = r4.load_cache_vloc(); Z1 = r4.load_cache_mwr()
    same_R = bool(np.array_equal(C1["Rloc"], Rloc)) and bool(np.array_equal(Z1["R"], Rn)) and bool(np.array_equal(Z1["R_d"], Rd))
    dV = float(np.abs(V["tot"] - (N_CELLS * C1["V_L"] + C1["V_NL"])).max()); dVL = float(np.abs(V["L"] - N_CELLS * C1["V_L"]).max()); dVN = float(np.abs(V["NL"] - C1["V_NL"]).max())
    dMw = float(np.abs(Mwr["tot"] - (N_CELLS * Z1["Mwr_L"] + Z1["Mwr_NL"])).max())
    i0 = r4.cell_index(Rn, (0, 0, 0)); iNN = [r4.cell_index(Rn, R) for R in NN_CELLS]; j0 = r4.cell_index(Rloc, (0, 0, 0))
    el = {}
    for tag in ("tot", "L", "NL"):
        Mw = Mwr[tag]
        el[tag] = dict(pz_vac_onsite=float(Mw[WF_PZ_A, i0, WF_PZ_A, i0].real), pz_B_onsite_R0=float(Mw[WF_PZ_B, i0, WF_PZ_B, i0].real),
                       pz_vac_pz_nn=[float(abs(Mw[WF_PZ_A, i0, WF_PZ_B, i])) for i in iNN], pz_nn_diag=[float(Mw[WF_PZ_B, i, WF_PZ_B, i].real) for i in iNN],
                       sp2_diag=[float(Mw[w, i0, w, i0].real) for w in WF_SIGMA], onsite_norm=float(np.linalg.norm(Mw[:, i0, :, i0])),
                       Vloc_max_abs=float(np.abs(V[tag]).max()), Vloc_herm_residual=float(herm[tag]),
                       Vloc_offblock_sigma_pi=float(max(np.abs(V[tag][np.ix_(isg, ip)]).max(), np.abs(V[tag][np.ix_(ip, isg)]).max())))
    v1 = {tag: dict(pz_vac_onsite=float(Z1[f"Mwr_{tag}"][WF_PZ_A, i0, WF_PZ_A, i0].real), onsite_norm=float(np.linalg.norm(Z1[f"Mwr_{tag}"][:, i0, :, i0]))) for tag in ("tot", "L", "NL")}
    out.update(R_cut=rc, nL=nL, dim=nL * NW, same_R_as_R4=same_R, Vloc_M2_vs_81VL1_plus_VNL1_max_abs_eV=dV, VL_M2_vs_81VL1=dVL, VNL_M2_vs_VNL1=dVN, Mwr_M2_vs_combo_v1=dMw,
               elements_eV=el, v1_elements_eV=v1, Vloc_M2_pz_vac=float(V["tot"][j0 * NW + WF_PZ_A, j0 * NW + WF_PZ_A].real), Vloc_v1_pz_vac=float(C1["V_tot"][j0 * NW + WF_PZ_A, j0 * NW + WF_PZ_A].real))
    log(f"[prep] R_d {Rd.tolist()}, R_loc {nL} mailles (dim {nL*NW}), mêmes R que R4 : {same_R} ; linéarité M_W {lin:.1e} ; "
        f"V_loc(M2) - (81 V_L + V_NL)(v1) : {dV:.2e} eV (L {dVL:.1e}, NL {dVN:.1e}) ; M_W idem {dMw:.2e}")
    log(f"[prep] p_z-p_z lacune (eV) : M2 tot {el['tot']['pz_vac_onsite']:+.4f} (L {el['L']['pz_vac_onsite']:+.4f}, NL {el['NL']['pz_vac_onsite']:+.4f}) ; v1 tot {v1['tot']['pz_vac_onsite']:+.4f} (L {v1['L']['pz_vac_onsite']:+.4f}) ; "
        f"||M_W(0,0)|| M2 {el['tot']['onsite_norm']:.4f} / v1 {v1['tot']['onsite_norm']:.4f} ; hors bloc sigma-pi de V_loc(M2) {el['tot']['Vloc_offblock_sigma_pi']:.2e} eV")
    save_json(os.path.join(d, "prep_results.json"), out)


# ----------------------------------------------------------------------------------------------- d4 : escalier
def cmd_d4(a):
    d = ensure("d4"); out = dict(gates={}, variants={}); t_start = time.time()
    B = r5.load_bases(with_dense=True); ng = B["ng"]; A_A = B["A_A"]
    Q = r5.qe_states_9x9(); hl = r5.highlighted(Q); mask2 = qg.inplane_disc_mask(ng, A_A, Q["s_vac"], r4.R_W2)
    b16 = B["bases"][16]; b20 = B["bases"][20]; idx81 = B["idx81"]; k81 = B["k81"]; nk = 81; nw = NW
    Ms = {16: {p: load_v2(FILES_COARSE[p]) for p in ("tot", "L", "NL")}, 20: {p: dense_restricted(p, idx81) for p in ("tot", "L", "NL")}}
    out["M2_linearity_eV"] = {str(nb): float(np.abs(Ms[nb]["tot"] - Ms[nb]["L"] - Ms[nb]["NL"]).max()) for nb in Ms}
    out["E_D"] = dict(uc16_at_K=b16["E_D"], uc20dense_at_K=b20["E_D"], SC_P=Q["E_D_SC"])
    V = out["variants"]
    # (a1) 16 bandes, tot / L / NL
    for part in ("tot", "L", "NL"):
        H = sf.bloch_folded_hamiltonian(b16["eps"], Ms[16][part], N_CELLS); bl, cpl, herm = r4.spectrum_blocks(H, b16["par"].reshape(-1), b16["E_D"])
        W = r5.window_weights(bl, b16["C"], b16["nG"], b16["flat"], 16, ng, mask2); r5.record_variant(V, f"a1_{part}", bl, W, b16["E_D"], Q["thr"], dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
        del H
    # (a2) 16 bandes du dense restreint, (a3) 20 bandes
    H = sf.bloch_folded_hamiltonian(b20["eps"][:16], Ms[20]["tot"][:16, :, :16, :], N_CELLS); bl, cpl, herm = r4.spectrum_blocks(H, b20["par"][:16].reshape(-1), b20["E_D"])
    W = r5.window_weights(bl, b20["C"][:16], b20["nG"], b20["flat"], 16, ng, mask2); r5.record_variant(V, "a2_16", bl, W, b20["E_D"], Q["thr"], dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
    H = sf.bloch_folded_hamiltonian(b20["eps"], Ms[20]["tot"], N_CELLS); bl, cpl, herm = r4.spectrum_blocks(H, b20["par"].reshape(-1), b20["E_D"])
    W = r5.window_weights(bl, b20["C"], b20["nG"], b20["flat"], 20, ng, mask2); r5.record_variant(V, "a3_20", bl, W, b20["E_D"], Q["thr"], dict(coupling_even_odd=cpl, herm=herm, dim=int(H.shape[0])))
    spectra = {"a1_tot": None, "a2_16": None, "a3_20": None}
    for tag, (eps, M, nb) in {"a1_tot": (b16["eps"], Ms[16]["tot"], 16), "a2_16": (b20["eps"][:16], Ms[20]["tot"][:16, :, :16, :], 16), "a3_20": (b20["eps"], Ms[20]["tot"], 20)}.items():
        Hh = sf.bloch_folded_hamiltonian(eps, M, N_CELLS); spectra[tag] = np.linalg.eigvalsh(0.5 * (Hh + Hh.conj().T)) - (b16["E_D"] if tag == "a1_tot" else b20["E_D"]); del Hh
    # (b) 5 WF : V^dag eps V (+) M_W(k,k')/81
    S = r4.setup(verbose=False); Hwr, Rw, nd, U, Ud = S["Hwr"], S["Rw"], S["nd"], S["U"], S["Ud"]
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
    r5.record_variant(V, "b_5wf", blocks_b, w_from_k(blocks_b), ED_b, Q["thr"], dict(coupling_even_odd=cpl, herm=herm, dim=int(Hb.shape[0]), E_D_def="paire centrale de V^dag eps V à K"))
    e_b = np.linalg.eigvalsh(0.5 * (Hb + Hb.conj().T)); spectra["b_5wf"] = e_b - ED_b
    # (c) repliement de M_W(M2) et de H(R) ; portes 1 et 2 de R4
    Z = np.load(os.path.join(WORK, "cache", "Mwr_M2_9x9.npz")); Mwr = Z["Mwr_tot"]; Rn = Z["R"]; R_d = Z["R_d"]; MP = tuple(int(x) for x in Z["MP"])
    C2 = np.load(os.path.join(WORK, "cache", "Vloc_M2_9x9.npz")); Rloc = C2["Rloc"]
    dev1, herm1, ok1 = sf.wannier_gate(Hwr, Rw, nd, N_SC); out["gates"]["gate1_fold_HR"] = dict(max_dev_eV=float(dev1), herm=float(herm1), ok=bool(ok1), threshold=1e-8)
    HS = sf.fold_hwr_to_supercell(Hwr, Rw, nd, N_SC); Hk1r = sf.kblocks_to_rbasis(Hk, k81, N_SC)
    Mall = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC); M3 = sf.fold_mwr_to_supercell(Mwr, Rn, N_SC, R_local=Rloc)
    e_c1 = np.linalg.eigvalsh(0.5 * ((Hk1r + Mall) + (Hk1r + Mall).conj().T)); dev2 = float(np.abs(e_b - e_c1).max())
    k27 = S["k27"]; Mfull = load_v2(FILES_DENSE["tot"]); Mwk_full = Mbk_to_Mwk(Mfull, U, Ud); del Mfull
    Mwr_snap, R_snap = Mwk_to_Mwr(Mwk_full, np.rint(27 * k27) / 27.0, MP); del Mwk_full
    Rn_snap, Rd_snap = lt.recenter_mwr(Mwr_snap, R_snap, MP); Mall_snap = sf.fold_mwr_to_supercell(Mwr_snap, Rn_snap, N_SC)
    e_c1s = np.linalg.eigvalsh(0.5 * ((Hk1r + Mall_snap) + (Hk1r + Mall_snap).conj().T)); dev2s = float(np.abs(e_b - e_c1s).max())
    Vsnap, _ = lt.extract_V_loc(Mwr_snap, Rn_snap, Rloc); j0 = r4.cell_index(Rloc, (0, 0, 0)) * nw + WF_PZ_A
    out["gates"]["gate2_fold_Mwr"] = dict(max_dev_eV_production_k=dev2, max_dev_eV_snapped_k=dev2s, ok=bool(dev2s < 1e-9), literal_ok=bool(dev2 < 1e-9), threshold=1e-9, R_d_snap=Rd_snap.tolist(),
                                          Vloc_prod_vs_snap_max_abs_eV=float(np.abs(C2["V_tot"] - Vsnap).max()), onsite_pz_prod=float(C2["V_tot"][j0, j0].real), onsite_pz_snap=float(Vsnap[j0, j0].real),
                                          Mall_prod_vs_snap_max_abs_eV=float(np.abs(Mall - Mall_snap).max()))
    del Mwr_snap, Mall_snap, Vsnap
    Hwk81, _, _ = Hwr_to_Hwk(Hwr, Rw, k81, ndegen=nd)
    out["gates"]["one_body_residual"] = dict(max_abs_Hwk_minus_VepsV_eV=float(np.abs(Hwk81 - Hk).max()))
    log(f"[D4] porte 1 (H(R) replié) : {dev1:.2e} eV -> {'OK' if ok1 else 'ÉCHEC'} ; porte 2 (repliement de M_W(M2)) : {dev2:.2e} eV (k du XML), {dev2s:.2e} eV (k = m/27) -> {'OK' if dev2s < 1e-9 else 'ÉCHEC'} ; "
        f"V_loc production vs k exacts {out['gates']['gate2_fold_Mwr']['Vloc_prod_vs_snap_max_abs_eV']:.2e} eV ; résidu H(R) vs V^+epsV {out['gates']['one_body_residual']['max_abs_Hwk_minus_VepsV_eV']:.2e} eV")
    _, EK, _ = Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); E_DW = float(0.5 * (EK[0, 3] + EK[0, 4])); out["E_D"]["wannier_at_K"] = E_DW
    phase_Rd = np.exp(-2j * np.pi * (k81 @ np.asarray(R_d, float)))

    def w_from_r(bl):
        Wr = {}
        for lab in bl:
            x = bl[lab]["x"]; sel = np.where((x >= WIN[0]) & (x <= WIN[1]))[0]; w = np.zeros(len(sel)); wsite = np.zeros(len(sel))
            for jj, j in enumerate(sel):
                cr = np.zeros(N_CELLS * nw, complex); cr[bl[lab]["idx"]] = bl[lab]["v"][:, j]
                ck = sf.rvec_to_kvec(cr, k81, N_SC, nw).reshape(nk, nw) * phase_Rd[:, None]
                dnk = np.einsum("kbw,kw->bk", V81, ck)
                w[jj] = sf.folded_density_2d(dnk, b20["C"], b20["nG"], b20["flat"], ng, workers=NT)[mask2].sum()
                p = np.abs(cr.reshape(N_CELLS, nw)) ** 2; c0 = sf.cell_of(np.array([0, 0, 0]), N_SC)
                wsite[jj] = p[c0, WF_PZ_A] + sum(p[sf.cell_of(np.array(R), N_SC), WF_PZ_B] for R in NN_CELLS) + p[c0, :3].sum()
            Wr[lab] = dict(sel=sel, w2=w, wsite=wsite)
        return Wr
    for tag, Msc in (("c_all", Mall), ("c_3", M3)):
        Hc = HS + Msc; bl, cpl, herm = r4.spectrum_blocks(Hc, parW, E_DW); Wr = w_from_r(bl)
        rec = r5.record_variant(V, tag, bl, Wr, E_DW, Q["thr"], dict(coupling_even_odd=cpl, herm=herm, dim=int(Hc.shape[0])))
        for lab in ("even", "odd"):
            rec["blocks"][lab]["w_site"] = [float(v) for v in Wr[lab]["wsite"]]
        spectra[tag] = np.linalg.eigvalsh(0.5 * (Hc + Hc.conj().T)) - E_DW
    # références : QE (D1), R4 (v1), R5 (variante M^L x 81)
    d4 = json.load(open(os.path.join(R4DIR, "d4", "d4_results.json"))); V["QE_D1"] = d4["variants"]["QE_D1"]
    a5 = json.load(open(os.path.join(R5DIR, "a", "a_results.json")))
    out["reference_R5_Lx81"] = {t: a5["A4"]["variants"][t]["blocks"] for t in ("a1_tot_Lx81", "a3_20_Lx81") if t in a5["A4"]["variants"]}
    out["reference_R4_v1"] = {t: {lab: d4["variants"][t]["blocks"][lab]["localized"] for lab in ("even", "odd")} for t in d4["variants"] if t != "QE_D1"}
    out["reference_R5_v1_corrected_w2"] = {t: {lab: a5["A4"]["variants"][t]["blocks"][lab]["localized"] for lab in ("even", "odd")} for t in a5["A4"]["variants"]}
    # décalage rigide résiduel (états QE < E_D - 4)
    eqe = qe_io.get_eigenvalues(r4.SC_D)[:, 0] * HA2EV - Q["shift"] - Q["E_D_SC"]; fits = {}
    for tag, spv in spectra.items():
        s, resid, n = al.rigid_shift_fit(eqe, spv, -4.0); fits[tag] = dict(shift_eV=float(s), max_residual_eV=float(resid), n_states=int(n))
    out["rigid_shift_fit"] = fits
    # marches
    ladder = ["QE_D1", "a1_tot", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"]; steps = {}

    def match(loc_a, loc_b):
        rows = []
        for xa, wa in loc_a:
            if loc_b:
                j = int(np.argmin([abs(xb - xa) for xb, _ in loc_b])); xb, wb = loc_b[j]; rows.append(dict(x_from=xa, w_from=wa, x_to=xb, w_to=wb, dx=xb - xa))
            else:
                rows.append(dict(x_from=xa, w_from=wa, x_to=None, w_to=None, dx=None))
        return rows
    for i in range(len(ladder) - 1):
        A_, B_ = ladder[i], ladder[i + 1]; steps[f"{A_}->{B_}"] = {}
        for lab in ("even", "odd"):
            la = V[A_]["blocks"][lab]["localized"]; lb = V[B_]["blocks"][lab]["localized"]
            xa = np.array(V[A_]["blocks"][lab]["x"]); xb = np.array(V[B_]["blocks"][lab]["x"])
            steps[f"{A_}->{B_}"][lab] = dict(n_from=len(xa), n_to=len(xb), localized_match=match(la, lb),
                                             median_abs_shift_sorted_eV=(float(np.median(np.abs(np.sort(xa) - np.sort(xb)))) if len(xa) == len(xb) else None))
    for tag in ("a1_tot", "a3_20", "b_5wf", "c_all", "c_3"):
        steps[f"QE_D1->{tag}"] = {lab: dict(localized_match=match(V["QE_D1"]["blocks"][lab]["localized"], V[tag]["blocks"][lab]["localized"])) for lab in ("even", "odd")}
    out["steps"] = steps; out["timing_s"] = time.time() - t_start
    np.savez(os.path.join(d, "d4_spectra.npz"), **spectra, qe_x_all=eqe)
    save_json(os.path.join(d, "d4_results.json"), out)
    write_tables_d4(out, Q); fig_d4(out)


def write_tables_d4(out, Q):
    V = out["variants"]; L = [f"# R6 étape 2 — escalier D4 avec M2 (9x9 ; ε − E_D en eV ; w₂ disque 2 Å, seuil {Q['thr']:.4f}) — {time.strftime('%Y-%m-%d %H:%M')}\n"]
    g = out["gates"]
    L.append(f"Porte 1 (H(R) replié) : {g['gate1_fold_HR']['max_dev_eV']:.2e} eV (seuil 1e-8) → {'OK' if g['gate1_fold_HR']['ok'] else 'ÉCHEC'}. Porte 2 (repliement de M_W(M2)) : {g['gate2_fold_Mwr']['max_dev_eV_production_k']:.2e} eV avec les k du XML, "
             f"{g['gate2_fold_Mwr']['max_dev_eV_snapped_k']:.2e} eV avec k = m/27 (seuil 1e-9) → {'OK' if g['gate2_fold_Mwr']['ok'] else 'ÉCHEC'} ; V_loc production contre k exacts {g['gate2_fold_Mwr']['Vloc_prod_vs_snap_max_abs_eV']:.2e} eV "
             f"(p_z–p_z {g['gate2_fold_Mwr']['onsite_pz_prod']:.6f} / {g['gate2_fold_Mwr']['onsite_pz_snap']:.6f}) ; résidu H(R) vs V†εV {g['one_body_residual']['max_abs_Hwk_minus_VepsV_eV']:.2e} eV.\n")
    L.append(f"E_D : maille 16 b {out['E_D']['uc16_at_K']:.6f}, dense 20 b {out['E_D']['uc20dense_at_K']:.6f}, Wannier {out['E_D']['wannier_at_K']:.6f}, parfaite {out['E_D']['SC_P']:.6f} eV ; M2 linéarité {out['M2_linearity_eV']}.\n")
    L.append("\n| variante | dim | couplage pair–impair (eV) | pairs (σ) localisés (ε − E_D ; w₂) | impairs (π) localisés | décalage rigide résiduel (meV) |")
    L.append("|---|---|---|---|---|---|")
    fits = out["rigid_shift_fit"]
    for tag in ("QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"):
        r = V[tag]; ev = "; ".join(f"{x:+.3f} ({w:.3f})" for x, w in r["blocks"]["even"]["localized"]) or "aucun"; od = "; ".join(f"{x:+.3f} ({w:.3f})" for x, w in r["blocks"]["odd"]["localized"]) or "aucun"
        if tag in ("c_all", "c_3"):
            ws = r["blocks"]["even"].get("w_site", []); x = r["blocks"]["even"]["x"]; ev += " ; poids WF site " + ", ".join(f"{w:.3f}" for w, xx in zip(ws, x) if any(abs(xx - xl) < 1e-6 for xl, _ in r["blocks"]["even"]["localized"]))
        L.append(f"| {tag} | {r.get('dim', '—')} | {r.get('coupling_even_odd', float('nan')):.1e} | {ev} | {od} | {fits[tag]['shift_eV']*1e3:+.1f} |" if tag in fits else f"| {tag} | — | — | {ev} | {od} | — |")
    L.append("\nRéférences : R5 (variante M^L × 81, attendue identique) : " + " ; ".join(f"{t} pairs {[(round(x,3), round(w,3)) for x, w in b['even']['localized']]}, impairs {[(round(x,3), round(w,3)) for x, w in b['odd']['localized']]}" for t, b in out["reference_R5_Lx81"].items()))
    L.append("\nR4 (v1, w₂ corrigés par R5) : " + " ; ".join(f"{t} pairs {[(round(x,3), round(w,3)) for x, w in b['even']]}, impairs {[(round(x,3), round(w,3)) for x, w in b['odd']]}" for t, b in out["reference_R5_v1_corrected_w2"].items() if t in ("a1_tot", "a3_20", "b_5wf", "c_all", "c_3")))
    L.append("\n\n| marche | parité | correspondances (ε_i ; w₂ → ε_{i+1} ; w₂ ; Δε) | médiane \\|Δε\\| spectres triés |")
    L.append("|---|---|---|---|")
    for k, st in out["steps"].items():
        for lab in ("even", "odd"):
            m = st[lab]["localized_match"]; txt = "; ".join((f"{r['x_from']:+.3f} ({r['w_from']:.3f}) → {r['x_to']:+.3f} ({r['w_to']:.3f}) ; {r['dx']:+.3f}" if r["x_to"] is not None else f"{r['x_from']:+.3f} → absent") for r in m) or "—"
            med = st[lab].get("median_abs_shift_sorted_eV"); L.append(f"| {k} | {lab} | {txt} | {('%.4f' % med) if med is not None else '—'} |")
    L.append("\n\nTous les états de la fenêtre par variante (ε − E_D ; w₂) :\n")
    for tag in ("QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"):
        for lab in ("even", "odd"):
            b = V[tag]["blocks"][lab]; L.append(f"- {tag}, {lab} : " + ", ".join(f"{x:+.3f} ({w:.3f})" for x, w in zip(b["x"], b["w2"])))
    with open(os.path.join(WORK, "d4", "d4_tables.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    log("saved d4/d4_tables.md")


def fig_d4(out):
    plt, pal = fig_style(); V = out["variants"]
    order = ["QE_D1", "a1_tot", "a1_L", "a1_NL", "a2_16", "a3_20", "b_5wf", "c_all", "c_3"]
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 4.2), sharey=True)
    for ax, lab, ttl in zip(axs, ("odd", "even"), ("(a) impairs (π)", "(b) pairs (σ)")):
        for i, tag in enumerate(order):
            if tag not in V: continue
            b = V[tag]["blocks"][lab]; x = np.array(b["x"]); w = np.array(b["w2"])
            ax.scatter(np.full(len(x), i), x, s=2 + 60 * w, color=pal.NAVY if tag != "QE_D1" else pal.ORANGE, alpha=0.7, lw=0)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=60, fontsize=6); ax.set_title(ttl); ax.set_ylim(-3, 1); ax.axhline(0, color=pal.MUTED, lw=0.5)
    axs[0].set_ylabel(r"$\varepsilon - E_D$ (eV)"); fig.suptitle("M2 (R6)", fontsize=8); fig.tight_layout(); savefig(fig, "d4_ladder_M2")


# ----------------------------------------------------------------------------------------------- d3 : critère de pôle avec V_loc(M2)
def cmd_d3(a):
    S = r4.setup(verbose=False); d = ensure("d3"); out = dict(nk={}); cfg = S["cfg"]; eta, npe = cfg["eta_eV"], cfg["ne_per_eta"]
    Hwr, Rw, nd, E_D = S["Hwr"], S["Rw"], S["nd"], S["E_D"]
    C2 = np.load(os.path.join(WORK, "cache", "Vloc_M2_9x9.npz")); V = dict(tot=C2["V_tot"], L=C2["V_L"], NL=C2["V_NL"]); Rloc = C2["Rloc"]; nL = len(Rloc)
    C1 = r4.load_cache_vloc(); assert np.array_equal(C1["Rloc"], Rloc), "R_loc différent de R4 : caches g0 inutilisables"
    ip, isg = C2["idx_pi"], C2["idx_sigma"]; grp = r4.index_groups(Rloc); blocks = dict(full=None, pi=ip, sigma=isg)
    _, EK, UK = Hwr_to_Hwk(Hwr, Rw, K_RED[None], ndegen=nd); ph = lt._phase(K_RED[None], Rloc)
    phi_K = np.einsum("kL,kwn->knLw", ph, UK, optimize=True).reshape(1, NW, nL * NW)[0][3:5]
    out["K"] = dict(E_pi_pistar=[float(EK[0, 3]), float(EK[0, 4])], E_D=E_D, R_cut=cfg["R_cut"], dim=nL * NW, eta=eta)
    curves = {}
    for nkint in (300, 600):
        kk = lt.mp_grid(nkint, nkint, 1); Hwk, _, _ = Hwr_to_Hwk(Hwr, Rw, kk, ndegen=nd)
        eg = r4.egrid(E_D, -3.0, 3.0 if nkint == 300 else 1.0, eta, npe); g0 = r4.g0_cached(f"g0_real_nk{nkint}", Hwk, kk, Rloc, eg, eta)
        x = eg - E_D; out["nk"][str(nkint)] = {}; curves[f"nk{nkint}_x"] = x
        for var, Vv in (("tot", V["tot"]), ("L_only", V["L"]), ("NL_only", V["NL"]), ("v1_tot", C1["V_tot"])):
            for bl, idx in blocks.items():
                t0 = time.time(); rr = pc.det_eig_criterion(Vv, g0, idx, vectors=True)
                jd = int(np.argmin(rr["logdet_rel"])); jl = int(np.argmin(rr["minlam"]))
                vfull = np.zeros(nL * NW, complex); vfull[idx if idx is not None else np.arange(nL * NW)] = rr["vec_min"][jl]
                ov = pc.branch_overlaps(rr["vec_min"]); ii, xc = pc.sign_changes(rr["lam_min"].real, x)
                rec = dict(min_det_rel=float(np.exp(rr["logdet_rel"][jd])), det_at=float(x[jd]), minima_det=[(float(x[i]), float(np.exp(rr["logdet_rel"][i])), float(rr["minlam"][i])) for i in pc.local_minima(rr["logdet_rel"], x)],
                           min_abs_lam=float(rr["minlam"][jl]), lam=complex(rr["lam_min"][jl]), lam_at=float(x[jl]), weights=pc.eigvec_weights(vfull, grp),
                           weight_pi=float((np.abs(vfull[ip]) ** 2).sum()), weight_sigma=float((np.abs(vfull[isg]) ** 2).sum()),
                           roots_Re_lam_min=[(float(e), float(ov[i])) for i, e in zip(ii, xc)], minima_lam=[(float(x[i]), float(rr["minlam"][i])) for i in pc.local_minima(rr["minlam"], x)])
                if bl in ("pi", "full"):
                    Vb = Vv if idx is None else Vv[np.ix_(idx, idx)]; gb = g0 if idx is None else g0[:, idx][:, :, idx]
                    tc = pc.local_t_cache(Vb, gb); PK = phi_K if idx is None else phi_K[:, idx]
                    _, tr = pc.tbar_pair(tc, PK); jt = int(np.argmax(-tr.imag))
                    rec["peak_minus_Im_Tbar"] = float(x[jt]); rec["max_minus_Im_Tbar"] = float(-tr.imag[jt]); rec["Tbar_at_ED"] = complex(np.interp(0.0, x, tr.real) + 1j * np.interp(0.0, x, tr.imag))
                    rec["Re_Tbar_zero_crossings"] = [float(v) for v in pc.sign_changes(tr.real, x)[1]]
                    curves[f"nk{nkint}_{var}_{bl}_ImTbar"] = -tr.imag; curves[f"nk{nkint}_{var}_{bl}_ReTbar"] = tr.real; del tc
                curves[f"nk{nkint}_{var}_{bl}_detrel"] = np.exp(rr["logdet_rel"]); curves[f"nk{nkint}_{var}_{bl}_minlam"] = rr["minlam"]
                out["nk"][str(nkint)][f"{var}_{bl}"] = rec
                log(f"[D3 M2] nk{nkint} {var} {bl}: min|det|/max {rec['min_det_rel']:.2e} à {rec['det_at']:+.3f} ; min|λ| {rec['min_abs_lam']:.4f} ({rec['lam']:.4f}) à {rec['lam_at']:+.3f} ; "
                    f"poids π {rec['weight_pi']:.3f} σ {rec['weight_sigma']:.3f} ; racines Re λ {[round(e, 3) for e, _ in rec['roots_Re_lam_min']]}"
                    + (f" ; pic −Im T̄ {rec['peak_minus_Im_Tbar']:+.3f} (max {rec['max_minus_Im_Tbar']:.2f}), Re T̄(E_D) {rec['Tbar_at_ED'].real:+.3f}" if "peak_minus_Im_Tbar" in rec else "") + f" ({time.time()-t0:.0f} s)")
    np.savez(os.path.join(d, "d3_curves_M2.npz"), **curves); save_json(os.path.join(d, "d3_results_M2.json"), out)
    # tableau + figure (300², M2 tot contre v1)
    L = [f"# R6 étape 2 — critère de pôle avec V_loc(M2), α = 1, R_cut 3 (dim {nL*NW}), η {eta} eV ({time.strftime('%Y-%m-%d %H:%M')})\n",
         "| grille | V | bloc | min \\|det\\|/max (ε − E_D) | min \\|λ\\| ; λ (ε − E_D) | poids π / σ du vecteur | racines de Re λ_min | pic −Im T̄(K) (max) | Re T̄(E_D) | minima locaux de \\|det\\| |", "|---|---|---|---|---|---|---|---|---|---|"]
    for nkint in ("300", "600"):
        for key, rec in out["nk"][nkint].items():
            L.append(f"| {nkint}² | {key.rsplit('_', 1)[0]} | {key.rsplit('_', 1)[1]} | {rec['min_det_rel']:.2e} ({rec['det_at']:+.3f}) | {rec['min_abs_lam']:.4f} ; {rec['lam'].real:+.4f}{rec['lam'].imag:+.4f}i ({rec['lam_at']:+.3f}) | {rec['weight_pi']:.3f} / {rec['weight_sigma']:.3f} | "
                     f"{', '.join(f'{e:+.3f}' for e, _ in rec['roots_Re_lam_min']) or '—'} | {(f'{rec[chr(112)+chr(101)+chr(97)+chr(107)+chr(95)+chr(109)+chr(105)+chr(110)+chr(117)+chr(115)+chr(95)+chr(73)+chr(109)+chr(95)+chr(84)+chr(98)+chr(97)+chr(114)]:+.3f} ({rec[chr(109)+chr(97)+chr(120)+chr(95)+chr(109)+chr(105)+chr(110)+chr(117)+chr(115)+chr(95)+chr(73)+chr(109)+chr(95)+chr(84)+chr(98)+chr(97)+chr(114)]:.2f})' if 'peak_minus_Im_Tbar' in rec else '—')} | "
                     f"{(f'{rec[chr(84)+chr(98)+chr(97)+chr(114)+chr(95)+chr(97)+chr(116)+chr(95)+chr(69)+chr(68)].real:+.3f}' if 'Tbar_at_ED' in rec else '—')} | {', '.join(f'{e:+.3f} ({dd:.1e})' for e, dd, _ in rec['minima_det'][:6])} |")
    with open(os.path.join(d, "d3_tables_M2.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    plt, pal = fig_style(); x = curves["nk300_x"]
    fig, ax = plt.subplots(1, 3, figsize=(6.5, 3.0))
    for bl, c in (("full", pal.NAVY), ("pi", pal.ORANGE), ("sigma", pal.CYCLE[2])):
        ax[0].plot(x, curves[f"nk300_tot_{bl}_detrel"], color=c, lw=0.8, label={"full": "complète", "pi": "bloc π", "sigma": "bloc σ"}[bl])
        ax[0].plot(x, curves[f"nk300_v1_tot_{bl}_detrel"], color=c, lw=0.6, ls="--")
        ax[1].plot(x, curves[f"nk300_tot_{bl}_minlam"], color=c, lw=0.8); ax[1].plot(x, curves[f"nk300_v1_tot_{bl}_minlam"], color=c, lw=0.6, ls="--")
        if bl != "sigma":
            ax[2].plot(x, curves[f"nk300_tot_{bl}_ImTbar"], color=c, lw=0.8); ax[2].plot(x, curves[f"nk300_v1_tot_{bl}_ImTbar"], color=c, lw=0.6, ls="--")
    ax[0].set_yscale("log"); ax[0].set_ylabel(r"$|\det|/\max$"); ax[1].set_yscale("log"); ax[1].set_ylabel(r"$\min_i |\lambda_i|$"); ax[2].set_ylabel(r"$-\mathrm{Im}\,\bar T(K)$ (eV)")
    for a_ in ax:
        a_.set_xlabel(r"$\varepsilon - E_D$ (eV)"); a_.axvline(-0.737, color=pal.MUTED, lw=0.5)
    ax[0].legend(fontsize=5); ax[0].set_title("M2 (plein) et v1 (tirets), $300^2$", fontsize=8); fig.tight_layout(); savefig(fig, "d3_pole_M2")


# ----------------------------------------------------------------------------------------------- b3 : (a1) à n bandes, M2 nb128
def cmd_b3(a):
    d = ensure("b3"); out = dict(n_list=NLIST); t_start = time.time()
    B = r5.load_bases(with_dense=False, nb128=True); ng = B["ng"]; A_A = B["A_A"]; b = B["bases"][128]
    Q = r5.qe_states_9x9(); ns = len(Q["e"]); hl = r5.highlighted(Q); mask2 = qg.inplane_disc_mask(ng, A_A, Q["s_vac"], r4.R_W2)
    Mt = load_v2(FILES_128["tot"]); out["E_D_128"] = b["E_D"]; out["max_abs_M2_128_eV"] = float(np.abs(Mt).max())
    t0 = time.time(); c128 = np.zeros((ns, 128, 81), complex)
    for s in range(ns):
        A = sp.sc_state_grid(Q["Cw"][s], Q["mill"], ng, gamma_only=Q["gamma_only"]); c128[s] = sp.bloch_overlaps(A, b["C"], b["nG"], b["flat"])
    log(f"[B3 M2] recouvrements 128 bandes en {time.time()-t0:.0f} s")
    d4 = json.load(open(os.path.join(R4DIR, "d4", "d4_results.json"))); qe_rec = d4["variants"]["QE_D1"]
    b3_r5 = json.load(open(os.path.join(R5DIR, "b", "b3_results.json"))) if os.path.exists(os.path.join(R5DIR, "b", "b3_results.json")) else None
    per_n = {}; prev = None
    for n in NLIST:
        t0 = time.time(); Mn = Mt[:n, :, :n, :]; eps_n = b["eps"][:n]; par_n = b["par"][:n].reshape(-1)
        H = sf.bloch_folded_hamiltonian(eps_n, Mn, N_CELLS); blocks, cpl, herm = r4.spectrum_blocks(H, par_n, b["E_D"])
        W = r5.window_weights(blocks, b["C"][:n], b["nG"], b["flat"], n, ng, mask2); Vv = {}; rec = r5.record_variant(Vv, f"n{n}", blocks, W, b["E_D"], Q["thr"], dict(coupling_even_odd=cpl, dim=int(H.shape[0])))
        odd = rec["blocks"]["odd"]["localized"]; even = rec["blocks"]["even"]["localized"]
        pi = max(odd, key=lambda t: t[1]) if odd else None; sig_near = [t for t in even if abs(t[0] - 0.101) < 0.3]
        deltas = {lab: float(1 - np.sum(np.abs(c128[s, :n, :]) ** 2)) for lab, s in hl.items()}
        r = dict(n=n, dim=int(H.shape[0]), coupling=cpl, blocks=rec["blocks"], pi_state=pi, sigma_localized=even, sigma_near_0p101=sig_near, delta=deltas,
                 n_window={lab: rec["blocks"][lab]["n_window"] for lab in ("even", "odd")},
                 median_dx_vs_qe={lab: (float(np.median(np.abs(np.sort(np.array(rec["blocks"][lab]["x"])) - np.sort(np.array(qe_rec["blocks"][lab]["x"]))))) if len(rec["blocks"][lab]["x"]) == len(qe_rec["blocks"][lab]["x"]) else None) for lab in ("even", "odd")})
        if prev is not None:
            r["dpi_vs_prev"] = (pi[0] - prev["pi_state"][0]) if (pi and prev["pi_state"]) else None
        if b3_r5 is not None and str(n) in b3_r5["per_n"]:
            r["R5_Lx81"] = dict(pi_state=b3_r5["per_n"][str(n)]["Lx81"]["pi_state"], sigma_localized=b3_r5["per_n"][str(n)]["Lx81"]["sigma_localized"])
        per_n[str(n)] = r; prev = r
        log(f"[B3 M2] n={n}: dim {H.shape[0]}, fenêtre σ {r['n_window']['even']} / π {r['n_window']['odd']}, π localisé {pi and (round(pi[0],3), round(pi[1],3))}, σ localisés {[(round(x,3), round(w,3)) for x, w in even]}, "
            f"σ à < 0,3 eV de +0,101 {[(round(x,3), round(w,3)) for x, w in sig_near]} ; δ(π_320) {deltas['pi_320']:.4f}, δ(σ_324) {deltas['sigma_324']:.4f} ; "
            f"R5 Lx81 π {r.get('R5_Lx81', {}).get('pi_state')} ({time.time()-t0:.0f} s)")
        del H, blocks
    out["per_n"] = per_n; out["first_n_sigma_near_0p101"] = next((n for n in NLIST if per_n[str(n)]["sigma_near_0p101"]), None)
    out["pi_vs_n"] = [(n, per_n[str(n)]["pi_state"]) for n in NLIST]; out["timing_s"] = time.time() - t_start
    save_json(os.path.join(d, "b3_results_M2.json"), out)
    L = [f"# R6 étape 2 — (a1) à n bandes avec le M2 à 128 bandes (ε − E_D ; w₂ ; seuil {Q['thr']:.4f}) — {time.strftime('%Y-%m-%d %H:%M')}\n",
         "| n | dim | fenêtre σ / π | π localisé (ε − E_D ; w₂) | Δε_π vs n précédent (meV) | σ localisés | σ à < 0,3 eV de +0,101 | δ(π_320) / δ(σ_324) | médiane \\|Δε\\| vs QE σ / π | R5 (M^L × 81) π |", "|---|---|---|---|---|---|---|---|---|---|"]
    for n in NLIST:
        r = per_n[str(n)]; pi = r["pi_state"]; m = r["median_dx_vs_qe"]; r5p = r.get("R5_Lx81", {}).get("pi_state")
        L.append(f"| {n} | {r['dim']} | {r['n_window']['even']} / {r['n_window']['odd']} | {(f'{pi[0]:+.3f} ({pi[1]:.3f})' if pi else 'aucun')} | {(f'{r[chr(100)+chr(112)+chr(105)+chr(95)+chr(118)+chr(115)+chr(95)+chr(112)+chr(114)+chr(101)+chr(118)]*1e3:+.1f}' if r.get('dpi_vs_prev') is not None else '—')} | "
                 f"{'; '.join(f'{x:+.3f} ({w:.3f})' for x, w in r['sigma_localized']) or 'aucun'} | {'; '.join(f'{x:+.3f} ({w:.3f})' for x, w in r['sigma_near_0p101']) or '—'} | {r['delta']['pi_320']:.4f} / {r['delta']['sigma_324']:.4f} | "
                 f"{(f'{m[chr(101)+chr(118)+chr(101)+chr(110)]:.3f}' if m['even'] is not None else '—')} / {(f'{m[chr(111)+chr(100)+chr(100)]:.3f}' if m['odd'] is not None else '—')} | {(f'{r5p[0]:+.3f} ({r5p[1]:.3f})' if r5p else '—')} |")
    with open(os.path.join(d, "b3_tables_M2.md"), "w") as f:
        f.write("\n".join(L) + "\n")
    plt, pal = fig_style(); fig, ax = plt.subplots(figsize=(3.4, 3.0))
    ax.plot(NLIST, [per_n[str(n)]["pi_state"][0] if per_n[str(n)]["pi_state"] else np.nan for n in NLIST], "o-", color=pal.NAVY, ms=3, label=r"$\pi$ localisé (M2)")
    sig = [min((t[0] for t in per_n[str(n)]["sigma_localized"] if t[0] > -1.5), default=np.nan) for n in NLIST]
    ax.plot(NLIST, sig, "s-", color=pal.ORANGE, ms=3, label=r"$\sigma$ localisé le plus haut (M2)")
    ax.axhline(-0.737, color=pal.NAVY, lw=0.5, ls="--"); ax.axhline(0.101, color=pal.ORANGE, lw=0.5, ls="--")
    ax.set_xlabel("nombre de bandes $n$"); ax.set_ylabel(r"$\varepsilon - E_D$ (eV)"); ax.legend(fontsize=6); fig.tight_layout(); savefig(fig, "b3_ladder_M2")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmds", nargs="+", choices=["prep", "d4", "d3", "b3"]); a = ap.parse_args()
    for c in a.cmds:
        log(f"===== r6_d4 {c} ====="); {"prep": cmd_prep, "d4": cmd_d4, "d3": cmd_d3, "b3": cmd_b3}[c](a)


if __name__ == "__main__":
    main()
