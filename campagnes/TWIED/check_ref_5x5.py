"""
check_ref_5x5.py
    Replays, on the two reference files written by make_ref_5x5.py (twied_dir: ref_5x5_inputs.h5, ref_5x5_chain.h5), the checks
    made while building them: structure and units, pairing of the two provenances, and each step of the chain recomputed from
    the inputs stored in the chain file. One PASS/FAIL line per check; exit code 0 if every check passes.

        --sources   also re-hashes the source files and compares /input with the QE / Wannier90 readers (local data needed)
        --full      also recomputes g0 (cluster_green_batch) and Gamma (scattering_rate_fast), about one minute

Usage (repository root): .venv/bin/python campagnes/TWIED/check_ref_5x5.py [--sources] [--full]
"""
import argparse
import os
import subprocess
import sys

import h5py
import numpy as np

from graphene_raman.config import ROOT, HA2EV, load_production, twied_dir, alignment_C
from graphene_raman.wannier.wannier_interpolation import Mbk_to_Mwk, Mwk_to_Mwr, _match_kpoint_order
from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk, dirac_point
from graphene_raman.defects.many_body import cluster_tmatrix as ct

RTOL = 1e-12            # same code on the same inputs: exact, or round-off at most
VARIANTS = ("unaligned", "kumagai_oba")
RESULTS = []

CHAIN_DATASETS = (
    ["M_coarse/" + n for n in ("M_L", "M_NL", "M")]
    + ["M_W/" + n for n in ("U", "U_dis", "k_coarse", "Mwk", "Mwr_raw", "R")]
    + ["G0/" + n for n in ("H_R", "R_w", "ndegen", "E_D", "gap", "R_cluster", "egrid", "g0")]
    + ["Gamma_grid/" + n for n in ("E_out", "e_window")]
    + [f"variants/{v}/{n}" for v in VARIANTS for n in ("C_N", "Mwr", "in_box", "Rn", "R_d", "dist", "wt", "R_cluster",
                                                         "M_cluster", "R_cut_1/R_cluster", "R_cut_1/M_cluster", "t", "Gamma")]
)
INPUTS_DATASETS = (
    ["input/unit_cell/" + n for n in ("C_nk", "nG", "G_red", "k_red", "eps", "A_cols", "Omega", "x_red", "ecut")]
    + [f"input/supercell_{s}/{n}" for s in ("p", "d") for n in ("A_cols", "Omega", "x_red", "V")]
    + ["input/pseudo/" + n for n in ("ekb_li", "fr_li", "rgrid", "V_L")]
    + ["input/wannier/" + n for n in ("U", "U_dis", "k_w90", "perm", "H_R", "R", "ndegen", "r_R", "lattice")]
    + ["M_coarse/prep/" + n for n in ("Ved", "Omega_sc")]
)


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ""), flush=True)


def same(a, b):
    """(ok, detail): exact equality, or max|a - b| <= RTOL max|b| for floats (NaN at the same places)."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False, f"shapes {a.shape} vs {b.shape}"
    inexact = np.issubdtype(a.dtype, np.inexact) or np.issubdtype(b.dtype, np.inexact)
    if np.array_equal(a, b, equal_nan=inexact):
        return True, "exact"
    if not inexact:
        return False, "integer/bool arrays differ"
    if not np.array_equal(np.isnan(a), np.isnan(b)):
        return False, "NaN at different places"
    dev = np.nanmax(np.abs(a - b)) / max(np.nanmax(np.abs(b)), 1e-300)
    return dev <= RTOL, f"rel. dev. {dev:.1e}"


def herm_residual(M):
    nb, nk = M.shape[0], M.shape[1]
    Mf = M.reshape(nb * nk, nb * nk)
    return float(np.abs(Mf - Mf.conj().T).max())


def check_structure(fi, fc, path_chain):
    for f, expected in ((fi, INPUTS_DATASETS), (fc, CHAIN_DATASETS)):
        missing = [n for n in expected if n not in f]
        check(f"{os.path.basename(f.filename)}: expected datasets present", not missing, f"missing {missing}" if missing else "")
        no_units = []
        f.visititems(lambda n, o: no_units.append(n) if isinstance(o, h5py.Dataset) and "units" not in o.attrs else None)
        check(f"{os.path.basename(f.filename)}: every dataset has 'units'", not no_units, f"{no_units}" if no_units else "")
    size = os.path.getsize(path_chain) / 2**20
    check("chain file under 50 MiB (GitHub warning threshold)", size < 50, f"{size:.1f} MiB")


def check_provenance(fi, fc):
    pi, pc = dict(fi["provenance"].attrs), dict(fc["provenance"].attrs)
    check("roles are 'inputs' and 'chain'", (pi["role"], pc["role"]) == ("inputs", "chain"))
    for key in ("run_id", "git_commit", "script"):
        check(f"same provenance '{key}' in both files", pi[key] == pc[key], str(pc[key]))
    check("same source files and md5 in both files",
          dict(fi["provenance/source_files"].attrs) == dict(fc["provenance/source_files"].attrs),
          f"{len(fc['provenance/source_files'].attrs)} files")
    check("produced from a clean working tree (strict mode)", bool(pc["git_clean"]))
    known = subprocess.run(["git", "cat-file", "-e", f"{pc['git_commit']}^{{commit}}"], cwd=ROOT,
                           capture_output=True).returncode == 0
    check("git_commit exists in this repository", known)
    check("HA2EV equals config.HA2EV", pc["HA2EV"] == HA2EV)


def check_M(fi, fc):
    g = fc["M_coarse"]
    M_L, M_NL, M = g["M_L"][()], g["M_NL"][()], g["M"][()]
    check("M = M_L + M_NL", *same(M, M_L + M_NL))
    for name, X in (("M_L", M_L), ("M_NL", M_NL), ("M", M)):
        r = herm_residual(X)
        check(f"{name} Hermitian, residual equals stored attribute", r < 1e-10 and r == g[name].attrs["herm_residual"], f"{r:.1e} Ha")
    Ved = fi["M_coarse/prep/Ved"][()]
    check("Ved = V_d - V_p (inputs file)", *same(Ved, fi["input/supercell_d/V"][()] - fi["input/supercell_p/V"][()]))


def check_MW(fi, fc):
    g = fc["M_W"]
    M_eV = fc["M_coarse/M"][()] * fc["provenance"].attrs["HA2EV"]
    U, U_dis, k, MP = g["U"][()], g["U_dis"][()], g["k_coarse"][()], tuple(g.attrs["MP"])
    perm = fi["input/wannier/perm"][()]
    check("M_W/U = input/wannier/U[perm]", *same(U, fi["input/wannier/U"][()][perm]))
    check("M_W/U_dis = input/wannier/U_dis[perm]", *same(U_dis, fi["input/wannier/U_dis"][()][perm]))
    check("M_W/k_coarse = input/unit_cell/k_red", *same(k, fi["input/unit_cell/k_red"][()]))
    Mwk = Mbk_to_Mwk(M_eV, U, U_dis)
    check("Mwk = Mbk_to_Mwk(M x HA2EV, U, U_dis)", *same(Mwk, g["Mwk"][()]))
    Mwr, R = Mwk_to_Mwr(g["Mwk"][()], k, MP)
    check("Mwr_raw = Mwk_to_Mwr(Mwk)", *same(Mwr, g["Mwr_raw"][()]))
    check("R = labels of Mwk_to_Mwr", *same(R, g["R"][()]))


def check_variants(fc, cfg):
    g = fc["M_W"]
    M_eV = fc["M_coarse/M"][()] * fc["provenance"].attrs["HA2EV"]
    U, U_dis, k, MP = g["U"][()], g["U_dis"][()], g["k_coarse"][()], tuple(g.attrs["MP"])
    n_box, Mwr_raw = int(g.attrs["n_box"]), g["Mwr_raw"][()]
    expected_C_N = {"unaligned": 0.0, "kumagai_oba": alignment_C(cfg, "5x5")}
    nw, nR = Mwr_raw.shape[0], Mwr_raw.shape[1]
    for v in VARIANTS:
        gv = fc[f"variants/{v}"]
        C_N, Mwr = float(gv["C_N"][()]), gv["Mwr"][()]
        check(f"[{v}] C_N matches the config", C_N == expected_C_N[v], f"{C_N:+.6f} eV")
        d = ct.defect_mwr(M_eV, U, U_dis, k, MP, n_box, C_N)
        for key in ("Mwr", "Rn", "R_d", "in_box"):
            check(f"[{v}] {key} = defect_mwr(...)['{key}']", *same(d[key], gv[key][()]))
        in_box = gv["in_box"][()]
        expected = Mwr_raw.copy()                                                   # the alignment, redone: same rounding
        for w in range(nw):
            expected[w, in_box, w, in_box] -= C_N
        check(f"[{v}] Mwr = Mwr_raw - C_N on the in_box diagonal, unchanged elsewhere", *same(Mwr, expected))
        dist, wt = ct.mwr_locality(Mwr, gv["Rn"][()])
        check(f"[{v}] dist, wt = mwr_locality(Mwr, Rn)", same(dist, gv["dist"][()])[0] and same(wt, gv["wt"][()])[0],
              f"on-site {wt[0]:.3f} eV")
        for sub, R_cut in (("", int(cfg["R_cut"])), ("R_cut_1/", 1)):
            Rc = ct.cluster_cells(gv["Rn"][()], R_cut)
            Mc, herm = ct.cluster_potential(Mwr, gv["Rn"][()], Rc)
            check(f"[{v}] {sub or 'R_cut=' + str(R_cut) + ' '}R_cluster, M_cluster recomputed",
                  same(Rc, gv[sub + "R_cluster"][()])[0] and same(Mc, gv[sub + "M_cluster"][()])[0],
                  f"{len(Rc)} cells, herm {herm:.1e} eV")
            check(f"[{v}] {sub}R_cluster attribute R_cut = {R_cut}", int(gv[sub + "R_cluster"].attrs["R_cut"]) == R_cut)
        check(f"[{v}] R_cluster = G0/R_cluster (shared g0)", *same(gv["R_cluster"][()], fc["G0/R_cluster"][()]))


def check_G0_t(fc, full):
    g = fc["G0"]
    H_R, R_w, nd = g["H_R"][()], g["R_w"][()], g["ndegen"][()]
    E = Hwr_to_Hwk(H_R, R_w, ct.mp_grid(90, 90, 1), ndegen=nd)[1]
    E_D, gap = dirac_point(E)
    check("E_D, gap = dirac_point(H on 90x90)", E_D == g["E_D"][()] and gap == g["gap"][()], f"E_D {E_D:.4f} eV, gap {gap*1e3:.1f} meV")
    egrid = g["egrid"][()]
    check("egrid centred on E_D", egrid[len(egrid) // 2] == E_D, f"{len(egrid)} energies")
    g0 = g["g0"][()]
    if full:
        k_int = ct.mp_grid(int(g.attrs["k_int"]))
        Hk = Hwr_to_Hwk(H_R, R_w, k_int, ndegen=nd)[0]
        check("[full] g0 = cluster_green_batch(...)", *same(ct.cluster_green_batch(Hk, k_int, g["R_cluster"][()], egrid,
                                                                                 float(g.attrs["eta"])), g0))
    for v in VARIANTS:
        gv = fc[f"variants/{v}"]
        Mc, t = gv["M_cluster"][()], gv["t"][()]
        check(f"[{v}] t = cluster_t(M_cluster, g0) at every energy",
              all(same(ct.cluster_t(Mc, g0[e]), t[e])[0] for e in range(len(egrid))))
        ls = max(np.abs(t[e] - (Mc + Mc @ g0[e] @ t[e])).max() for e in range(len(egrid)))
        check(f"[{v}] Lippmann-Schwinger t = M + M g0 t at every energy", ls < 1e-10, f"max residual {ls:.1e} eV")


def check_Gamma(fc, full):
    gg, g = fc["Gamma_grid"], fc["G0"]
    E_out, (lo, hi) = gg["E_out"][()], gg["e_window"][()]
    check("e_window = E_D -/+ half_width", np.isclose(lo, g["E_D"][()] - gg.attrs["half_width"], rtol=0, atol=1e-12)
          and np.isclose(hi, g["E_D"][()] + gg.attrs["half_width"], rtol=0, atol=1e-12), f"({lo:.4f}, {hi:.4f}) eV")
    inside = ((E_out >= lo) & (E_out <= hi)).T                                   # (band, k_out), the layout of Gamma
    for v in VARIANTS:
        G = fc[f"variants/{v}/Gamma"][()]
        check(f"[{v}] Gamma finite exactly inside e_window", np.array_equal(np.isfinite(G), inside), f"{inside.sum()} states")
        check(f"[{v}] Gamma >= 0", np.nanmin(G) >= -1e-8, f"min {np.nanmin(G)*1e3:.1f} meV, median {np.nanmedian(G)*1e3:.1f} meV")
    if full:
        k_out, k_int = ct.mp_grid(int(gg.attrs["N"])), ct.mp_grid(int(g.attrs["k_int"]))
        for v in VARIANTS:
            gv = fc[f"variants/{v}"]
            G = ct.scattering_rate_fast(g["H_R"][()], g["R_w"][()], g["ndegen"][()], gv["M_cluster"][()], gv["R_cluster"][()],
                                        k_out, float(g.attrs["eta"]), k_int=k_int, e_window=(lo, hi),
                                        ne_per_eta=int(gg.attrs["ne_per_eta"]))
            check(f"[full][{v}] Gamma = scattering_rate_fast(...)", *same(G, gv["Gamma"][()]))


def check_sources(fi):
    from graphene_raman.io import qe_io, pseudo_io, wannier_io
    from make_ref_5x5 import INPUTS, POT_PATH, PSEUDO_PATH, md5_of_file          # the paths of the producing script
    src = dict(fi["provenance/source_files"].attrs)
    bad = [p for p, m in src.items() if not os.path.isfile(os.path.join(ROOT, p)) or md5_of_file(os.path.join(ROOT, p)) != m]
    check("[sources] md5 of every source file unchanged", not bad, f"{len(src)} files" + (f", differ: {bad}" if bad else ""))
    uc, g = INPUTS["unit_cell"], fi["input/unit_cell"]
    C, nG = qe_io.get_C_nk(uc)
    A, Om = qe_io.get_A_volume(uc)
    ref = {"C_nk": C, "nG": nG, "G_red": qe_io.get_G_red(uc), "k_red": qe_io.get_k_red(uc), "eps": qe_io.get_eigenvalues(uc),
           "A_cols": A, "Omega": Om, "x_red": qe_io.get_x_red(uc), "ecut": qe_io.get_ecut(uc)}
    check("[sources] input/unit_cell = qe_io readers", all(same(v, g[n][()])[0] for n, v in ref.items()))
    for s in ("p", "d"):
        g = fi[f"input/supercell_{s}"]
        A, Om = qe_io.get_A_volume(INPUTS[s])
        ok = (same(A, g["A_cols"][()])[0] and same(Om, g["Omega"][()])[0]
              and same(qe_io.get_x_red(INPUTS[s]), g["x_red"][()])[0] and same(qe_io.read_filplot(POT_PATH[s])["V"], g["V"][()])[0])
        check(f"[sources] input/supercell_{s} = qe_io readers", ok)
    ekb, fr, r, lmax, imax, V_L = pseudo_io.read_upf(PSEUDO_PATH)
    g = fi["input/pseudo"]
    check("[sources] input/pseudo = read_upf", all(same(v, g[n][()])[0] for n, v in
                                                  (("ekb_li", ekb), ("fr_li", fr), ("rgrid", r), ("V_L", V_L)))
          and (g.attrs["lmax"], g.attrs["imax"]) == (lmax, imax))
    W, g = INPUTS["wannier"], fi["input/wannier"]
    U, k_w90 = wannier_io.read_w90_mat(f"{W}/wannier_u.mat")
    U_dis, _ = wannier_io.read_w90_mat(f"{W}/wannier_u_dis.mat")
    H_R, R, nd, r_R, lat = wannier_io.read_w90_tb(f"{W}/wannier_tb.dat")
    perm = _match_kpoint_order(k_w90, qe_io.get_k_red(uc))
    ref = {"U": U, "U_dis": U_dis, "k_w90": k_w90, "perm": perm, "H_R": H_R, "R": R, "ndegen": nd, "r_R": r_R, "lattice": lat}
    check("[sources] input/wannier = wannier_io readers", all(same(v, g[n][()])[0] for n, v in ref.items()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[2].strip())
    ap.add_argument("--sources", action="store_true", help="re-hash the source files and compare /input with the readers")
    ap.add_argument("--full", action="store_true", help="also recompute g0 and Gamma (about one minute)")
    a = ap.parse_args()
    cfg = load_production(verbose=False)
    path_inputs = os.path.join(twied_dir(cfg), "ref_5x5_inputs.h5")
    path_chain = os.path.join(twied_dir(cfg), "ref_5x5_chain.h5")
    with h5py.File(path_inputs, "r") as fi, h5py.File(path_chain, "r") as fc:
        print(f"checking {os.path.relpath(path_inputs, ROOT)} and {os.path.relpath(path_chain, ROOT)} "
              f"(run_id {fc['provenance'].attrs['run_id']}, commit {fc['provenance'].attrs['git_commit'][:7]})", flush=True)
        check_structure(fi, fc, path_chain)
        check_provenance(fi, fc)
        check_M(fi, fc)
        check_MW(fi, fc)
        check_variants(fc, cfg)
        check_G0_t(fc, a.full)
        check_Gamma(fc, a.full)
        if a.sources:
            check_sources(fi)
    n_fail = RESULTS.count(False)
    print(f"RESULT: {'PASS' if n_fail == 0 else 'FAIL'} ({len(RESULTS) - n_fail}/{len(RESULTS)} checks)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
