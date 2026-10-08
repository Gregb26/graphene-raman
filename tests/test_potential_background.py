"""R11 (2026-10-08): qe_io.read_filplot and the background measures of alignment.py (kumagai_oba_constant, plane_radial_profile,
z_profile, potential_background). Synthetic checks without data; on the local 5x5 pair (needs_data) the header route reproduces the XML
geometry, get_pot, C_5 of config/production.json, the vacuum offset of tab:alignement and the masked profile of results/M2_plateau.
Run: .venv/bin/python -m pytest tests/test_potential_background.py"""
import csv
import json
import os
import subprocess
import sys

import numpy as np
import pytest
from conftest import REPO, local_data, _paths
from graphene_raman.config import HA2EV, load_production, alignment_C, results_dir
from graphene_raman.io import qe_io
from graphene_raman.defects import alignment as al
from test_r9_functions import _graphene_supercell

BOHR = 0.529177210903
POT = (_paths.pot_d("5x5"), _paths.pot_p("5x5"), _paths.sc_d("5x5"), _paths.sc_p("5x5"))
has_5x5 = local_data(*POT)


def _write_filplot(path, V_ry, A_cols_bohr, x_red, alat, plot_num=1):
    """Minimal pp.x filplot (iflag=3, ibrav=0) with the layout read by read_filplot / skipped by get_pot. V_ry [ix, iy, iz]."""
    n1, n2, n3 = V_ry.shape
    at = A_cols_bohr.T / alat
    tau = (x_red @ A_cols_bohr.T) / alat
    L = ["title", f"{n1} {n2} {n3} {n1} {n2} {n3} {len(x_red)} 1", f"0 {alat:.10f} 0 0 0 0 0"]
    L += [" ".join(f"{v:.17g}" for v in row) for row in at]
    L += [f"100.0 4.0 50.0 {plot_num}", "1 C 4.00"]
    L += [f"{i + 1} " + " ".join(f"{v:.12f}" for v in t) + " 1" for i, t in enumerate(tau)]
    flat = V_ry.transpose(2, 1, 0).ravel()                                   # x fastest, as pp.x
    L += [" ".join(f"{v:.9E}" for v in flat[i:i + 5]) for i in range(0, flat.size, 5)]
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def test_read_filplot_roundtrip(tmp_path):
    rng = np.random.default_rng(3)
    A_A, x = _graphene_supercell(2, c=10.0)
    A_b = A_A / BOHR
    alat = float(np.linalg.norm(A_b[:, 0]))
    V = rng.standard_normal((6, 5, 8))
    p = str(tmp_path / "Vks_test")
    _write_filplot(p, V, A_b, x, alat, plot_num=11)
    h = qe_io.read_filplot(p)
    V_pot, ng = qe_io.get_pot(p, subtract_mean=False, to_hartree=True)
    assert h["ngfft"] == (6, 5, 8) == ng and h["plot_num"] == 11 and h["species"] == ["C"] and np.all(h["ityp"] == 1)
    assert np.array_equal(h["V"], V_pot.transpose(2, 1, 0))                  # same values, same layout as get_pot
    assert np.allclose(h["V"], 0.5 * V, atol=1e-8) and h["V"].flags["C_CONTIGUOUS"]
    assert np.allclose(h["A_cols"], A_b, atol=1e-9) and abs(h["alat"] - alat) < 1e-9
    assert np.allclose(h["x_red"], x, atol=1e-9)
    assert np.allclose(qe_io.read_filplot(p, to_hartree=False)["V"], V, atol=1e-8)


def test_constant_offset_is_seen_by_the_three_measures():
    """V_d = V_p + c0: C_N, in-plane background and vacuum all equal c0 with zero spread; the vacancy is found."""
    rng = np.random.default_rng(0)
    A, x_p = _graphene_supercell(4, c=16.0)
    iv = 11
    x_d = np.delete(x_p, iv, axis=0)
    n = (48, 48, 64)
    g = np.stack(np.meshgrid(*[np.arange(m) / m for m in n], indexing="ij"), -1)
    V_p = np.cos(2 * np.pi * g[..., 0]) + 0.3 * np.sin(2 * np.pi * (g[..., 1] + g[..., 2])) + 0.1 * rng.standard_normal(n)
    c0 = -0.025
    r = al.potential_background(V_p + c0, V_p, x_d, x_p, A)
    assert np.allclose(r["s_vac"], np.mod(x_p[iv], 1.0))
    assert abs(r["C_N"] - c0) < 1e-13 and r["C_N_rms"] < 1e-13 and 0 < r["C_N_n"] < len(x_d)
    assert abs(r["plane"] - c0) < 1e-13 and r["plane_rms"] < 1e-13 and all(abs(v - c0) < 1e-13 for v in r["plane_at"].values())
    assert abs(r["vacuum"] - c0) < 1e-13 and r["vacuum_std"] < 1e-13 and abs(r["mean3d"] - c0) < 1e-13
    rp = r["radial"]
    assert 0 < rp["masked_fraction"] < 1 and abs(rp["r_max"] - 0.5 * np.linalg.norm(A[:, 0])) < 1e-12
    assert np.all(rp["rad"][np.isfinite(rp["rad"])] - c0 < 1e-13) and rp["cnt"].sum() > 0
    zp = r["z"]
    assert zp["n_vacuum"] > 0 and zp["c"] == 16.0 and np.abs(zp["z"]).max() <= 8.0 + 1e-12
    # the plateau rule and the sub-functions agree with the assembled result
    ko = al.kumagai_oba_constant(V_p + c0, V_p, x_d, x_p, A)
    assert ko["C_N"] == r["C_N"] and ko["n"] == r["C_N_n"] and ko["d_min"] == 0.75 * ko["r_max"]
    assert np.all(ko["shifts"]["dist"][ko["mask"]] >= ko["d_min"])


def test_z_profile_requires_a3_along_z_and_some_vacuum():
    A, x = _graphene_supercell(2, c=6.0)
    dV = np.zeros((8, 8, 12))
    with pytest.raises(ValueError):
        al.z_profile(dV, A, x, z_vac=5.0)                                     # c / 2 = 3 A < 5 A: no vacuum plane
    A_tilt = A.copy()
    A_tilt[0, 2] = 1.0
    with pytest.raises(ValueError):
        al.z_profile(dV, A_tilt, x, z_vac=1.0)


@pytest.mark.needs_data
@pytest.mark.skipif(not has_5x5, reason="5x5 supercell pair and Vks absent")
def test_header_geometry_matches_xml_and_get_pot(repo_cwd):
    for pot, save in ((POT[0], POT[2]), (POT[1], POT[3])):
        h = qe_io.read_filplot(pot)
        A_xml, _ = qe_io.get_A_volume(save)
        x_xml = qe_io.get_x_red(save)
        assert np.abs(h["A_cols"] - A_xml).max() < 1e-6
        assert np.abs((h["x_red"] - x_xml + 0.5) % 1.0 - 0.5).max() < 1e-8
        assert h["ngfft"] == tuple(int(v) for v in qe_io.get_ngfft(save)) and h["plot_num"] == 1
        V_pot, _ = qe_io.get_pot(pot, subtract_mean=False, to_hartree=True)
        assert np.array_equal(h["V"], V_pot.transpose(2, 1, 0))


@pytest.mark.needs_data
@pytest.mark.skipif(not has_5x5, reason="5x5 supercell pair and Vks absent")
def test_5x5_reproduces_config_thesis_table_and_ved_analysis(repo_cwd):
    """C_5 of config/production.json (bit), vacuum of tab:alignement (E), masked profile of results/M2_plateau/ved_analysis.npz
    and sampling_table.csv (r = 2, 3 A)."""
    cfg = load_production(verbose=False)
    hd = qe_io.read_filplot(POT[0])
    hp = qe_io.read_filplot(POT[1])
    r = al.potential_background(hd["V"] * HA2EV, hp["V"] * HA2EV, hd["x_red"], hp["x_red"], hd["A_cols"] * BOHR)
    assert abs(r["C_N"] - alignment_C(cfg, "5x5")) < 1e-9                   # measured 4e-15 eV
    assert abs(r["C_N_rms"] * 1e3 - 94.1) < 0.05 and r["C_N_n"] == 13
    assert abs(r["vacuum"] * 1e3 - (-10.9)) < 0.05 and abs(r["vacuum_std"] * 1e3 - 9.8) < 0.05
    RES = results_dir(cfg)
    npz = os.path.join(RES, "ved_analysis.npz")
    if os.path.exists(npz):
        V = np.load(npz)
        rp = r["radial"]
        assert np.allclose(rp["rc"], V["5x5_rc_masked"]) and np.allclose(rp["rad"], V["5x5_rad_masked"], atol=1e-9, equal_nan=True)
    table = os.path.join(RES, "sampling_table.csv")
    if os.path.exists(table):
        row = next(x for x in csv.DictReader(open(table)) if x["N"] == "5")
        assert abs(r["plane_at"][2.0] * 1e3 - float(row["Ved_radial_masked_2A_meV"])) < 0.05
        assert abs(r["plane_at"][3.0] * 1e3 - float(row["Ved_radial_masked_3A_meV"])) < 0.05


@pytest.mark.needs_data
@pytest.mark.skipif(not has_5x5, reason="5x5 supercell pair and Vks absent")
def test_script_gate_5x5(repo_cwd, tmp_path):
    out = tmp_path / "bg.json"
    cmd = [sys.executable, os.path.join(REPO, "scripts", "m", "potential_background.py"), POT[0], POT[1], "--gate", "5x5", "--json", str(out),
           "--npz", str(tmp_path / "bg.npz")]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO)
    assert p.returncode == 0, p.stdout + p.stderr
    s = json.load(open(out))
    assert s["gate"]["pass_"] and s["units"] == "eV" and s["n_atoms"] == {"d": 49, "p": 50}
    assert abs(s["C_N"]["value"] - alignment_C(load_production(verbose=False), "5x5")) < 1e-9
    z = np.load(tmp_path / "bg.npz")
    assert z["rad_eV"].shape == z["rc_A"].shape and z["dist_A"].shape == (49,) and z["plateau_mask"].sum() == 13
