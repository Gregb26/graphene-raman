"""
EM2, step A.4: direct DFT (QE 7.5 pw.x 'bands' + bands.x lp) against the Wannier velocity, k by k.

Inputs (read only):
- em2_kpoints_table.txt (em2_kpoints.py), in the order of the K_POINTS list;
- QE eigenvalues from the XML of the bands run (full precision, Hartree), outdir of bands.in;
- em2_p_avg.dat (bands.x, filp): |<c|p_a|v>|^2 for a = x, y, z (Cartesian axes of the QE cell), c = bands
  nbnd_occ+1..nbnd (rows), v = 1..nbnd_occ (columns), in (hbar/bohr)^2. p = (i/2)[H, x] in Rydberg units, i.e.
  m_e v with the commutator [V_NL, r] included (PP/src/compute_ppsi.f90, PW/src/commutator_Hx_psi.f90), so
  hbar v = (hbar^2/m_e) p = 2 Ry bohr p = P_TO_HV p.
- Wannier: compute_velocity on wannier/27x27/wannier_tb.dat, three variants (full, centres_only, no_berry).

pi and pi* are chosen by energy around E_D (v: highest band below E_D, c: lowest above), in QE and in Wannier.
Compares eps_pi, eps_pi* (meV) and the in-plane |hbar v_cv| = sqrt(|hbar v^x_cv|^2 + |hbar v^y_cv|^2) (gauge
invariant, no node). Dirac circle (q = 0.005 1/A): hbar v_F^DFT = <(eps_pi* - eps_pi)/(2 q)>_theta against
<|hbar v_cv|>_theta from p (units check independent of Wannier) and against 5.469 eV A (Wannier).

Writes em2_A_table.txt (per k), em2_A_summary.txt, em2_A.npz. Run: python em2_A_compare.py
"""

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from electron_defect_interaction.electron_photon import (centres_only, compute_velocity, fermi_velocity,
                                                         make_grid_tb, make_wannier_tb)

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
TB = REPO / "wannier" / "27x27" / "wannier_tb.dat"
E_D = -4.238895                     # eV
HA_EV = 27.211386245988             # QE 7.5 (CODATA 2018)
RY_EV = HA_EV / 2
BOHR_A = 0.529177210903
P_TO_HV = 2 * RY_EV * BOHR_A        # eV A per (hbar/bohr): hbar^2/m_e = 2 Ry bohr^2
Q_DIRAC = 0.005                     # 1/A
HV_F_W90 = 5.469187853771915        # eV A, campagnes/EM/M4_sigma/em_scalars.json (inter)


def outdir_from(infile):
    m = re.search(r"outdir\s*=\s*'([^']+)'", infile.read_text())
    p = re.search(r"prefix\s*=\s*'([^']+)'", infile.read_text())
    return Path(m.group(1)) / f"{p.group(1)}.save"


def read_qe_xml(save):
    """k (Cartesian, 2 pi/alat), eigenvalues (eV), alat (bohr), a_i (bohr, rows), Fermi energy (eV)."""
    root = ET.parse(save / "data-file-schema.xml").getroot()
    out = root.find("output")
    alat = float(out.find("atomic_structure").get("alat"))
    cell = np.array([[float(x) for x in out.find(f"atomic_structure/cell/a{i}").text.split()] for i in (1, 2, 3)])
    bs = out.find("band_structure")
    k, eps = [], []
    for ks in bs.findall("ks_energies"):
        k.append([float(x) for x in ks.find("k_point").text.split()])
        eps.append([float(x) for x in ks.find("eigenvalues").text.split()])
    ef = float(bs.find("fermi_energy").text) * HA_EV
    return np.array(k), np.array(eps) * HA_EV, alat, cell, ef


def read_p_avg(path, nk_expected):
    """xk (2 pi/alat), nbnd_occ, |p|^2 as (nk, 3, nbnd, nbnd) with entry [c, v] (0-based), NaN elsewhere."""
    lines = path.read_text().splitlines()
    head = re.match(r"\s*&p_mat nbnd=\s*(\d+), nks=\s*(\d+)", lines[0])
    nbnd, nks = int(head.group(1)), int(head.group(2))
    assert nks == nk_expected, (nks, nk_expected)
    xk = np.zeros((nks, 3)); nocc = np.zeros(nks, int)
    P2 = np.full((nks, 3, nbnd, nbnd), np.nan)
    i = 1
    for ik in range(nks):
        parts = lines[i].split(); xk[ik] = [float(x) for x in parts[:3]]; nocc[ik] = int(parts[3]); i += 1
        n_o = nocc[ik]
        for ipol in range(3):
            assert int(lines[i]) == ipol + 1; i += 1
            for c in range(n_o, nbnd):
                vals = []
                while len(vals) < n_o:            # 5 values per line (5f15.8)
                    vals += [float(x) for x in lines[i].split()]; i += 1
                P2[ik, ipol, c, :n_o] = vals
    assert i == len(lines), (i, len(lines))
    return xk, nocc, P2


def pi_pair(eps):
    """Indices (v, c) of pi and pi* by energy around E_D, per k."""
    n_below = np.sum(eps < E_D, axis=1)
    return n_below - 1, n_below


# --- k list ---
tab = np.loadtxt(HERE / "em2_kpoints_table.txt", comments="#", dtype=str)
sets = tab[:, 1]; hw = tab[:, 2].astype(float); theta = np.radians(tab[:, 3].astype(float))
k_red = tab[:, 4:7].astype(float)
nk = k_red.shape[0]

# --- QE ---
save = outdir_from(HERE / "bands.in")
xk_xml, eps_qe, alat, cell, ef = read_qe_xml(save)
xk_p, nocc, P2 = read_p_avg(HERE / "em2_p_avg.dat", nk)
bg = np.linalg.inv(cell).T * alat                       # rows b_j in 2 pi/alat
xk_list = k_red @ bg
assert np.abs(xk_xml - xk_list).max() < 1e-9, np.abs(xk_xml - xk_list).max()
assert np.abs(xk_p - xk_list).max() < 1e-6, np.abs(xk_p - xk_list).max()   # f10.6 in the p file
v_qe, c_qe = pi_pair(eps_qe)
assert np.all(v_qe == 3) and np.all(c_qe == 4), "pi/pi* are not QE bands 4/5 at every k"
assert np.all(nocc == 4), "nbnd_occ != 4: pi/pi* not the last valence/first conduction band of the p file"
p2_cv = P2[np.arange(nk), :, c_qe, v_qe]                 # (nk, 3), (hbar/bohr)^2
hv_qe = P_TO_HV * np.sqrt(p2_cv[:, 0] + p2_cv[:, 1])     # eV A, in-plane norm
hvz_qe = P_TO_HV * np.sqrt(p2_cv[:, 2])

# --- Wannier ---
tb = make_wannier_tb(TB)
K = make_grid_tb(tb, 3).K
B = 2 * np.pi * np.linalg.inv(tb.lattice).T              # columns b_j, 1/A
k_cart = k_red @ B.T
variants = {"full": (tb, "berry"), "centres_only": (centres_only(tb), "berry"), "no_berry": (tb, "no_berry")}
hv_w, dP, eps_w = {}, {}, None
for name, (model, mode) in variants.items():
    _, eps, _, hv = compute_velocity(model, k_cart, mode)
    v_w, c_w = pi_pair(eps)
    hv_cv = hv[np.arange(nk), :, c_w, v_w]               # (nk, 3)
    hv_w[name] = np.sqrt(np.abs(hv_cv[:, 0])**2 + np.abs(hv_cv[:, 1])**2)
    P_w = np.abs(hv_cv[:, :2])**2                         # (nk, 2) |hbar v^x_cv|^2, |hbar v^y_cv|^2, each gauge invariant
    dP[name] = (P_w - P_TO_HV**2 * p2_cv[:, :2]) / HV_F_W90**2
    if name == "full":
        eps_w = eps; vw, cw = v_w, c_w
de_v = (eps_qe[np.arange(nk), v_qe] - eps_w[np.arange(nk), vw]) * 1e3   # meV
de_c = (eps_qe[np.arange(nk), c_qe] - eps_w[np.arange(nk), cw]) * 1e3
gap_qe = eps_qe[np.arange(nk), c_qe] - eps_qe[np.arange(nk), v_qe]

# --- Dirac circle ---
d = sets == "dirac"
q_qe = np.linalg.norm((k_red[d] @ bg - np.array([2 / 3, 1 / 3, 0]) @ bg) * 2 * np.pi / (alat * BOHR_A), axis=1)
hvF_dft = np.mean(gap_qe[d] / (2 * q_qe))
hvF_p = np.mean(hv_qe[d])
fv = fermi_velocity(tb, K, E_D, q=Q_DIRAC, ntheta=360)

# --- outputs ---
rel = {n: hv_w[n] / hv_qe - 1 for n in hv_w}
with open(HERE / "em2_A_table.txt", "w") as f:
    f.write(f"# EM2 A: QE (bands + bands.x lp) vs Wannier, per k. P_TO_HV = {P_TO_HV:.10f} eV A per hbar/bohr\n")
    f.write("# ik set hw theta(deg) eps_pi^QE-E_D eps_pi*^QE-E_D (eV) d_eps_pi d_eps_pi* (meV, QE-W) "
            "|hv_cv|_QE |hv_cv|_full |hv_cv|_centres |hv_cv|_noberry (eV A) |hv^z_cv|_QE rel_full rel_centres rel_noberry\n")
    for i in range(nk):
        f.write(f"{i + 1:3d} {sets[i]:5s} {hw[i]:4.2f} {np.degrees(theta[i]):8.3f} "
                f"{eps_qe[i, v_qe[i]] - E_D:10.6f} {eps_qe[i, c_qe[i]] - E_D:10.6f} {de_v[i]:8.3f} {de_c[i]:8.3f} "
                f"{hv_qe[i]:9.5f} {hv_w['full'][i]:9.5f} {hv_w['centres_only'][i]:9.5f} {hv_w['no_berry'][i]:9.5f} "
                f"{hvz_qe[i]:.1e} {rel['full'][i]:+.5f} {rel['centres_only'][i]:+.5f} {rel['no_berry'][i]:+.5f}\n")

lines = [f"EM2 A — QE 7.5 (bands, bands.x lp) vs Wannier 27x27, {nk} k; E_D = {E_D} eV; QE ef (XML) = {ef:.6f} eV",
         f"p -> hbar v: x {P_TO_HV:.8f} eV A per hbar/bohr (2 Ry bohr, CODATA 2018 as in QE 7.5)",
         f"pi/pi* = QE bands 4/5 at all k (by energy), nbnd_occ = 4 at all k; max |p^z_cv| (as hbar v) = {hvz_qe.max():.1e} eV A",
         ""]
groups = [("ring", 1.96), ("ring", 2.33), ("ring", 2.54), ("dirac", 0.0)]
lines.append("set      hw  n | eps_pi QE-W (meV) min/mean/max | eps_pi* QE-W (meV) min/mean/max | |hv_cv|_QE (eV A) min-max")
for s, w in groups:
    m = (sets == s) & (hw == w)
    lines.append(f"{s:5s} {w:5.2f} {m.sum():2d} | {de_v[m].min():+7.2f} {de_v[m].mean():+7.2f} {de_v[m].max():+7.2f} | "
                 f"{de_c[m].min():+7.2f} {de_c[m].mean():+7.2f} {de_c[m].max():+7.2f} | {hv_qe[m].min():.4f}-{hv_qe[m].max():.4f}")
lines.append("")
lines.append("|hv_cv|_W / |hv_cv|_QE - 1 (%), min / mean / max")
lines.append("set      hw  n | full                     | centres only             | no Berry")
for s, w in groups:
    m = (sets == s) & (hw == w)
    cols = " | ".join(f"{100 * rel[n][m].min():+7.3f} {100 * rel[n][m].mean():+7.3f} {100 * rel[n][m].max():+7.3f}"
                      for n in ("full", "centres_only", "no_berry"))
    lines.append(f"{s:5s} {w:5.2f} {m.sum():2d} | {cols}")
lines.append("")
lines.append("Components (each gauge invariant): max over the set of | |hv^a_cv|^2_W - |hv^a_cv|^2_QE | / (hbar v_F)^2, a = x ; y")
lines.append("set      hw  n | full            | centres only    | no Berry")
for s, w in groups:
    m = (sets == s) & (hw == w)
    cols = " | ".join(f"{np.abs(dP[n][m, 0]).max():.1e} ; {np.abs(dP[n][m, 1]).max():.1e}" for n in ("full", "centres_only", "no_berry"))
    lines.append(f"{s:5s} {w:5.2f} {m.sum():2d} | {cols}")
lines.append("")
lines.append(f"Dirac circle q = {Q_DIRAC} 1/A (|k-K| QE frame {q_qe.min():.9f}-{q_qe.max():.9f}), 12 k:")
lines.append(f"  hbar v_F^DFT = <(eps_pi* - eps_pi)/(2q)> = {hvF_dft:.5f} eV A (spread {np.ptp(gap_qe[d] / (2 * q_qe)):.1e})")
lines.append(f"  <|hbar v_cv|>_p (bands.x lp)          = {hvF_p:.5f} eV A (spread {np.ptp(hv_qe[d]):.1e}); ratio p/DFT - 1 = {hvF_p / hvF_dft - 1:+.2e}")
lines.append(f"  Wannier on the same circle: fermi_velocity(q = {Q_DIRAC}) inter {fv['inter_avg']:.5f}, intra {fv['intra_avg']:.5f} eV A; "
             f"<|hv_cv|>_full at the 12 k = {hv_w['full'][d].mean():.5f}")
lines.append(f"  hbar v_F^DFT / 5.469188 (M4, q = 1e-3) - 1 = {hvF_dft / HV_F_W90 - 1:+.2e}")
(HERE / "em2_A_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))

np.savez(HERE / "em2_A.npz", k_red=k_red, sets=sets, hw=hw, theta=theta, eps_qe=eps_qe, p2_cv=p2_cv, hv_qe=hv_qe,
         eps_w=eps_w, **{f"hv_{n}": hv_w[n] for n in hv_w}, de_v=de_v, de_c=de_c, E_D=E_D, P_TO_HV=P_TO_HV,
         hvF_dft=hvF_dft, hvF_p=hvF_p)
