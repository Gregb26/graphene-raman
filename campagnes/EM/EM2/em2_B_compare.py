"""
EM2, step B.4-B.5: postw90 3.1.0 `kubo` against the production sigma(omega) of M4 (campagnes/EM/M4_sigma/).

postw90 writes seedname-kubo_S_ab.dat: omega (eV), Re sigma^H_ab (absorptive, symmetric part), Im sigma^AH_ab, in S/cm
for the 3D cell (berry.F90: fac = 1e8 e^2 / (hbar V_c), V_c in A^3). Conversion:
    sigma_2D = sigma_3D * 100 * c   (S/cm -> S/m, times c in m),   sigma / sigma_0 with sigma_0 = e^2 / (4 hbar),
with e and hbar of the postw90 build (CODATA 2006: its bohr is 0.52917720859 A, see em2_kpoints.log), so that
sigma/sigma_0 does not depend on the constants. postw90 has NO spin factor (Kubo sum over spinless Wannier states,
berry.tex eq. sig-H, berry.F90): the thesis convention (M4, kubo.py) includes g_s = 2, applied here once.

Usage: python em2_B_compare.py DIR [DIR ...]   (postw90/DIR/, e.g. k301 k301_nows k1201)
The first DIR is the reference of the pairwise differences; em2_postw90_sigma.npz (hw, sigma (nw, 3, 3) in
sigma/sigma_0, format of scripts/fig/make_figures_em.py) is written from k1201 (default postw90 settings) only
when k1201 is among the DIRs. Also compares each DIR with em2_B_ours.npz (same grid and prefactor) when it exists.
"""

import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
M4 = REPO / "campagnes" / "EM" / "M4_sigma"
E_SI, HBAR_SI = 1.602176487e-19, 1.054571628e-34      # w90_constants, CODATA 2006 block
BOHR_W90 = 0.52917720859                               # A, same block
C_BOHR = 30.0                                          # |a3| in bohr (unit_cell_cart of wannier.win)
G_S = 2                                                # spin degeneracy, absent from postw90
SIGMA0 = E_SI**2 / (4 * HBAR_SI)                       # S
LASERS = [1.96, 2.33, 2.54]
PAIRS = {"xx": (0, 0), "yy": (1, 1), "zz": (2, 2), "xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}


def read_postw90(d):
    """hw (nw,), sigma/sigma_0 (nw, 3, 3) from the symmetric kubo_S files of postw90/d, spin included."""
    sig = None
    for ab, (i, j) in PAIRS.items():
        data = np.loadtxt(HERE / "postw90" / d / f"wannier-kubo_S_{ab}.dat")
        if sig is None:
            hw = data[:, 0]; sig = np.zeros((len(hw), 3, 3))
        assert np.array_equal(data[:, 0], hw)
        s2d = data[:, 1] * 100 * C_BOHR * BOHR_W90 * 1e-10 * G_S / SIGMA0
        sig[:, i, j] = s2d; sig[:, j, i] = s2d
    return hw, sig


def at(hw, e):
    return int(np.argmin(np.abs(hw - e)))


def main():
    m4 = np.load(M4 / "em_sigma_full_N1200_eta0.04.npz")
    hw_m4, sig_m4 = m4["hw"], m4["sigma"]
    assert float(m4["eta"]) == 0.04 and int(m4["N"]) == 1200

    dirs = sys.argv[1:]
    res = {d: read_postw90(d) for d in dirs}
    out = [f"EM2 B — postw90 kubo vs M4 (full, N = 1200 shifted, eta = 0.04 eV); sigma_0 = e^2/4hbar = {SIGMA0:.6e} S "
           f"(postw90 constants), c = {C_BOHR * BOHR_W90:.7f} A, g_s = {G_S} applied to postw90", ""]
    for d, (hw, sig) in res.items():
        assert np.abs(hw - hw_m4).max() < 1e-6, f"{d}: frequency grids differ"
        diff = sig - sig_m4
        iL = [at(hw, e) for e in LASERS]; i02 = at(hw, 0.2)
        ivh, ivh4 = int(np.argmax(sig[:, 0, 0])), int(np.argmax(sig_m4[:, 0, 0]))
        out.append(f"[{d}]")
        out.append(f"  sigma_xx/sigma_0 at 0.2 eV: postw90 {sig[i02, 0, 0]:.5f}, M4 {sig_m4[i02, 0, 0]:.5f}  "
                   f"(postw90 raw, no spin: {sig[i02, 0, 0] / G_S:.5f})")
        out.append("  hw (eV) | xx postw90  M4       diff      | yy postw90  M4       diff      | xy postw90  M4")
        for e, i in zip(LASERS, iL):
            out.append(f"  {e:5.2f}   | {sig[i, 0, 0]:.5f}  {sig_m4[i, 0, 0]:.5f}  {diff[i, 0, 0]:+.2e} | "
                       f"{sig[i, 1, 1]:.5f}  {sig_m4[i, 1, 1]:.5f}  {diff[i, 1, 1]:+.2e} | {sig[i, 0, 1]:+.1e}  {sig_m4[i, 0, 1]:+.1e}")
        out.append(f"  van Hove peak (xx): postw90 {hw[ivh]:.2f} eV, {sig[ivh, 0, 0]:.4f}; M4 {hw_m4[ivh4]:.2f} eV, {sig_m4[ivh4, 0, 0]:.4f}")
        for ab in ("xx", "yy", "xy"):
            i, j = PAIRS[ab]
            a = np.abs(diff[:, i, j]); k = int(np.argmax(a))
            rel = np.abs(diff[:, i, j]) / np.abs(sig_m4[:, i, j]).clip(min=1e-3) if ab != "xy" else None
            s = f"  whole curve {ab}: max |diff| {a.max():.2e} at {hw[k]:.2f} eV, mean |diff| {a.mean():.2e}"
            if rel is not None:
                kr = int(np.argmax(rel)); s += f", max rel {rel.max():.2e} at {hw[kr]:.2f} eV"
            out.append(s)
        out.append(f"  postw90 |yy - xx| max {np.abs(sig[:, 1, 1] - sig[:, 0, 0]).max():.1e} (M4 {np.abs(sig_m4[:, 1, 1] - sig_m4[:, 0, 0]).max():.1e}); "
                   f"|xy| max {np.abs(sig[:, 0, 1]).max():.1e}; zz max {np.abs(sig[:, 2, 2]).max():.1e}")
        out.append("")
    ours = np.load(HERE / "em2_B_ours.npz") if (HERE / "em2_B_ours.npz").exists() else None
    if ours is not None:
        out.append("Against our sigma with the postw90 conventions (em2_B_ours.py: prefactor 1/Delta, Gamma-centred N x N grid;")
        out.append("left: use_ws_distance and r(R) rebuilt by postw90 from the .mmn):")
        for d, (hw, sig) in res.items():
            N = int(open(HERE / "postw90" / d / "wannier.win").read().split("berry_kmesh =")[1].split()[0])
            for key in (f"gam{N}", "shift_delta"):
                if key not in ours.files:
                    continue
                dd = sig - ours[key]
                k = int(np.argmax(np.abs(dd[:, 0, 0])))
                out.append(f"  {d} - ours[{key}]: max |diff| xx {np.abs(dd[:, 0, 0]).max():.2e} at {hw[k]:.2f} eV, "
                           f"yy {np.abs(dd[:, 1, 1]).max():.2e}, xy {np.abs(dd[:, 0, 1]).max():.1e}; mean |diff| xx "
                           f"{np.abs(dd[:, 0, 0]).mean():.1e}; at 0.2 eV and lasers xx "
                           + ", ".join(f"{dd[at(hw, e), 0, 0]:+.1e}" for e in [0.2] + LASERS))
        out.append("")
    ref = dirs[0]
    for d in dirs[1:]:
        dd = res[d][1] - res[ref][1]
        out.append(f"[{d} - {ref}] max |diff| xx {np.abs(dd[:, 0, 0]).max():.2e}, yy {np.abs(dd[:, 1, 1]).max():.2e}, "
                   f"xy {np.abs(dd[:, 0, 1]).max():.2e}; at the lasers xx "
                   + ", ".join(f"{dd[at(res[d][0], e), 0, 0]:+.1e}" for e in LASERS))
    text = "\n".join(out) + "\n"
    print(text)
    tag = "_".join(dirs)
    (HERE / f"em2_B_compare_{tag}.txt").write_text(text)

    if "k1201" in dirs:          # the postw90 curve of the figure: default settings, production grid
        hw, sig = res["k1201"]
        np.savez(HERE / "em2_postw90_sigma.npz", hw=hw, sigma=sig, berry_kmesh="k1201", eta=0.04, width_w90=0.0565685425,
                 mu=-4.238895, g_s=G_S, sigma0_SI=SIGMA0, c_A=C_BOHR * BOHR_W90, source="postw90/k1201/wannier-kubo_S_*.dat")
        print("wrote em2_postw90_sigma.npz from k1201")


if __name__ == "__main__":
    main()
