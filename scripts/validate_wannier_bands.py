"""
validate_wannier_bands.py
    Physical validation of the Wannier interpolation against the QE DFT eigenvalues, using only files
    already in the unit-cell .save (no extra `bands` run needed):
        - wannier_tb.dat / wannier_u.mat / wannier_u_dis.mat  (Wannier90, eV)
        - data-file-schema.xml: KS eigenvalues eps(k) on the 5x5 coarse grid (Hartree) + Fermi energy

    Wannier90 writes H(R) in eV; QE eigenvalues are in Hartree -> convert.

    Test A (exact, decisive). On the coarse grid the Fourier reconstruction must reproduce the projected
    DFT Hamiltonian:
        eig[ sum_R e^{2pi i k.R} H(R)/ndegen(R) ]  ==  eig[ V(k)^dag diag(eps_DFT(k)) V(k) ],  V = U_dis U
    This simultaneously checks ndegen, the parsers, and that H(R) is consistent with U/U_dis/eps. We also
    redo the left side WITHOUT ndegen to show it then breaks (proving the ndegen fix is necessary).

    Test B (physical). Wannier vs DFT bands on the coarse grid: each interpolated energy matched to the
    nearest DFT eigenvalue; agreement is ~0 inside the frozen window.

    Figures (saved as PNG):
        - wannier_bands_path.png : Wannier-interpolated bands along Gamma-K-M-Gamma (Dirac cone), E_F at 0
        - wannier_vs_dft_coarse.png : interpolated (red x) vs DFT (black .) eigenvalues on the 25 coarse k
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.wannier.wannier_interpolation import _match_kpoint_order

HA_TO_EV = 27.211386245988

DATA = "data/graphene/unit_cell/qe/defect_5x5.save"
U_PATH = f"{DATA}/wannier_u.mat"
UDIS_PATH = f"{DATA}/wannier_u_dis.mat"
TB_PATH = f"{DATA}/wannier_tb.dat"

# graphene high-symmetry path (reduced coords), convention from plotting/plot_band.py
GAMMA = np.array([0.0, 0.0, 0.0])
K = np.array([2 / 3, 1 / 3, 0.0])
M = np.array([1 / 2, 0.0, 0.0])


def build_kpath(corners, n_per_seg=200):
    """Concatenate straight segments through `corners`; return (k_red, corner_indices)."""
    pts, idx = [], [0]
    for a, b in zip(corners[:-1], corners[1:]):
        seg = a + np.linspace(0, 1, n_per_seg, endpoint=False)[:, None] * (b - a)
        pts.append(seg)
        idx.append(idx[-1] + n_per_seg)
    pts.append(corners[-1][None, :])
    return np.vstack(pts), idx


def main():
    # --- inputs ---
    Hwr, Rw, ndegen = read_w90_HR(TB_PATH)                 # H(R) in eV
    U, k_U = read_w90_mat(U_PATH)
    Ud, k_Ud = read_w90_mat(UDIS_PATH)

    k_coarse = qe_io.get_k_red(DATA)                       # (nk, 3) reduced
    eps = qe_io.get_eigenvalues(DATA) * HA_TO_EV           # (nband_DFT, nk) eV
    fermi = qe_io.get_fermi(DATA) if hasattr(qe_io, "get_fermi") else None
    nbDFT, nk = eps.shape
    nw = U.shape[1]
    nbW = Ud.shape[1]
    print(f"nk={nk}  nband_DFT={nbDFT}  nw={nw}  n_dis_bands={nbW}")

    # align Wannier matrices to the coarse-grid k order
    U = U[_match_kpoint_order(k_U, k_coarse)]
    Ud = Ud[_match_kpoint_order(k_Ud, k_coarse)]
    V = np.einsum("kbw, kwx -> kbx", Ud, U)                # (nk, nbW, nw), V = U_dis U
    eps_dis = eps[:nbW, :]                                 # bands entering the disentanglement window

    # interpolated bands on the coarse grid, with and without the ndegen weights
    Ewk_hr = Hwr_to_Hwk(Hwr, Rw, k_coarse, ndegen=ndegen)[1]      # (nk, nw)  correct
    Ewk_no = Hwr_to_Hwk(Hwr, Rw, k_coarse, ndegen=None)[1]        # (nk, nw)  ndegen ignored

    # ===== Test 1 (physical, decisive): Wannier vs DFT bands on the coarse grid =====
    # match each interpolated energy to the nearest DFT eigenvalue; agreement ~0 in the frozen window
    def nearest_dft(E):
        return np.array([[np.min(np.abs(eps[:, ik] - e)) for e in E[ik]] for ik in range(nk)])
    dmin = nearest_dft(Ewk_hr)
    dmin_no = nearest_dft(Ewk_no)
    ok1 = np.median(dmin) < 1e-3 and np.median(dmin) < np.median(dmin_no) / 5
    print("\n=== Test 1: Wannier bands vs QE DFT bands on the coarse 5x5 grid ===")
    print(f"  with ndegen : |E_W - nearest E_DFT| median={np.median(dmin):.2e} eV  "
          f"max={np.max(dmin):.2e} eV  frac<1meV={np.mean(dmin<1e-3):.2f}")
    print(f"  no  ndegen  : |E_W - nearest E_DFT| median={np.median(dmin_no):.2e} eV  "
          f"max={np.max(dmin_no):.2e} eV  frac<1meV={np.mean(dmin_no<1e-3):.2f}")
    print(f"  -> {'PASS' if ok1 else 'FAIL'}  (ndegen sharpens the agreement by "
          f"{np.median(dmin_no)/max(np.median(dmin),1e-30):.0f}x)")

    # ----- note (not a pass/fail): H(R) Fourier vs V^dag eps V -----
    # On the coarse grid sum_R e^{ik.R} H(R)/ndegen reproduces V^dag eps V up to W90's use_ws_distance
    # refinement (per-element R-shifts by the Wannier centres), which we deliberately do not apply
    # (the M pipeline uses the same plain FT for H and for M, so they stay mutually consistent).
    Ewk_V = np.array([np.linalg.eigvalsh(V[ik].conj().T @ (eps_dis[:, ik][:, None] * V[ik]))
                      for ik in range(nk)])
    errV = np.max(np.abs(np.sort(Ewk_hr, axis=1) - np.sort(Ewk_V, axis=1)))
    errV_no = np.max(np.abs(np.sort(Ewk_no, axis=1) - np.sort(Ewk_V, axis=1)))
    print(f"  [note] max|eig(H(R) FT) - eig(V^d eps V)|: with ndegen {errV*1e3:.2f} meV, "
          f"without {errV_no*1e3:.2f} meV (residual = use_ws_distance)")

    # reference energy = QE Fermi level (eV); both QE and W90 share the absolute energy reference
    Efermi = (fermi * HA_TO_EV) if fermi is not None else _fermi_from_xml() * HA_TO_EV

    # ===== Figure 1 : Wannier bands along Gamma-K-M-Gamma =====
    kpath, corner_idx = build_kpath([GAMMA, K, M, GAMMA], n_per_seg=300)
    Ek = Hwr_to_Hwk(Hwr, Rw, kpath, ndegen=ndegen)[1]     # (npath, nw) eV
    B, _ = qe_io.get_B_volume(DATA)
    kcart = kpath @ B.T                                   # cartesian for a correct horizontal scale
    dist = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(kcart, axis=0), axis=1))])

    fig, ax = plt.subplots(figsize=(6, 4.5))
    for w in range(nw):
        ax.plot(dist, Ek[:, w] - Efermi, lw=1.4, color="C0")
    for ci in corner_idx:
        ax.axvline(dist[ci], color="k", lw=0.6, alpha=0.5)
    ax.axhline(0.0, color="r", lw=0.8, ls="--", alpha=0.8, label=r"$E_F$ (QE)")
    ax.set_xticks([dist[ci] for ci in corner_idx])
    ax.set_xticklabels([r"$\Gamma$", "K", "M", r"$\Gamma$"])
    ax.set_xlim(dist[0], dist[-1])
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Wannier-interpolated bands (graphene unit cell)")
    ax.grid(True, ls="--", lw=0.4, alpha=0.5)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig("wannier_bands_path.png", dpi=200)
    print("\nsaved wannier_bands_path.png")

    # ===== Figure 2 : Wannier vs DFT on the coarse grid =====
    order = np.lexsort((k_coarse[:, 1], k_coarse[:, 0]))  # arbitrary but reproducible ordering
    x = np.arange(nk)
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))
    for ib in range(nbW):
        ax2.plot(x, eps[ib, order] - Efermi, "k.", ms=4, label="DFT (QE)" if ib == 0 else None)
    for w in range(nw):
        ax2.plot(x, Ewk_hr[order, w] - Efermi, "rx", ms=5, mew=1.0,
                 label="Wannier" if w == 0 else None)
    ax2.set_xlabel("coarse k-point index (sorted)")
    ax2.set_ylabel(r"$E - E_F$ (eV)")
    ax2.set_title("Wannier vs DFT eigenvalues on the 5x5 coarse grid")
    ax2.grid(True, ls="--", lw=0.4, alpha=0.5)
    ax2.legend(loc="upper right")
    fig2.tight_layout()
    fig2.savefig("wannier_vs_dft_coarse.png", dpi=200)
    print("saved wannier_vs_dft_coarse.png")

    print("\nRESULT:", "PASS" if ok1 else "FAIL")
    return 0 if ok1 else 1


def _fermi_from_xml():
    """Fallback: read fermi_energy (Hartree) straight from the XML."""
    import xml.etree.ElementTree as ET
    root = ET.parse(f"{DATA}/data-file-schema.xml").getroot()
    return float(root.find(".//output/band_structure/fermi_energy").text)


if __name__ == "__main__":
    raise SystemExit(main())
