"""
run.py
    Serial driver: compute the electron-defect scattering matrix M = M^L + M^NL from Quantum ESPRESSO
    outputs and save it to disk.

    Inputs are QE `prefix.save/` directories plus pp.x local-potential dumps (plot_num=1 filplot) and the
    UPF pseudopotential.

    Note on memory: compute_ML_R builds the unit-cell wavefunctions folded onto the full supercell real-space
    grid, which is large for big supercells. Restrict `BANDS` accordingly, or use the MPI reciprocal-space
    driver run_mpi.py for the local part on many bands.
"""

import numpy as np

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.pseudo_io import read_upf
from electron_defect_interaction.defects.local_R import compute_ML_R
from electron_defect_interaction.defects.non_local import compute_M_NL


def main():
    # --- Paths (QE) ---
    uc_save   = "data/graphene/unit_cell/qe/defect_5x5.save"        # unit-cell wavefunctions
    sc_p_save = "data/graphene/supercell/qe/defect_5x5_p.save"      # pristine supercell
    sc_d_save = "data/graphene/supercell/qe/defect_5x5_d.save"      # defective supercell
    pot_p     = f"{sc_p_save}/Vks_5x5_p"                            # pristine local potential (pp.x plot_num=1)
    pot_d     = f"{sc_d_save}/Vks_5x5_d"                            # defective local potential
    upf       = f"{uc_save}/C.upf"                                  # pseudopotential

    # Band subset to compute (keep small for the real-space local part; None = all bands)
    BANDS = [2, 3, 4, 5]

    # --- Local part M^L (real space) ---
    print("Computing local part M^L ...")
    M_L = compute_ML_R(
        uc_save, sc_p_save, pot_p, pot_d,
        subtract_mean=False, bands=BANDS,
        io=qe_io,
    )

    # --- Non-local part M^NL (all bands, then restrict to the same subset) ---
    print("Computing non-local part M^NL ...")
    M_NL_full = compute_M_NL(
        uc_save, sc_p_save, sc_d_save, upf,
        io=qe_io, pseudo_reader=read_upf,
    )
    M_NL = M_NL_full[BANDS][:, :, BANDS]   # (nbands_sel, nk, nbands_sel, nk)

    # --- Total matrix ---
    M = M_L + M_NL

    np.save("M_ed.npy", M)
    print(f"Done. M shape={M.shape}, saved to M_ed.npy")


if __name__ == "__main__":
    main()
