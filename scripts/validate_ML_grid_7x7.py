"""Reference M^L (7x7 coarse, bands 0-3) with the serial kernel compute_ML_R (psi folded in G space on the supercell grid,
no commensurability assumption). Then compare with the MPI-shared kernel (fixed: Fourier-resampled V_ed) and with the old file."""
import sys, numpy as np
from electron_defect_interaction.io import qe_io
from electron_defect_interaction.defects.local_R import compute_ML_R
D = "data/graphene"; S = "7x7"
if sys.argv[1] == "ref":
    M = compute_ML_R(f"{D}/unit_cell/qe/defect_{S}.save", f"{D}/supercell/qe/defect_{S}_p.save", f"{D}/supercell/qe/defect_{S}_p.save/Vks_{S}_p",
                     f"{D}/supercell/qe/defect_{S}_d.save/Vks_{S}_d", subtract_mean=False, bands=[0, 1, 2, 3], io=qe_io)
    np.save("results/M/_ML_7x7_ref_b4.npy", M); print("ref saved", M.shape, "max", np.abs(M).max())
else:
    ref = np.load("results/M/_ML_7x7_ref_b4.npy"); new = np.load("results/M/M_L_dense_7x7_coarsecheck.npy"); old = np.load("results/M/M_L_7x7.npy")[:4, :, :4, :]
    for name, X in (("fixed MPI kernel (resampled 224)", new), ("OLD M_L_7x7 (ix % 30 on a 216 grid)", old)):
        print(f"{name}: max|X - ref| / max|ref| = {np.abs(X - ref).max() / np.abs(ref).max():.3e}; diag rel {np.abs(np.einsum('nknk->nk', X) - np.einsum('nknk->nk', ref)).max() / np.abs(np.einsum('nknk->nk', ref)).max():.3e}")
