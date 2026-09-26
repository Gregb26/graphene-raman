import numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.config import load_production, results_dir
RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)
print(f"{'N':>4} {'N_cells':>8} {'max|M|':>10} {'||M||_F':>10} {'max|M|*Nc':>10} {'||M||_F/Nc':>11}")
for N in ["5x5","6x6","7x7","8x8"]:
    matrix_io.check_manifest(f"{RES}/M_ed_{N}.npy", require_normalization=matrix_io.M_NORM_V2); M=np.load(f"{RES}/M_ed_{N}.npy"); n=int(N.split("x")[0]); Nc=n*n
    print(f"{N:>4} {Nc:>8} {np.max(np.abs(M)):>10.4e} {np.linalg.norm(M):>10.4e} {np.max(np.abs(M))*Nc:>10.4e} {np.linalg.norm(M)/Nc:>11.4e}")
print("=== assert structurel: aligned_eigenvalues ===")
for N in ["5x5","11x11"]:
    uc=f"data/graphene/unit_cell/qe/defect_{N}.save"; matrix_io.check_manifest(f"{RES}/M_ed_{N}.npy", require_normalization=matrix_io.M_NORM_V2); M=np.load(f"{RES}/M_ed_{N}.npy")
    try:
        e=qe_io.aligned_eigenvalues(uc, nk_expected=M.shape[1], shift_Fermi=True)
        print(f"{N}: OK eps{e.shape} aligne sur M (nk={M.shape[1]})")
    except ValueError as ex:
        print(f"{N}: ASSERT DECLENCHE -> {ex}")
