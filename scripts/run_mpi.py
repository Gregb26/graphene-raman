import os
import numpy as np

from electron_defect_interaction.defects.local_G import prep_reciprocal_inputs, compute_ML_G_mpi
from mpi4py import MPI

def main():

    wfk_uc = "wfk_uc.nc"
    wfk_sc = "wfk_p_sc.nc"
    pot_p  = "pot_p_sc.nc"
    pot_d  = "pot_d_sc.nc"

    prep = prep_reciprocal_inputs(wfk_uc, wfk_sc, pot_p, pot_d, subtract_mean=False, pristine=False)

    ML_G = compute_ML_G_mpi(prep, block_size=128, show_tqdm=True)

    if MPI.COMM_WORLD.Get_rank() == 0:
        np.save("MLG.npy", ML_G)

if __name__ == "__main__":

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

    main()