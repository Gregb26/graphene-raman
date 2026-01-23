import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
from mpi4py import MPI
from electron_defect_interaction.defects.local_G import (
    prep_reciprocal_inputs,
    compute_ML_G_mpi,
)
from electron_defect_interaction.defects.non_local import compute_M_NL


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Local part: all ranks participate, only rank 0 gets the full ML_G
    prep = None
    if rank == 0:
        wfk_uc = "wfk_uc.nc"
        wfk_p  = "wfk_p.nc"   
        wfk_d  = "wfk_d.nc"   
        pot_p  = "pot_p.nc"
        pot_d  = "pot_d.nc"
        psp8   = "C.psp8"
        prep = prep_reciprocal_inputs(
        wfk_uc, wfk_d, pot_p, pot_d,
        subtract_mean=False, pristine=False
    )
    prep = comm.bcast(prep, root=0)

    ML_G = compute_ML_G_mpi(prep, block_size=128, show_tqdm=(rank == 0))

    # Non-local part: compute only on rank 0
    if rank == 0:
        M_NL = compute_M_NL(wfk_uc, wfk_d, psp8) - compute_M_NL(wfk_uc, wfk_p, psp8)
        M = ML_G + M_NL
        np.save("M_test.npy", M)

    comm.Barrier()


if __name__ == "__main__":
    main()