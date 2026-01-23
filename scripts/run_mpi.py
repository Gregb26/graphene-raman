import os

# Set threads BEFORE importing numpy / scipy to avoid oversubscription.
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

    wfk_uc = "wfk_uc.nc"
    wfk_p  = "wfk_p.nc"   # pristine supercell WFK (or whatever your prep expects)
    wfk_d  = "wfk_d.nc"   # defect supercell WFK
    pot_p  = "pot_p.nc"
    pot_d  = "pot_d.nc"
    psp8   = "C.psp8"

    # Local part: all ranks participate, only rank 0 gets the full ML_G
    prep = prep_reciprocal_inputs(
        wfk_uc, wfk_d, pot_p, pot_d,
        subtract_mean=False, pristine=False
    )
    ML_G = compute_ML_G_mpi(prep, block_size=128, show_tqdm=(rank == 0))

    # Non-local part: compute ONLY on rank 0
    if rank == 0:
        # Replace with the correct difference for your use-case.
        # Example: defect - pristine
        M_NL = compute_M_NL(wfk_uc, wfk_d, psp8) - compute_M_NL(wfk_uc, wfk_p, psp8)

        M = ML_G + M_NL
        np.save("M_test.npy", M)

    # Ensure clean exit for non-root ranks
    comm.Barrier()


if __name__ == "__main__":
    main()