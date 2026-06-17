"""
run_mpi.py
    MPI driver: compute the electron-defect scattering matrix M = M^L + M^NL from Quantum ESPRESSO outputs.

    The local part M^L is computed in reciprocal space with local_G (memory-light, no folded real-space
    grid) and parallelised over (k', k) blocks across MPI ranks. The non-local part M^NL is cheap and is
    computed on rank 0. Run with e.g.:

        mpirun -n 64 python scripts/run_mpi.py

    For ABINIT inputs, drop the io=qe_io / pseudo_reader=read_upf arguments (defaults use the ABINIT
    backend) and pass the corresponding .nc / .psp8 paths.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
from mpi4py import MPI

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.pseudo_io import read_upf
from electron_defect_interaction.defects.local_G import (
    prep_reciprocal_inputs,
    compute_ML_G_mpi,
)
from electron_defect_interaction.defects.non_local import compute_M_NL


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # --- Paths (QE) ---
    uc_save   = "data/graphene/unit_cell/qe/defect_5x5.save"
    sc_p_save = "data/graphene/supercell/qe/defect_5x5_p.save"
    sc_d_save = "data/graphene/supercell/qe/defect_5x5_d.save"
    pot_p     = f"{sc_p_save}/Vks_5x5_p"
    pot_d     = f"{sc_d_save}/Vks_5x5_d"
    upf       = f"{uc_save}/C.upf"

    BANDS = None   # None = all bands (reciprocal-space path is memory-light)

    # --- Local part: rank 0 builds the inputs (reads files), then broadcasts ---
    prep = None
    if rank == 0:
        prep = prep_reciprocal_inputs(
            uc_save, sc_p_save, pot_p, pot_d,
            subtract_mean=False, bands=BANDS, io=qe_io,
        )
    prep = comm.bcast(prep, root=0)

    # All ranks cooperate; only rank 0 receives the assembled M^L
    M_L = compute_ML_G_mpi(prep, block_size=128, show_tqdm=(rank == 0))

    # --- Non-local part: compute on rank 0 only ---
    if rank == 0:
        M_NL = compute_M_NL(
            uc_save, sc_p_save, sc_d_save, upf,
            io=qe_io, pseudo_reader=read_upf,
        )
        if BANDS is not None:
            M_NL = M_NL[BANDS][:, :, BANDS]

        M = M_L + M_NL
        np.save("M_ed.npy", M)
        print(f"Done. M shape={M.shape}, saved to M_ed.npy")

    comm.Barrier()


if __name__ == "__main__":
    main()
