#!/usr/bin/env python
"""
compute_M_cluster.py
    Cluster driver: compute the full electron-defect scattering matrix M = M^L + M^NL from
    Quantum ESPRESSO outputs, using the MPI reciprocal-space local part.

    The local part M^L is parallelised over (k', k) blocks (compute_ML_G_mpi); the non-local part
    M^NL is computed on rank 0. Rank 0 writes the assembled M to --out.

    Inputs (all from a finished QE run):
        --uc      unit-cell prefix.save dir (wavefunctions on the full k-grid: use nosym+noinv)
        --sc-p    pristine supercell prefix.save dir
        --sc-d    defective supercell prefix.save dir
        --pot-p   pristine supercell local potential (pp.x plot_num=1 filplot)
        --pot-d   defective supercell local potential
        --upf     UPF pseudopotential

    Example SLURM launch (mpi4py built against MPICH -> use the MPICH mpirun/srun):

        srun -n 256 python scripts/compute_M_cluster.py \
            --uc   graphene/uc.save \
            --sc-p graphene/sc_p.save  --pot-p graphene/sc_p.save/Vks_p \
            --sc-d graphene/sc_d.save  --pot-d graphene/sc_d.save/Vks_d \
            --upf  graphene/uc.save/C.upf \
            --out  M_ed.npy

    For ABINIT inputs use --backend abinit and pass .nc / .psp8 paths instead.
"""

import os

# Pin BLAS to one thread per rank: the parallelism is over MPI ranks, not threads.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import numpy as np
from mpi4py import MPI

from electron_defect_interaction.defects.local_G import prep_reciprocal_inputs, compute_ML_G_mpi
from electron_defect_interaction.defects.local_R import prep_realspace_inputs, compute_ML_R_mpi
from electron_defect_interaction.defects.non_local import compute_M_NL_mpi


def parse_args():
    p = argparse.ArgumentParser(description="Compute M = M^L + M^NL (MPI) from QE/ABINIT outputs.")
    p.add_argument("--uc", required=True, help="unit-cell wavefunctions (.save dir for QE, WFK.nc for ABINIT)")
    p.add_argument("--sc-p", required=True, help="pristine supercell")
    p.add_argument("--sc-d", required=True, help="defective supercell")
    p.add_argument("--pot-p", required=True, help="pristine supercell local potential")
    p.add_argument("--pot-d", required=True, help="defective supercell local potential")
    p.add_argument("--upf", required=True, help="pseudopotential (.upf for QE, .psp8 for ABINIT)")
    p.add_argument("--out", default="M_ed.npy", help="output .npy file for M (default: M_ed.npy)")
    p.add_argument("--backend", choices=["qe", "abinit"], default="qe")
    p.add_argument("--bands", default="all", help="comma-separated band indices, or 'all' (default)")
    p.add_argument("--method", choices=["real", "reciprocal"], default="real",
                   help="M^L method: 'real' (grid-distributed real space, fast for moderate supercells) "
                        "or 'reciprocal' (better asymptotic scaling for very large supercells)")
    p.add_argument("--block-size", type=int, default=512,
                   help="reciprocal: G' block size; real: grid block size per rank")
    return p.parse_args()


def main():
    args = parse_args()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Backend selection
    if args.backend == "qe":
        from electron_defect_interaction.io import qe_io as io
        from electron_defect_interaction.io.pseudo_io import read_upf as pseudo_reader
    else:
        from electron_defect_interaction.io import abinit_io as io
        from electron_defect_interaction.io.pseudo_io import read_psp8 as pseudo_reader

    bands = None if args.bands == "all" else [int(b) for b in args.bands.split(",")]

    # --- Local part: rank 0 reads files and builds inputs, then broadcasts ---
    prep = None
    if rank == 0:
        print(f"[rank0] backend={args.backend}, method={args.method}, bands={args.bands}, "
              f"nranks={comm.Get_size()}", flush=True)
        prep_fn = prep_realspace_inputs if args.method == "real" else prep_reciprocal_inputs
        prep = prep_fn(args.uc, args.sc_p, args.pot_p, args.pot_d,
                       subtract_mean=False, bands=bands, io=io)
    prep = comm.bcast(prep, root=0)

    if args.method == "real":
        # grid-distributed real space; Allreduce -> M^L on every rank
        M_L = compute_ML_R_mpi(prep, grid_block=args.block_size if args.block_size > 1000 else 200_000)
    else:
        # reciprocal over (k', k) blocks; assembled M^L on rank 0 only
        M_L = compute_ML_G_mpi(prep, block_size=args.block_size, show_tqdm=(rank == 0))

    # --- Non-local part: distributed over k' across ranks (returns on every rank) ---
    if rank == 0:
        print("[rank0] computing M^NL (MPI) ...", flush=True)
    M_NL = compute_M_NL_mpi(
        args.uc, args.sc_p, args.sc_d, args.upf,
        io=io, pseudo_reader=pseudo_reader, bands=bands,
    )

    # --- Assemble and save on rank 0 ---
    if rank == 0:
        M = M_L + M_NL
        np.save(args.out, M)

        # quick sanity: M must be Hermitian as an operator over (band, k)
        nb, nk, _, _ = M.shape
        Op = M.reshape(nb * nk, nb * nk)
        herm = np.max(np.abs(Op - Op.conj().T))
        print(f"[rank0] saved {args.out}  shape={M.shape}  max|M-M^dag|={herm:.2e}", flush=True)

    comm.Barrier()


if __name__ == "__main__":
    main()
