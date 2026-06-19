#!/usr/bin/env python
"""
compute_M.py
    Single driver for the electron-defect scattering matrix M = M^L + M^NL from Quantum
    ESPRESSO (or ABINIT) outputs. Works for ANY supercell size via a two-stage split:

        stage 'ml'  (MPI, run under srun): compute the local part M^L (real-space,
                    grid-distributed over ranks) and save it to --ml-out.
        stage 'nl'  (serial, plain python, FRESH process): load M^L, compute the non-local
                    part M^NL serially, add them, save the full M to --out.

    The two stages MUST run as separate processes. compute_ML_R_mpi leaves a latent heap
    corruption for nk >= 81 (k-grid >= 9x9) that aborts the *next* heavy allocation in the
    same process (observed as "munmap_chunk(): invalid pointer" / "double free" inside
    compute_M_NL). A fresh process for stage 'nl' starts from a clean heap, so the split
    works identically for 5x5 and 20x20. See submit_M.sh for the canonical launch.

    Inputs (all from a finished QE run; use nosym+noinv for the unit-cell full k-grid):
        --uc      unit-cell prefix.save dir
        --sc-p    pristine supercell prefix.save dir
        --sc-d    defective supercell prefix.save dir   (stage nl only)
        --pot-p   pristine supercell local potential, pp.x plot_num=1 filplot (stage ml only)
        --pot-d   defective supercell local potential                          (stage ml only)
        --upf     UPF pseudopotential                                          (stage nl only)

    Examples:
        # stage 1 (MPI):
        srun -n 32 python scripts/compute_M.py --stage ml \
            --uc uc.save --sc-p sc_p.save --pot-p sc_p.save/Vks_p --pot-d sc_d.save/Vks_d \
            --ml-out results/M/M_L_9x9.npy
        # stage 2 (serial, fresh process):
        python scripts/compute_M.py --stage nl \
            --uc uc.save --sc-p sc_p.save --sc-d sc_d.save --upf uc.save/C.upf \
            --ml results/M/M_L_9x9.npy --out results/M/M_ed_9x9.npy

    For ABINIT inputs pass --backend abinit and .nc / .psp8 paths instead.
"""

import os

# Pin every threading layer to one thread BEFORE python imports BLAS: parallelism is over
# MPI ranks, not threads. Oversubscribing FlexiBLAS/OpenBLAS corrupts the heap.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("FLEXIBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Compute M = M^L + M^NL from QE/ABINIT outputs (two-stage).")
    p.add_argument("--stage", choices=["ml", "nl"], required=True,
                   help="'ml': MPI local part M^L (run under srun); "
                        "'nl': serial non-local part M^NL + combine (fresh process)")
    p.add_argument("--uc", required=True, help="unit-cell wavefunctions (.save for QE, WFK.nc for ABINIT)")
    p.add_argument("--sc-p", required=True, help="pristine supercell")
    p.add_argument("--sc-d", help="defective supercell (stage nl)")
    p.add_argument("--pot-p", help="pristine supercell local potential (stage ml)")
    p.add_argument("--pot-d", help="defective supercell local potential (stage ml)")
    p.add_argument("--upf", help="pseudopotential, .upf for QE / .psp8 for ABINIT (stage nl)")
    p.add_argument("--ml", help="precomputed M^L .npy to load (stage nl)")
    p.add_argument("--ml-out", help="output .npy for M^L (stage ml)")
    p.add_argument("--out", help="output .npy for the full M (stage nl)")
    p.add_argument("--backend", choices=["qe", "abinit"], default="qe")
    p.add_argument("--bands", default="all", help="comma-separated band indices, or 'all' (default)")
    p.add_argument("--block-size", type=int, default=50_000, help="real-space grid block per rank (stage ml)")
    return p.parse_args()


def get_io(backend):
    if backend == "qe":
        from electron_defect_interaction.io import qe_io as io
        from electron_defect_interaction.io.pseudo_io import read_upf as pseudo_reader
    else:
        from electron_defect_interaction.io import abinit_io as io
        from electron_defect_interaction.io.pseudo_io import read_psp8 as pseudo_reader
    return io, pseudo_reader


def hermiticity(M):
    nb, nk, _, _ = M.shape
    Op = M.reshape(nb * nk, nb * nk)
    return np.max(np.abs(Op - Op.conj().T))


def stage_ml(args):
    """MPI: build inputs on rank 0, broadcast, compute grid-distributed M^L, save on rank 0."""
    if not args.ml_out:
        raise SystemExit("--ml-out is required for --stage ml")
    if not (args.pot_p and args.pot_d):
        raise SystemExit("--pot-p and --pot-d are required for --stage ml")

    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import prep_realspace_inputs, compute_ML_R_mpi

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    io, _ = get_io(args.backend)
    bands = None if args.bands == "all" else [int(b) for b in args.bands.split(",")]

    prep = None
    if rank == 0:
        print(f"[rank0] stage=ml backend={args.backend} bands={args.bands} nranks={comm.Get_size()}", flush=True)
        prep = prep_realspace_inputs(args.uc, args.sc_p, args.pot_p, args.pot_d,
                                     subtract_mean=False, bands=bands, io=io)
    prep = comm.bcast(prep, root=0)

    M_L = compute_ML_R_mpi(prep, grid_block=args.block_size)

    if rank == 0:
        np.save(args.ml_out, M_L)
        print(f"[rank0] saved {args.ml_out}  shape={M_L.shape}", flush=True)
    comm.Barrier()


def stage_nl(args):
    """Serial (fresh process): load M^L, compute M^NL, assemble M = M^L + M^NL, save."""
    from electron_defect_interaction.defects.non_local import compute_M_NL

    if not (args.ml and args.out):
        raise SystemExit("--ml and --out are required for --stage nl")
    if not (args.sc_d and args.upf):
        raise SystemExit("--sc-d and --upf are required for --stage nl")

    io, pseudo_reader = get_io(args.backend)
    bands = None if args.bands == "all" else [int(b) for b in args.bands.split(",")]

    M_L = np.load(args.ml)
    print(f"loaded M^L {M_L.shape} from {args.ml}", flush=True)

    M_NL = compute_M_NL(args.uc, args.sc_p, args.sc_d, args.upf,
                        io=io, pseudo_reader=pseudo_reader, bands=bands)
    print(f"computed M^NL {M_NL.shape}", flush=True)

    M = M_L + M_NL
    np.save(args.out, M)
    print(f"saved {args.out}  shape={M.shape}  max|M-M^dag|={hermiticity(M):.2e}", flush=True)


def main():
    args = parse_args()
    if args.stage == "ml":
        stage_ml(args)
    else:
        stage_nl(args)


if __name__ == "__main__":
    main()
