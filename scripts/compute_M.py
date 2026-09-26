#!/usr/bin/env python
"""
compute_M.py
    Single driver for the electron-defect scattering matrix M = M^L + M^NL from Quantum
    ESPRESSO outputs. Works for ANY supercell size. Three orthogonal stages let
    you compute the local part alone, the non-local part alone, or assemble the full matrix:

        stage 'ml'       (MPI, run under srun): local part M^L (real-space, grid-distributed
                         over ranks) -> --out
        stage 'nl'       (serial, plain python, FRESH process): non-local part M^NL -> --out
        stage 'combine'  (serial, trivial): M = M^L + M^NL from --ml + --nl -> --out

    Why ml and nl must be SEPARATE processes: compute_ML_R_mpi leaves a latent heap
    corruption for nk >= 81 (k-grid >= 9x9) that aborts the *next* heavy allocation in the
    same process (observed as "munmap_chunk(): invalid pointer" / "double free" inside
    compute_M_NL). A fresh process for stage 'nl' starts from a clean heap, so the pipeline
    works identically for 5x5 and 20x20. 'combine' is a cheap numpy add and can run anywhere.
    See submit_M.sh for the canonical launch chaining all three.

    Inputs (all from a finished QE run; use nosym+noinv for the unit-cell full k-grid):
        --uc      unit-cell prefix.save dir                 (ml, nl)
        --sc-p    pristine supercell prefix.save dir        (ml, nl)
        --sc-d    defective supercell prefix.save dir       (nl)
        --pot-p   pristine supercell local potential, pp.x plot_num=1 filplot   (ml)
        --pot-d   defective supercell local potential                           (ml)
        --upf     UPF pseudopotential                                           (nl)
        --ml      precomputed M^L .npy   (combine)
        --nl      precomputed M^NL .npy  (combine)

    Examples:
        # local part only (MPI):
        srun -n 32 python scripts/compute_M.py --stage ml \
            --uc uc.save --sc-p sc_p.save --pot-p sc_p.save/Vks_p --pot-d sc_d.save/Vks_d \
            --out results/M/M_L_9x9.npy
        # non-local part only (serial, fresh process):
        python scripts/compute_M.py --stage nl \
            --uc uc.save --sc-p sc_p.save --sc-d sc_d.save --upf uc.save/C.upf \
            --out results/M/M_NL_9x9.npy
        # full matrix from the two parts:
        python scripts/compute_M.py --stage combine \
            --ml results/M/M_L_9x9.npy --nl results/M/M_NL_9x9.npy \
            --out results/M/M_ed_9x9.npy
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
import time
import numpy as np

from electron_defect_interaction.io import matrix_io


def parse_args():
    p = argparse.ArgumentParser(description="Compute M = M^L + M^NL from Quantum ESPRESSO outputs.")
    p.add_argument("--stage", choices=["ml", "nl", "combine"], required=True,
                   help="'ml': MPI local part M^L (run under srun); "
                        "'nl': serial non-local part M^NL (fresh process); "
                        "'combine': M = M^L + M^NL from precomputed parts")
    p.add_argument("--uc", help="unit-cell wavefunctions (prefix.save dir)")
    p.add_argument("--sc-p", help="pristine supercell")
    p.add_argument("--sc-d", help="defective supercell (stage nl)")
    p.add_argument("--pot-p", help="pristine supercell local potential (stage ml)")
    p.add_argument("--pot-d", help="defective supercell local potential (stage ml)")
    p.add_argument("--upf", help="UPF pseudopotential (stage nl)")
    p.add_argument("--ml", help="precomputed M^L .npy (stage combine)")
    p.add_argument("--nl", help="precomputed M^NL .npy (stage combine)")
    p.add_argument("--out", required=True, help="output .npy for this stage's result")
    p.add_argument("--bands", default="all", help="comma-separated band indices, or 'all' (default)")
    p.add_argument("--block-size", type=int, default=50_000, help="real-space grid block per rank (stage ml)")
    return p.parse_args()


def require(args, names, stage):
    missing = [f"--{n.replace('_', '-')}" for n in names if not getattr(args, n)]
    if missing:
        raise SystemExit(f"--stage {stage} requires {', '.join(missing)}")


def get_io():
    from electron_defect_interaction.io import qe_io as io
    from electron_defect_interaction.io.pseudo_io import read_upf as pseudo_reader
    return io, pseudo_reader


def hermiticity(M):
    nb, nk, _, _ = M.shape
    Op = M.reshape(nb * nk, nb * nk)
    return np.max(np.abs(Op - Op.conj().T))


def stage_ml(args):
    """MPI: build inputs on rank 0, broadcast, compute grid-distributed M^L, save on rank 0."""
    require(args, ["uc", "sc_p", "pot_p", "pot_d"], "ml")

    from mpi4py import MPI
    from electron_defect_interaction.defects.local_R import prep_realspace_inputs, compute_ML_R_mpi

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    io, _ = get_io()
    bands = None if args.bands == "all" else [int(b) for b in args.bands.split(",")]

    prep = None
    if rank == 0:
        print(f"[rank0] stage=ml bands={args.bands} nranks={comm.Get_size()}", flush=True)
        prep = prep_realspace_inputs(args.uc, args.sc_p, args.pot_p, args.pot_d,
                                     subtract_mean=False, bands=bands, io=io)
    prep = comm.bcast(prep, root=0)

    M_L = compute_ML_R_mpi(prep, grid_block=args.block_size)

    if rank == 0:
        # raw part (unit-cell-normalized Bloch states); the 1/N_cells factor is applied at combine.
        matrix_io.save_M(args.out, M_L, matrix_io.UNIT_CELL, part="M_L", M_normalization="v2 : L et NL en norme unit_cell, 2026-09-25")
        print(f"[rank0] saved M^L {args.out}  shape={M_L.shape} (bloch_norm=unit_cell)", flush=True)
    comm.Barrier()


def stage_nl(args):
    """Serial (fresh process): compute the non-local part M^NL and save it."""
    require(args, ["uc", "sc_p", "sc_d", "upf"], "nl")

    from electron_defect_interaction.defects.non_local import compute_M_NL

    io, pseudo_reader = get_io()
    bands = None if args.bands == "all" else [int(b) for b in args.bands.split(",")]

    M_NL = compute_M_NL(args.uc, args.sc_p, args.sc_d, args.upf,
                        io=io, pseudo_reader=pseudo_reader, bands=bands)
    # raw part (unit-cell-normalized Bloch states); the 1/N_cells factor is applied at combine.
    matrix_io.save_M(args.out, M_NL, matrix_io.UNIT_CELL, part="M_NL", M_normalization="v2 : L et NL en norme unit_cell, 2026-09-25")
    print(f"saved M^NL {args.out}  shape={M_NL.shape} (bloch_norm=unit_cell)", flush=True)


def stage_combine(args):
    """Assemble the full matrix M = M^L + M^NL from two precomputed parts."""
    require(args, ["ml", "nl"], "combine")

    # both parts in unit-cell Bloch norm (R6 kernels) -> M = M^L + M^NL, unit_cell, same convention as the
    # dense chain (compute_M_dense_stages.py combine) and as every production consumer; the 1/N_cells factor
    # is applied downstream (bloch_folded_hamiltonian, compute_T), never here.
    M_L = matrix_io.load_M_checked(args.ml, require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    M_NL = matrix_io.load_M_checked(args.nl, require_bloch_norm=matrix_io.UNIT_CELL, units=matrix_io.HARTREE)
    if M_L.shape != M_NL.shape:
        raise SystemExit(f"shape mismatch: M^L {M_L.shape} vs M^NL {M_NL.shape}")
    N_cells = M_L.shape[1]
    M = M_L + M_NL
    matrix_io.save_M(args.out, M, matrix_io.UNIT_CELL, N_cells=int(N_cells), M_normalization="v2 : L et NL en norme unit_cell, 2026-09-25")
    print(f"saved M {args.out}  shape={M.shape}  N_cells={N_cells}  bloch_norm=unit_cell  M_normalization=v2  "
          f"max|M-M^dag|={hermiticity(M):.2e}", flush=True)


def main():
    args = parse_args()
    t0 = time.perf_counter()
    {"ml": stage_ml, "nl": stage_nl, "combine": stage_combine}[args.stage](args)
    dt = time.perf_counter() - t0
    # One timing line per stage; under MPI only rank 0 reports.
    rank = 0
    if args.stage == "ml":
        from mpi4py import MPI
        rank = MPI.COMM_WORLD.Get_rank()
    if rank == 0:
        print(f"[timing] stage={args.stage} elapsed={dt:.1f}s", flush=True)


if __name__ == "__main__":
    main()
