#!/usr/bin/env python
"""Diagnostic: run compute_M_NL_mpi under MPI with faulthandler so any heap-corruption
abort prints a Python traceback pinpointing the offending line (einsum loop vs Allreduce).
Usage: srun python _diag_mnl_mpi.py UC SCP SCD UPF OUT"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("FLEXIBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import faulthandler
faulthandler.enable()

import sys
import numpy as np
from mpi4py import MPI

from electron_defect_interaction.io import qe_io
from electron_defect_interaction.io.pseudo_io import read_upf
from electron_defect_interaction.defects.non_local import compute_M_NL_mpi

uc, scp, scd, upf, out = sys.argv[1:6]
rank = MPI.COMM_WORLD.Get_rank()

if rank == 0:
    print(f"[rank0] compute_M_NL_mpi nranks={MPI.COMM_WORLD.Get_size()}", flush=True)
M = compute_M_NL_mpi(uc, scp, scd, upf, io=qe_io, pseudo_reader=read_upf, bands=None)
if rank == 0:
    np.save(out, M)
    print(f"[rank0] MPI M^NL saved {out} shape={M.shape}", flush=True)
MPI.COMM_WORLD.Barrier()
