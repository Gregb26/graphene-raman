#!/usr/bin/env python
"""
compute_M_dense_stages.py
    Dense electron-defect matrix M on the (p*N)x(p*N) primitive k-grid by EXACT zero-padding of the
    supercell defect potential (local_G.compute_ML_G_dense_mpi), in three isolated stages -- NOT the
    monolithic compute_M_dense, which runs M^L (MPI) and M^NL in one process (the pattern that
    corrupts the heap for nk >= 81):

        --stage ml      (MPI, srun):   M^L_dense  -> --out   [bloch_norm=unit_cell]
        --stage nl      (serial, fresh process): M^NL on the dense primitive .save -> --out
        --stage combine (serial): M = M^L + M^NL -> --out (unit_cell) and --out-norm = M/N_kd
                        (supercell; N_kd = (p*N)^2 = cells of the padded supercell = one defect)

    Sizes: 5x5->p=5 (25x25), 7x7->p=4 (28x28), 8x8->p=4 (32x32), 9x9->p=3 (27x27).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("FLEXIBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import argparse
import numpy as np
from electron_defect_interaction.io import qe_io, matrix_io

PFAC = {"5x5": 5, "7x7": 4, "8x8": 4, "9x9": 3}
SCRATCH = "/home/gregb26/links/scratch/qe_tmp"


def paths(size):
    n = int(size.split("x")[0]); p = PFAC[size]; D = n * p
    return dict(p=p, D=D,
                uc=f"data/graphene/unit_cell/qe/defect_{size}.save",
                uc_dense=f"{SCRATCH}/defect_uc_dense_{D}/defect_uc_dense_{D}.save",
                scp=f"data/graphene/supercell/qe/defect_{size}_p.save",
                scd=f"data/graphene/supercell/qe/defect_{size}_d.save",
                pot_p=f"data/graphene/supercell/qe/defect_{size}_p.save/Vks_{size}_p",
                pot_d=f"data/graphene/supercell/qe/defect_{size}_d.save/Vks_{size}_d",
                upf=f"data/graphene/unit_cell/qe/defect_{size}.save/C.upf")


def stage_ml(a):
    from mpi4py import MPI
    from electron_defect_interaction.defects.local_G import prep_reciprocal_inputs, compute_ML_G_dense_mpi
    P = paths(a.size); comm = MPI.COMM_WORLD; rank = comm.Get_rank()
    bands = None if a.bands == "all" else [int(b) for b in a.bands.split(",")]
    prep = None
    if rank == 0:
        print(f"[rank0] dense ml size={a.size} p={P['p']} D={P['D']} nranks={comm.Get_size()}", flush=True)
        prep = prep_reciprocal_inputs(P["uc"], P["scp"], P["pot_p"], P["pot_d"],
                                      subtract_mean=False, bands=bands, io=qe_io)
    prep = comm.bcast(prep, root=0)
    C_d, nG_d = qe_io.get_C_nk(P["uc_dense"]); G_d = qe_io.get_G_red(P["uc_dense"]); k_d = qe_io.get_k_red(P["uc_dense"])
    if bands is not None:
        C_d = C_d[list(bands), ...]
    assert len(k_d) == P["D"] ** 2, f"dense save has {len(k_d)} k, expected {P['D']**2}"
    M_L = compute_ML_G_dense_mpi(prep, P["p"], k_d, C_d, G_d, nG_d, block_size=a.block_size, show_tqdm=(rank == 0))
    if rank == 0:
        matrix_io.save_M(a.out, M_L, matrix_io.UNIT_CELL, part="M_L_dense", p=P["p"], D=P["D"])
        print(f"[rank0] saved {a.out} shape={M_L.shape}", flush=True)
    comm.Barrier()


def stage_nl(a):
    from electron_defect_interaction.io.pseudo_io import read_upf
    from electron_defect_interaction.defects.non_local import compute_M_NL
    P = paths(a.size); bands = None if a.bands == "all" else [int(b) for b in a.bands.split(",")]
    M_NL = compute_M_NL(P["uc_dense"], P["scp"], P["scd"], P["upf"], io=qe_io, pseudo_reader=read_upf, bands=bands)
    matrix_io.save_M(a.out, M_NL, matrix_io.UNIT_CELL, part="M_NL_dense", p=P["p"], D=P["D"])
    print(f"saved {a.out} shape={M_NL.shape}", flush=True)


def stage_combine(a):
    P = paths(a.size)
    M_L = matrix_io.load_M_checked(a.ml, require_bloch_norm=matrix_io.UNIT_CELL)
    M_NL = matrix_io.load_M_checked(a.nl, require_bloch_norm=matrix_io.UNIT_CELL)
    assert M_L.shape == M_NL.shape, (M_L.shape, M_NL.shape)
    M = M_L + M_NL; N_kd = M.shape[1]
    matrix_io.save_M(a.out, M, matrix_io.UNIT_CELL, p=P["p"], D=P["D"], N_kd=int(N_kd))
    if a.out_norm:
        matrix_io.save_M(a.out_norm, M / N_kd, matrix_io.SUPERCELL, N_cells=int(N_kd), p=P["p"], D=P["D"])
    nb, nk = M.shape[:2]; Op = M.reshape(nb * nk, nb * nk)
    print(f"saved {a.out} (+norm) shape={M.shape} N_kd={N_kd} max|M-M^dag|={np.max(np.abs(Op-Op.conj().T)):.2e}", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--stage", choices=["ml", "nl", "combine"], required=True)
    p.add_argument("--size", required=True, choices=list(PFAC))
    p.add_argument("--out", required=True); p.add_argument("--out-norm", default=None)
    p.add_argument("--ml"); p.add_argument("--nl")
    p.add_argument("--bands", default="all"); p.add_argument("--block-size", type=int, default=128)
    a = p.parse_args()
    {"ml": stage_ml, "nl": stage_nl, "combine": stage_combine}[a.stage](a)


if __name__ == "__main__":
    main()
