"""
Diagnostics of a WannierTB (M1): numbers that describe a model without changing it.
"""

import numpy as np
from electron_defect_interaction.electron_photon.velocity_operator import dagger, fourier
from electron_defect_interaction.electron_photon.tb_model import centres_only

def hermiticity_report(tb, blocks, k):
    """
    Hermiticity defect of the raw r, block by block, in R and in k (F10, M1).

    D(R) = r(R) - r(-R)^dagger and A(k) - A(k)^dagger (A not hermitized); frob_R = |D| / |S| in
    Frobenius norm, S = r(R) + r(-R)^dagger WITHOUT the centres (they depend on the origin).
    frob_R is nan when S vanishes in the block (e.g. the M0 model, which has only centres).

    Inputs:
        tb     : WannierTB
        blocks : dict name -> (rows, cols), lists of Wannier indices, e.g. 'pz': ([3,4], [3,4])
        k      : (Nk, 3) float, 1/Angstrom, Cartesian k points (e.g. make_grid_tb(tb, N).k_cart)
    Returns:
        dict name -> {'max_R': (3,) Angstrom, 'frob_R': float, 'max_k': (3,) Angstrom}, x y z
    """

    # compute defect in R
    Dr = tb.r_R - dagger(tb.r_R[tb.minus])

    r_no_centers = tb.r_R - centres_only(tb).r_R
    S = r_no_centers + dagger(r_no_centers[tb.minus])

    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, k)
    DA = A_k - dagger(A_k)

    report = {}
    for name, block in blocks.items():

        Dr_b = Dr[..., block[0], :][..., block[1]]
        S_b  = S[..., block[0], :][..., block[1]]
        DA_b = DA[..., block[0], :][..., block[1]]

        max_R = np.abs(Dr_b).max(axis=(0,2,3)) #(3,)
        frob_R = np.linalg.norm(Dr_b) / np.linalg.norm(S_b) # float
        max_k = np.abs(DA_b).max(axis=(0,2,3)) # (3,)

        report[name] = {'max_R': max_R, 'frob_R': frob_R, 'max_k': max_k}

    return report

