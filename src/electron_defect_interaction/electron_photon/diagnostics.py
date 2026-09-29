"""
Diagnostics of a WannierTB (M1): numbers that describe a model without changing it
(hermiticity_report, symmetry_report).
"""

import numpy as np
from electron_defect_interaction.electron_photon.velocity_operator import dagger, fourier
from electron_defect_interaction.electron_photon.tb_model import centres_only

def extract_block(X, rows, columns):
    """Sub-block rows x columns of the last two axes of X (X[..., rows, columns] would pair them)."""
    return X[..., rows, :][..., columns]

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

    r_R_nc = tb.r_R - centres_only(tb).r_R
    S = r_R_nc + dagger(r_R_nc[tb.minus])

    A_k = fourier(tb.r_R, tb.R_cart, tb.ndegen, k)
    DA = A_k - dagger(A_k)

    report = {}
    for name, block in blocks.items():

        Dr_b = extract_block(Dr, block[0], block[1]) 
        S_b  = extract_block(S, block[0], block[1])
        DA_b = extract_block(DA, block[0], block[1])
        
        max_R = np.abs(Dr_b).max(axis=(0,2,3)) #(3,)
        frob_R = np.linalg.norm(Dr_b) / np.linalg.norm(S_b) # float
        max_k = np.abs(DA_b).max(axis=(0,2,3)) # (3,)

        report[name] = {'max_R': max_R, 'frob_R': frob_R, 'max_k': max_k}

    return report

def symmetry_report(tb, even, odd):
    """
    Selection rules of the mirror sigma_h of the sheet (F11, M1). The sigma bonds are even, the p_z
    odd; H, x and y are even, z is odd. So H and r^{x,y} vanish between even and odd functions, and
    r^z within each parity, measured without the centres (they carry the height z0 of the sheet).

    Inputs:
        tb   : WannierTB
        even : list of indices of the even Wannier functions (sigma)
        odd  : list of indices of the odd Wannier functions (p_z)
    Returns:
        dict {'H_mixed': float eV, 'r_mixed': (3,) Angstrom (x y forbidden, z allowed),
              'rz_same': float Angstrom}
    """

    H_mixed = np.max(np.abs(extract_block(tb.H_R, even, odd)))
    r_mixed = np.max(np.abs(extract_block(tb.r_R, even, odd)), axis=(0,2,3))

    # for rz work with r_R without centers
    r_R_nc = tb.r_R - centres_only(tb).r_R

    rz_even = np.max(np.abs(extract_block(r_R_nc[:,2], even, even)))
    rz_odd  = np.max(np.abs(extract_block(r_R_nc[:,2], odd, odd)))
    rz_same = np.max([rz_even, rz_odd])

    symmetries = {'H_mixed': H_mixed,  'r_mixed': r_mixed, 'rz_same': rz_same}

    return symmetries