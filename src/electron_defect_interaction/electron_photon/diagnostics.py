"""
Diagnostics of a WannierTB: numbers that describe a model without changing it
(hermiticity_report, symmetry_report in M1; frozen_window_limit in M4).
"""

import numpy as np
from electron_defect_interaction.electron_photon.velocity_operator import dagger, fourier
from electron_defect_interaction.electron_photon.tb_model import centres_only, extract_block
from electron_defect_interaction.electron_photon.kgrid import make_grid_tb


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

def frozen_window_limit(tb, N, froz_max, mu, chunk=int(1e5)):
    """
    hbar omega_froz = min of eps_c - eps_v over the k whose eps_c lies above the frozen window (c, v on
    either side of mu). Above it, some transitions end on states that the interpolation no longer
    reproduces exactly at the coarse k points: the vertical line of the sigma(omega) figure. A minimum
    over grid points, so it converges from above with N.

    Inputs:
        tb       : WannierTB
        N        : int, N x N grid (make_grid_tb)
        froz_max : float, eV, top of the frozen window (dis_froz_max of the .win)
        mu       : float, eV, chemical potential (E_D for the real data)
        chunk    : int, k points per block (memory: see `fourier`)
    Returns:
        float, eV; inf if eps_c never exceeds froz_max
    """

    grid = make_grid_tb(tb, N)

    mins = [] # store mins per block
    # loop over kpoints in blocks
    for block in range(0, grid.nk, chunk):
        # only keep k's in this block
        k = grid.k_cart[block:block+chunk]

        # find eigenvalues for this block
        H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, k) # (chunk, nW, nW)
        eps = np.linalg.eigvalsh(H_k) # (chunk, nW)

        # find valence and conduction band indices
        n_occ = np.sum(eps < mu, axis=1) # (chunk,)
        n = n_occ[0]
        assert np.allclose(n_occ, n)
        v = n - 1 # valence band
        c = n     # conduction band

        # keep the k whose conduction state lies above the frozen window
        mask = eps[:, c] > froz_max
        # find min of eps_c - eps_v on these k's
        mins.append(np.min(eps[mask, c] - eps[mask, v], initial=np.inf)) # inf for a block without such k

    # global minimum over the blocks
    return np.min(mins)