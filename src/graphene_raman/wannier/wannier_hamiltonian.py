"""
wannier_hamiltonian.py
    Wannier Hamiltonian at arbitrary k-points (Hwr_to_Hwk) and the Dirac energy of graphene from band energies
    (dirac_point).
"""

import numpy as np

def Hwr_to_Hwk(Hwr, Rw, k, ndegen=None):
    """
    Wannier Hamiltonian at any k-points, H(k) = sum_R e^{2 pi i k.R} H(R) / ndegen(R), and its eigenpairs.
    Inputs:
        Hwr:    (nRw, nw, nw) complex, H(R) from Wannier90 (eV).
        Rw:     (nRw, 3) ints, its R vectors (reduced coordinates).
        k:      (nk, 3) floats, k-points (reduced coordinates).
        ndegen: (nRw,) ints or None, Wigner-Seitz degeneracies (None = 1; omitting them is wrong for Wannier90 data).
    Returns:
        Hwk: (nk, nw, nw) complex, Wannier gauge.
        Ewk: (nk, nw) floats, band energies, ascending.
        Uwk: (nk, nw, nw) complex, eigenvectors (columns): rotation to the smooth Bloch gauge.
    """

    # divide out the Wigner-Seitz degeneracies: H(R) -> H(R)/ndegen(R)
    if ndegen is not None:
        Hwr = Hwr / ndegen[:, None, None]

    # precompute phase
    phase = np.exp(2j*np.pi * (k @ Rw.T)) # (nkf, nrpts)
    nk, nR = phase.shape; nW = Hwr.shape[-1]

    # Sum over Rw
    Hwk = phase @ Hwr.reshape(nR, -1) # (nk, nW*nW)
    Hwk = Hwk.reshape(nk, nW, nW)
    
    # Diagonalize to get Wannier interpolated eigenvalues and eigenvectors
    Ewk, Uwk = np.linalg.eigh(Hwk) # (nk, nW), (nk, nW, nW)

    return Hwk, Ewk, Uwk

def dirac_point(E, bands=(3, 4)):
    """
    Dirac energy of graphene from band energies: midpoint of the smallest gap between the pi and pi* bands.
    Inputs:
        E:     (nk, nw) floats, band energies sorted per k (e.g. Ewk of Hwr_to_Hwk), on a grid or path through K.
        bands: (2,) ints, indices of pi and pi* (3, 4 in the 5-function wannierization).
    Returns:
        E_D:   float, (E[i, pi] + E[i, pi*]) / 2 at the k of the smallest gap.
        gap:   float, that smallest gap (same unit as E).
    """
    lo, hi = bands
    gap = E[:, hi] - E[:, lo]
    i = int(np.argmin(gap))
    return float(0.5 * (E[i, lo] + E[i, hi])), float(gap[i])
