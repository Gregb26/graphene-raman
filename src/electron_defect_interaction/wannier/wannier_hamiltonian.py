"""
Python module containing functions to transform the Hamiltonian in Wannier (real space) basis to reciprocal space and interpolate eigenvalues and eigenvectors on a fine kpoint grid.
"""

import numpy as np

def Hwr_to_Hwk(Hwr, Rw, k, ndegen=None):
    """
    Transform Hamiltonian in Wannier (real space) basis to reciprocal space in Wannier gauge on an arbitrary kpoint grid.
    Inputs:
        Hwr: (nrpts, nw, nw) array of complex
            Hamiltonian in Wannier basis computed by Wannier90. nrpts is the number of Rw vectors used internally by Wannier90 to compute the Hamiltonian and nw is the number of Wannier functions.
        Rw: (nrpts, 3) array of ints
            Real space lattice vectors in reduced coordinates used internally by Wannier90 to compute the Hamiltonian.
        k: (nkf, 3) array of floats
            k vectors in reduced coordinate on which to transform Hwr. This can be any grid.
        ndegen: (nrpts,) array of ints, optional
            Wigner-Seitz degeneracy of each Rw (from the Wannier90 _hr.dat/_tb.dat file). The correct
            interpolation is H(k) = sum_R e^{2pi i k.R} H(R)/ndegen(R). If None, ndegen is taken to be 1.
    Returns:
        Hwk: (nkf, nw, nw) array of complex
            Hamiltonian in reciprocal space in Wannier gauge on an arbitrary grid k.
        Ewk: (nkf, nw) array of floats
            Wannier-interpolated band energies (eigenvalues of Hwk).
        Uwk: (nkf, nw, nw) array of complex
            Eigenvectors of Hwk; the per-k unitary rotation from the Wannier gauge to the smooth Bloch gauge.
    """

    # divide out the Wigner-Seitz degeneracies: H(R) -> H(R)/ndegen(R)
    if ndegen is not None:
        Hwr = Hwr / ndegen[:, None, None]

    # precompute phase
    phase = np.exp(2j*np.pi * (k @ Rw.T)) # (nkf, nrpts)

    # Sum over Rw
    Hwk = np.einsum("kr, rwW -> kwW", phase, Hwr, optimize=True)

    # Diagonalize to get Wannier interpolated eigenvalues and eigenvectors
    Ewk, Uwk = np.linalg.eigh(Hwk)

    return Hwk, Ewk, Uwk