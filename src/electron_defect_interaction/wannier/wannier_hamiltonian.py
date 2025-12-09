"""
Python module containing functions to transform the Hamiltonian in Wannier (real space) basis to reciprocal space and interpolate eigenvalues and eigenvectors on a fine kpoint grid.
"""

import numpy as np

def Hwr_to_Hwk(Hwr, Rw, k):
    """
    Transform Hamiltonian in Wannier (real space) basis to reciprocal space in Wannier gauge on an arbitrary kpoint grid.
    Inputs:
        Hwr: (nrpts, nw, nw) array of complex
            Hamiltonian in Wannier basis computed by Wannier90. nrpts is the number of Rw vectors used internally by Wannier90 to compute the Hamiltonian and nw is the number of Wannier functions.
        Rw: (nrpts, 3) array of ints
            Real space lattice vectors in reduced coordinates used internally by Wannier90 to compute the Hamiltonian. 
        k: (nkf, 3) array of floats
            k vectors in reduced coordinate on which to transform Hwr. This can be any grid.
    Returns:
        Hwk: (nkf, nw, nw) array of complex
            Hamiltonian in reciprocal space in Wannier gauge on an arbitrary grid k.
    """

    # precompute phase
    phase = np.exp(2j*np.pi * (k @ Rw.T)) # (nw, nrpts)

    # Sum over Rw
    Hwk = np.einsum("kr, rwW -> kwW", phase, Hwr, optimize=True)

    # TODO Check hermiticity

    # Diagonalize to get Wannier interpolated eigenvalues and eigenvectors
    Ewk, Uwk = np.linalg.eigh(Hwk)
    
    return Hwk, Ewk, Uwk