"""
wannier_interpolation.py
    Python module containing functions to interpolate objects computed on a coarse kpoint grid onto a fine kpoint grid
    via Maximally Localized Wannier Functions. 
"""

import numpy as np
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk

def Mbk_to_Mwk(Mbk, U, U_dis=None):
    """
    Rotates an object Mbk in Bloch gauge into Wannier gauge using the Wannier gauge matrix U (or V=U_dis @ U in the case of entangled bands)
    Inputs:
        Mbk: (nb, nk, nb, nk) array of complex
            Matrix in Bloch gauge to rotate. This can be e.g. a Hamiltonian or a scattering matrix. nk is the number of kpoints, nb is the number
            of Bloch bands and nw is the number of Wannier functions.
        U: (nk, nb, nw) or (nk, nw, nw) array of complex
            Wannier gauge matrix. First shape if no entangled bands (U_dis is None). Second shape if entangled bands (U_dis is not None).
        U_dis: (nk, nb, nw) array of complex
            Disentanglement matrix if the number of Wannier functions is less than the number of Bloch bands.
    Returns:
        Mwk: (nw, nk, nw, nk)
            Object in the Wannier gauge: Mwk = V^dag Mbk V
    """

    nk, nw, _ = U.shape
    # Entangled case, rotation matrix is V = U_dis @ U
    if U_dis is not None:
        nb = U_dis.shape[1]
        V = np.zeros((nk, nb, nw), dtype=complex)
        for ik in range(nk):

            V[ik] = U_dis[ik] @ U[ik] # (nk, nb, nw)

            # testing
            with np.errstate(all='ignore'):
                P = V[ik] @ V[ik].conj().T
                assert np.allclose(P, P.conj().T, atol=1e-10)     
                assert np.allclose(P @ P, P, atol=1e-8)           
                assert np.allclose(np.trace(P), V.shape[-1], 1e-10)

    else:
        V = U # (nk, nb, nw) with nw = nb
    
    # Rotate Mb from Bloch gauge to Wannier gauge Mw
    Mwk = np.zeros((nw, nk, nw, nk), dtype=complex)
    for ik in range(nk):
        Vk_h = V[ik].conj().T # Hermitian conjugate
        for ikp in range(nk):

            Mwk[:, ik, :, ikp] = Vk_h @ Mbk[:, ik, :, ikp] @ V[ikp] # (nw, nk, nw, nk)
    
    return Mwk

def Mwk_to_Mwr(Mwk, k_red, MP_grid):
    """
    Transforms the object Mwk from reciprocal space to real space using a double Fourier transform.
    Inputs:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space to transform to real space, nk is the number of kpoints and nw the number of Wannier functions. This is e.g.
            a scattering matrix in reciprocal space in the Wannier gauge.
        MP_grid: tuple of ints
            Monkhorst-Pack grid used to define the kpoint grid. Assuming no symmetry reduction has been applied! This is important! 
    Returns
        Mwr: (nw, nr, nw, nr) array of complex
            Fourier transform of Mk, nr is the number of R vectors, inferred from the Monkhorst-Pack grid.
        R:  (nr, 3) array of ints
            R vectors in real space used to compute the double FT. Is dual to the kpoint grid.
    """

    nk = Mwk.shape[1]

    # Build the grid in real space based on the kpoint grid
    N1, N2, N3 = MP_grid
    r1 = np.arange(N1); r2 = np.arange(N2); r3 = np.arange(N3)
    rr1, rr2, rr3 = np.meshgrid(r1, r2, r3)
    R = np.stack((rr1, rr2, rr3), axis=-1) # (N1, N2, N3, 3)
    R = R.reshape(N1*N2*N3, 3) # (nr, 3)

    # Compute phase
    phase_kp = np.exp(2j*np.pi * (k_red @ R.T)) # (nk, nr)
    phase_k = phase_kp.conj() # (nk, nr)

    # Double sum over k and k'
    Mwr = np.einsum('kr, nkNK, KR -> nrNR', phase_kp, Mwk, phase_k, optimize=True)  / nk **2 # (nw, nr, nk, nr)

    return Mwr, R

def Mwr_to_Mwk(Mwr, R, k):
    """
    Transforms the object in Mwr in real space in Wannier gauge to Mwk in reciprocal space in Wannier gauge. Here, the kpoint grid k can be any arbitrarily dense grid hence this function can be used to interpolate Mwk on a finer grid.
    Inputs:
        Mwr: (nw, nr, nk, nr) array of complex
            Object in real space in Wannier gauge to transforms to reciprocal space in Wannier gauge. nr is the number of R vectors in real sapce and nw is the number of Wannier functions.
        R: (nr, 3) array of ints
            Real space lattice vectors in reduced coordinates for which Mwr is defined.
        k: (nw, 3) array of floats
            k vectors in reduced coordinates for which to transform Mwr.
    Returns:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space in Wannier gauge.
    """
    
    # precompute phases
    phase_kp = np.exp(-2j*np.pi * (k @ R.T)) # (nk, nr)
    phase_k = phase_kp.conj() # (nk, nr)

    # sum over R 
    Mwk = np.einsum("kr, wrWR, KR -> wkWK", phase_kp, Mwr, phase_k, optimize=True)

    return Mwk

def Mwk_to_Mbk(Mwk, Hwr, Rw, k):
    """
    Transforms the object Mwk in reciprocal space in Wannier gauge to Bloch gauge using the Wannier Hamiltonian.
    Inputs:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space in Wannier gauge to transform to Bloch gauge. nw is the number of Wannier functions and nk is the number of kpoints.
        Hwr: (nrpts, nw, nw) array of complex
            Hamtilontian in Wannier (real space) basis. This Hamiltonian is computed on the kpoint grid k via a Fourier transformed and diagonalized to obtain the unitary gauge rotation matrix U. nrw is the number of R vectors used internally by Wannier90 to compute the Hamiltonian in real space.
        Rw: (nrpts, 3) array of ints:
            Real space lattice vectors used internally by Wannier90.
        k: (nk, 3) array of floats
            kpoint grid on which to transform Hwr. Can arbitrarily dense and must match the one on which Mwk is defined.
    Returns:
        Mbk: (nw, nk, nw, nk) array of complex
            Object transformed to Bloch basis. 
    """

    # compute rotation matrix from Wannier Hamiltonian
    _, _, Uwk = Hwr_to_Hwk(Hwr, Rw, k) # (nk, nw, nw)

    # Rotate Mwk to Bloch gauge
    Mbk = np.einsum('kwb, wkWK, KWB -> bkBK', Uwk.conj(), Mwk, Uwk)

    return Mbk

