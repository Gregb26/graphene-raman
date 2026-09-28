"""
Wannier interpolation of H, dH/dk and the Berry connection A (single Fourier routine `fourier`),
and the velocity matrix elements in the band basis (`velocity`),

    hbar v(k) = V(k)^dagger [ dH(k)/dk + i [H(k), A(k)] ] V(k).
"""

import numpy as np

def dagger(X):
    """
    Hermitian conjugate over the last two axes; the leading axes (k, mu) are untouched.
    """

    return np.conj(np.swapaxes(X, -2, -1))

def hermitize(X):
    """
    Hermitian part (X + X^dagger)/2 over the last two axes.

    For A(k): the r(R) of Wannier90 are Hermitian only approximately (finite differences). Applied
    by the caller after `fourier`, so that the raw non-Hermiticity stays visible to the tests.
    """

    return (X + dagger(X))/2

def fourier(X_R, R_cart, ndegen, k, deriv=None):
    """
    Lattice-gauge Fourier transform to arbitrary k points, the single one of the module:

        X(k)        = sum_R e^{+ik.R} X(R) / ndegen(R)
        dX(k)/dk_mu = sum_R i R_mu e^{+ik.R} X(R) / ndegen(R)      (deriv=True)

    Used for H(k) (Nk, nW, nW) in eV, dH/dk (Nk, 3, nW, nW) in eV*Angstrom and A(k) (Nk, 3, nW, nW)
    in Angstrom. Same phase as Hwr_to_Hwk (k_cart . R_cart = 2 pi k_red . R_int). The phase matrix
    takes 16 Nk nR bytes (38 GB for 1800^2 k and 741 R): large grids go through in blocks of
    1e4-1e5 k points.

    Inputs:
        X_R    : (nR, ...) complex; only (nR, nW, nW) can be differentiated
        R_cart : (nR, 3) float, Angstrom
        ndegen : (nR,) int
        k      : (Nk, 3) float, 1/Angstrom (a single k point as k[None, :])
        deriv  : bool, return dX/dk instead of X(k)
    Returns:
        X_k : (Nk, ...) complex, or (Nk, 3, nW, nW) with deriv
    """

    # shapes
    nR = X_R.shape[0]; Nk = k.shape[0]

    # compute phase scaled by ndegen: phase[k, R] = e^{ik.R} / ndegen(R)
    phase = np.exp(1j * k @ R_cart.T) / ndegen[None, :] # (nK, nR)

    # Fourier transform
    if deriv:
        # a 4D r_R would broadcast silently into a wrong shape, hence the check
        assert X_R.ndim == 3 # differentiate hamiltonian only and not position operator
        # (nR, 1, nW, nW) * (nR, 3, 1, 1) -> (nR, 3, nW, nW): element [R, mu] = i R_mu X(R)
        X_R = X_R[:, None, :, :] * 1j*R_cart[:,:, None, None]

    # trailing shape, read AFTER the derivative (which adds the Cartesian axis)
    tail = X_R.shape[1:]

    # flatten the trailing axes, sum over R as one matrix product, restore the trailing axes (C order)
    X_R = X_R.reshape(nR, -1) # (nR, M)
    X_k = phase @ X_R # (nK, M)
    X_k = X_k.reshape(Nk, *tail)

    return X_k

def velocity(H_k, dH_k, A_k=None):
    """
    Band energies, eigenvectors and hbar v_mn(k) in the band basis,

        hbar v^mu = V^dagger [ dH/dk_mu + i [H, A_mu] ] V      (Wang et al., PRB 74, 195118 (2006)).

    The Berry term i[H, A] accounts for where the Wannier functions actually are. It vanishes on the
    diagonal (hbar v_nn = d eps_n/dk either way) and makes the interband elements independent of the
    home cell of each Wannier function. No consistency check inside (they would cost a
    diagonalization per block; they live in the tests). Interband elements between degenerate bands
    (at K) are ill-defined: grids and rings avoid K.

    Inputs:
        H_k  : (Nk, nW, nW) complex, eV
        dH_k : (Nk, 3, nW, nW) complex, eV*Angstrom (not modified)
        A_k  : (Nk, 3, nW, nW) complex, Angstrom, Hermitian; None drops the Berry term
               (eq. (2.5.7) of the thesis as written)
    Returns:
        eps : (Nk, nW) float, eV, ascending at every k
        V   : (Nk, nW, nW) complex, eigenvectors in columns (Wannier -> band rotation); the phase of
              each column is arbitrary at every k, so only |hbar v_mn|^2 is comparable between calls
        hv  : (Nk, 3, nW, nW) complex, eV*Angstrom, Hermitian, mu on axis 1
              (hbar v_F = 5.751 eV*Angstrom <-> v_F = 8.7e5 m/s)
    """
    nk, nW, _ = H_k.shape

    # diagonalize H(k) at every k at once: eps ascending, eigenvectors in the columns of V
    eps, V = np.linalg.eigh(H_k) # (nk, nW) floats, (nk, nW, nW) complex

    # velocity operator (times hbar) in the Wannier basis, W = dH + i [H, A];
    # copy(): W = dH_k would alias the caller's array and the += below would overwrite it
    W = dH_k.copy() # (nk, 3, nW, nW) complex
    if A_k is not None:
        assert A_k.shape == (nk, 3, nW, nW), 'wrong shape for A_k'
        # H_k[:, None]: same H(k) for the three components A_mu; i makes the commutator Hermitian
        W += 1j*(H_k[:, None, ...] @ A_k - A_k @ H_k[:, None, ...]) # add berry connection term

    assert W.shape == (nk, 3, nW, nW), 'wrong shape for W'

    # rotation to the band basis, with the same V(k) for the three components mu
    hv = dagger(V)[:, None, ...] @ W @ V[:, None, ...] # (nk, 3, nW, nW) complex
    assert hv.shape == (nk, 3, nW, nW), 'wrong shape for hv'

    return eps, V, hv
