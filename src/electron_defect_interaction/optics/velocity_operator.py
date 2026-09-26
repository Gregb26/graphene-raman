"""
Functions to compute Wannier-interpolated velocity operator from DFT outputs
"""

import numpy as np
from dataclasses import dataclass

@dataclass
class WannierTB:
    """Class for storing tight binding parameters"""
    lattice: np.ndarray
    R_int: np.ndarray
    ndegen: np.ndarray
    H_R: np.ndarray
    r_R: np.ndarray
    t: float

def make_graphene_tb(t=2.7, a_cc=1.42, c=15.0, shift_B=(0,0,0)):
    """ Construct graphene's tight-binding Hamiltonian in real space for a p_z orbital per atom """

    # build lattice vectors in cartesian coords
    a = np.sqrt(3)*a_cc # graphene lattice constant
    a1 = a*np.array([1,0,0]); a2 = a*np.array([1,np.sqrt(3),0])/2; a3=np.array([0,0,c])
    lattice = np.column_stack((a1, a2, a3)) # lattice[:,i] = a_i

    # cartesian vectors that connect A to its three B neighbours
    delta1 = a_cc*np.array([np.sqrt(3)/2, 0.5, 0]); delta2 = a_cc*np.array([-np.sqrt(3)/2, 0.5, 0]); delta3 = a_cc*np.array([0, -1.0, 0])
    delta = np.column_stack((delta1, delta2, delta3)) # delta[:,i] = delta_i

    # atom positions in the unit cell cartesian coords
    tau_A = np.zeros(3)
    tau_B = tau_A + delta[:,0] + lattice @ shift_B

    # Hopping : one hop per bond
    hops = []
    for d in range(3):
        R_red = np.linalg.solve(lattice, tau_A + delta[:,d] - tau_B)
        R = np.rint(R_red).astype(int)
        assert np.allclose(R_red, R, atol=1e-8), 'bond does not end on a B site'
        hops.append(R)

    # R list : A-B hops, their Hermitian partners at B-A at -R and R=0 for centers
    R_int = np.unique(np.array(hops + [-R for R in hops]+[np.zeros(3, int)]), axis=0)
    index = {tuple(R): iR for iR, R in enumerate(R_int)}
    nR = len(R_int)
    nW = 2 # two p_z orbitals in a unit cell

    H_R = np.zeros((nR, nW, nW), complex)
    r_R= np.zeros((nR, 3, nW, nW), complex)
    ndegen = np.ones(nR, dtype=int)

    A,B = 0,1
    for R in hops:
        H_R[index[tuple(R)], A, B] = -t # <0A|H|RB>
        H_R[index[tuple(-R)], B, A] = - t # <0B|H|-RA>

    i0 = index[(0,0,0)]
    r_R[i0, :, A, A] = tau_A
    r_R[i0,:, B, B] = tau_B

    return WannierTB(lattice=lattice, R_int=R_int, ndegen=ndegen, H_R=H_R, r_R=r_R, t=t)

def make_graphene_tb_analytic(t, k, A):
    """
    Inputs:
        -t : float, hopping parameter
        k  : (Nk, 3) ndarray of floats, kpoint grid in cartesian coordinates
        A  : (3, 3)  ndarray of floats, primitive lattice vectors stored in columns A[:,i] = a_i
    Returns:
        TB Hamiltonian : H_AB(k) = -t(1+e^{-ik.a1}+e^{-ik.a2})
    """

    return -t*(1+np.exp(-1j*k @ A[:,0].T) + np.exp(-1j*k @ A[:,1].T))

def reciprocal(lattice):
    """ compute primitive reciprocal lattice vectors from primitive lattice vectors """

    B = 2*np.pi*np.linalg.inv(lattice.T) # B[:,i] = b_i

    assert np.allclose(lattice @ B.T, 2*np.pi*np.eye(3)), 'not lattice vectors'

    return B

def k_grid(B, N, shift=0.5):
    """ build kpoint grid in reduced and cartesian coords """

    u = (np.arange(N) + shift)/N

    u1, u2 = np.meshgrid(u,u, indexing='ij')
    u3 = np.zeros(N**2)

    k_red = np.column_stack((u1.ravel(), u2.ravel(), u3))
    k_cart = k_red @ B.T

    K = (2*B[:,0] + B[:,1])/3 # Dirac point

    return k_red, k_cart, K

def hermitize(X):
    """ Compute hermitian conjugate of X """
    return np.conj(np.swapaxes(X, -2, -1))

def fourier(X_R, R_cart, ndegen, k, deriv=None):
    """
    Inputs:
        X_R : H_R (nR, nW, nW) deriv=None or H_R (nR, nW, nW) deriv=True or r_R (nR, 3, nW, nW) deriv=None
        R_cart : (nR, 3)
        ndegen : (nR,)
        k : (Nk, 3)
    Outputs : sum_R e^{ik.R} X(R)/ndegen(R)
    """

    # shapes
    nR = X_R.shape[0]; Nk = k.shape[0]

    # compute phase scaled by ndegen
    phase = np.exp(1j * k @ R_cart.T) / ndegen[None, :] # (nK, nR)

    # Fourier transform 
    if deriv:
        X_R = X_R[:, None, :, :] * 1j*R_cart[:,:, None, None]
        
    tail = X_R.shape[1:]

    X_R = X_R.reshape(nR, -1) # (nR, M)
    X_k = phase @ X_R # (nK, M)
    X_k = X_k.reshape(Nk, *tail)

    # hermitian test
    X_k_dag = hermitize(X_k)
    assert np.allclose(np.abs(X_k - X_k_dag), 0, atol=10e-11), 'object is not hermitian'

    return X_k
