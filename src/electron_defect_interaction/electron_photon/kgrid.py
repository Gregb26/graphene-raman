"""
Reciprocal lattice and k grids: `reciprocal`, `k_grid`, and the `GridTB` container of the callers.
"""

import numpy as np
from dataclasses import dataclass

@dataclass
class GridTB:
    """
    k-grid quantities for the callers (tests, driver), built by `make_grid_tb`. The computing
    functions take plain k arrays instead, since they also run on blocks, rings and paths.

    Fields:
        B      : (3, 3) float, 1/Angstrom, B[:, j] = b_j
        k_red  : (N^2, 3) float, reduced k points
        k_cart : (N^2, 3) float, 1/Angstrom
        K      : (3,) float, 1/Angstrom, Dirac point (2 b1 + b2)/3
        nk     : int, number of k points, N^2
    """
    B: np.ndarray
    k_red: np.ndarray
    k_cart: np.ndarray
    K: np.ndarray
    nk: int

def reciprocal(lattice):
    """
    Reciprocal vectors from a_i . b_j = 2 pi delta_ij, i.e. B = 2 pi (lattice^T)^{-1}.

    Inputs:
        lattice : (3, 3) float, Angstrom, lattice[:, i] = a_i
    Returns:
        B : (3, 3) float, 1/Angstrom, B[:, j] = b_j
    """

    B = 2*np.pi*np.linalg.inv(lattice.T) # B[:,i] = b_i; inv, not pinv: a singular lattice must fail

    return B

def k_grid(B, N, shift=0.5):
    """
    Shifted N x N grid over one reciprocal cell, k = ((i + s)/N) b1 + ((j + s)/N) b2.

    Any complete cell is equivalent for Brillouin-zone sums (H(k + G) = H(k) in the lattice gauge).
    The half-step shift keeps K and K' off the grid when N is a multiple of 3: the bands are
    degenerate there and the occupation would be ambiguous.

    Inputs:
        B     : (3, 3) float, 1/Angstrom, B[:, j] = b_j
        N     : int, number of points along b1 and along b2
        shift : float, offset in units of the grid step
    Returns:
        k_red  : (N^2, 3) float, reduced coordinates, j (along b2) varying fastest
        k_cart : (N^2, 3) float, 1/Angstrom, k_red @ B.T
        K      : (3,) float, 1/Angstrom, Dirac point (2 b1 + b2)/3
    """

    u = (np.arange(N) + shift)/N # reduced coordinates along one axis, in (0, 1)

    # all (u1, u2) pairs: two (N, N) arrays, flattened to N^2 aligned entries
    u1, u2 = np.meshgrid(u,u, indexing='ij')
    u3 = np.zeros(N**2)

    k_red = np.column_stack((u1.ravel(), u2.ravel(), u3))
    k_cart = k_red @ B.T # rows of reduced coords -> rows of cartesian vectors

    K = (2*B[:,0] + B[:,1])/3 # Dirac point

    return k_red, k_cart, K

def make_grid_tb(tb, N):
    """
    GridTB of `tb` (only its lattice is used): reciprocal vectors, shifted N x N grid (`k_grid`)
    and Dirac point.
    """
    B = reciprocal(tb.lattice)
    k_red, k_cart, K = k_grid(B,N)
    nk = k_red.shape[0]

    return GridTB(B=B, k_red=k_red, k_cart=k_cart, K=K, nk=nk)
