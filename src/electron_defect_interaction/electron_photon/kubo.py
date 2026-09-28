"""
Absorptive optical conductivity sigma(omega)/sigma_0 by the Kubo formula, accumulated block by
block (`kubo_accumulate`, `kubo_normalize`), and the driver `sigma_on_grid`.
"""

import numpy as np

from .kgrid import make_grid_tb
from .velocity_operator import compute_velocity

def gaussian_eta(x, eta):
    """
    Normalized Gaussian of standard deviation eta (not the FWHM = 2.355 eta), in 1/eV: the broadened
    delta(hw - gap) of the Kubo sum.
    """
    return np.exp(-(x**2)/(2*eta**2)) / (eta*np.sqrt(2*np.pi))

def kubo_accumulate(eps, hv, hw, mu=0, eta=0.04):
    """
    Partial Kubo sum over one block of k points,

        S_ab(w) = sum_k sum_{eps_v < mu < eps_c} Re[(hbar v^a_cv)* hbar v^b_cv] g_eta(hw - eps_c + eps_v).

    Transitions are selected by energy (c empty, v occupied at each k), all pairs at once, then
    flattened to a list of n_t transitions so that the sum over them is one matrix product
    (nw, n_t) @ (n_t, 9). The weights do not depend on the arbitrary eigh phases. Blocks add up;
    normalize once at the end with `kubo_normalize`.

    Inputs:
        eps : (nk_block, nW) float, eV, from `velocity`
        hv  : (nk_block, 3, nW, nW) complex, eV*Angstrom, from `velocity`
        hw  : (nw,) float or scalar, eV, photon energies
        mu  : float, eV, chemical potential
        eta : float, eV, Gaussian broadening (standard deviation)
    Returns:
        S : (nw, 3, 3) float, eV*Angstrom^2 (block contribution, not normalized)
    """
    hw = np.atleast_1d(hw) 
    # allowed vertical transitions: pair[k, c, v] = c empty and v occupied
    occ = eps < mu # (nk_block, nW) bool, gives true if the state is occupied and false otherwise
    pair = (~occ)[:, :, None] & occ[:, None, :] # (nk_block, nW, nW) bool
    ik, ic, iv = np.nonzero(pair) # (nt,) each: k point, conduction band and valence band of each transition
    nt = ik.shape[0]; nw = len(hw)

    gaps = eps[ik, ic] - eps[ik, iv] # (nt,) optical gaps at each k
    hv_cv = hv[ik, :, ic, iv] # (nt, 3) velocity operator for allowed transitions

    W = np.real(np.conj(hv_cv[:, :, None]) * hv_cv[:, None, :]) # (nt, 3, 3) # Re[(ħv^α_cv)* ħv^β_cv] 
    
    G = gaussian_eta(hw[:, None] - gaps[None, :], eta=eta) # (nw, nt)

    W_ = W.reshape(nt, 9) # (nt, 9); 9 rather than -1, so that a block without transition (nt = 0) works
    S  = G @ W_ # (nw, 9)

    return S.reshape(nw, 3, 3)

def kubo_normalize(S, hw, nk, A_cell):
    """
    Absorptive conductivity in units of sigma_0 = e^2 / (4 hbar), spin included:

        sigma_ab / sigma_0 = 8 pi S_ab / (hw N_k A_cell).

    Kubo gives Re sigma = (pi e^2 g_s / (omega N_k A_cell)) sum |v_cv|^2 delta(hw - gap), with
    v = (hbar v)/hbar and g_s = 2; dividing by sigma_0 leaves 8 pi. A perfect Dirac cone gives
    exactly 1 (universal absorption pi alpha).

    Inputs:
        S      : (nw, 3, 3) float, eV*Angstrom^2, sum of `kubo_accumulate` over all blocks
        hw     : (nw,) float, eV
        nk     : int, number of k points of the WHOLE grid (not of a block)
        A_cell : float, Angstrom^2, in-plane cell area |a1 x a2|
    Returns:
        (nw, 3, 3) float, sigma / sigma_0
    """
    return 8 * np.pi * S / (hw[:, None, None] * nk * A_cell)

def sigma_on_grid(tb, N, hw, mu=0, eta=0.04, mode='berry', chunk=int(1e5)):
    """
    Optical conductivity sigma(omega)/sigma_0 on the shifted N x N grid: the driver of M0 and M4.

    Block by block: compute_velocity -> kubo_accumulate, then one kubo_normalize with N_k = N^2.
    The result does not depend on `chunk`. Resolving the resonant ring needs a grid step below
    ~eta / (2 hbar v_F), i.e. N >~ 800 for eta = 0.04 eV (N = 300 is already within 3e-4 in M0).
    Cost in M0: ~3 s per call at N = 1800.

    Inputs:
        tb    : WannierTB
        N     : int, grid of N x N points over the whole Brillouin zone (both valleys)
        hw    : (nw,) float or scalar, eV, photon energies
        mu    : float, eV, chemical potential (0 in M0, E_D in M1)
        eta   : float, eV, Gaussian broadening (standard deviation)
        mode  : 'berry' or 'no_berry', passed to `compute_velocity`
        chunk : int, number of k points per block (memory: see `fourier`)
    Returns:
        (nw, 3, 3) float, sigma_ab / sigma_0
    """

    hw = np.atleast_1d(hw)

    grid = make_grid_tb(tb, N)  

    # loop over blocks and accumulate
    S = 0
    for block in range(0, grid.nk, chunk):

        k = grid.k_cart[block:block+chunk] # only keep kpoints in this block (the last one may be shorter)

        _, eps, _, hv = compute_velocity(tb, k, mode) # ValueError for an unknown mode

        S += kubo_accumulate(eps, hv, hw, mu, eta) # (nw, 3, 3)

    # normalize once, with the size of the whole grid and the in-plane cell area
    A_cell = np.linalg.norm(np.cross(tb.lattice[:,0], tb.lattice[:,1]))

    return kubo_normalize(S, hw, grid.nk, A_cell)
