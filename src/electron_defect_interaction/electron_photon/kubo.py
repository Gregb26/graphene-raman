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

def kubo_accumulate(eps, hv, hw, mu=0, eta=0.04, kT=0):
    """
    Partial Kubo sum over one block of k points,

        S_ab(w) = sum_k sum_{eps_v < mu < eps_c} Re[(hbar v^a_cv)* hbar v^b_cv] g_eta(hw - eps_c + eps_v),

    over all pairs (c empty, v occupied), selected by energy. Blocks add up; normalize once at the end
    with `kubo_normalize`.

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

    # compute occupations
    if kT == 0:
        f = np.zeros_like(eps) # (nk_block, nW)
        f[np.where(eps < mu)] = 1 # step function for T=0

    else:
        f = 0.5 * (1 - np.tanh((eps - mu) /(2*kT) )) # Fermi-Dirac for T>0, (nk_block, nW)

    # compute difference in occupation
    w = f[:, None, :] - f[..., None] # (nk_block, nW, nW)

    # get indices of allowed transitions
    ik, ic, iv = np.nonzero(w > 1e-14) # 3* (nt, )
    nt = ik.shape[0]; nw = len(hw)

    # get relevant quantities for each allowed transition 
    gaps = eps[ik, ic] - eps[ik, iv] # (nt,) optical gaps at each k
    hv_cv = hv[ik, :, ic, iv] # (nt, 3) velocity operator for allowed transitions
    w_cv  = w[ik, ic, iv] # (nt, ) transition weights

    W = np.real(np.conj(hv_cv[:, :, None]) * hv_cv[:, None, :]) # (nt, 3, 3) # Re[(ħv^α_cv)* ħv^β_cv] 
    
    G = gaussian_eta(hw[:, None] - gaps[None, :], eta=eta) # (nw, nt)

    W_ = W.reshape(nt, 9) # (nt, 9)
    W_weighted = W_ * w_cv[:, None]
    
    S  = G @ W_weighted # (nw, 9)

    return S.reshape(nw, 3, 3)

def kubo_normalize(S, hw, nk, A_cell):
    """
    Absorptive conductivity in units of sigma_0 = e^2 / (4 hbar), spin included:
    sigma_ab / sigma_0 = 8 pi S_ab / (hw N_k A_cell), 1 for a perfect Dirac cone (derivation: EM.md, M0).

    Inputs:
        S      : (nw, 3, 3) float, eV*Angstrom^2, sum of `kubo_accumulate` over all blocks
        hw     : (nw,) float, eV
        nk     : int, number of k points of the WHOLE grid (not of a block)
        A_cell : float, Angstrom^2, in-plane cell area |a1 x a2|
    Returns:
        (nw, 3, 3) float, sigma / sigma_0
    """
    return 8 * np.pi * S / (hw[:, None, None] * nk * A_cell)

def sigma_on_grid(tb, N, hw, mu=0, eta=0.04, mode='berry', kT=0, chunk=int(1e5)):
    """
    Optical conductivity sigma(omega)/sigma_0 on the shifted N x N grid (driver of M0 and M4): block by
    block compute_velocity -> kubo_accumulate, then kubo_normalize with N_k = N^2; independent of
    `chunk`. Converges in N much faster than the grid step suggests (EM.md, M4 pilot).

    Inputs:
        tb    : WannierTB
        N     : int, grid of N x N points over the whole Brillouin zone (both valleys)
        hw    : (nw,) float or scalar, eV, photon energies
        mu    : float, eV, chemical potential (E_D for the real data)
        eta   : float, eV, Gaussian broadening (standard deviation)
        mode  : 'berry' or 'no_berry', passed to `compute_velocity`
        chunk : int, number of k points per block (memory: see `fourier`)
    Returns:
        (nw, 3, 3) float, sigma_ab / sigma_0
    """
    assert kT >= 0, 'k_BT must be positive or zero !'
    hw = np.atleast_1d(hw)

    grid = make_grid_tb(tb, N)  

    # loop over blocks and accumulate
    S = 0
    for block in range(0, grid.nk, chunk):

        k = grid.k_cart[block:block+chunk] # only keep kpoints in this block (the last one may be shorter)

        _, eps, _, hv = compute_velocity(tb, k, mode) # ValueError for an unknown mode

        S += kubo_accumulate(eps, hv, hw, mu, eta, kT=kT) # (nw, 3, 3)

    # normalize once, with the size of the whole grid and the in-plane cell area
    A_cell = np.linalg.norm(np.cross(tb.lattice[:,0], tb.lattice[:,1]))

    return kubo_normalize(S, hw, grid.nk, A_cell)


def kubo_doped_finite_T_analytical(hw, mu, kT):
    """ mu measured wrt to E_D """

    assert kT > 0, 'k_BT must be positive to use this formula ! '
    sigma = 0.5 * (np.tanh((hw + 2*mu) / (4*kT)) + np.tanh((hw - 2*mu) / (4*kT)))
    return sigma