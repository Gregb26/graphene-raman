"""
Optical conductivity sigma(omega)/sigma_0 by the Kubo formula, accumulated block by block (`kubo_accumulate`,
`kubo_normalize`), and the driver `sigma_on_grid`: absorptive part by default (Gaussian, 1/hw prefactor, M4
production), full complex sigma with a complex kernel (`gaussian_complex`, `lorentzian_complex`; perspective B).
Doping and temperature through mu and kT, with the Dirac-cone references `kubo_doped_finite_T_analytical` and
`kubo_doped_complex_analytical`.
"""

import numpy as np

from .kgrid import make_grid_tb
from .velocity_operator import compute_velocity
from scipy.special import wofz

def gaussian_eta(x, eta):
    """
    Normalized Gaussian of standard deviation eta (not the FWHM = 2.355 eta), in 1/eV: the broadened
    delta(hw - gap) of the Kubo sum.
    """
    return np.exp(-(x**2)/(2*eta**2)) / (eta*np.sqrt(2*np.pi))

def lorentzian_complex(x, eta):
    """
    Complex Lorentzian kernel K(x) = 1/(x + i eta), in 1/eV: Re[iK] = eta/(x^2 + eta^2) (pi times the normalized
    Lorentzian of half-width eta), Im[iK] = x/(x^2 + eta^2) -> 1/x. Its 1/x^2 tails bias Re sigma at O(eta);
    it gives the closed-form Dirac-cone reference. K(-x) = -K(x)*.
    """
    return 1 / (x + 1j*eta)

def gaussian_complex(x, eta):
    """
    Complex Gaussian kernel K(x) = -i sqrt(pi)/(sqrt(2) eta) w(x/(sqrt(2) eta)), w the Faddeeva function (scipy
    wofz), in 1/eV: Re[iK] = pi gaussian_eta(x, eta) exactly (eta = standard deviation), Im[iK] its Hilbert
    transform (Dawson function), -> 1/x for |x| >> eta. K(-x) = -K(x)*.
    """
    t = x / (np.sqrt(2)*eta)
    return -1j*(np.sqrt(np.pi) / (np.sqrt(2)*eta)) * wofz(t)

def kubo_accumulate(eps, hv, hw, mu=0, eta=0.04, kT=0, kernel=None):
    """
    Partial Kubo sum over one block of k points, over the pairs m < n with w = f(eps_m) - f(eps_n) > 0
    (step at kT = 0, Fermi-Dirac above), Delta = eps_n - eps_m > 0 and W_ab = Re[(hbar v^a_nm)* hbar v^b_nm]:

        kernel None : S_ab(w) = sum_k sum_{m<n} w W_ab g_eta(hw - Delta)                        (absorptive)
        kernel K    : S_ab(w) = sum_k sum_{m<n} (w/Delta) W_ab [K(hw - Delta) + K(hw + Delta)]  (complex)

    The complex form needs 1/Delta (1/hw would give a spurious 1/omega in Im sigma) and the antiresonant
    term K(hw + Delta) (Im sigma, sigma(-w) = sigma(w)*). Interband only (no Drude term). Blocks add up;
    normalize once at the end with `kubo_normalize` and the same kernel.

    Inputs:
        eps    : (nk_block, nW) float, eV, from `velocity`
        hv     : (nk_block, 3, nW, nW) complex, eV*Angstrom, from `velocity`
        hw     : (nw,) float or scalar, eV, photon energies
        mu     : float, eV, chemical potential
        eta    : float, eV, broadening (standard deviation of the Gaussians, half-width of the Lorentzian)
        kT     : float, eV, k_B T (0.02585 eV at 300 K); 0 gives the step
        kernel : None, or a complex kernel K(x, eta): gaussian_complex, lorentzian_complex
    Returns:
        S : (nw, 3, 3), eV*Angstrom^2, float (kernel None) or complex (block contribution, not normalized)
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
    deltas = eps[ik, ic] - eps[ik, iv] # (nt,) optical gaps at each k
    hv_cv = hv[ik, :, ic, iv] # (nt, 3) velocity operator for allowed transitions
    w_cv  = w[ik, ic, iv] # (nt, ) transition weights

    W = np.real(np.conj(hv_cv[:, :, None]) * hv_cv[:, None, :]) # (nt, 3, 3) # Re[(ħv^α_cv)* ħv^β_cv] 
    W_ = W.reshape(nt, 9) # (nt, 9)

    # compute kernel
    if kernel is None:  
        K = gaussian_eta(hw[:, None] - deltas[None, :], eta=eta) # (nw, nt)
        W_weighted = W_ * w_cv[:, None]

    else:
        K = kernel(hw[:, None] - deltas[None, :], eta=eta) + kernel(hw[:, None] + deltas[None, :], eta=eta) # (nw, nt)
        W_weighted = W_ * w_cv[:, None] / deltas[:, None]

    S  = K @ W_weighted # (nw, 9)

    return S.reshape(nw, 3, 3)

def kubo_normalize(S, hw, nk, A_cell, kernel=None):
    """
    sigma_ab / sigma_0, sigma_0 = e^2 / (4 hbar), spin included, from the summed S:

        kernel None : 8 pi S_ab / (hw N_k A_cell), absorptive, 1 for a perfect Dirac cone (EM.md, M0)
        kernel K    : 8 i S_ab / (N_k A_cell), complex (hw unused)

    Inputs:
        S      : (nw, 3, 3), eV*Angstrom^2, sum of `kubo_accumulate` over all blocks (same kernel)
        hw     : (nw,) float, eV
        nk     : int, number of k points of the WHOLE grid (not of a block)
        A_cell : float, Angstrom^2, in-plane cell area |a1 x a2|
        kernel : None or the complex kernel given to `kubo_accumulate`
    Returns:
        (nw, 3, 3) float or complex, sigma / sigma_0
    """
    if kernel is None:
        return 8 * np.pi * S / (hw[:, None, None] * nk * A_cell)

    else:
        return 8j * S / (nk * A_cell)
    
def sigma_on_grid(tb, N, hw, mu=0, eta=0.04, mode='berry', kT=0, chunk=int(1e5), kernel=None):
    """
    Optical conductivity sigma(omega)/sigma_0 on the shifted N x N grid (driver of M0 and M4): block by
    block compute_velocity -> kubo_accumulate, then kubo_normalize with N_k = N^2; independent of
    `chunk`. Converges in N much faster than the grid step suggests (EM.md, M4 pilot). Doping and
    temperature through mu and kT (Pauli blocking below hw = 2|mu - E_D|). With a complex kernel, undoped
    at kT = 0, the 1/Delta weight is singular at the Dirac point: Im sigma converges only as 1/N (kT > 0
    or doping removes it; EM.md, perspective B).

    Inputs:
        tb     : WannierTB
        N      : int, grid of N x N points over the whole Brillouin zone (both valleys)
        hw     : (nw,) float or scalar, eV, photon energies
        mu     : float, eV, chemical potential (E_D for the real data)
        eta    : float, eV, broadening (standard deviation; half-width for the Lorentzian)
        mode   : 'berry' or 'no_berry', passed to `compute_velocity`
        kT     : float, eV, k_B T >= 0, passed to `kubo_accumulate` (0: T = 0)
        chunk  : int, number of k points per block (memory: see `fourier`; complex kernels double it)
        kernel : None (absorptive part, M4 production) or a complex kernel (full sigma, perspective B)
    Returns:
        (nw, 3, 3) float or complex, sigma_ab / sigma_0
    """
    assert kT >= 0, 'k_BT must be positive or zero !'
    hw = np.atleast_1d(hw)

    grid = make_grid_tb(tb, N)  

    # loop over blocks and accumulate
    S = 0
    for block in range(0, grid.nk, chunk):

        k = grid.k_cart[block:block+chunk] # only keep kpoints in this block (the last one may be shorter)

        _, eps, _, hv = compute_velocity(tb, k, mode) # ValueError for an unknown mode

        S += kubo_accumulate(eps, hv, hw, mu, eta, kT=kT, kernel=kernel) # (nw, 3, 3)

    # normalize once, with the size of the whole grid and the in-plane cell area
    A_cell = np.linalg.norm(np.cross(tb.lattice[:,0], tb.lattice[:,1]))

    return kubo_normalize(S, hw, grid.nk, A_cell, kernel=kernel)

def kubo_doped_finite_T_analytical(hw, mu, kT):
    """
    Interband sigma/sigma_0 of a perfect Dirac cone at chemical potential mu and temperature kT,

        0.5 [tanh((hw + 2 mu) / 4kT) + tanh((hw - 2 mu) / 4kT)],

    i.e. f(-hw/2) - f(hw/2) at the resonance: Pauli blocking below hw = 2|mu|. Reference of the
    finite-T tests (EM.md, perspective A).

    Inputs:
        hw : (nw,) float or scalar, eV
        mu : float, eV, measured from the Dirac point
        kT : float, eV, > 0
    Returns:
        (nw,) float, sigma / sigma_0
    """

    assert kT > 0, 'k_BT must be positive to use this formula ! '
    sigma = 0.5 * (np.tanh((hw + 2*mu) / (4*kT)) + np.tanh((hw - 2*mu) / (4*kT)))
    return sigma

def kubo_doped_complex_analytical(hw, mu, eta):
    """
    Complex interband sigma/sigma_0 of a perfect Dirac cone at T = 0 with the Lorentzian kernel
    (`lorentzian_complex`, half-width eta), mu measured from the Dirac point:

        1 - (i/pi) [ln(hw + 2|mu| + i eta) - ln(hw - 2|mu| + i eta)]

    As eta -> 0: Re -> Pauli step at hw = 2|mu| (1/2 + eta/(4 pi |mu|) on the edge), Im -> -(1/pi)
    ln|(hw + 2mu)/(hw - 2mu)|; exactly 1 at mu = 0. Reference of the complex tests (EM.md, perspective B).

    Inputs:
        hw  : (nw,) float or scalar, eV
        mu  : float, eV, measured from the Dirac point
        eta : float, eV, > 0
    Returns:
        (nw,) complex, sigma / sigma_0
    """
    return 1 - (1j/np.pi) * (np.log(hw + 2*np.abs(mu) +1j*eta) - np.log(hw - 2*np.abs(mu) + 1j*eta))
