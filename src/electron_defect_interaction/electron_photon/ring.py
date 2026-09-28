"""
Resonant ring around a Dirac point: the k points where eps_c - eps_v = hbar omega (`ring`).
"""

import numpy as np

from .velocity_operator import fourier

def _gap(tb, K, q, e, mu):
    """
    Gap eps_c - eps_v at k = K + q e along each ray, for trial distances q: the function whose
    root `ring` looks for. c and v are the bands just above and below mu, found by energy (fixed
    indices would fail with the 5-band data). Never called at q = 0, where the bands are degenerate
    at mu.

    Inputs:
        tb : WannierTB
        K  : (3,) float, 1/Angstrom
        q  : (ntheta,) float, 1/Angstrom, trial distance per ray
        e  : (ntheta, 3) float, unit vectors of the rays
        mu : float, eV
    Returns:
        gap : (ntheta,) float, eV
    """

    # one k point per ray
    k = K[None, :] + q[:, None] * e # (ntheta, 3)
    H_k = fourier(tb.H_R, tb.R_cart, tb.ndegen, k) # (ntheta, nW, nW)
    eps = np.linalg.eigvalsh(H_k) # (ntheta, nW), ascending
    n_occ = np.sum(eps < mu, axis=1) # (ntheta,), number of occupied states (with energy less than chemical potential) per k
    v = n_occ - 1 # index of valence band
    c = n_occ  # index of conduction band

    # eps[i, v[i]] row by row; [:, 0] turns (ntheta, 1) into (ntheta,)
    eps_v = np.take_along_axis(eps, v[:, None], axis=1)[:,0]
    eps_c = np.take_along_axis(eps, c[:, None], axis=1)[:,0]

    gap = eps_c - eps_v # (ntheta, )

    return gap

def ring(tb, K, hw, mu=0.0, ntheta=720, tol=1e-12, maxsteps=100):
    """
    Resonant ring around K: one point per direction theta where eps_c - eps_v = hw,

        k(theta) = K + q(theta) (cos theta, sin theta, 0).

    These states absorb light of energy hw exactly (no broadening); they serve for the angular
    checks of the Berry term, the ring averages of M3, the figure of section 2.5 and the k list of
    EM2. A circle of radius q0 = hw / (2 hbar v_F) for a perfect cone, a rounded triangle with
    trigonal warping. Vectorized bisection on [0, 2 q0], valid while hw stays well below the gap 2t
    at M (checked by the bracket assert): ~39 steps, gap exact to ~1e-11 eV.

    Inputs:
        tb       : WannierTB (t and a_cc only set q0)
        K        : (3,) float, 1/Angstrom, centre of the ring (Dirac point)
        hw       : float, eV, photon energy
        mu       : float, eV, chemical potential (0 in M0, E_D in M1)
        ntheta   : int, multiple of 6 (theta + pi and theta + 120 degrees stay on the grid)
        tol      : float, 1/Angstrom, final bisection width
        maxsteps : int, cap on the number of bisection steps
    Returns:
        k     : (ntheta, 3) float, 1/Angstrom
        theta : (ntheta,) float, rad, equally spaced, endpoint excluded (plain means are angular averages)
        q     : (ntheta,) float, 1/Angstrom, |k - K|
        q0    : float, 1/Angstrom, Dirac estimate hw / (2 hbar v_F)
    """
    assert ntheta % 3 == 0 and ntheta % 2 == 0, 'number of angles must be a multiple of 2 and 3'

    # directions: unit vectors (cos theta, sin theta, 0), one per row
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False) # (ntheta,)
    e = np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))) # (ntheta, 3)

    # bracket: q = 0 (gap 0 < hw, known, not evaluated) and twice the Dirac estimate
    hv_F = 3*tb.t * tb.a_cc / 2
    q0 = hw / (2*hv_F) # float
    lo = np.zeros(ntheta); hi = 2*q0*np.ones(ntheta) # (ntheta, )
    assert np.all(_gap(tb, K, hi, e, mu) > hw ), 'hw too large: gap below hw at 2 q0 on some ray (ring not closed?)'

    # bisection on all rays at once: gap(mid) < hw -> crossing further out (lo <- mid), else hi <- mid
    nsteps = 0
    while (hi - lo).max() > tol and nsteps < maxsteps:
        mid = (hi + lo) / 2
        gap_mid = _gap(tb, K, mid, e, mu)
        cond = gap_mid < hw # (ntheta,) bool

        lo = np.where(cond, mid, lo); hi = np.where(cond, hi, mid)
        nsteps += 1
    assert nsteps < maxsteps, 'bisection did not converge within maxsteps'

    q = (hi + lo) / 2
    k = K[None, :] + q[:, None]*e # (ntheta, 3)

    assert k.shape == (ntheta, 3), 'wrong shape for k'
    assert theta.shape == (ntheta,), 'wrong shape for theta'
    assert q.shape == (ntheta,), 'wrong shape for q'

    return k, theta, q, q0
