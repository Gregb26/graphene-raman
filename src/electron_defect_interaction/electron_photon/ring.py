"""
Around a Dirac point: the resonant ring, k points where eps_c - eps_v = hbar omega (`ring`), the
Fermi velocity on a small circle (`fermi_velocity`), the ring statistics of the three velocity
variants (`ring_stats`, M3), the ring points in crystal coordinates for EM2 (`ring_kpoints_crystal`)
and the maps around K of the figure (`map_around_K`, M4).
"""

import numpy as np

from .velocity_operator import fourier, compute_velocity
from .tb_model import centres_only

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

def fermi_velocity(tb, K, mu, q=1e-3, ntheta=360, mode='berry'):
    """
    Fermi velocity on a circle of radius q around the Dirac point K (F13, M2), two ways. Intraband:
    |hbar v_nn| of pi and pi*, averaged over directions, then over the two bands (the electron-hole
    asymmetry, linear in q, cancels). Interband: |hbar v_cv|, equal to hbar v_F for a Dirac cone.
    v and c are the bands just below and above mu, the same all around the circle (asserted).

    Inputs:
        tb     : WannierTB
        K      : (3,) float, 1/Angstrom, Dirac point
        mu     : float, eV, between pi and pi* (E_D for Wannier90 data)
        q      : float, 1/Angstrom, radius of the circle
        ntheta : int, number of directions
        mode   : 'berry' or 'no_berry', passed to compute_velocity
    Returns:
        dict {'intra_avg', 'inter_avg', 'pi', 'pi_star'} of floats, eV*Angstrom
    """

    # compute circle of radius q around Dirac point K
    theta = np.linspace(0, 2*np.pi, ntheta, endpoint=False) # (nt,)
    k = K[None, :] + q*np.column_stack((np.cos(theta), np.sin(theta), np.zeros_like(theta))) # (nt, 3)

    # compute eigenvalues and velocity operator on the circle
    _, eps, _, hv = compute_velocity(tb, k, mode=mode) # (nt, nW), (nt, 3, nW, nW)
    # compute number of occupied bands per point on the circle
    n_occ = np.sum(eps < mu, axis=1) # (nt, )
    assert np.all(n_occ == n_occ[0])

    n = n_occ[0]
    v = n - 1 # index of valence band on circle
    c = n # conduction band

    # intraband pi-pi
    hv_xvv = hv[:, 0, v, v].real; hv_yvv = hv[:,1, v, v].real # (nt,)
    hv_vv = np.sqrt(hv_xvv**2 + hv_yvv**2) # (nt, )
    pi = np.mean(hv_vv)

    # intraband pi_star-pi_star
    hv_xcc = hv[:, 0, c, c].real; hv_ycc = hv[:,1, c, c].real # (nt,)
    hv_cc = np.sqrt(hv_xcc**2 + hv_ycc**2) # (nt, )
    pi_star = np.mean(hv_cc)

    # intraband average
    intra_avg = np.mean([pi, pi_star])

    # interband pi-pi_star (same as pi_star-pi)
    hv_xcv = hv[:, 0, c, v]; hv_ycv = hv[:,1, c, v] # (nt,)
    hv_cv = np.sqrt(np.abs(hv_xcv)**2 + np.abs(hv_ycv)**2)
    inter_avg = np.mean(hv_cv)

    return {"intra_avg": intra_avg, "inter_avg": inter_avg, "pi": pi, "pi_star": pi_star}

def ring_stats(tb, K, mu, hws, ntheta=720):
    """
    In-plane interband velocity on the resonant rings, for the three variants of the thesis: 'full'
    (tb, berry), 'centres_only' (centres_only(tb), berry) and 'no_berry'. One ring per hw, shared by
    the three (same k); c and v chosen by energy around mu. A perfect cone gives avg = 1, node = 0.

    Inputs:
        tb     : WannierTB
        K      : (3,) float, 1/Angstrom, Dirac point
        mu     : float, eV, chemical potential (E_D for the real data)
        hws    : float or list of float, eV, photon energies
        ntheta : int, angles per ring (multiple of 6)
    Returns:
        dict hw -> variant -> {
            'avg'  : (2,) <|hbar v^x_cv|^2>, <|hbar v^y_cv|^2> over theta, in units of (hbar v_F)^2 / 2
            'node' : (2,) degrees, angle of the zero of |e.hbar v_cv|^2 minus that of e, for e_x, e_y;
                     searched in the half ring facing +e, resolution 360/ntheta
            'ratio': (min, max) over theta of |hbar v_cv| / |hbar v_cv|_full, in-plane norms }
    """
    hws = np.atleast_1d(hws)
    # compute Fermi velocity on circle of radius q around Dirac point
    hv_F = fermi_velocity(tb, K, mu, q=1e-3, ntheta=ntheta)['intra_avg']
    hv_F2 = (hv_F)**2 / 2
    # construct each variant
    variants = [('full', tb, 'berry'), ('centres_only', centres_only(tb), 'berry'), ('no_berry', tb, 'no_berry')]

    stats_dict = {}
    # loop over laser energies
    for hw in hws:
        # construct resonant ring per energy
        k, theta, _, _ = ring(tb, K, hw, mu, ntheta=ntheta) # (ntheta, 3)

        name_dict = {}
        # loop over variants
        for name, model, mode in variants:
            # compute eigenvalues and velocity operator on ring, for each model
            _, eps, _, hv = compute_velocity(model, k, mode=mode) # (ntheta, nW), (ntheta, 3, nW, nW)
            # compute number of occupied bands per point on the circle
            n_occ = np.sum(eps < mu, axis=1) # (nt, )
            n = n_occ[0]
            assert np.all(n_occ == n)

            v = n - 1 # valence band index
            c = n     # conduction band index

            # Compute P = |\hbar v^{x,y}_cv|^2 of shape (ntheta, 2)
            P = np.abs(hv[:, :2, c, v])**2 # (ntheta, 2)
            P_avg = np.mean(P, axis=0) / hv_F2 # (2,)

            if  name == 'full':
                P_full = P # store average for ratio below

            # loop over polarisations (angles phi)
            phis = [0, np.pi/2]
            node = np.zeros(len(phis))
            for i, phi in enumerate(phis):
                # compute deviation and bring back into (-pi, pi]
                d = (theta - phi - np.pi) % (2*np.pi) - np.pi

                mask = np.abs(d) < np.pi / 2 # keep only |d| < pi/2
                id = np.argmin(P[mask][:,i]) 

                delta = np.degrees(d[mask][id])
                node[i] = delta

            ratios = np.sqrt(P[:,0] + P[:,1]) / np.sqrt(P_full[:,0] + P_full[:,1])
            ratio = (min(ratios), max(ratios))

            name_dict[name] = {
                'avg': P_avg, # (2,)
                'node': node, # (2,)
                'ratio': ratio # (min, max)
            }
        stats_dict[hw] = name_dict

    return stats_dict

def ring_kpoints_crystal(tb, K, mu, npoints):
    """
    k points of the resonant rings in crystal coordinates, k_red_i = a_i.k / 2 pi (tb.lattice), for
    the direct DFT check of EM2 (QE `K_POINTS crystal`, same cell as Wannier90).

    Inputs:
        tb      : WannierTB
        K       : (3,) float, 1/Angstrom, Dirac point
        mu      : float, eV, chemical potential (E_D for the real data)
        npoints : dict hw (eV) -> number of points on that ring (multiple of 6), e.g. {2.33: 48}
    Returns:
        dict hw -> (k_red (n, 3), theta (n,) rad), the points of `ring` with ntheta = n
    """
    kpoints = {}
    for w, nk in npoints.items():
        k, theta, _, _ = ring(tb, K, w, mu, ntheta=nk) # (nk, 3), (nk,)
        k_red = k @ tb.lattice / (2*np.pi) # (nk, 3)

        kpoints[w] = (k_red, theta)

    return kpoints

def map_around_K(tb, K, mu, half_width, nq, mode='berry'):
    """
    Maps around K for panel (a) of the figure: Delta eps = eps_c - eps_v and |hbar v^{x,y}_cv|^2 on the
    square K + (qx, qy, 0), qx and qy = np.linspace(-half_width, half_width, nq). The iso-lines
    Delta eps = hw are the rings of `ring`; |hbar v^x_cv|^2 vanishes on a dark line through K (the node
    of `ring_stats`). One call to compute_velocity: ~1 GB at nq = 300 (see `fourier`).

    Inputs:
        tb         : WannierTB
        K          : (3,) float, 1/Angstrom, Dirac point
        mu         : float, eV, chemical potential (E_D for the real data)
        half_width : float, 1/Angstrom
        nq         : int, points per axis, even (K itself, degenerate, stays off the grid)
        mode       : 'berry' or 'no_berry' (centres only: pass centres_only(tb))
    Returns:
        deps : (nq, nq) float, eV, indexed [iy, ix] (plot with contour(q, q, deps))
        P    : (nq, nq, 2) float, eV^2 Angstrom^2, |hbar v^x_cv|^2 and |hbar v^y_cv|^2, same indexing
    """
    assert nq % 2 == 0, 'nq must be even'

    h = half_width
    q = np.linspace(-h, h, nq)
    Qx, Qy = np.meshgrid(q,q, indexing='xy')

    Q = np.column_stack((Qx.ravel(), Qy.ravel(), np.zeros_like(Qx.ravel())))
    k = K + Q

    _, eps, _, hv = compute_velocity(tb, k, mode=mode) # (nq^2, nW), (nq^2, 3, nW, nW)

    # find valence and conduction bands indices
    n_occ = np.sum(eps < mu, axis=1)
    n = n_occ[0]
    assert np.allclose(n_occ, n)
    v = n - 1 # valence
    c = n     # conduction

    deps = eps[:, c] - eps[:, v] # (nq^2, )
    deps = deps.reshape(nq, nq)
    
    P = np.abs(hv[:, :2, c, v])**2 # (nq^2, 2)
    P = P.reshape(nq, nq, -1)

    return deps, P 
    