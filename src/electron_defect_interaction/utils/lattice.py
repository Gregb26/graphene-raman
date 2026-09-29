"""
lattice_utils.py
    Python module containing helper functions to perform operations on lattice quantities
"""

import numpy as np

def red_to_cart(x_red, X):
    """
    Converts the vector x in reduced coordinates to cartesian coordinates using the scaling matrix X:

    Inputs:
    -------
        x_red: (.., nd) array:
            Vector in reduced coordinates to convert to cartesian coordinates
        
        X: (nd, nd) array:
            Primitive vectors to use to convert to cartesian coordinates. For real space use A where A[:,i] is a_i. For 
            reciprocal space use B where B[:, i] is b_i

    Returns:
    --------
        x: (nd, ):
            Vector in cartesian coordinates
    """

    x_red = np.asarray(x_red, dtype=np.float64)   # (..., nd)
    X     = np.asarray(X,     dtype=np.float64)   # (nd, nd)

    x = np.einsum('...j,ij->...i', x_red, X, optimize=True)

    return x

def monkhorst_pack_grid(ngkpt, signed=True):
    """
    Generates a uniform (no symmetry reduction) Monkhorst-Pack kpoint grid in reduced coordinates from the tuple ngkpt.
    With signed=True the grid is folded into [-0.5, 0.5).
    Inputs:
        ngkpt: tuple (N1, N2, N3) of ints
            Number of kpoints in each direction
        signed: Bool:
            If True, fold [0,1) -> [-0.5, 0.5) (signed convention). Default is True.
    Returns:
        k_grid: (nk, 3) array of floats
            kpoint grid in reduced coordinates
    TODO: implement shift
    """

    N1, N2, N3 = ngkpt
    N = np.prod(ngkpt)
    k_grid = np.zeros((N,3), dtype=float)
    
    ik = 0
    for i3 in range(N3):
        k3 = i3 / N3
        for i2 in range(N2):
            k2 = i2 / N2
            for i1 in range(N1):
                k1 = i1 / N1
                k_grid[ik, :] = (k1, k2, k3)
                ik += 1
    if signed:
        k_grid = ((k_grid + 0.5) % 1.0) - 0.5 # fold [0, 1) -> [-0.5, 0.5)

    return k_grid

def build_k_path(high_sym_points, nk, B):
    """
    Builds a path in k-space through the high-symmetry points, with a number of points on each segment proportional to its
    Cartesian length (R8, 2026-09-28; the former version ignored nk and forced 100 points per segment whatever its length).

    The nk - 1 intervals are shared between the segments in proportion to their lengths |B (k_{i+1} - k_i)| (largest-remainder
    rounding, at least one interval per segment); every corner appears exactly once, at the index returned in idx.

    Inputs:
        high_sym_points: list of (label, k) pairs, k (3,) in reduced coordinates; the path joins them in order.
        nk: int, total number of k-points on the path (corners included), nk >= number of corners.
        B: (3, 3) array, reciprocal primitive vectors as columns (B[:, i] = b_i, the red_to_cart convention), any length unit.
    Returns:
        k: (nk, 3) floats, path in reduced coordinates.
        labels: list of the corner labels.
        idx: list of ints, index in k of each corner (idx[0] = 0, idx[-1] = nk - 1).
        s: (nk,) floats, cumulative Cartesian distance along the path (unit of B).
    """
    labels = [lab for lab, _ in high_sym_points]
    points = [np.asarray(p, dtype=np.float64) for _, p in high_sym_points]
    nseg = len(points) - 1
    if nk < nseg + 1:
        raise ValueError(f"nk = {nk} < number of corners ({nseg + 1})")

    lengths = np.array([np.linalg.norm(red_to_cart(points[i + 1] - points[i], B)) for i in range(nseg)])
    share = (nk - 1) * lengths / lengths.sum()
    n = np.maximum(np.floor(share).astype(int), 1)            # intervals per segment
    while n.sum() < nk - 1:                                   # largest remainders first
        n[np.argmax(share - n)] += 1
    while n.sum() > nk - 1:                                   # only if the minimum of one interval overshot
        n[np.argmax(np.where(n > 1, n - share, -np.inf))] -= 1

    ks, idx = [points[0][None]], [0]
    for i in range(nseg):
        t = np.arange(1, n[i] + 1) / n[i]                    # excludes the start corner, ends exactly on the next one
        ks.append((1 - t)[:, None] * points[i] + t[:, None] * points[i + 1])
        idx.append(idx[-1] + int(n[i]))
    k = np.vstack(ks)
    for i, j in enumerate(idx):                               # corners exactly as given (no rounding from the interpolation)
        k[j] = points[i]
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(red_to_cart(np.diff(k, axis=0), B), axis=1))])

    return k, labels, idx, s

import numpy as np

def generate_mp_grid(nk1, nk2, nk3, fold=True):
    """
    Generates a uniform, Gamma-centered, Monkhorst-Pack kpoint grid from the number of kpoints
    in each direction nki.
    """
    kpoints = []
    for i in range(nk1):
        for j in range(nk2):
            for k in range(nk3):
                k1 = i / nk1
                k2 = j / nk2
                k3 = k / nk3 if nk3 > 1 else 0.0 # handles 2d case
                if fold:
                    k1 -= round(k1)
                    k2 -= round(k2)
                    k3 -= round(k3)
                kpoints.append([k1, k2, k3])

    return np.array(kpoints)

def write_kpoints(nk1, nk2, nk3, print_weights=True, fold=True):
    """
    Write the kpoint list in a format suitable for Quantum-ESPRESSO input.
    """
    kpts = generate_mp_grid(nk1, nk2, nk3, fold)
    w = 1.0 / (nk1 * nk2 * nk3) # kpoint weight, constant because uniform sampling

    print("K_POINTS crystal")
    print(len(kpts))

    if print_weights:
        for k1, k2, k3 in kpts:
            print(f"{k1:.8f} {k2:.8f} {k3:.8f} {w:.8f}")
    else:
        for k1, k2, k3 in kpts:
            print(f"{k1:.8f} {k2:.8f} {k3:.8f}")


