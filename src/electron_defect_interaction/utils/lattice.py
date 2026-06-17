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
    Generates a uniform (no symmetry reduction) Monkhorst-Pack kpoint grid in reduced coordinates from the tuple ngkpt. Returns the same grid
    as Abinit if signed=True
    Inputs:
        ngkpt: tuple (N1, N2, N3) of ints
            Number of kpoints in each direction
        signed: Bool:
            If True, fold [0,1) -> [-0.5, 0.5). Signed convention Abinit uses. Default is True.
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

def build_k_path(high_sym_points, nk):
    """
    Builds a path in kspace joining the points specified in high_sym_points
    Inputs:
        high_sym_points: list of np (label, k) arrays, list of high symmetry points in kspace. The path joins this points.
        nk: int, number of kpoints to use between the high symmetry points
    Ouputs:
        kpath: (np*nk, 3), path in kspace
    """

    labels = [lab for lab, _ in high_sym_points]
    points = [pts for _, pts in high_sym_points]

    nk=100
    ks = []
    idx = []
    count=0
    t = np.linspace(0, 1, nk)
    for i in range(len(points)-1):
        seg = (1 -t)[:, None]*points[i] + t[:, None]*points[i+1]
        if i > 0:
            seg=seg[1:] # avoid duplicate at junction
            count-= 1
        ks.append(seg)
        idx.append(count)
        count+=len(seg)
        
    k = np.vstack(ks)
    idx.append(len(k)-1)

    return k, labels, idx

import numpy as np

def generate_mp_grid(nk1, nk2, nk3):
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
                kpoints.append([k1, k2, k3])

    return np.array(kpoints)

def write_kpoints(nk1, nk2, nk3, print_weights=True):
    """
    Write the kpoint list in a format suitable for Quantum-ESPRESSO input.
    """
    kpts = generate_mp_grid(nk1, nk2, nk3)
    w = 1.0 / (nk1 * nk2 * nk3) # kpoint weight, constant because uniform sampling

    print("K_POINTS crystal")
    print(len(kpts))

    if print_weights:
        for k1, k2, k3 in kpts:
            print(f"{k1:.8f} {k2:.8f} {k3:.8f} {w:.8f}")
    else:
        for k1, k2, k3 in kpts:
            print(f"{k1:.8f} {k2:.8f} {k3:.8f}")


