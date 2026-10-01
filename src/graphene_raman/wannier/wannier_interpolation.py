"""
wannier_interpolation.py
    Python module containing functions to interpolate objects computed on a coarse kpoint grid onto a fine kpoint grid
    via Maximally Localized Wannier Functions. 
"""

import numpy as np
from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk
from graphene_raman.io.wannier_io import read_w90_mat, read_w90_tb

def _rotate(M, X):
    """
    Gauge rotation at both k-points, M'(k', k) = X(k')^dag M(k', k) X(k).
    Inputs:
        M: (nb, nk, nb, nk) complex, index [bra band, k', ket band, k].
        X: (nk, nb, nw) complex, rotation matrix at each k.
    Returns:
        (nw, nk, nw, nk) complex, contiguous.
    """
    X_dag = X.transpose(0, 2, 1).conj()

    M_ = X_dag[:, None, ...] @ M.transpose(1, 3, 0, 2) @ X[None, ...]

    return np.ascontiguousarray(M_.transpose(2, 0, 3, 1))
    

def Mbk_to_Mwk(Mbk, U, U_dis=None):
    """
    Bloch gauge -> Wannier gauge, Mwk(k', k) = V(k')^dag Mbk(k', k) V(k) with V = U_dis U (or U alone); checks that
    V V^dag is a Hermitian projector of rank nw.
    Inputs:
        Mbk:   (nb, nk, nb, nk) complex, index [bra band, k', ket band, k].
        U:     (nk, nw, nw) complex with U_dis, or (nk, nb, nw) without, Wannier90 gauge matrix.
        U_dis: (nk, nb, nw) complex or None, disentanglement matrix.
    Returns:
        Mwk: (nw, nk, nw, nk) complex.
    """

    # Entangled case, rotation matrix is V = U_dis @ U
    if U_dis is not None:

        V = U_dis @ U # (nk, nb, nw)

        # V must be an isometry: P = V V^dag is a Hermitian projector of rank nw
        P = V @ V.conj().transpose(0, 2, 1) # (nk, nb, nb)

        assert np.allclose(P, P.conj().transpose(0, 2, 1), atol=1e-10)     
        assert np.allclose(P @ P, P, atol=1e-8)           
        assert np.allclose(np.trace(P, axis1=1, axis2=2), V.shape[-1], 1e-10)

    else:
        V = U # (nk, nb, nw) 
    
    Mwk = _rotate(Mbk, V)

    return Mwk

def Mwk_to_Mwr(Mwk, k_red, MP_grid):
    """
    Double Fourier transform k -> R on a full MP grid, Mwr(R, R') = (1/nk^2) sum_{k', k} e^{+2 pi i k'.R} Mwk(k', k)
    e^{-2 pi i k.R'}, on the box arange(N) - N//2 per axis (labels defined modulo N; no ndegen, no Wigner-Seitz).
    Inputs:
        Mwk:     (nw, nk, nw, nk) complex, Wannier gauge.
        k_red:   (nk, 3) floats, the full unshifted MP grid (reduced coordinates), any order.
        MP_grid: (3,) ints, (N1, N2, N3) with N1 N2 N3 = nk (no symmetry reduction).
    Returns:
        Mwr: (nw, nR, nw, nR) complex, nR = nk.
        R:   (nR, 3) ints, the R labels ('ij' order).
    """

    nk = Mwk.shape[1]

    # Build the grid in real space based on the kpoint grid
    N1, N2, N3 = MP_grid
    def centered(N):
        return np.arange(N) - (N//2)
    r1, r2, r3 = centered(N1), centered(N2), centered(N3)
    rr1, rr2, rr3 = np.meshgrid(r1, r2, r3, indexing='ij')
    R = np.stack((rr1, rr2, rr3), axis=-1).reshape(-1, 3)

    # Compute phase
    bra = np.exp(2j*np.pi * (k_red @ R.T)) # (nk, nr)
    ket = bra.conj() # (nk, nr)

    # M is (nw', nk', nw, nk), want (nr', nk') @ (nw', nw, nk', nk) @ (nk, nr)
    Mwr = bra.T @ Mwk.transpose(0, 2, 1, 3) @ ket / (nk**2) # (nw', nw, nr', nr)

    # want (nw', nr', nw, nr)
    return np.ascontiguousarray(Mwr.transpose(0, 2, 1, 3)), R

def Mwr_to_Mwk(Mwr, R, k, ws=None):
    """
    Inverse transform R -> k at any k-points (square case of Mwr_to_Mwk_pairs), Mwk(k', k) = sum_{R, R'} e^{-2 pi i k'.R}
    Mwr(R, R') e^{+2 pi i k.R'} (no 1/nk^2).
    Inputs:
        Mwr: (nw, nR, nw, nR) complex.
        R:   (nR, 3) ints, labels of Mwr.
        k:   (nk, 3) floats, k-points (reduced coordinates), on or off the MP grid.
        ws:  dict or None, Wigner-Seitz images of R (ws_images), for off-grid k.
    Returns:
        Mwk: (nw, nk, nw, nk) complex.
    """

    return Mwr_to_Mwk_pairs(Mwr, R, k, k, ws=ws)

def Mwr_to_Mwk_pairs(Mwr, R, k_bra, k_ket, ws=None):
    """
    Mwr_to_Mwk with independent bra and ket k-lists, Mwk[w, k', W, k] = sum_{R, R'} e^{-2 pi i k'.R} Mwr[w, R, W, R']
    e^{+2 pi i k.R'}.
    Inputs:
        Mwr:   (nw, nR, nw, nR) complex.
        R:     (nR, 3) ints, labels of Mwr.
        k_bra: (nk', 3) floats, bra k-points k'.
        k_ket: (nk, 3) floats, ket k-points k.
        ws:    dict or None, Wigner-Seitz images of R (ws_images); phases ws_phase(k, ws, nR, -/+1).
    Returns:
        Mwk: (nw, nk', nw, nk) complex.
    """

    R = np.asarray(R, dtype=float)
    nR = R.shape[0]

    if ws is not None:
        assert len(ws["dist"]) == nR, "ws was built on other R labels than Mwr (len(ws['dist']) != len(R))"
        bra = ws_phase(k_bra, ws, nR, -1)
        ket = ws_phase(k_ket, ws, nR, +1)

    else:
        bra = np.exp(-2j*np.pi * (np.asarray(k_bra, dtype=float) @ R.T)) # (nk', nr)
        ket = np.exp(+2j*np.pi * (np.asarray(k_ket, dtype=float) @ R.T)) # (nk, nr)

    # M is (nw', nr', nw, nr), want (nk', nr') @ (nw', nw ,nr', nr) @ (nr, nk)
    Mwk = bra @ Mwr.transpose(0, 2, 1, 3) @ ket.T # (nw', nw, nk', nk)

    # want (nw', nk', nw, nk)
    return np.ascontiguousarray(Mwk.transpose(0, 2, 1, 3))


def Mwk_to_Mbk(Mwk, Hwr, Rw, k, ndegen=None):
    """
    Wannier gauge -> smooth Bloch gauge (eigenbasis of the interpolated H(k)), Mbk(k', k) = U(k')^dag Mwk(k', k) U(k).
    Inputs:
        Mwk:    (nw, nk, nw, nk) complex.
        Hwr:    (nRw, nw, nw) complex, H(R) of the Wannier model (eV).
        Rw:     (nRw, 3) ints, its R vectors.
        k:      (nk, 3) floats, the k-points of Mwk.
        ndegen: (nRw,) ints or None, Wigner-Seitz degeneracies.
    Returns:
        Mbk: (nw, nk, nw, nk) complex, band basis.
    """

    # compute rotation matrix from Wannier Hamiltonian
    _, _, Uwk = Hwr_to_Hwk(Hwr, Rw, k, ndegen=ndegen) # (nk, nw, nw)

    Mbk = _rotate(Mwk, Uwk)

    return Mbk

def _match_kpoint_order(k_from, k_to, tol=1e-5):
    """
    Permutation `perm` such that k_from[perm] == k_to (compared modulo 1, with periodic distance).
    Used to bring Wannier90's k-ordering (in the .mat files) into the ordering of the coarse grid
    on which M was computed.
    """
    kf = np.mod(k_from, 1.0)
    kt = np.mod(k_to, 1.0)
    perm = np.empty(len(kt), dtype=int)
    used = np.zeros(len(kf), dtype=bool)
    for i, kk in enumerate(kt):
        d = np.abs(kf - kk)
        d = np.minimum(d, 1.0 - d)          # periodic distance per component
        dist = d.sum(axis=1)
        dist[used] = np.inf                 # one-to-one matching
        j = int(np.argmin(dist))
        if dist[j] > tol:
            raise ValueError(f"k-point {kk} from the target grid has no match in the U-matrix grid")
        perm[i] = j
        used[j] = True
    return perm

def _infer_mp_grid(k):
    """
    Infer the (N1, N2, N3) Monkhorst-Pack size from an unshifted Gamma-centered grid by counting the
    distinct reduced coordinates along each axis. Raises if N1*N2*N3 != nk (i.e. not a full MP grid).
    """
    grid = []
    for ax in range(3):
        v = np.mod(np.round(k[:, ax], 6), 1.0)   # fold to [0,1), killing -1e-15 -> 0.999... noise
        v[v > 1.0 - 1e-4] = 0.0                   # treat ~1 as 0 (periodic wrap)
        v = np.sort(v)
        uniq = [v[0]]
        for x in v[1:]:
            if abs(x - uniq[-1]) > 1e-4:
                uniq.append(x)
        grid.append(len(uniq))
    grid = tuple(grid)
    if np.prod(grid) != len(k):
        raise ValueError(f"inferred MP grid {grid} ({int(np.prod(grid))} pts) != nk={len(k)}; "
                         "k_coarse must be a full, unshifted, Gamma-centered MP grid")
    return grid

def wannier_interpolate(M, k_coarse, k_fine, wannier_tb, u_path, u_dis_path=None):
    """
    Interpolates a matrix M computed on a coarse kpoint grid onto a fine kpoint grid (or k-path) using
    Maximally Localized Wannier Function interpolation. The chain is

        M(b,k) --[V^dag . V : Mbk_to_Mwk]--> M_wk --[double FT: Mwk_to_Mwr]--> M_wr
              --[inverse FT to fine grid: Mwr_to_Mwk]--> M_wk(fine)
              --[U(fine) from H(R): Mwk_to_Mbk]--> M_bk(fine)

    where V = U_dis @ U is the (possibly disentangling) Wannier gauge matrix and U(fine) is obtained by
    Fourier-interpolating and diagonalising the tight-binding Hamiltonian H(R).

    Inputs:
        M:          (nb, nk, nb, nk) array of complex, object on the coarse grid to interpolate.
        k_coarse:   (nk, 3) array of floats, coarse kpoint grid in reduced coords. Must be an UNSHIFTED
                    Gamma-centered MP grid (it is the dual of the R grid used for the double FT).
        k_fine:     (nkf, 3) array of floats, fine kpoint grid (or path) in reduced coords. Any grid.
        wannier_tb: str, path to the Wannier90 tight-binding Hamiltonian file (seedname_tb.dat).
        u_path:     str, path to the Wannier90 U matrix (seedname_u.mat).
        u_dis_path: str or None, path to the disentanglement matrix (seedname_u_dis.mat), if used.
    Returns:
        M_bk_fine:  (nw, nkf, nw, nkf) array of complex, M interpolated onto k_fine in the smooth Bloch
                    (Hamiltonian eigenstate) gauge.
        E_fine:     (nkf, nw) array of floats, Wannier-interpolated band energies on k_fine.
    """
    # 1. Wannier gauge matrices, reordered from Wannier90's k-order to the coarse-grid order
    U, k_U = read_w90_mat(u_path)
    U = U[_match_kpoint_order(k_U, k_coarse)]
    U_dis = None
    if u_dis_path is not None:
        U_dis, k_Ud = read_w90_mat(u_dis_path)
        U_dis = U_dis[_match_kpoint_order(k_Ud, k_coarse)]

    # 2. tight-binding Hamiltonian H(R) + Wigner-Seitz degeneracies
    Hwr, Rw, ndegen, _, _ = read_w90_tb(wannier_tb)

    # 3. coarse Bloch -> Wannier gauge -> real space -> fine grid -> smooth Bloch gauge
    MP_grid = _infer_mp_grid(k_coarse)
    Mwk = Mbk_to_Mwk(M, U, U_dis)                                  # (nw, nk, nw, nk)
    Mwr, R = Mwk_to_Mwr(Mwk, k_coarse, MP_grid)                    # (nw, nr, nw, nr)
    Mwk_fine = Mwr_to_Mwk(Mwr, R, k_fine)                          # (nw, nkf, nw, nkf)
    M_bk_fine = Mwk_to_Mbk(Mwk_fine, Hwr, Rw, k_fine, ndegen=ndegen)

    # Wannier-interpolated band energies on the fine grid (useful for lifetimes / band plots)
    _, E_fine, _ = Hwr_to_Hwk(Hwr, Rw, k_fine, ndegen=ndegen)

    return M_bk_fine, E_fine


def _as_int(X, name):
    """
    X as an int array, rounded to the nearest integer (np.rint, never truncation). Raises ValueError naming `name` if X is not
    integer to 1e-6 (e.g. a k-point passed instead of an R label).
    """
    X = np.asarray(X, float)
    if np.allclose(X, np.rint(X), atol=1e-6, rtol=0):
        return np.rint(X).astype(int)
    else:
        raise ValueError(f'integer rounding error in {name}')


def ws_images(R, R_center, MP, A_cols, tol=1e-6):
    """
    Wigner-Seitz images of the M_W labels (defined modulo the MP superlattice D s), centred on the defect: among
    R - R_center + D s, s in {-1, 0, 1}^2, keep those nearest in CARTESIAN distance, weight 1/n_tie each (Wannier90's ndegen).
    Inputs:
        R:        (nR, 3) ints, labels of Mwk_to_Mwr (raw) or of recenter_mwr (recentred).
        R_center: (3,) ints, R_d for raw labels, 0 for recentred ones.
        MP:       (3,) ints, MP grid (D1, D2, 1) of M, i.e. the period of the labels.
        A_cols:   (3, 3) floats, unit-cell vectors in columns.
        tol:      float, absolute tie tolerance, in the length unit of A_cols.
    Returns: dict with
        R_img: (n_img, 3) ints, kept images (absolute labels, = R[idx] modulo D).
        idx:   (n_img,) ints, row of R each image belongs to.
        w:     (n_img,) floats, 1/n_tie[idx].
        dist:  (nR,) floats, true distance to R_center, in the length unit of A_cols.
        n_tie: (nR,) ints, number of images kept per label.
    """

    A_cols = np.asarray(A_cols, float)

    # convert arrays to integer
    R = _as_int(R, 'R') # (nR, 3)
    R_center = _as_int(R_center, 'R_center') # (3,)
    D = _as_int(MP, 'MP') # (3,)

    # compute relative vector
    d = R - R_center # (nR, 3) int

    # bring d into [-D/2, D/2)
    d = (d + D//2 ) % D - D//2

    # construct table of the supercell's nine translations
    s = np.array([-1, 0, 1])
    T = np.array([(i,j,0) for i in s for j in s])
    T = D * T # (9, 3)

    # compute candidates
    c = d[:, None, :] + T[None, :, :] # (nR, 9, 3)
    c_cart = c @ A_cols.T

    L = np.linalg.norm(c_cart, axis=-1) # (nR, 9), lengths
    dist = np.min(L, axis=-1) # (nR, ), minimum distance

    mask = L <= dist[:, None] + tol # (nR, 9), bool
    n_tie = np.sum(mask, axis=1) # (nR,)

    idx, j = np.nonzero(mask) # 2*(n_img, )

    R_img = R_center[None, :] + c[idx, j, :] # (n_img, 3) int labels
    w = np.asarray(1 / n_tie[idx], float) # (n_img)

    return {'R_img': R_img, 'idx': idx, 'w': w, 'dist': dist, 'n_tie': n_tie}

def ws_phase(k, ws, nR, sign):
    """
    Phase matrix of the Wigner-Seitz images, Phi(k, r) = sum_{a : idx_a = r} w_a exp(sign 2 pi i k.R_img_a); on the MP grid
    it equals exp(sign 2 pi i k.R_r).
    Inputs:
        k:    (nk, 3) floats, k-points in reduced coordinates.
        ws:   dict returned by ws_images (uses R_img, idx, w).
        nR:   int, number of labels R (rows of the R passed to ws_images).
        sign: +1 (ket side) or -1 (bra side).
    Returns:
        Phi:  (nk, nR) complex.
    """

    k = np.atleast_2d(np.asarray(k, float)
                      )
    R_img = ws["R_img"] # (n_img, 3), ints
    dot = k @ R_img.T # (nk, n_img)
    phase = sign * 2j*np.pi * dot

    w = ws["w"] # (n_img)
    idx = ws['idx'] # (n_img)
    n_img = len(w)


    I = np.zeros((len(w), nR)) # (n_img, nR)
    I[np.arange(n_img), idx] = w 

    Phi = np.exp(phase) @ I # (nk, nR)

    return Phi