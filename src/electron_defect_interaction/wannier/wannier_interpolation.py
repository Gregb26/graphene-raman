"""
wannier_interpolation.py
    Python module containing functions to interpolate objects computed on a coarse kpoint grid onto a fine kpoint grid
    via Maximally Localized Wannier Functions. 
"""

import numpy as np
from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.io.wannier_io import read_w90_mat, read_w90_HR

def Mbk_to_Mwk(Mbk, U, U_dis=None):
    """
    Rotates an object Mbk in Bloch gauge into Wannier gauge using the Wannier gauge matrix U (or V=U_dis @ U in the case of entangled bands)
    Inputs:
        Mbk: (nb, nk, nb, nk) array of complex
            Matrix in Bloch gauge to rotate. This can be e.g. a Hamiltonian or a scattering matrix. nk is the number of kpoints, nb is the number
            of Bloch bands and nw is the number of Wannier functions.
        U: (nk, nb, nw) or (nk, nw, nw) array of complex
            Wannier gauge matrix. First shape if no entangled bands (U_dis is None). Second shape if entangled bands (U_dis is not None).
        U_dis: (nk, nb, nw) array of complex
            Disentanglement matrix if the number of Wannier functions is less than the number of Bloch bands.
    Returns:
        Mwk: (nw, nk, nw, nk)
            Object in the Wannier gauge: Mwk = V^dag Mbk V
    """

    nk, nw, _ = U.shape
    # Entangled case, rotation matrix is V = U_dis @ U
    if U_dis is not None:
        nb = U_dis.shape[1]
        V = np.zeros((nk, nb, nw), dtype=complex)
        for ik in range(nk):

            V[ik] = U_dis[ik] @ U[ik] # (nk, nb, nw)

            # testing
            with np.errstate(all='ignore'):
                P = V[ik] @ V[ik].conj().T
                assert np.allclose(P, P.conj().T, atol=1e-10)     
                assert np.allclose(P @ P, P, atol=1e-8)           
                assert np.allclose(np.trace(P), V.shape[-1], 1e-10)

    else:
        V = U # (nk, nb, nw) with nw = nb
    
    # Rotate Mb from Bloch gauge to Wannier gauge Mw
    Mwk = np.zeros((nw, nk, nw, nk), dtype=complex)
    for ik in range(nk):
        Vk_h = V[ik].conj().T # Hermitian conjugate
        for ikp in range(nk):

            Mwk[:, ik, :, ikp] = Vk_h @ Mbk[:, ik, :, ikp] @ V[ikp] # (nw, nk, nw, nk)
    
    return Mwk

def Mwk_to_Mwr(Mwk, k_red, MP_grid):
    """
    Transforms the object Mwk from reciprocal space to real space using a double Fourier transform.
    Inputs:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space to transform to real space, nk is the number of kpoints and nw the number of Wannier functions. This is e.g.
            a scattering matrix in reciprocal space in the Wannier gauge.
        MP_grid: tuple of ints
            Monkhorst-Pack grid used to define the kpoint grid. Assuming no symmetry reduction has been applied! This is important! 
    Returns
        Mwr: (nw, nr, nw, nr) array of complex
            Double Fourier transform of Mwk; nr is the number of R vectors, inferred from the Monkhorst-Pack grid.
        R:  (nr, 3) array of ints
            R vectors in real space used to compute the double FT. Is dual to the kpoint grid.
    """

    nk = Mwk.shape[1]

    # Build the grid in real space based on the kpoint grid
    N1, N2, N3 = MP_grid
    r1 = np.arange(N1); r2 = np.arange(N2); r3 = np.arange(N3)
    rr1, rr2, rr3 = np.meshgrid(r1, r2, r3)
    R = np.stack((rr1, rr2, rr3), axis=-1) # (N1, N2, N3, 3)
    R = R.reshape(N1*N2*N3, 3) # (nr, 3)

    # Compute phase
    phase_kp = np.exp(2j*np.pi * (k_red @ R.T)) # (nk, nr)
    phase_k = phase_kp.conj() # (nk, nr)

    # Double sum over k and k'
    Mwr = np.einsum('kr, nkNK, KR -> nrNR', phase_kp, Mwk, phase_k, optimize=True)  / nk **2 # (nw, nr, nk, nr)

    return Mwr, R

def Mwr_to_Mwk(Mwr, R, k):
    """
    Transforms the object in Mwr in real space in Wannier gauge to Mwk in reciprocal space in Wannier gauge. Here, the kpoint grid k can be any arbitrarily dense grid hence this function can be used to interpolate Mwk on a finer grid.
    Inputs:
        Mwr: (nw, nr, nw, nr) array of complex
            Object in real space in Wannier gauge to transform to reciprocal space in Wannier gauge. nr is the number of R vectors in real space and nw is the number of Wannier functions.
        R: (nr, 3) array of ints
            Real space lattice vectors in reduced coordinates for which Mwr is defined.
        k: (nk, 3) array of floats
            k vectors in reduced coordinates for which to transform Mwr.
    Returns:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space in Wannier gauge.
    """
    
    # precompute phases
    phase_kp = np.exp(-2j*np.pi * (k @ R.T)) # (nk, nr)
    phase_k = phase_kp.conj() # (nk, nr)

    # sum over R 
    Mwk = np.einsum("kr, wrWR, KR -> wkWK", phase_kp, Mwr, phase_k, optimize=True)

    return Mwk

def Mwk_to_Mbk(Mwk, Hwr, Rw, k, ndegen=None):
    """
    Transforms the object Mwk in reciprocal space in Wannier gauge to Bloch gauge using the Wannier Hamiltonian.
    Inputs:
        Mwk: (nw, nk, nw, nk) array of complex
            Object in reciprocal space in Wannier gauge to transform to Bloch gauge. nw is the number of Wannier functions and nk is the number of kpoints.
        Hwr: (nrpts, nw, nw) array of complex
            Hamtilontian in Wannier (real space) basis. This Hamiltonian is computed on the kpoint grid k via a Fourier transformed and diagonalized to obtain the unitary gauge rotation matrix U. nrw is the number of R vectors used internally by Wannier90 to compute the Hamiltonian in real space.
        Rw: (nrpts, 3) array of ints:
            Real space lattice vectors used internally by Wannier90.
        k: (nk, 3) array of floats
            kpoint grid on which to transform Hwr. Can arbitrarily dense and must match the one on which Mwk is defined.
        ndegen: (nrpts,) array of ints, optional
            Wigner-Seitz degeneracies passed through to Hwr_to_Hwk (H(k)=sum_R e^{ik.R} H(R)/ndegen(R)).
    Returns:
        Mbk: (nw, nk, nw, nk) array of complex
            Object transformed to the smooth Bloch gauge (eigenbasis of the interpolated Hamiltonian).
    """

    # compute rotation matrix from Wannier Hamiltonian
    _, _, Uwk = Hwr_to_Hwk(Hwr, Rw, k, ndegen=ndegen) # (nk, nw, nw)

    # Rotate Mwk to Bloch gauge: per-k unitary similarity Mbk[:,k,:,K] = Uwk[k]^dag Mwk[:,k,:,K] Uwk[K]
    Mbk = np.einsum('kwb, wkWK, KWB -> bkBK', Uwk.conj(), Mwk, Uwk, optimize=True)

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
    Hwr, Rw, ndegen = read_w90_HR(wannier_tb)

    # 3. coarse Bloch -> Wannier gauge -> real space -> fine grid -> smooth Bloch gauge
    MP_grid = _infer_mp_grid(k_coarse)
    Mwk = Mbk_to_Mwk(M, U, U_dis)                                  # (nw, nk, nw, nk)
    Mwr, R = Mwk_to_Mwr(Mwk, k_coarse, MP_grid)                    # (nw, nr, nw, nr)
    Mwk_fine = Mwr_to_Mwk(Mwr, R, k_fine)                          # (nw, nkf, nw, nkf)
    M_bk_fine = Mwk_to_Mbk(Mwk_fine, Hwr, Rw, k_fine, ndegen=ndegen)

    # Wannier-interpolated band energies on the fine grid (useful for lifetimes / band plots)
    _, E_fine, _ = Hwr_to_Hwk(Hwr, Rw, k_fine, ndegen=ndegen)

    return M_bk_fine, E_fine

