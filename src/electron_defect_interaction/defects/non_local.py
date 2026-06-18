"""
nonlocal.py
    Python modules containing the functions to compute the non-local part of the electron-defect scattering matrix.

TODO: 
    - Generalize to supercells with different atoms (right now assumes one atom type)
    - Generalize to atoms with different l quantum number (right now assums a single l since one atom type)
    - Generalize to atoms with different number of projector channels per l (now assums 2 for l=0,1)
"""

import numpy as np

from electron_defect_interaction.utils.lattice import red_to_cart
from electron_defect_interaction.utils.planewaves import mask_invalid_G
from electron_defect_interaction.io.pseudo_io import read_psp8, read_upf, fq_from_fr


# NO MPI

def build_K_vectors(k_red, G_red, keep, B_uc):
    """
    Compute the K=k+G vectors, their norm and their unit vectors.

    Inputs:
        k_red: (nkpt, 3) array of ints:
            kpoints in reduced coordinates.
        G_red: (nkpt, nG_max, 3) array of ints:
            Reciprocal lattice vectors in reduced coordinates
        mask: (nkpt, nG_max) array of boolean:
            Boolean mask that selects only active G vectors for a given k 
            i.e. those that satisfy 0.5*(2pi)^2(k+G)^2 < ecut
        B_uc: (3, 3) array of floats:
            Columns are the primitive reciprocal lattice vectors B[:,i] = B_i
        
    Returns
        K: (nkpt, nG_max, 3) array of floats:
            K=k+G vectors in cartesian coords
        K_norm: (nkpt, nG_max) array of floats:
            Norm of K=k+G vectors
        K_hat: (nkpt, nG_max, 3):
            Unit vectors in the direction of the K=k+G vectors.
    """

    # Compute K=k+G in cartesian coordinates
    K_red = k_red[:, np.newaxis, :] + G_red
    Ks = red_to_cart(K_red, B_uc)
    K = np.where(keep[..., np.newaxis], Ks, 0.0)

    # Compute norm of K=k+G vectors
    norms = np.linalg.norm(K, axis=2) # (nkpt, nG_max)
    K_norm = np.where(keep, norms, 0.0) # (nkpt, nG_max)

    # Compute unit vectors 
    valid_Ks = keep & (K_norm > 0) # selects Ks with valig G and non-zero norm (to avoid division by zero)
    K_hat = np.zeros_like(K) # (nkpt, nG_max, 3)
    K_hat = np.divide(K, K_norm[:, :, None], out=K_hat, where=valid_Ks[:,:,None]) # (nkpt, nG_max, 3)

    return K, K_norm, K_hat

def compute_phase(K, tau_as):
    """
    Compute phase exp(i tau_{a}\\cdot K)

    Inputs:
        K: (nkpt, nG_max, 3):
            K=k+G vectors in cartesian coordinates
        tau_as: (natom, 3):
            Atomic positions of atom s of type a in cartesian coordinates

    Returns:
        phase_ksg: (nkpt, natom, nG_max)      
    """

    dot = np.einsum("kgd, sd -> ksg", K, tau_as, optimize = True)
    phase_ksg = np.exp(-1j * dot)

    return phase_ksg

def compute_angular_part(K_hat, lmax):
    """
    Compute the spherical harmonics Y_l^m(K_hat)

    Inputs:
        K_hat: (nkpt, nG_max) array of floats:
            Norm of the K=k+G vectors
        lmax: float
            Maximum orbital angular momentum quantum number
        
    Returns:
        Y_kglm (nkgpt, nG_max, lmax, 2lmax+1) array of complex
            Spherical harmonics at l, m, phi and theta, where phi and theta are the spherical angles of the unit vectors K_hat
    """

    nkpt, nG_max, _ = K_hat.shape
    # compute angles in K
    theta = np.arccos(np.clip(K_hat[...,2], -1.0, 1.0)) # theta = arccos(Kz) (0, pi)
    phi = np.arctan2(K_hat[..., 1], K_hat[..., 0]) + np.pi # phi = arctan(Ky/Kx) (0, 2pi)

    Y_kglm = np.zeros((nkpt, nG_max, lmax+1, 2*lmax+1), dtype=complex)

    from scipy.special import sph_harm_y
    for l in range(lmax+1):
        for m in range(-l, l+1):
            mi = m + lmax
            Y_kglm[..., l, mi] = sph_harm_y(l, m, theta, phi) # (nkpt, nG_max, lmax, 2lmax+1)

    return Y_kglm

def compute_M_NL(uc_wfk_path, sc_p_wfk_path, sc_d_wfk_path, psp8_path, io=None, pseudo_reader=None, bands=None):
    """
    Computes the non local part of the electron-defect interaction matrix.

    Inputs:
        uc_wfk_path: str
            Path to the unit-cell wavefunctions (ABINIT WFK.nc / QE prefix.save dir).
        sc_p_wfk_path: str
            Path to the pristine supercell (its atomic positions and geometry are used).
        sc_d_wfk_path: str
            Path to the defective supercell (its atomic positions are used).
        psp8_path: str
            Path to the pseudopotential (ABINIT .psp8 / QE .upf).
        io: module, optional
            I/O backend (defaults to qe_io for Quantum ESPRESSO; pass abinit_io for ABINIT).
        pseudo_reader: callable, optional
            Pseudopotential reader returning (ekb_li, fr_li, rgrid, lmax, imax, V_L).
            Defaults to read_upf (QE); pass read_psp8 for ABINIT.
    Returns:
        M_NL: (nband, nkpt, nband, nkpt) array of complex:
            Non-local part of the electron-defect interaction matrix
    """
    from scipy.interpolate import CubicSpline

    if io is None:
        from electron_defect_interaction.io import qe_io as io
    if pseudo_reader is None:
        pseudo_reader = read_upf

    # Get non local part of the pseudopotentials (ekb energies, r*beta projectors, radial grid, ...)
    ekb_li, fr_li, rgrid, lmax, imax, _ = pseudo_reader(psp8_path)

    # Get necessary unit cells quantities
    C_nkg, nG = io.get_C_nk(uc_wfk_path) # planewave coefficients (nband, nkpt, nG_max) and number of active G per k (nkpt)
    if bands is not None:
        C_nkg = C_nkg[list(bands), ...]  # restrict to the requested bands (the M sub-block is independent)
    G_red = io.get_G_red(uc_wfk_path) # reciprocal lattice vectors in reduced coords of the unit cell (nkpt, nG_max, 3)
    k_red = io.get_k_red(uc_wfk_path) # kpoints in reduced coords of the unit cell (nkpt, 3)
    B_uc, _ = io.get_B_volume(uc_wfk_path) # primitive reciprocal lattice vectors of the unit cell B[:,i] = b_i
    A_uc, Omega_uc = io.get_A_volume(uc_wfk_path) # unit-cell volume sets the M^NL prefactor (4pi)^2/Omega_uc
    ecut = io.get_ecut(uc_wfk_path) # planewave energy cutoff used in the calculation

    assert np.isclose(ecut, io.get_ecut(sc_p_wfk_path)), 'Must use the same ecut in both unit cell and supercell calculations!'

    # Get necessary super cell quantities
    A_sc, Omega_sc = io.get_A_volume(sc_p_wfk_path) # primitive lattice vectors and cell volume of the supercell

    x_red_p = io.get_x_red(sc_p_wfk_path) # atomic positions in reduced coords of the supercell (natom, 3)
    tau_s_p = red_to_cart(x_red_p, A_sc) # atomoic positions of atoms in the supercell in cartesian coords

    x_red_d = io.get_x_red(sc_d_wfk_path) # atomic positions in reduced coords of the supercell (natom, 3)
    tau_s_d = red_to_cart(x_red_d, A_sc) # atomoic positions of atoms in the supercell in cartesian coords

    # Compute boolean mask that selects only the active recripocal lattice vector G for each k-point and mask invalid C's (pad with zeros)
    keep = mask_invalid_G(nG)
    C_nkg = np.where(keep, C_nkg, 0.0) # (nband, nkpt, nG_max)

    # Compute K=k+G vectors
    K, K_norm, K_hat = build_K_vectors(k_red, G_red, keep, B_uc) # (nkpt, nG_max, 3), (nkpt. nG_max), (nkpt, nG_max, 3)

    # Transform the radial from factor to q space
    qmax = 2*np.sqrt(2*ecut) # choose a qmax two times the recommended for safety ...
    q = np.linspace(0, qmax, 2000)
    fq_li = fq_from_fr(rgrid, fr_li, q)

    # Interpolate radial form factors to be able to evaluate them at K=|k+G| vectors
    Fq_li = CubicSpline(q, fq_li, axis=-1, extrapolate=False)   
    assert np.max(K_norm) < qmax, 'Maximum |K|=|k+G| must be within the qgrid'
    F_likg = Fq_li(K_norm) # (lmax+1, imax+1, nkpt, nG_max)

    # Compute phases
    phase_ksg_p = compute_phase(K, tau_s_p) # (nkpt, natom, nG_max)
    phase_ksg_d = compute_phase(K, tau_s_d)
    # Compute angular part
    Y_kglm = compute_angular_part(K_hat, lmax) # (nkpt, nG_max, lmax+1, 2lmax+1)

    # Compute overlaps by summing over all G vectors.
    # Prefactor uses the UNIT-CELL volume: M^NL = (4pi)^2/Omega_uc * sum ... (equation convention).
    pref = (4*np.pi / np.sqrt(Omega_uc))
    B_p = pref * np.einsum("nkg,likg,kglm,ksg->nkslim",
                        np.conj(C_nkg), F_likg, Y_kglm, phase_ksg_p, optimize=True)
    B_d = pref * np.einsum("nkg,likg,kglm,ksg->nkslim",
                        np.conj(C_nkg), F_likg, Y_kglm, phase_ksg_d, optimize=True)

    M_p = np.einsum("li,nkslim,jpslim->nkjp", ekb_li, B_p, np.conj(B_p), optimize=True)
    M_d = np.einsum("li,nkslim,jpslim->nkjp", ekb_li, B_d, np.conj(B_d), optimize=True)

    M_NL = M_d - M_p

    return M_NL 

# MPI

from mpi4py import MPI

def split_counts(n, size):
    """
    For n total objects and size MPI ranks, distribute the objects across the ranks. counts[r] is how many objects rank r should process
    and displs[r] is the starting index of rank r inside the global array
    """
    counts = np.array([n // size + (1 if r < (n % size) else 0) for r in range(size)], dtype=np.int64)
    displs = np.zeros(size, dtype=np.int64)
    displs[1:] = np.cumsum(counts[:-1])

    return counts, displs

def local_slice(n, comm):
    """
    Slice for current MPI rank
    """
    size = comm.Get_size()
    rank = comm.Get_rank()
    counts, displs = split_counts(n, size)
    start = displs[rank]
    stop = start + counts[rank]

    return slice(start, stop), counts, displs



def compute_M_NL_mpi(uc_wfk_path, sc_p_wfk_path, sc_d_wfk_path, psp8_path,
                     io=None, pseudo_reader=None, bands=None):
    """
    MPI version of compute_M_NL. The separable factor B_nk(atom,l,i,m_l) is built on every rank
    (it is small, especially with `bands` restricted), and the O(nk^2) contraction that forms
    M_mn(k',k) = (4pi)^2/Omega_uc sum E_li B_mk'* B_nk is distributed over the bra k-index k'.
    Each rank fills its k' slices of M and the ranks Allreduce the full result.

    Returns M[bra_band, k', ket_band, k] (Hartree) on every rank.
    """
    from scipy.interpolate import CubicSpline

    if io is None:
        from electron_defect_interaction.io import qe_io as io
    if pseudo_reader is None:
        pseudo_reader = read_upf

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank(); size = comm.Get_size()

    ekb_li, fr_li, rgrid, lmax, imax, _ = pseudo_reader(psp8_path)

    C_nkg, nG = io.get_C_nk(uc_wfk_path)
    if bands is not None:
        C_nkg = C_nkg[list(bands), ...]
    G_red = io.get_G_red(uc_wfk_path)
    k_red = io.get_k_red(uc_wfk_path)
    B_uc, _ = io.get_B_volume(uc_wfk_path)
    _, Omega_uc = io.get_A_volume(uc_wfk_path)
    ecut = io.get_ecut(uc_wfk_path)

    A_sc, _ = io.get_A_volume(sc_p_wfk_path)
    tau_p = red_to_cart(io.get_x_red(sc_p_wfk_path), A_sc)
    tau_d = red_to_cart(io.get_x_red(sc_d_wfk_path), A_sc)

    nb, nk, _ = C_nkg.shape
    keep = mask_invalid_G(nG)
    C = np.where(keep, C_nkg, 0.0)

    K, K_norm, K_hat = build_K_vectors(k_red, G_red, keep, B_uc)
    q = np.linspace(0, 2 * np.sqrt(2 * ecut), 2000)
    Fq = CubicSpline(q, fq_from_fr(rgrid, fr_li, q), axis=-1, extrapolate=False)
    F_likg = Fq(K_norm)                       # (l, i, nk, nG)
    Y_kglm = compute_angular_part(K_hat, lmax)  # (nk, nG, l, m)
    pref = 4 * np.pi / np.sqrt(Omega_uc)

    def build_B(tau):
        # B[n,k,s,l,i,m] built one k at a time to avoid the (nk,natom,nG) phase array
        natom = tau.shape[0]
        B = np.zeros((nb, nk, natom, lmax + 1, imax, 2 * lmax + 1), dtype=np.complex128)
        for ik in range(nk):
            phase = np.exp(-1j * (K[ik] @ tau.T)).T   # (natom, nG)
            B[:, ik] = pref * np.einsum("ng,lig,glm,sg->nslim",
                                        np.conj(C[:, ik]), F_likg[:, :, ik], Y_kglm[ik], phase,
                                        optimize=True)
        return B

    B_p = build_B(tau_p)
    B_d = build_B(tau_d)

    # Distribute the bra k-index (k') over ranks
    counts = [nk // size + (1 if r < (nk % size) else 0) for r in range(size)]
    displs = [sum(counts[:r]) for r in range(size)]
    my_kp = range(displs[rank], displs[rank] + counts[rank])

    M_local = np.zeros((nb, nk, nb, nk), dtype=np.complex128)
    for ikp in my_kp:
        Mp = np.einsum("li,nslia,jpslia->njp", ekb_li, B_p[:, ikp], np.conj(B_p), optimize=True)
        Md = np.einsum("li,nslia,jpslia->njp", ekb_li, B_d[:, ikp], np.conj(B_d), optimize=True)
        M_local[:, ikp, :, :] = Md - Mp

    M = np.zeros((nb, nk, nb, nk), dtype=np.complex128)
    comm.Allreduce(M_local, M, op=MPI.SUM)
    return M
