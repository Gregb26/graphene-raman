"""
local_G.py
    Python module used to compute the local part of the electron-defect scattering matrix in reciprocal space. 
    Quite a bit slow due to the number of FFTs to do ...
"""

import numpy as np

from electron_defect_interaction.utils.fft_utils import map_G_to_fft_grid

from tqdm import tqdm
from mpi4py import MPI

def prep_reciprocal_inputs(
    uc_wfk_path,
    sc_wfk_path,
    sc_p_pot_path,
    sc_d_pot_path,
    subtract_mean=True,
    pristine=False,
    bands=None,
    io=None,
):
    if io is None:
        from electron_defect_interaction.io import qe_io as io

    # get unit cell quantities
    C_nkg, nG = io.get_C_nk(uc_wfk_path)   # PW coefficients (ng,nk,nG_max) and number of active PWs per k (nk,)
    G_red = io.get_G_red(uc_wfk_path)      # reciprocal lattice vectors of unit cell in reduced coords (nk, nG_max, 3)
    k_red = io.get_k_red(uc_wfk_path)      # kpoints of unit cell in reduced coords (nk, 3)
    A_uc, _ = io.get_A_volume(uc_wfk_path) # primitive lattice vectors of unit cell A[:,i] = (a_i) (3,3)

    # Optionally restrict to a subset of bands
    if bands is not None:
        C_nkg = C_nkg[list(bands), ...]

    nb, nk, _ = C_nkg.shape

    # supercell geometry
    A_sc, Omega_sc = io.get_A_volume(sc_wfk_path)    # primtive lattice vectors of supercell (3,3) and supercell volume float

    # Local potential of pristine and defective supercell
    Vp, _ = io.get_pot(sc_p_pot_path, subtract_mean); Vp = Vp.transpose(2, 1, 0) # (Nx, Ny, Nz)
    Vd, _ = io.get_pot(sc_d_pot_path, subtract_mean); Vd = Vd.transpose(2, 1, 0) # (Nx, Ny, Nz)
    ngfft = Vp.shape  # FFT grid size

    if pristine:
        Ved = Vp
    else:
        Ved = Vd - Vp # defect potential

    #dDefect potential in reciprocal space
    Ved_G = np.fft.fftn(Ved) * (1 / np.prod(ngfft))          

    # supercell scaling factor, to fold unit cell quantities into supercell
    Ndiag = np.diag(np.rint(A_sc @ np.linalg.pinv(A_uc))).astype(int)  # (3,)

    # precompute mapping: reduced coords -> FFT index
    map_dict_x, map_dict_y, map_dict_z = map_G_to_fft_grid(ngfft)

    # cache reciprocal lattice vectors, FFT indises and PW coeffs per l
    G_sc_int = []   # list length nkpt, each (nG_k, 3) int
    idx_sc = []     # list length nkpt, each tuple(jx,jy,jz) int arrays length nG_k
    Ck_list = []    # list length nkpt, each (nband_sel, nG_k) complex

    for ik in range(nk):
        nG_k = int(nG[ik]) # active PWs for this k
        # fold unit cell Gs to supercell
        Gk_sc = G_red[ik, :nG_k, :] * Ndiag[np.newaxis, :]            # float but should be integer-ish
        Gk_sc_int = np.rint(Gk_sc).astype(int)                        # enforce integer triplets

        # FFT-grid indices for these G vectors
        jx = np.fromiter((map_dict_x[int(gx)] for gx in Gk_sc_int[:, 0]), dtype=np.int64, count=nG_k)
        jy = np.fromiter((map_dict_y[int(gy)] for gy in Gk_sc_int[:, 1]), dtype=np.int64, count=nG_k)
        jz = np.fromiter((map_dict_z[int(gz)] for gz in Gk_sc_int[:, 2]), dtype=np.int64, count=nG_k)

        G_sc_int.append(Gk_sc_int)
        idx_sc.append((jx, jy, jz))
        Ck_list.append(np.ascontiguousarray(C_nkg[:, ik, :nG_k]))      # (nb, nG_k)

    # return dict
    return {
        "Ved_G": Ved_G,                 # (Nx,Ny,Nz)
        "ngfft": ngfft,                 # (Nx,Ny,Nz)
        "Ndiag": Ndiag,                 # (3,)
        "Omega_sc": Omega_sc,
        "k_red": k_red,                 # (nkpt,3)
        "G_sc_int": G_sc_int,           # list of (nG_k,3) ints
        "idx_sc": idx_sc,               # list of (jx,jy,jz)
        "Ck": Ck_list,                  # list of (nband_sel,nG_k)
        "nb": nb,
        "nk": nk,
    }

########################################
# COMMENSURATE GRIDS: NO INTERPOLATION #
########################################

# NO MPI

def compute_ML_G(prep, block_size=512):

    # get inputs from prep function
    Ved_G   = prep["Ved_G"]
    ngfft   = prep["ngfft"]
    Ndiag   = prep["Ndiag"]
    k_red   = prep["k_red"]
    G_sc    = prep["G_sc_int"]
    Ck      = prep["Ck"]
    nb      = prep["nb"]
    nk      = prep["nk"]

    Nx, Ny, Nz = ngfft

    M = np.zeros((nb, nk, nb, nk), dtype=np.complex128)
    with tqdm(total=nk * nk, desc="(k',k) blocks") as pbar:
        for ikp in range(nk):
            Gp = G_sc[ikp]          # (nG_kp, 3)
            Ckp = Ck[ikp]           # (nband, nG_kp)
            nG_kp = Gp.shape[0]   

            for ik in range(nk):
                G = G_sc[ik]        # (nG_k, 3)
                Ck_ = Ck[ik]        # (nband, nG_k)

                # k'-k in supercell coordinates, round to integer and do direct indexing
                delta_k_sc = np.rint((k_red[ikp] - k_red[ik]) * Ndiag).astype(int)

                # compute M per block
                M_block = np.zeros((nb, nb), dtype=np.complex128)

                # split the sum over G' into managable blocks of size block_size
                for start in range(0, nG_kp, block_size):
                    stop = min(start + block_size, nG_kp)

                    # G' block
                    Gp_blk = Gp[start:stop]                    # (B,3)
                    Ckp_blk = Ckp[:, start:stop]               # (nband,B)

                    # compute q = k' - k + G' - G, wrapped into the FFT grid (periodic)
                    qx = (delta_k_sc[0] + Gp_blk[:, None, 0] - G[None, :, 0]) % Nx
                    qy = (delta_k_sc[1] + Gp_blk[:, None, 1] - G[None, :, 1]) % Ny
                    qz = (delta_k_sc[2] + Gp_blk[:, None, 2] - G[None, :, 2]) % Nz

                    # get V(q) for this block
                    V_block = Ved_G[qx, qy, qz]                   # (B, nG_k)

                    # contract: M += C_kp_block.conj() @ (V_block @ Ck0.T)
                    tmp = V_block @ Ck_.T
                    M_block += Ckp_blk.conj() @ tmp

                M[:, ikp, :, ik] = M_block
                pbar.update(1)

    return M

# MPI

def compute_block_M(Ved_G, ngfft, delta_k_sc, Gp, G, Ckp, Ck0, block_size):
    Nx, Ny, Nz = ngfft
    nband = Ckp.shape[0]
    nG_kp = Gp.shape[0]

    M_block = np.zeros((nband, nband), dtype=np.complex128)

    for start in range(0, nG_kp, block_size):
        stop = min(start + block_size, nG_kp)

        Gp_blk  = Gp[start:stop]              # (B,3)
        Ckp_blk = Ckp[:, start:stop]          # (nband,B)

        qx = (delta_k_sc[0] + Gp_blk[:, None, 0] - G[None, :, 0]) % Nx
        qy = (delta_k_sc[1] + Gp_blk[:, None, 1] - G[None, :, 1]) % Ny
        qz = (delta_k_sc[2] + Gp_blk[:, None, 2] - G[None, :, 2]) % Nz

        Vblk = Ved_G[qx, qy, qz]              # (B, nG_k)
        tmp  = Vblk @ Ck0.T                   # (B, nband)
        M_block += Ckp_blk.conj() @ tmp       # (nband, nband)

    return M_block

def compute_ML_G_mpi(prep, block_size=128, show_tqdm=True):
    """
    MPI over (k',k) blocks.
    Each rank computes a subset of pairs and rank 0 assembles the full M.
    """

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    Ved_G   = prep["Ved_G"]
    ngfft   = prep["ngfft"]
    Ndiag   = prep["Ndiag"]
    k_red   = prep["k_red"]
    G_sc    = prep["G_sc_int"]
    Ck      = prep["Ck"]
    nb   = prep["nb"]
    nk    = prep["nk"]

    # Build all (ikp, ik) pairs in a deterministic order
    pairs = [(ikp, ik) for ikp in range(nk) for ik in range(nk)]
    npairs = len(pairs)

    # Partition pairs across ranks (contiguous chunks)
    counts = [npairs // size + (1 if r < (npairs % size) else 0) for r in range(size)]
    displs = [sum(counts[:r]) for r in range(size)]
    local_pairs = pairs[displs[rank] : displs[rank] + counts[rank]]
    nlocal = len(local_pairs)

    # Local storage: blocks in the same order as local_pairs
    local_blocks = np.zeros((nlocal, nb, nb), dtype=np.complex128)

    it = range(nlocal)
    if show_tqdm and rank == 0:
        # only rank 0 shows progress for entire job
        it = tqdm(it, desc="(k',k) blocks (rank0 view)", total=nlocal, dynamic_ncols=True, mininterval=0.5)

    for i in it:
        ikp, ik = local_pairs[i]

        Gp  = G_sc[ikp]
        G   = G_sc[ik]
        Ckp = Ck[ikp]
        Ck0 = Ck[ik]

        delta_k_sc = np.rint((k_red[ikp] - k_red[ik]) * Ndiag).astype(int)

        local_blocks[i] = compute_block_M(
            Ved_G=Ved_G, ngfft=ngfft, delta_k_sc=delta_k_sc,
            Gp=Gp, G=G, Ckp=Ckp, Ck0=Ck0,
            block_size=block_size
        )
  
    # Gather all blocks to rank 0 (Gatherv on flattened complex128)
    sendbuf = local_blocks.reshape(-1)
    sendcount = sendbuf.size

    recvcounts = comm.gather(sendcount, root=0)

    if rank == 0:
        rdispls = np.zeros(size, dtype=np.int64)
        rdispls[1:] = np.cumsum(recvcounts[:-1], dtype=np.int64)
        recvbuf = np.empty(sum(recvcounts), dtype=np.complex128)
    else:
        rdispls = None
        recvbuf = None

    comm.Gatherv(
        sendbuf,
        [recvbuf, recvcounts, rdispls, MPI.DOUBLE_COMPLEX],
        root=0
    )

    if rank != 0:
        return None

    # Assemble full M on rank 0
    M = np.zeros((nb, nk, nb, nk), dtype=np.complex128)

    # Reconstruct the global ordering of blocks from counts/displs
    # recvbuf contains blocks for rank 0 chunk, then rank 1 chunk, etc.
    offset = 0
    for r in range(size):
        nblocks_r = counts[r]
        nvals_r = nblocks_r * nb * nb
        blocks_r = recvbuf[offset : offset + nvals_r].reshape(nblocks_r, nb, nb)
        offset += nvals_r

        pairs_r = pairs[displs[r] : displs[r] + counts[r]]
        for j, (ikp, ik) in enumerate(pairs_r):
            M[:, ikp, :, ik] = blocks_r[j]

    return M


###############################################################################
# DENSE GRIDS BY EXACT ZERO-PADDING OF THE DEFECT POTENTIAL                    #
#                                                                             #
# Densifying the primitive k-grid by an integer factor p means querying the   #
# defect potential on a p-times FINER reciprocal lattice (spacing b_sc/p).    #
# Because the defect potential has compact support inside the supercell, its  #
# continuous Fourier transform is fully determined by the existing samples    #
# (Whittaker-Shannon); the finer samples are obtained by sinc interpolation   #
# in reciprocal space == zero-padding in REAL space. Padding a real array     #
# keeps the spectrum Hermitian, so the padded potential is real with no       #
# manual Nyquist bookkeeping. We normalise so that the coincident coefficients#
# are preserved exactly (Ved_G_pad[::p, ::p, :] == Ved_G), which makes the    #
# dense matrix reduce to the original one at the shared k-points to machine    #
# precision and lets us reuse compute_block_M verbatim.                       #
###############################################################################


def zero_pad_potential(Ved, p):
    """
    Zero-pad a real defect potential by an integer factor p in the in-plane directions (x, y; z is
    left unchanged for 2D systems such as graphene) and return its Fourier transform on the resulting
    (p*Nx, p*Ny, Nz) reciprocal grid.

    The padding is done in REAL space (embedding the supercell potential, which has compact support,
    into a p-times larger box) which is the exact sinc interpolation of the reciprocal samples onto
    the p-times finer reciprocal lattice. The normalisation is chosen so that the coefficients at the
    coincident reciprocal points are preserved exactly:

        Ved_G_pad[::p, ::p, :] == np.fft.fftn(Ved) / Ved.size            (machine precision)

    and, since the padded real array is real, ifftn(Ved_G_pad) is real to machine precision with

        ifftn(Ved_G_pad).real[:Nx, :Ny, :] / p**2 == Ved                 (the padding region is 0).

    Inputs:
        Ved: (Nx, Ny, Nz) real array
            Defect potential V_d - V_p in real space (same layout as prep["Ved_G"] before its FFT).
        p: int
            In-plane densification factor (>= 1).
    Returns:
        Ved_G_pad: (p*Nx, p*Ny, Nz) complex array
            Defect potential on the finer reciprocal grid, same convention as prep["Ved_G"].
    """
    p = int(p)
    if p < 1:
        raise ValueError(f"densification factor p must be >= 1, got {p}")
    Ved = np.asarray(Ved)
    Nx, Ny, Nz = Ved.shape

    if p == 1:
        return np.fft.fftn(Ved) / Ved.size

    # real-space zero-pad in the plane; support stays in the corner so that frequency p*m maps to m
    Ved_pad = np.zeros((p * Nx, p * Ny, Nz), dtype=np.result_type(Ved.dtype, np.float64))
    Ved_pad[:Nx, :Ny, :] = Ved.real

    # normalise by the ORIGINAL grid size -> coincident coefficients preserved (drop-in for Ved_G)
    Ved_G_pad = np.fft.fftn(Ved_pad) / (Nx * Ny * Nz)
    return Ved_G_pad


def _reconstruct_Ved(prep):
    """Real-space defect potential from prep["Ved_G"] (Ved_G = fftn(Ved)/prod -> Ved = ifftn*prod)."""
    return (np.fft.ifftn(prep["Ved_G"]) * np.prod(prep["ngfft"])).real


def _prep_dense_caches(prep, p, C_dense, G_dense, nG_dense):
    """
    Build the per-k caches needed by the dense convolution: the padded potential on the finer
    reciprocal grid, the effective supercell folding/grid, and (G_sc, Ck) lists folded with Ndiag*p.

    C_dense, G_dense, nG_dense are the raw QE quantities of a PRIMITIVE calculation on the dense grid
    (same ecut as the coarse one): C_dense (nb, nk, nG_max), G_dense (nk, nG_max, 3) reduced ints,
    nG_dense (nk,) active plane waves per k. Returns a dict mirroring the fields compute_ML_G uses.
    """
    Ved_G_pad = zero_pad_potential(_reconstruct_Ved(prep), p)
    ngfft_eff = Ved_G_pad.shape                                  # (p*Nx, p*Ny, Nz)
    Ndiag_eff = prep["Ndiag"] * np.array([p, p, 1], dtype=int)   # finer supercell folding

    nb = C_dense.shape[0]
    nk = len(nG_dense)

    G_sc, Ck = [], []
    for ik in range(nk):
        nGk = int(nG_dense[ik])
        # fold primitive G onto the effective (p-times finer) supercell reciprocal lattice
        Gk_sc = np.rint(G_dense[ik, :nGk, :] * Ndiag_eff[np.newaxis, :]).astype(int)
        G_sc.append(Gk_sc)
        Ck.append(np.ascontiguousarray(C_dense[:, ik, :nGk]))

    return {
        "Ved_G": Ved_G_pad, "ngfft": ngfft_eff, "Ndiag": Ndiag_eff,
        "G_sc_int": G_sc, "Ck": Ck, "nb": nb, "nk": nk,
    }


# NO MPI

def compute_ML_G_dense(prep, p, k_dense, C_dense, G_dense, nG_dense, block_size=512):
    """
    Local part M^L on a primitive k-grid p-times denser than the supercell grid, by exact zero-padding
    of the defect potential (no supercell recomputation, no interpolation).

    The dense grid must be the commensurate (p*N1) x (p*N2) x N3 grid so that delta_k_sc =
    rint((k'-k) * Ndiag*p) is integer and indexes Ved_G_pad exactly.

    Inputs:
        prep: dict from prep_reciprocal_inputs (provides Ved_G, ngfft, Ndiag of the coarse problem).
        p: int, in-plane densification factor.
        k_dense: (nk_d, 3) reduced coords of the dense grid (order must match C_dense/G_dense, i.e.
            qe_io.get_k_red(dense_uc.save)).
        C_dense, G_dense, nG_dense: raw QE quantities of the primitive calculation on the dense grid
            (same ecut). If prep was built with a band subset, restrict C_dense to the same bands.
        block_size: int, G'-blocking for the contraction (as in compute_ML_G).
    Returns:
        M: (nb, nk_d, nb, nk_d) complex, M^L[bra_band, k', ket_band, k] on the dense grid.
    """
    d = _prep_dense_caches(prep, p, C_dense, G_dense, nG_dense)
    Ved_G, ngfft, Ndiag = d["Ved_G"], d["ngfft"], d["Ndiag"]
    G_sc, Ck, nb, nk = d["G_sc_int"], d["Ck"], d["nb"], d["nk"]
    k_dense = np.asarray(k_dense)

    M = np.zeros((nb, nk, nb, nk), dtype=np.complex128)
    with tqdm(total=nk * nk, desc="(k',k) dense blocks") as pbar:
        for ikp in range(nk):
            for ik in range(nk):
                delta_k_sc = np.rint((k_dense[ikp] - k_dense[ik]) * Ndiag).astype(int)
                M[:, ikp, :, ik] = compute_block_M(
                    Ved_G=Ved_G, ngfft=ngfft, delta_k_sc=delta_k_sc,
                    Gp=G_sc[ikp], G=G_sc[ik], Ckp=Ck[ikp], Ck0=Ck[ik],
                    block_size=block_size,
                )
                pbar.update(1)
    return M


# MPI

def compute_ML_G_dense_mpi(prep, p, k_dense, C_dense, G_dense, nG_dense, block_size=128, show_tqdm=True):
    """
    MPI version of compute_ML_G_dense: same dense zero-padding, with the (k',k) pairs distributed over
    ranks and gathered on rank 0 (mirrors compute_ML_G_mpi). Returns M on rank 0, None elsewhere.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    d = _prep_dense_caches(prep, p, C_dense, G_dense, nG_dense)
    Ved_G, ngfft, Ndiag = d["Ved_G"], d["ngfft"], d["Ndiag"]
    G_sc, Ck, nb, nk = d["G_sc_int"], d["Ck"], d["nb"], d["nk"]
    k_dense = np.asarray(k_dense)

    pairs = [(ikp, ik) for ikp in range(nk) for ik in range(nk)]
    npairs = len(pairs)

    counts = [npairs // size + (1 if r < (npairs % size) else 0) for r in range(size)]
    displs = [sum(counts[:r]) for r in range(size)]
    local_pairs = pairs[displs[rank] : displs[rank] + counts[rank]]
    nlocal = len(local_pairs)

    local_blocks = np.zeros((nlocal, nb, nb), dtype=np.complex128)

    it = range(nlocal)
    if show_tqdm and rank == 0:
        it = tqdm(it, desc="(k',k) dense blocks (rank0 view)", total=nlocal, dynamic_ncols=True, mininterval=0.5)

    for i in it:
        ikp, ik = local_pairs[i]
        delta_k_sc = np.rint((k_dense[ikp] - k_dense[ik]) * Ndiag).astype(int)
        local_blocks[i] = compute_block_M(
            Ved_G=Ved_G, ngfft=ngfft, delta_k_sc=delta_k_sc,
            Gp=G_sc[ikp], G=G_sc[ik], Ckp=Ck[ikp], Ck0=Ck[ik],
            block_size=block_size,
        )

    sendbuf = local_blocks.reshape(-1)
    recvcounts = comm.gather(sendbuf.size, root=0)

    if rank == 0:
        rdispls = np.zeros(size, dtype=np.int64)
        rdispls[1:] = np.cumsum(recvcounts[:-1], dtype=np.int64)
        recvbuf = np.empty(sum(recvcounts), dtype=np.complex128)
    else:
        rdispls = None
        recvbuf = None

    comm.Gatherv(sendbuf, [recvbuf, recvcounts, rdispls, MPI.DOUBLE_COMPLEX], root=0)

    if rank != 0:
        return None

    M = np.zeros((nb, nk, nb, nk), dtype=np.complex128)
    offset = 0
    for r in range(size):
        nvals_r = counts[r] * nb * nb
        blocks_r = recvbuf[offset : offset + nvals_r].reshape(counts[r], nb, nb)
        offset += nvals_r
        for j, (ikp, ik) in enumerate(pairs[displs[r] : displs[r] + counts[r]]):
            M[:, ikp, :, ik] = blocks_r[j]
    return M


def compute_M_dense(prep, p, uc_dense_path, sc_p_wfk_path, sc_d_wfk_path, pseudo_path,
                    k_dense, C_dense, G_dense, nG_dense, io=None, pseudo_reader=None,
                    bands=None, block_size=512, use_mpi=False):
    """
    Full dense electron-defect matrix M = M^L_dense + M^NL_dense on the p-times denser grid.

    M^L is the zero-padded reciprocal evaluation (compute_ML_G_dense[_mpi]); M^NL is analytic in k and
    needs no densification machinery -- compute_M_NL is simply called on the dense primitive .save
    (uc_dense_path), which already carries the dense k-grid. The two share the [bra_band, k', ket_band,
    k] convention; pass the SAME `bands` to both and make sure C_dense/k_dense are in the order of
    qe_io.get_k_red(uc_dense_path) so the band/k axes line up before summing.

    Returns M[bra_band, k', ket_band, k] on the dense grid (Hartree).
    """
    from electron_defect_interaction.defects.non_local import compute_M_NL

    if use_mpi:
        M_L = compute_ML_G_dense_mpi(prep, p, k_dense, C_dense, G_dense, nG_dense, block_size=block_size)
    else:
        M_L = compute_ML_G_dense(prep, p, k_dense, C_dense, G_dense, nG_dense, block_size=block_size)

    M_NL = compute_M_NL(uc_dense_path, sc_p_wfk_path, sc_d_wfk_path, pseudo_path,
                        io=io, pseudo_reader=pseudo_reader, bands=bands)

    if M_L is None:          # non-root rank in the MPI path
        return None
    return M_L + M_NL
