"""
local_G.py
    Python module used to compute the local part of the electron-defect scattering matrix in reciprocal space. 
    Quite a bit slow due to the number of FFTs to do ...
"""

from electron_defect_interaction.io.abinit_io import *
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
):

    # get unit cell quantities
    C_nkg, nG = get_C_nk(uc_wfk_path)   # PW coefficients (ng,nk,nG_max) and number of active PWs per k (nk,)
    G_red = get_G_red(uc_wfk_path)      # reciprocal lattice vectors of unit cell in reduced coords (nk, nG_max, 3)
    k_red = get_k_red(uc_wfk_path)      # kpoints of unit cell in reduced coords (nk, 3)
    A_uc, _ = get_A_volume(uc_wfk_path) # primitive lattice vectors of unit cell A[:,i] = (a_i) (3,3)

    nb, nk, _ = C_nkg.shape

    # supercell geometry
    A_sc, Omega_sc = get_A_volume(sc_wfk_path)    # primtive lattice vectors of supercell (3,3) and supercell volume float

    # Local potential of pristine and defective supercell
    Vp, _ = get_pot(sc_p_pot_path, subtract_mean); Vp = Vp.transpose(2, 1, 0) # (Nx, Ny, Nz)
    Vd, _ = get_pot(sc_d_pot_path, subtract_mean); Vd = Vd.transpose(2, 1, 0) # (Nx, Ny, Nz)
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

# NO MPI

def compute_ML_G_blocked(prep, block_size=512):

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

                # k'-k in supercell coordinates
                delta_k_sc = np.rint(
                    (k_red[ikp] - k_red[ik]) * Ndiag
                ).astype(int)

                # compute M per block
                M_block = np.zeros((nb, nb), dtype=np.complex128)

                # split the sum over G' into managable blocks of size block_size
                for start in range(0, nG_kp, block_size):
                    stop = min(start + block_size, nG_kp)

                    # G' block
                    Gp_blk = Gp[start:stop]                    # (B,3)
                    Ckp_blk = Ckp[:, start:stop]               # (nband,B)

                    # compute q = k' - k + G' - G 
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

    Ved_G   = prep["Ved_G"].copy()
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
    report_every = 5  # report every N completed local blocks (tune: 1..20)


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
