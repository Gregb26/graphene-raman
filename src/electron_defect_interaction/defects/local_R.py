"""
local_R.py
    Python module to compute the local part of the electron-defect scattering matrix by directly evaluating the 
    integral in real space. Surprisingly fast.
"""

import numpy as np

from electron_defect_interaction.wavefunctions.fold_wfk_to_sc import compute_psi_nk_fold_sc
from electron_defect_interaction.utils.fft_utils import map_G_to_fft_grid
from tqdm import tqdm

# NO MPI 

def compute_ML_R(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path, subtract_mean=True, pristine=False, bands=None, io=None):
    """
    Computes the local part of the electron-defect interaction matrix in real space.

    Inputs:
        uc_wfk_path: str
            Path to the unit-cell wavefunctions (the prefix.save dir).
        sc_wfk_path:
            Path to the pristine supercell wavefunctions (only its geometry/volume is used here).
        sc_p_pot_path: str
            Path to the local potential of the pristine supercell (pp.x plot_num=1 filplot).
        sc_d_pot_path: str
            Path to the local potential of the defective supercell.
        bands: list of ints, optional
            Band indices to compute. If None, all bands are used (warning: the folded real-space grid
            for the supercell can be very large, so restrict the bands for big supercells).
        io: module, optional
            I/O backend exposing get_C_nk/get_G_red/get_k_red/get_A_volume/get_pot. Defaults to qe_io.
    """
    if io is None:
        from electron_defect_interaction.io import qe_io as io

    # Get necessary unit cells quantities
    C_nkg, nG = io.get_C_nk(uc_wfk_path) # planewave coeffs (nband, nkpt, nG_max) and number of active G per k (nkpt, )
    G_red = io.get_G_red(uc_wfk_path) # reciprocal lattice vectors in reduced coords of unit cell (nkpt, nG_max, 3)
    k_red = io.get_k_red(uc_wfk_path) # kpoints in reduced coords of unit cell (nkpt, 3)
    A_uc, _ = io.get_A_volume(uc_wfk_path) # primitive lattice vectors of the unit cell A[:, i]=a_i

    # Optionally restrict to a subset of bands (keeps the folded real-space arrays manageable)
    nband_all, nkpt, _ = C_nkg.shape
    bands = list(range(nband_all)) if bands is None else list(bands)
    nband = len(bands)

    # Get necessary super cell quantities
    A_sc, Omega_sc = io.get_A_volume(sc_wfk_path) # primitive lattice vectors and cell volume of the supercell
    Vp, _ = io.get_pot(sc_p_pot_path, subtract_mean); Vp = Vp.transpose(2,1,0) # pristine supercell local potential
    Vd, _ = io.get_pot(sc_d_pot_path, subtract_mean); Vd = Vd.transpose(2,1,0) # defective supercell local potential
    ngfft = Vp.shape # FFT grid shape

    # Compute defect potential
    if pristine:
        Ved = Vp
    else:
        Ved = Vd - Vp

    # Supercell scaling factor
    Ndiag = tuple(np.diag(np.rint(A_sc @ np.linalg.pinv(A_uc))))

    # Compute unict wavefunctions unfolded onto supercell from unit cell planewave coefficients
    print('Computing wavefunctions')
    psi = compute_psi_nk_fold_sc(C_nkg, nG, G_red, k_red, Omega_sc, Ndiag, ngfft, bands=bands)

    R = np.prod(ngfft)

    psi_r = psi.reshape(nband, nkpt, R) # (nband, nkpt, R)
    Ved_r = Ved.reshape(R) # (R,)
    Bk = nband * nkpt
    dV = Omega_sc / np.prod(ngfft)
    Psi = np.ascontiguousarray(psi_r.reshape(Bk, R)) # (Bk, R)
 
    def accumulate_M(Psi, Ved_r, dV, chunk=100_000):
        # Psi: (BK, R) complex128; Ved_r: (R,) real/complex
        BK, R = Psi.shape
        M = np.zeros((BK, BK), dtype=np.complex128)

        with tqdm(total=R, unit="r", desc="Accumulating M") as pbar:
            for start in range(0, R, chunk):
                sl = slice(start, min(start + chunk, R))
                PsiB = Psi[:, sl]                              # (BK, r)
                M += dV * ((PsiB.conj() * Ved_r[sl]) @ PsiB.T) # (BK, BK)
                pbar.update(sl.stop - sl.start)
        return M
    with np.errstate(all='ignore'):
        M_L = accumulate_M(Psi, Ved_r, dV)
    print('Done!')

    # Index convention [bra_band, k', ket_band, k], matching compute_ML_G and compute_M_NL.
    return M_L.reshape(nband, nkpt, nband, nkpt)

# MPI (grid-distributed real-space)

from mpi4py import MPI


def fourier_resample(V, new_shape):
    """Exact Fourier interpolation of a periodic real grid function onto a finer grid (zero padding in G)."""
    Vg = np.fft.fftshift(np.fft.fftn(V))
    old = np.asarray(V.shape); new = np.asarray(new_shape)
    if np.any(new < old):
        raise ValueError(f"fourier_resample: new grid {new_shape} must be >= old {V.shape} (no truncation)")
    out = np.zeros(tuple(new), dtype=complex)
    lo = (new - old) // 2
    out[lo[0]:lo[0] + old[0], lo[1]:lo[1] + old[1], lo[2]:lo[2] + old[2]] = Vg
    return np.fft.ifftn(np.fft.ifftshift(out)).real * (np.prod(new) / np.prod(old))


def prep_realspace_inputs(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path,
                          subtract_mean=False, pristine=False, bands=None, io=None):
    """
    Build the inputs for compute_ML_R_mpi. Reads the unit-cell wavefunctions, supercell geometry and
    the pristine/defective local potentials. Returns the unit-cell plane-wave coefficients (small) plus
    the defect potential on the supercell grid; each rank rebuilds the wavefunctions locally on its grid
    slab, so this dict stays small (no full real-space psi is ever materialised).
    """
    if io is None:
        from electron_defect_interaction.io import qe_io as io

    C_nkg, nG = io.get_C_nk(uc_wfk_path)
    G_red = io.get_G_red(uc_wfk_path)
    k_red = io.get_k_red(uc_wfk_path)
    A_uc, _ = io.get_A_volume(uc_wfk_path)
    A_sc, Omega_sc = io.get_A_volume(sc_wfk_path)

    Vp, _ = io.get_pot(sc_p_pot_path, subtract_mean); Vp = Vp.transpose(2, 1, 0)
    Vd, _ = io.get_pot(sc_d_pot_path, subtract_mean); Vd = Vd.transpose(2, 1, 0)
    Ved = Vp if pristine else (Vd - Vp)
    ngfft = Vp.shape

    Ndiag = np.rint(np.diag(A_sc @ np.linalg.pinv(A_uc))).astype(int)

    # The MPI kernels sample u_nk on the unit-cell grid ngfft // Ndiag and use ix % nxu: this REQUIRES the
    # supercell FFT grid to be Ndiag x (unit-cell grid). QE may pick an FFT-friendly size that is not
    # (graphene 7x7: 216 != 7 x 30). In that case V_ed is Fourier-resampled (zero padding in G space, exact for
    # a band-limited potential) onto the next commensurate grid; dV = Omega_sc / prod(ngfft) follows.
    ng = np.asarray(ngfft, int)
    if np.any(ng % Ndiag):
        ng_new = ((ng + Ndiag - 1) // Ndiag) * Ndiag
        Ved = fourier_resample(Ved, tuple(int(x) for x in ng_new)); ngfft = Ved.shape
        print(f"[prep_realspace_inputs] supercell FFT grid {tuple(int(x) for x in ng)} not a multiple of Ndiag={tuple(int(x) for x in Ndiag)}: "
              f"V_ed Fourier-resampled onto {ngfft}", flush=True)

    if bands is not None:
        C_nkg = C_nkg[list(bands), ...]

    return {
        "C_nkg": C_nkg, "nG": nG, "G_red": G_red, "k_red": k_red,
        "Ved": Ved, "ngfft": ngfft, "Ndiag": Ndiag, "Omega_sc": Omega_sc,
    }


def _build_u_uc(C_nkg, nG, G_red, ngfft_uc):
    """Bloch-periodic part u_nk(p) on the unit-cell grid (no k-phase, no 1/sqrt(Omega))."""
    nb, nk, _ = C_nkg.shape
    nxu, nyu, nzu = ngfft_uc
    Nuc = nxu * nyu * nzu
    mx, my, mz = map_G_to_fft_grid(ngfft_uc)
    u = np.zeros((nb, nk, nxu, nyu, nzu), dtype=np.complex128)
    for ik in range(nk):
        nGk = int(nG[ik]); Gk = G_red[ik, :nGk, :]
        ix = np.array([mx[int(g)] for g in Gk[:, 0]])
        iy = np.array([my[int(g)] for g in Gk[:, 1]])
        iz = np.array([mz[int(g)] for g in Gk[:, 2]])
        Cg = np.zeros((nb, nxu, nyu, nzu), dtype=np.complex128)
        Cg[:, ix, iy, iz] = C_nkg[:, ik, :nGk]
        u[:, ik] = np.fft.ifftn(Cg, axes=(1, 2, 3)) * Nuc
    return u.reshape(nb, nk, Nuc)


def compute_ML_R_mpi(prep, grid_block=200_000):
    """
    Real-space M^L parallelised by distributing the supercell real-space grid across MPI ranks.

    M_ab = dV * sum_r psi*_a(r) Ved(r) psi_b(r), a,b over (band, k). Each rank owns a slab of grid
    points, rebuilds psi on that slab from the stored unit-cell Bloch parts u_nk (psi_a(r) =
    u_nk(p(r)) e^{2pi i k.(i/N_uc)} / sqrt(Omega_sc)), forms the partial M with a BLAS matmul, and the
    ranks Allreduce the (Bk, Bk) result. Memory per rank ~ u (unit-cell sized) + one grid block of psi,
    independent of the number of ranks. Returns M[bra_band, k', ket_band, k] on every rank.
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank(); size = comm.Get_size()

    C_nkg = prep["C_nkg"]; nG = prep["nG"]; G_red = prep["G_red"]; k_red = prep["k_red"]
    Ved = prep["Ved"]; ngfft = prep["ngfft"]; Ndiag = prep["Ndiag"]; Omega_sc = float(prep["Omega_sc"])
    nb, nk, _ = C_nkg.shape
    Nx, Ny, Nz = ngfft
    nxu, nyu, nzu = Nx // Ndiag[0], Ny // Ndiag[1], Nz // Ndiag[2]
    if (Nx % Ndiag[0]) or (Ny % Ndiag[1]) or (Nz % Ndiag[2]):
        raise ValueError(f"compute_ML_R_mpi: supercell grid {ngfft} not commensurate with Ndiag={tuple(Ndiag)}; use prep_realspace_inputs (resamples)")

    N_cells = np.prod(Ndiag)

    u = _build_u_uc(C_nkg, nG, G_red, (nxu, nyu, nzu))   # (nb, nk, Nuc)
    Bk = nb * nk
    Ved_flat = Ved.reshape(-1)                            # C-order [ix,iy,iz]
    Ntot = Nx * Ny * Nz
    dV = Omega_sc / Ntot
    inv_sqrtO = 1.0 / np.sqrt(Omega_sc/N_cells)
    twopi = 2.0 * np.pi

    # Contiguous slab of grid points for this rank
    counts = [Ntot // size + (1 if r < (Ntot % size) else 0) for r in range(size)]
    displs = [sum(counts[:r]) for r in range(size)]
    g0, g1 = displs[rank], displs[rank] + counts[rank]

    M_local = np.zeros((Bk, Bk), dtype=np.complex128)
    for start in range(g0, g1, grid_block):
        stop = min(start + grid_block, g1)
        idx = np.arange(start, stop)
        ix = idx // (Ny * Nz); iy = (idx // Nz) % Ny; iz = idx % Nz

        # unit-cell linear index p(r) and the Bloch phase e^{2pi i k . (i/N_uc)} (per k)
        pu = (ix % nxu) * (nyu * nzu) + (iy % nyu) * nzu + (iz % nzu)
        fx = ix / nxu; fy = iy / nyu; fz = iz / nzu
        arg = twopi * (k_red[:, 0:1] * fx[None, :] + k_red[:, 1:2] * fy[None, :] + k_red[:, 2:3] * fz[None, :])
        phase = np.exp(1j * arg)                          # (nk, nblock)

        # psi on this block: (nb, nk, nblock) -> (Bk, nblock)
        Psi = (u[:, :, pu] * phase[None, :, :]) * inv_sqrtO
        Psi = Psi.reshape(Bk, -1)

        Vb = Ved_flat[start:stop]
        M_local += dV * (Psi.conj() * Vb[None, :]) @ Psi.T

    # Reduce in sub-buffer chunks: OpenMPI's large-message Allreduce on DOUBLE_COMPLEX
    # corrupts the heap (latent "double free") once (nb*nk)^2 crosses ~1.3M elements
    # (i.e. nk >= 81). Chunking each Allreduce well under that threshold avoids it.
    M = np.zeros((Bk, Bk), dtype=np.complex128)
    flat_local = M_local.reshape(-1)
    flat = M.reshape(-1)
    chunk = 1_000_000
    for s in range(0, flat.size, chunk):
        e = min(s + chunk, flat.size)
        comm.Allreduce(flat_local[s:e].copy(), flat[s:e], op=MPI.SUM)
    return M.reshape(nb, nk, nb, nk)


# ---------------------------------------------------------------------------------------------- #
# Node-shared variant for LARGE k sets (dense primitive grids -> exact zero-padded dense M^L).     #
# ---------------------------------------------------------------------------------------------- #
def _node_shared(nodecomm, shape, dtype):
    """One physical copy per node of an array of `shape`, exposed to every rank of `nodecomm`."""
    itemsize = np.dtype(dtype).itemsize
    n = int(np.prod(shape))
    win = MPI.Win.Allocate_shared(n * itemsize if nodecomm.Get_rank() == 0 else 0, itemsize, comm=nodecomm)
    buf, _ = win.Shared_query(0)
    return np.ndarray(buffer=buf, dtype=dtype, shape=shape), win


def compute_ML_R_mpi_shared(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path,
                            subtract_mean=False, bands=None, io=None, grid_block=2000, verbose=True):
    """
    Same physics and conventions as compute_ML_R_mpi -- M[bra_band, k', ket_band, k] =
    dV sum_r psi*_{n'k'}(r) Ved(r) psi_{nk}(r) over the SUPERCELL grid, psi built pointwise from the
    unit-cell Bloch part u_nk and the phase e^{2 pi i k.r} -- but organised for LARGE k sets:

      * u_nk (nb, nk, N_uc) is stored ONCE per node in MPI shared memory (tens of GB on dense grids)
        instead of once per rank; the node-local root reads the inputs and fills it;
      * the supercell grid is distributed over ALL ranks (contiguous slabs), each rank accumulating
        its (Bk, Bk) partial with threaded BLAS on blocks of `grid_block` points;
      * the partial sums are Reduce'd to rank 0 in sub-buffer chunks (OpenMPI large-message bug).

    Because the Bloch phase is evaluated per grid point, k need NOT be commensurate with the
    supercell. Fed with the wavefunctions of a DENSE primitive calculation and the N x N supercell
    potential, this is EXACTLY the zero-padded dense M^L of local_G (Ved = 0 outside the supercell,
    so the padded-cell integral reduces to the supercell one), in O(N_r Bk^2) BLAS flops instead of
    the O(nk^2 nG^2) random gathers of compute_ML_G_dense_mpi. Returns M on rank 0, None elsewhere.
    """
    import time
    if io is None:
        from electron_defect_interaction.io import qe_io as io
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank(); size = comm.Get_size()
    node = comm.Split_type(MPI.COMM_TYPE_SHARED)
    nrank = node.Get_rank()

    meta = None
    prep = None
    if nrank == 0:
        prep = prep_realspace_inputs(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path,
                                     subtract_mean=subtract_mean, bands=bands, io=io)
        nb, nk, _ = prep["C_nkg"].shape
        Nx, Ny, Nz = prep["ngfft"]; Ndiag = prep["Ndiag"]
        meta = dict(nb=int(nb), nk=int(nk), ngfft=(int(Nx), int(Ny), int(Nz)),
                    ngfft_uc=(int(Nx // Ndiag[0]), int(Ny // Ndiag[1]), int(Nz // Ndiag[2])),
                    Omega_sc=float(prep["Omega_sc"]), k_red=np.asarray(prep["k_red"], dtype=float), Ndiag=tuple(int(x) for x in Ndiag))
    meta = node.bcast(meta, root=0)
    nb, nk = meta["nb"], meta["nk"]
    Nx, Ny, Nz = meta["ngfft"]; nxu, nyu, nzu = meta["ngfft_uc"]
    if (Nx % nxu) or (Ny % nyu) or (Nz % nzu) or (Nx // nxu, Ny // nyu, Nz // nzu) != tuple(int(x) for x in meta["Ndiag"]):
        raise ValueError(f"compute_ML_R_mpi_shared: supercell grid {meta['ngfft']} not commensurate with Ndiag={meta['Ndiag']}")
    Nuc = nxu * nyu * nzu; Ntot = Nx * Ny * Nz
    Omega_sc = meta["Omega_sc"]; k_red = meta["k_red"]

    u, win_u = _node_shared(node, (nb, nk, Nuc), np.complex128)
    Ved_flat, win_v = _node_shared(node, (Ntot,), np.float64)
    win_u.Sync(); win_v.Sync()
    if nrank == 0:
        t0 = time.time()
        C_nkg, nG, G_red = prep["C_nkg"], prep["nG"], prep["G_red"]
        mx, my, mz = map_G_to_fft_grid((nxu, nyu, nzu))
        for ik in range(nk):
            nGk = int(nG[ik]); Gk = G_red[ik, :nGk, :]
            ix = np.fromiter((mx[int(g)] for g in Gk[:, 0]), dtype=np.int64, count=nGk)
            iy = np.fromiter((my[int(g)] for g in Gk[:, 1]), dtype=np.int64, count=nGk)
            iz = np.fromiter((mz[int(g)] for g in Gk[:, 2]), dtype=np.int64, count=nGk)
            Cg = np.zeros((nb, nxu, nyu, nzu), dtype=np.complex128)
            Cg[:, ix, iy, iz] = C_nkg[:, ik, :nGk]
            u[:, ik, :] = (np.fft.ifftn(Cg, axes=(1, 2, 3)) * Nuc).reshape(nb, Nuc)
        Ved_flat[:] = np.ascontiguousarray(prep["Ved"]).reshape(-1)
        del prep, C_nkg, Cg
        if verbose and rank == 0:
            print(f"[ML_R shared] u_nk built: nb={nb} nk={nk} Nuc={Nuc} ({nb*nk*Nuc*16/1e9:.1f} GB/node) "
                  f"grid={Nx}x{Ny}x{Nz} Bk={nb*nk} in {time.time()-t0:.0f}s", flush=True)
    win_u.Sync(); win_v.Sync()
    node.Barrier()
    win_u.Sync(); win_v.Sync()

    N_cells = np.prod(meta["Ndiag"])
    Bk = nb * nk
    dV = Omega_sc / Ntot
    inv_sqrtO = 1.0 / np.sqrt(Omega_sc/N_cells)
    twopi = 2.0 * np.pi

    counts = [Ntot // size + (1 if r < (Ntot % size) else 0) for r in range(size)]
    displs = [sum(counts[:r]) for r in range(size)]
    g0, g1 = displs[rank], displs[rank] + counts[rank]

    M_local = np.zeros((Bk, Bk), dtype=np.complex128)
    nblocks = max(1, -(-(g1 - g0) // grid_block)); t0 = time.time(); ib = 0
    for start in range(g0, g1, grid_block):
        stop = min(start + grid_block, g1)
        idx = np.arange(start, stop)
        ix = idx // (Ny * Nz); iy = (idx // Nz) % Ny; iz = idx % Nz
        pu = (ix % nxu) * (nyu * nzu) + (iy % nyu) * nzu + (iz % nzu)
        fx = ix / nxu; fy = iy / nyu; fz = iz / nzu
        arg = twopi * (k_red[:, 0:1] * fx[None, :] + k_red[:, 1:2] * fy[None, :] + k_red[:, 2:3] * fz[None, :])
        phase = np.exp(1j * arg)                                     # (nk, nblock)
        Psi = (u[:, :, pu] * phase[None, :, :]) * inv_sqrtO          # (nb, nk, nblock)
        Psi = Psi.reshape(Bk, -1)
        Vb = Ved_flat[start:stop]
        M_local += dV * (Psi.conj() * Vb[None, :]) @ Psi.T
        ib += 1
        if verbose and rank == 0 and (ib % max(1, nblocks // 20) == 0 or ib == nblocks):
            el = time.time() - t0
            print(f"[ML_R shared] rank0 block {ib}/{nblocks}  elapsed {el/60:.1f} min  ETA {el/ib*(nblocks-ib)/60:.1f} min", flush=True)
    del Psi

    M = np.zeros((Bk, Bk), dtype=np.complex128) if rank == 0 else None
    flat_local = M_local.reshape(-1)
    flat = M.reshape(-1) if rank == 0 else None
    chunk = 1_000_000
    for s in range(0, Bk * Bk, chunk):
        e = min(s + chunk, Bk * Bk)
        comm.Reduce(flat_local[s:e], flat[s:e] if rank == 0 else None, op=MPI.SUM, root=0)
    del M_local, flat_local
    comm.Barrier()
    win_u.Free(); win_v.Free(); node.Free()
    return M.reshape(nb, nk, nb, nk) if rank == 0 else None
