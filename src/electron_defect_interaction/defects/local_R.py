"""
local_R.py
    Python module to compute the local part of the electron-defect scattering matrix by directly evaluating the 
    integral in real space. Surprisingly fast.
"""

import numpy as np

from electron_defect_interaction.io.abinit_io import *
from electron_defect_interaction.wavefunctions.wfk import compute_psi_nk
from electron_defect_interaction.wavefunctions.fold_wfk_to_sc import compute_psi_nk_fold_sc
from electron_defect_interaction.utils.fft_utils import map_G_to_fft_grid
from tqdm import tqdm

def compute_ML_R(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path, subtract_mean=True, pristine=False, bands=None, io=None):
    """
    Computes the local part of the electron-defect interaction matrix in real space.

    Inputs:
        uc_wfk_path: str
            Path to the unit-cell wavefunctions. For ABINIT a WFK.nc file; for QE the prefix.save dir.
        sc_wfk_path:
            Path to the pristine supercell wavefunctions (only its geometry/volume is used here).
        sc_p_pot_path: str
            Path to the local potential of the pristine supercell (ABINIT POT.nc / QE pp.x plot_file).
        sc_d_pot_path: str
            Path to the local potential of the defective supercell.
        bands: list of ints, optional
            Band indices to compute. If None, all bands are used (warning: the folded real-space grid
            for the supercell can be very large, so restrict the bands for big supercells).
        io: module, optional
            I/O backend exposing get_C_nk/get_G_red/get_k_red/get_A_volume/get_pot. Defaults to abinit_io;
            pass electron_defect_interaction.io.qe_io for Quantum ESPRESSO inputs.
    """
    if io is None:
        from electron_defect_interaction.io import abinit_io as io

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

def compute_ML_R_uc(wfk_uc_path, pot_uc_path, subtract_mean=False):
    """
    Computes the local part of the electron-defect interaction matrix in real space.

    Inputs:
        uc_wfk_path: str
            Path to the ABINIT WFK.nc output file for the wavefunctions of the unit cell
        sc_wfk_path:
            Path to the ABINIT WFK.nc output file for the wavefunctions of the pristine super cell
        sc_p_pot_path: str
            Path to the ABINIT POT.nc output file for the local potential of the pristine super cell
        sc_d_pot_path: str
            Path to the ABINIT POT.nc output file for the local potential of the defective super cell
        bands: list of ints:
            Band indices at which to compute the matrix elements
    """

    # Get necessary unit cells quantities
    C_nkg, nG = get_C_nk(wfk_uc_path) # planewave coeffs (nband, nkpt, nG_max) and number of active G per k (nkpt, )
    G_red = get_G_red(wfk_uc_path) # reciprocal lattice vectors in reduced coords of unit cell (nkpt, nG_max, 3)
    k_red = get_k_red(wfk_uc_path) # kpoints in reduced coords of unit cell (nkpt, 3)
    _, Omega = get_A_volume(wfk_uc_path) # primitive lattice vectors of the unit cell A[:, i]=a_i
    nband, nkpt, _ = C_nkg.shape

    V,_ = get_pot(pot_uc_path, subtract_mean); V = V.transpose(2,1,0) # potential
    ngfft = V.shape # FFT grid shape
    print("Potential FFT grid: ", ngfft)
    # Compute unict wavefunctions unfolded onto supercell from unit cell planewave coefficients
    print('Computing wavefunctions')
    psi, ngfft_wfk = compute_psi_nk(C_nkg, nG, G_red, k_red, Omega, check_normalize=True, ngfft=ngfft)
    print("Wfk FFT grid: ", ngfft_wfk)
    # V = np.ones(ngfft)
    print(V)
    nr = np.prod(ngfft) # number of grid points
    
    psi_r = psi.reshape(nband, nkpt, nr) # (nband, nkpt, nr)
    Vr = V.reshape(nr) # (nr,)
    Bk = nband * nkpt
    dV = Omega / np.prod(ngfft)
    Psi = np.ascontiguousarray(psi_r.reshape(Bk, nr)) # (Bk, R)
 
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
        M_L = accumulate_M(Psi, Vr, dV)
    print('Done!')

    # Index convention [bra_band, k', ket_band, k], matching compute_ML_G and compute_M_NL.
    return M_L.reshape(nband, nkpt, nband, nkpt)

# MPI (grid-distributed real-space)

from mpi4py import MPI


def prep_realspace_inputs(uc_wfk_path, sc_wfk_path, sc_p_pot_path, sc_d_pot_path,
                          subtract_mean=False, pristine=False, bands=None, io=None):
    """
    Build the inputs for compute_ML_R_mpi. Reads the unit-cell wavefunctions, supercell geometry and
    the pristine/defective local potentials. Returns the unit-cell plane-wave coefficients (small) plus
    the defect potential on the supercell grid; each rank rebuilds the wavefunctions locally on its grid
    slab, so this dict stays small (no full real-space psi is ever materialised).
    """
    if io is None:
        from electron_defect_interaction.io import abinit_io as io

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

    u = _build_u_uc(C_nkg, nG, G_red, (nxu, nyu, nzu))   # (nb, nk, Nuc)
    Bk = nb * nk
    Ved_flat = Ved.reshape(-1)                            # C-order [ix,iy,iz]
    Ntot = Nx * Ny * Nz
    dV = Omega_sc / Ntot
    inv_sqrtO = 1.0 / np.sqrt(Omega_sc)
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

    M = np.zeros((Bk, Bk), dtype=np.complex128)
    comm.Allreduce(M_local, M, op=MPI.SUM)
    return M.reshape(nb, nk, nb, nk)
