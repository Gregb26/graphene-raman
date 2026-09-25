"""
supercell_fold.py
    Gamma-point supercell Hamiltonians for the D4 reconstruction of R4 (2026-09-25).

    Folded Bloch basis |n k>, k on the N x N unshifted MP grid (which folds onto Gamma of the N x N supercell):
        bloch_folded_hamiltonian  -- H = diag(eps_nk) + M[n,k,n',k'] / N_cells   (M with unit-cell Bloch norm)
    Wannier bases on the same 81 k (|w k>) or on the N x N cells (|w r>):
        kblocks_to_rbasis         -- H_r = F^dagger H_k F, F_{k r} = e^{-2 pi i k.r}/N   (block-diagonal H_k)
        fold_hwr_to_supercell     -- H_SC[(r,w),(r',w')] = sum_{R = r'-r mod N} H(R)[w,w']/ndegen(R)
        fold_mwr_to_supercell     -- M_SC[(r,w),(r',w')] = sum_{R = r, R' = r' mod N} Mwr[w,R,w',R']
                                     (optionally restricted to R, R' in R_local before folding)
    Index layout of the r-basis: c*nw + w with c = r1*N + r2 (r = (r1, r2, 0), 0 <= ri < N), the same
    "cell-major" layout as local_tmatrix (L*nw + w). Real-space weights of eigenvectors:
        sc_planewave_index        -- supercell FFT index of each unit-cell plane wave k + G
        folded_density_2d         -- |Psi(r)|^2 summed over z for Psi = sum_nk d_nk psi_nk (unit-cell Bloch states)
"""
import numpy as np

from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk


def supercell_cells(N):
    """(N*N, 3) integer cells r = (r1, r2, 0) in cell-major order c = r1*N + r2."""
    r1, r2 = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    return np.stack([r1.ravel(), r2.ravel(), np.zeros(N * N, int)], axis=1)


def cell_of(R, N):
    """Cell index c = (R1 mod N)*N + (R2 mod N) of integer vectors R (..., 3)."""
    R = np.asarray(R, int)
    return (R[..., 0] % N) * N + (R[..., 1] % N)


def bloch_folded_hamiltonian(eps, M, N_cells):
    """
    eps (nb, nk), M (nb, nk, nb, nk) [bra_band, k', ket_band, k], unit-cell Bloch normalization.
    Returns H (nb*nk, nb*nk) with index n*nk + k: diag(eps) + M/N_cells.
    """
    nb, nk = eps.shape
    H = np.asarray(M).reshape(nb * nk, nb * nk) / float(N_cells)
    H = H.astype(complex)
    H[np.diag_indices(nb * nk)] += np.asarray(eps).reshape(-1)
    return H


def fourier_matrix(k_red, N):
    """F (nk, N*N) with F[k, c] = e^{-2 pi i k.r_c} / N (unitary when k is the full N x N grid)."""
    cells = supercell_cells(N)
    return np.exp(-2j * np.pi * (np.asarray(k_red, float) @ cells.T)) / N


def kblocks_to_rbasis(Hk, k_red, N):
    """H_r[(c,w),(c',w')] = sum_k conj(F[k,c]) Hk[k,w,w'] F[k,c']; layout c*nw + w. Hk (nk, nw, nw)."""
    F = fourier_matrix(k_red, N); nw = Hk.shape[1]
    Hr = np.einsum("kc,kwv,kd->cwdv", np.conj(F), Hk, F, optimize=True)
    return Hr.reshape(N * N * nw, N * N * nw)


def kbasis_matrix(Hk_blocks, Mwk, k_red, N_cells):
    """H_k[(k,w),(k',w')] = delta_kk' Hk_blocks[k] + Mwk[w,k,w',k']/N_cells ; layout k*nw + w."""
    nk, nw, _ = Hk_blocks.shape
    H = np.transpose(np.asarray(Mwk), (1, 0, 3, 2)).reshape(nk * nw, nk * nw) / float(N_cells)
    H = H.astype(complex)
    for k in range(nk):
        H[k * nw:(k + 1) * nw, k * nw:(k + 1) * nw] += Hk_blocks[k]
    return H


def kvec_to_rvec(c_k, k_red, N, nw):
    """Eigenvector in the k-basis (layout k*nw + w) -> r-basis (layout c*nw + w): c_r = F^dagger c_k."""
    F = fourier_matrix(k_red, N)
    ck = np.asarray(c_k).reshape(len(k_red), nw)
    return np.einsum("kc,kw->cw", np.conj(F), ck).reshape(-1)


def rvec_to_kvec(c_r, k_red, N, nw):
    """Inverse of kvec_to_rvec: c_k = F c_r."""
    F = fourier_matrix(k_red, N)
    cr = np.asarray(c_r).reshape(N * N, nw)
    return np.einsum("kc,cw->kw", F, cr).reshape(-1)


def fold_hwr_to_supercell(Hwr, Rw, ndegen, N):
    """H_SC (N*N*nw, N*N*nw) at Gamma from the Wannier90 H(R) (ndegen divided out), layout c*nw + w."""
    nw = Hwr.shape[1]
    Hd = np.zeros((N * N, nw, nw), dtype=complex)                  # class d = (R mod N)
    for r, R in enumerate(Rw):
        Hd[cell_of(R, N)] += Hwr[r] / ndegen[r]
    cells = supercell_cells(N)
    H = np.zeros((N * N, nw, N * N, nw), dtype=complex)
    for c in range(N * N):
        d = cell_of(cells - cells[c], N)                            # class of r' - r for every r'
        H[c, :, :, :] = np.transpose(Hd[d], (1, 0, 2))              # (nw, N*N, nw)
    return H.reshape(N * N * nw, N * N * nw)


def fold_mwr_to_supercell(Mwr, R_mwr, N, R_local=None):
    """
    Mwr (nw, nr, nw, nr) on the periodic R box R_mwr (nr, 3) -> M_SC (N*N*nw, N*N*nw), layout c*nw + w.
    R_local: if given, only R, R' in R_local are kept (V_loc = P Mwr P) before folding.
    """
    nw, nr = Mwr.shape[0], Mwr.shape[1]
    keep = np.ones(nr, bool)
    if R_local is not None:
        keep[:] = False
        for R in np.asarray(R_local, int):
            keep |= (np.asarray(R_mwr, int) == R).all(axis=1)
    P = np.zeros((nr, N * N))
    P[np.arange(nr)[keep], cell_of(np.asarray(R_mwr)[keep], N)] = 1.0
    M = np.einsum("rc,wrvs,sd->cwdv", P, Mwr, P, optimize=True)
    return M.reshape(N * N * nw, N * N * nw)


def wannier_gate(Hwr, Rw, ndegen, N, atol=1e-8):
    """Gate 1 of D4: eigenvalues of fold_hwr_to_supercell == union of the Wannier eigenvalues on the N x N grid."""
    from electron_defect_interaction.defects.many_body.local_tmatrix import mp_grid
    k = mp_grid(N, N, 1)
    _, E, _ = Hwr_to_Hwk(Hwr, Rw, k, ndegen=ndegen)
    HS = fold_hwr_to_supercell(Hwr, Rw, ndegen, N)
    herm = float(np.abs(HS - HS.conj().T).max())
    e = np.linalg.eigvalsh(0.5 * (HS + HS.conj().T))
    dev = float(np.abs(e - np.sort(E.reshape(-1))).max())
    return dev, herm, (dev < atol)


def sc_planewave_index(k_red, G_red, nG, Ndiag, ngfft_sc):
    """
    For each unit-cell k and its active G: supercell Miller index g_sc = Ndiag * (k + G) (integer for k on the
    N x N grid) and the flat index on the supercell FFT grid. Returns list of (nG_k,) int arrays.
    """
    n = np.asarray(ngfft_sc, int); Nd = np.asarray(Ndiag, float)
    out = []
    for ik in range(len(k_red)):
        g = np.rint((k_red[ik][None, :] + G_red[ik, :nG[ik]]) * Nd[None, :]).astype(int)
        out.append(np.ravel_multi_index((g[:, 0] % n[0], g[:, 1] % n[1], g[:, 2] % n[2]), tuple(n)))
    return out


def folded_density_2d(d_nk, C_nkg, nG, flat_idx, ngfft_sc, workers=8):
    """
    |Psi(r)|^2 summed over z on the supercell in-plane grid, Psi = sum_{n,k} d_nk psi_nk, psi_nk the unit-cell
    Bloch states (plane-wave coefficients C_nkg (nb, nk, nGmax)); d_nk (nb, nk). Normalized to unit sum.
    """
    from scipy import fft as sfft
    n1, n2, n3 = (int(x) for x in ngfft_sc)
    A = np.zeros(n1 * n2 * n3, dtype=complex)
    for ik in range(len(flat_idx)):
        A[flat_idx[ik]] += d_nk[:, ik] @ C_nkg[:, ik, :nG[ik]]
    psi = sfft.ifftn(A.reshape(n1, n2, n3), workers=workers)
    r2 = (psi.real ** 2 + psi.imag ** 2).sum(axis=2)
    return r2 / r2.sum()
