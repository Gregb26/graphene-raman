"""
r5_sc_projection.py  (R5, 2026-09-25 ; module de campagne, répertoire de travail, hors dépôt)
    Gamma-point supercell states (plane waves, gamma trick) against the folded unit-cell Bloch basis.

    A unit-cell Bloch state |nk> (coefficients C_nk(G), sum_G |C|^2 = 1 on the unit cell, k on the N x N grid) is, on the
    N x N supercell, the state with the SAME coefficients at the supercell Miller indices g_sc = N (k + G); normalised on
    the supercell it is |nk>_sc = psi_nk / sqrt(N_cells). A supercell Gamma state Psi is stored by QE as D(g) on half the
    sphere (gamma_only), D(-g) = D*(g), |D(0)|^2 + 2 sum_{g != 0} |D(g)|^2 = 1. Then
        c_nk = <nk|Psi>_sc = sum_G C*_nk(G) D(N(k + G)),     sum_nk |c_nk|^2 <= 1,
    and delta = 1 - sum |c|^2 is the weight outside the retained bands (the union of the folded unit-cell spheres is the
    supercell sphere for equal ecut: 2 * 385 710 - 1 = 771 419 = sum_k igwx for the 9x9).
    Everything is done on the supercell FFT grid: A[flat] = D(g) (sphere completed), flat = supercell_fold.sc_planewave_index
    evaluated with the k of the .save that provides C (see r5_alignment_ext.check_k_reference).
"""
import numpy as np


def sc_full_sphere(C, mill):
    """Half-sphere coefficients C (..., ngw) and Miller indices mill (ngw, 3) -> (D (..., nfull), mill_full (nfull, 3))."""
    mill = np.asarray(mill, np.int64)
    nonzero = ~(mill == 0).all(axis=1)
    mill_full = np.concatenate([mill, -mill[nonzero]], axis=0)
    C = np.asarray(C)
    D = np.concatenate([C, np.conj(C[..., nonzero])], axis=-1)
    return D, mill_full


def grid_flat_index(mill, ngfft):
    """Flat C-order index on the (n1, n2, n3) FFT grid of Miller triplets (negative indices wrapped)."""
    n = np.asarray(ngfft, int); m = np.asarray(mill, np.int64)
    return np.ravel_multi_index((m[:, 0] % n[0], m[:, 1] % n[1], m[:, 2] % n[2]), tuple(int(x) for x in n))


def sc_state_grid(C, mill, ngfft, gamma_only=True):
    """Coefficients of one supercell state (C (ngw,)) on the flat FFT grid; A[-g] = conj(A[g]) added when gamma_only."""
    n = np.asarray(ngfft, int); N = int(np.prod(n))
    mill = np.asarray(mill, np.int64)
    A = np.zeros(N, dtype=complex)
    A[grid_flat_index(mill, n)] = C
    if gamma_only:
        nonzero = ~(mill == 0).all(axis=1)
        A[grid_flat_index(-mill[nonzero], n)] = np.conj(C[nonzero])
    return A


def grid_to_real(A, ngfft, Omega_sc, workers=8):
    """Psi(r) on the (n1, n2, n3) grid: Psi = N_grid * ifftn(A) / sqrt(Omega_sc), so that sum_r |Psi|^2 (Omega_sc/N_grid) = sum_g |A|^2."""
    from scipy import fft as sfft
    n1, n2, n3 = (int(x) for x in ngfft); N = n1 * n2 * n3
    return sfft.ifftn(np.asarray(A).reshape(n1, n2, n3), workers=workers) * (N / np.sqrt(float(Omega_sc)))


def bloch_overlaps(A, C_nkg, nG, flat_idx):
    """c_nk = sum_G conj(C_nk(G)) A[flat_idx[k][G]] for all (n, k); returns (nb, nk) complex."""
    nb, nk, _ = C_nkg.shape
    c = np.zeros((nb, nk), dtype=complex)
    for ik in range(nk):
        c[:, ik] = np.conj(C_nkg[:, ik, :nG[ik]]) @ A[flat_idx[ik]]
    return c


def bloch_state_grid(C_nkg, nG, flat_idx, n, k, ngfft):
    """|nk>_sc on the flat supercell grid (coefficients C_nk(G) at g_sc = N(k + G))."""
    A = np.zeros(int(np.prod(ngfft)), dtype=complex)
    A[flat_idx[k]] = C_nkg[n, k, :nG[k]]
    return A


def bloch_superposition_grid(d_nk, C_nkg, nG, flat_idx, ngfft):
    """sum_nk d_nk |nk>_sc on the flat supercell grid; d_nk (nb, nk)."""
    A = np.zeros(int(np.prod(ngfft)), dtype=complex)
    for ik in range(len(flat_idx)):
        A[flat_idx[ik]] += d_nk[:, ik] @ C_nkg[:, ik, :nG[ik]]
    return A


def check_planewave_union(flat_idx, n_expected=None):
    """(n_distinct, n_total, ok): the folded unit-cell plane waves must be distinct on the FFT grid (and equal n_expected if given)."""
    allidx = np.concatenate([np.asarray(f) for f in flat_idx])
    nu = int(len(np.unique(allidx)))
    ok = (nu == len(allidx)) and (n_expected is None or nu == int(n_expected))
    return nu, int(len(allidx)), bool(ok)
