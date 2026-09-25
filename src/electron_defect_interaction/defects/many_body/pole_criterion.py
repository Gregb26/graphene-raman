"""
pole_criterion.py
    Pole criterion of the local t-matrix, as a function (R4, 2026-09-25). The production script
    scripts/resonance_criteria.py evaluates, inline (lines 39-54), for A(e) = 1 - V_loc g0(e):
        |det A| relative to its maximum on the window, the eigenvalue of A closest to zero, local minima.
    This module reproduces those lines (same products, same numpy calls) for an arbitrary sub-block of
    indices (parity blocks), adds the right eigenvector of the smallest eigenvalue, sign changes of Re lambda
    with a branch-identity check, and the K-point pi-pair trace of the t-matrix (resonance_metrics.py l. 86-90).
    Nothing in the production scripts is modified.
"""
import os

import numpy as np


def _workers():
    """Python-level threads for batched LAPACK calls (set R4_EIG_WORKERS or SLURM_CPUS_PER_TASK; default 1).
    Batched np.linalg.eig on ~150 x 150 matrices runs ~10x faster with 1 BLAS thread and several Python threads
    than with a multithreaded BLAS (measured 2026-09-25: 61 ms vs 4.7 ms per matrix on 16 cores)."""
    return max(1, int(os.environ.get("R4_EIG_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1"))))


def _batched(fn, A, workers=None):
    """Apply fn to chunks of the batch A (first axis) in Python threads (LAPACK releases the GIL); concatenate."""
    from concurrent.futures import ThreadPoolExecutor
    w = workers or _workers()
    if w <= 1 or A.shape[0] < 2 * w:
        return fn(A)
    chunks = np.array_split(A, w)
    with ThreadPoolExecutor(w) as ex:
        parts = list(ex.map(fn, chunks))
    if isinstance(parts[0], tuple):
        return tuple(np.concatenate([p[i] for p in parts]) for i in range(len(parts[0])))
    return np.concatenate(parts)


def block_indices(nL, nw, wfs):
    """Flat indices L*nw + w (the local_tmatrix layout) of the Wannier functions `wfs` on all nL cells."""
    wfs = list(wfs)
    return np.array([L * nw + w for L in range(nL) for w in wfs], dtype=int)


def det_eig_criterion(V, g0, idx=None, vectors=True):
    """
    For each energy j: A_j = I - V_b @ g0_j,b with V_b = V[idx][:, idx], g0_j,b = g0[j][idx][:, idx]
    (resonance_criteria.py l. 41: A = I - V g0). Returns dict:
        logabs      (nE,)   log|det A|            (slogdet)
        logdet_rel  (nE,)   logabs - max(logabs) over the energies given
        lam_all     (nE, n) eigenvalues of A      (np.linalg.eig, unsorted)
        minlam      (nE,)   min_i |lam_i|
        lam_min     (nE,)   the eigenvalue of smallest modulus
        vec_min     (nE, n) its right eigenvector (unit norm), if vectors
    """
    if idx is not None:
        idx = np.asarray(idx, int)
        Vb = V[np.ix_(idx, idx)]; gb = g0[:, idx][:, :, idx]
    else:
        Vb = V; gb = g0
    n = Vb.shape[0]; nE = gb.shape[0]
    I = np.eye(n)
    A = I[None] - Vb[None] @ gb
    sign, logabs = _batched(np.linalg.slogdet, A)
    if vectors:
        lam, vec = _batched(np.linalg.eig, A)
    else:
        lam = _batched(np.linalg.eigvals, A); vec = None
    ilam = np.abs(lam).argmin(1)
    lam_min = lam[np.arange(nE), ilam]
    out = dict(logabs=logabs, logdet_rel=logabs - logabs.max(), lam_all=lam, minlam=np.abs(lam).min(1),
               lam_min=lam_min, n=n)
    if vectors:
        vm = vec[np.arange(nE), :, ilam]
        out["vec_min"] = vm / np.linalg.norm(vm, axis=1)[:, None]
    return out


def nontrivial_eigenvalues(V, g0, support):
    """
    Eigenvalues of A = I - V g0 that differ from 1 when V is supported on the index set `support`:
    they are the eigenvalues of I_S - V_SS g0_SS (the other ones are exactly 1). Returns (nE, |S|) complex.
    """
    S = np.asarray(support, int)
    Vs = V[np.ix_(S, S)]; gs = g0[:, S][:, :, S]
    A = np.eye(len(S))[None] - Vs[None] @ gs
    return _batched(np.linalg.eigvals, A)


def local_minima(y, x, n=6):
    """resonance_criteria.py l. 45-47: indices of the n deepest local minima of y, sorted by x."""
    idx = [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] <= y[i + 1]]
    idx = sorted(idx, key=lambda i: y[i])[:n]
    return sorted(idx, key=lambda i: x[i])


def sign_changes(y, x):
    """Energies (linear interpolation) where the real array y changes sign between consecutive points."""
    y = np.asarray(y); x = np.asarray(x)
    s = np.sign(y); i = np.where((s[:-1] * s[1:]) < 0)[0]
    xc = x[i] - y[i] * (x[i + 1] - x[i]) / (y[i + 1] - y[i])
    return i, xc


def branch_overlaps(vec_min):
    """|<v_j | v_{j+1}>| between the smallest-eigenvalue eigenvectors of consecutive energies; (nE-1,)."""
    v = np.asarray(vec_min)
    return np.abs(np.einsum("ji,ji->j", np.conj(v[:-1]), v[1:]))


def eigvec_weights(vec, groups):
    """|v_i|^2 (unit-norm vector) summed over index groups: dict name -> array of flat indices."""
    p = np.abs(vec) ** 2
    p = p / p.sum()
    return {g: float(p[np.asarray(ix, int)].sum()) for g, ix in groups.items()}


def local_t_cache(V, g0):
    """t_j = V [1 - g0_j V]^{-1} for every energy (local_tmatrix.local_t, batched over energies)."""
    n = V.shape[0]; I = np.eye(n)
    A = I[None] - g0 @ V[None]
    X = _batched(lambda a: np.linalg.solve(a, np.broadcast_to(I, a.shape).copy()), A)
    return V[None] @ X


def tbar_pair(t_cache, phi_pair):
    """
    resonance_metrics.py l. 87-89: Tbar_j = P* t_j P^T for the pair phi_pair (2, dim); returns
    (Tbar (nE, 2, 2), tr = trace/2 (nE,)).
    """
    PK = np.asarray(phi_pair)
    Tbar = np.array([PK.conj() @ t_cache[j] @ PK.T for j in range(t_cache.shape[0])])
    return Tbar, 0.5 * np.trace(Tbar, axis1=1, axis2=2)
