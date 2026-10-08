"""
tb_models.py
    Defect potentials M_cluster for synthetic tight-binding benches (D6 bench of R4, 2026-09-25). Test material. The
    model itself is electron_photon.tb_model.make_graphene_tb: (tb.H_R, tb.R_int, tb.ndegen) plays (Hwr, Rw, ndegen).

    onsite_M_cluster   -- M_cluster with a single on-site element U on one Wannier function of one cell.
    removed_site_M_cluster -- M_cluster = -H on the bonds of a site (all hoppings from (R0, wf)), plus U on the site.
"""
import numpy as np


def cell_index(R_cluster, R):
    hit = np.where((np.asarray(R_cluster, int) == np.asarray(R, int)).all(axis=1))[0]
    if len(hit) == 0:
        raise ValueError(f"cell {R} not in R_cluster")
    return int(hit[0])


def onsite_M_cluster(R_cluster, nw, wf, U, R0=(0, 0, 0)):
    """M_cluster (nL*nw, nL*nw) with the single element U at (R0, wf)."""
    nL = len(R_cluster); M_cluster = np.zeros((nL * nw, nL * nw), dtype=complex)
    i = cell_index(R_cluster, R0) * nw + wf
    M_cluster[i, i] = U
    return M_cluster, [i]


def removed_site_M_cluster(Hwr, Rw, ndegen, R_cluster, nw, wf, U, R0=(0, 0, 0)):
    """
    M_cluster = -H on every hopping from (R0, wf): M[(R0,wf),(R0+R, w)] = -H(R)[wf, w]/ndegen(R) for w != wf or
    R != 0, plus the Hermitian counterpart, plus U on (R0, wf). Cells R0+R must be in R_cluster.
    Returns (V, support indices).
    """
    nL = len(R_cluster); M_cluster = np.zeros((nL * nw, nL * nw), dtype=complex)
    i0 = cell_index(R_cluster, R0) * nw + wf
    supp = {i0}
    for r, R in enumerate(Rw):
        for w in range(nw):
            h = Hwr[r, wf, w] / ndegen[r]
            if abs(h) == 0 or (w == wf and (np.asarray(R) == 0).all()):
                continue
            j = cell_index(R_cluster, np.asarray(R0) + np.asarray(R)) * nw + w
            M_cluster[i0, j] += -h; M_cluster[j, i0] += -np.conj(h); supp.add(j)
    M_cluster[i0, i0] += U
    return M_cluster, sorted(supp)
