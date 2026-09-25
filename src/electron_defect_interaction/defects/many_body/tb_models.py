"""
tb_models.py
    Synthetic tight-binding material for the D6 test bench of R4 (2026-09-25). Test material: the report
    says whether it stays in src/ or moves under tests/.

    graphene_pz_tb     -- nearest-neighbour pi model of graphene written in the 5-Wannier-function layout of
                          the production wannierization (0-2: sp2 on A, 3: pz on A, 4: pz on B); the sigma
                          functions are decoupled (on-site energy far from the window, no hopping). The three
                          A-B bonds are at the same R as in wannier_tb.dat (H(R)[3,4] != 0 for R = 0, (-1,0,0),
                          (0,-1,0)), with H(R)[m,n] = <m 0| H |n R>.
    pz_bands_analytic  -- +-t |1 + e^{-2 pi i k1} + e^{-2 pi i k2}| on reduced k.
    check_pz_model     -- bands on Gamma-K-M-Gamma vs analytic, van Hove at +-t, E_D = 0, on-site g0 from
                          local_green_batch vs a direct sum over the grid.
    onsite_vloc        -- V_loc with a single on-site element U on one Wannier function of one cell.
    removed_site_vloc  -- V_loc = -H on the bonds of a site (all hoppings from (R0, wf)), plus U on the site.
"""
import numpy as np

from electron_defect_interaction.wannier.wannier_hamiltonian import Hwr_to_Hwk
from electron_defect_interaction.defects.many_body import local_tmatrix as lt


def graphene_pz_tb(t=2.7, e_sigma=-100.0, e_pz=0.0):
    """Returns (Hwr (5, 5, 5), Rw (5, 3), ndegen (5,)) with hopping -t (production sign) on the pz A-B bonds."""
    Rw = np.array([[0, 0, 0], [-1, 0, 0], [0, -1, 0], [1, 0, 0], [0, 1, 0]], dtype=int)
    H = np.zeros((5, 5, 5), dtype=complex)
    for w in range(3):
        H[0, w, w] = e_sigma
    H[0, 3, 3] = e_pz; H[0, 4, 4] = e_pz
    H[0, 3, 4] = -t; H[0, 4, 3] = -t            # bond within the cell
    H[1, 3, 4] = -t                             # <3, 0| H |4, (-1,0,0)>
    H[2, 3, 4] = -t                             # <3, 0| H |4, (0,-1,0)>
    H[3, 4, 3] = -t                             # Hermiticity: H(-R) = H(R)^dagger
    H[4, 4, 3] = -t
    return H, Rw, np.ones(len(Rw), dtype=int)


def pz_bands_analytic(k_red, t):
    k = np.asarray(k_red, float)
    f = 1.0 + np.exp(-2j * np.pi * k[:, 0]) + np.exp(-2j * np.pi * k[:, 1])
    return np.abs(f) * t


def check_pz_model(Hwr, Rw, ndegen, t, eta=0.02, nk=300, eps_test=-1.0, npath=61):
    """
    Consistency checks of the synthetic model. Returns a dict of raw numbers:
      path_max_dev  : max |E_pi(k) - (-t|f|)| and |E_pi*(k) - (+t|f|)| on Gamma-K-M-Gamma (eV)
      E_K_pi, E_K_pistar, E_M_pi, E_M_pistar : energies at K and M (van Hove expected at -t, +t)
      gap_K_grid    : min pi/pi* gap on the nk x nk grid (E_D = midpoint)
      E_D           : midpoint of the minimal gap on the nk x nk grid
      g0_batch, g0_direct : on-site pz(A) g0 at eps_test (cluster {0}, local_green_batch) vs direct sum
    """
    G = np.array([0, 0, 0.]); K = np.array([2 / 3, 1 / 3, 0.]); M = np.array([0.5, 0, 0.])
    path = np.concatenate([np.linspace(G, K, npath), np.linspace(K, M, npath), np.linspace(M, G, npath)])
    _, E, _ = Hwr_to_Hwk(Hwr, Rw, path, ndegen=ndegen)
    Ea = pz_bands_analytic(path, t)
    dev = max(np.abs(E[:, 3] + Ea).max(), np.abs(E[:, 4] - Ea).max())      # bands 3, 4 = pi, pi* (sigma far below)
    _, EK, _ = Hwr_to_Hwk(Hwr, Rw, K[None], ndegen=ndegen); _, EM, _ = Hwr_to_Hwk(Hwr, Rw, M[None], ndegen=ndegen)
    kg = lt.mp_grid(nk, nk, 1)
    Hwk, Eg, _ = Hwr_to_Hwk(Hwr, Rw, kg, ndegen=ndegen)
    gap = Eg[:, 4] - Eg[:, 3]; iD = int(np.argmin(gap)); E_D = 0.5 * (Eg[iD, 3] + Eg[iD, 4])
    g0b = lt.local_green_batch(Hwk, kg, np.array([[0, 0, 0]]), np.array([eps_test]), eta)[0][3, 3]
    A = (eps_test + 1j * eta) * np.eye(5)[None] - Hwk
    g0d = np.linalg.inv(A)[:, 3, 3].mean()
    return dict(path_max_dev=float(dev), E_K_pi=float(EK[0, 3]), E_K_pistar=float(EK[0, 4]),
                E_M_pi=float(EM[0, 3]), E_M_pistar=float(EM[0, 4]), gap_K_grid=float(gap[iD]), E_D=float(E_D),
                g0_batch=complex(g0b), g0_direct=complex(g0d))


def cell_index(R_local, R):
    hit = np.where((np.asarray(R_local, int) == np.asarray(R, int)).all(axis=1))[0]
    if len(hit) == 0:
        raise ValueError(f"cell {R} not in R_local")
    return int(hit[0])


def onsite_vloc(R_local, nw, wf, U, R0=(0, 0, 0)):
    """V_loc (nL*nw, nL*nw) with the single element U at (R0, wf)."""
    nL = len(R_local); V = np.zeros((nL * nw, nL * nw), dtype=complex)
    i = cell_index(R_local, R0) * nw + wf
    V[i, i] = U
    return V, [i]


def removed_site_vloc(Hwr, Rw, ndegen, R_local, nw, wf, U, R0=(0, 0, 0)):
    """
    V_loc = -H on every hopping from (R0, wf): V[(R0,wf),(R0+R, w)] = -H(R)[wf, w]/ndegen(R) for w != wf or
    R != 0, plus the Hermitian counterpart, plus U on (R0, wf). Cells R0+R must be in R_local.
    Returns (V, support indices).
    """
    nL = len(R_local); V = np.zeros((nL * nw, nL * nw), dtype=complex)
    i0 = cell_index(R_local, R0) * nw + wf
    supp = {i0}
    for r, R in enumerate(Rw):
        for w in range(nw):
            h = Hwr[r, wf, w] / ndegen[r]
            if abs(h) == 0 or (w == wf and (np.asarray(R) == 0).all()):
                continue
            j = cell_index(R_local, np.asarray(R0) + np.asarray(R)) * nw + w
            V[i0, j] += -h; V[j, i0] += -np.conj(h); supp.add(j)
    V[i0, i0] += U
    return V, sorted(supp)
