"""dirac_point (2026-09-30): replaces the inline E_D of six scripts (midpoint of the smallest pi/pi* gap, bands 3 and 4)."""
import numpy as np

from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk, dirac_point
from graphene_raman.defects.many_body import cluster_tmatrix as ct
from graphene_raman.defects.many_body.tb_models import graphene_pz_tb


def test_dirac_point_matches_inline_expression():
    """Same value, bit for bit, as the expression it replaces in the scripts."""
    E = np.sort(np.random.default_rng(3).normal(size=(500, 5)), axis=1)
    gap = E[:, 4] - E[:, 3]; iD = int(np.argmin(gap))
    assert dirac_point(E) == (float(0.5 * (E[iD, 3] + E[iD, 4])), float(gap[iD]))


def test_dirac_point_graphene_model():
    """Nearest-neighbour model with e_pz = 0.3 eV: E_D = e_pz exactly, gap 0 at K (on the 90 x 90 grid since 90 = 3 x 30)."""
    Hwr, Rw, nd = graphene_pz_tb(e_pz=0.3)
    E = Hwr_to_Hwk(Hwr, Rw, ct.mp_grid(90), ndegen=nd)[1]
    E_D, gap = dirac_point(E)
    assert abs(E_D - 0.3) < 1e-12 and abs(gap) < 1e-12


def test_dirac_point_other_bands():
    """bands= selects the pair (here a two-band array, as for interpolated pi / pi* on a path)."""
    E = np.array([[-1.0, 2.0], [-0.2, 0.4], [-0.5, 0.9]])
    assert dirac_point(E, bands=(0, 1)) == (0.1, 0.6000000000000001)
