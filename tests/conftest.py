"""
Shared fixtures of the electron_photon tests (loaded by pytest for every test file).

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture: gauge-independent properties must hold in both. Reference values and the
derivations behind the tolerances: memoire/EM/EM.md (M0, implementation notes). The real data of M1
come from the 27 x 27 wannierisation tracked in wannier/27x27/ (fixtures w90_dir, tb_w90, eig_w90).

Run with:  .venv/bin/python -m pytest tests -v
"""

from pathlib import Path

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import make_graphene_tb, make_grid_tb, make_wannier_tb
from electron_defect_interaction.io.wannier_io import read_w90_mat

N = 100    # the k grid has N x N = 1e4 points: fast, yet covers the whole Brillouin zone

HW = 2.33  # eV, 532 nm laser

W90_DIR = Path(__file__).resolve().parents[1] / "wannier" / "27x27" # repo root from this file

@pytest.fixture
def tb(request):
    """
    Graphene model, shift_B = (0,0,0) by default. Other gauges through
    @pytest.mark.parametrize("tb", [...], indirect=True), which puts each value in request.param.
    """
    shift = getattr(request, "param", (0,0,0))
    return make_graphene_tb(shift_B=shift)


@pytest.fixture
def grid(tb):
    """
    N x N GridTB of `tb`, rebuilt for every gauge of a parametrized test.
    """
    return make_grid_tb(tb, N)


@pytest.fixture(scope="session")
def w90_dir():
    """
    Directory of the tracked 27 x 27 wannierisation: wannier_tb.dat, .wout, .eig, wannier_u.mat.
    """
    return W90_DIR


@pytest.fixture(scope="session")
def tb_w90(w90_dir):
    """
    WannierTB of the real 27 x 27 data (M1), read once per session.
    """
    return make_wannier_tb(w90_dir / "wannier_tb.dat")


@pytest.fixture(scope="session")
def eig_w90(w90_dir):
    """
    DFT eigenvalues of the coarse grid: k_red (729, 3), in the order of the .eig (the k list of u.mat),
    and E (729, 20) in eV.
    """
    eig = np.loadtxt(w90_dir / "wannier.eig") # lines: band, k, energy (eV), band fastest
    nb, nk = int(eig[:, 0].max()), int(eig[:, 1].max())
    _, k_red = read_w90_mat(w90_dir / "wannier_u.mat")
    return k_red, eig[:, 2].reshape(nk, nb)
