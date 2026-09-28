"""
Shared fixtures of the electron_photon tests (loaded by pytest for every test file).

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture: gauge-independent properties must hold in both. Reference values and the
derivations behind the tolerances: memoire/EM/EM.md (M0, implementation notes).

Run with:  .venv/bin/python -m pytest tests -v
"""

import pytest
from electron_defect_interaction.electron_photon import make_graphene_tb, make_grid_tb

N = 100    # the k grid has N x N = 1e4 points: fast, yet covers the whole Brillouin zone

HW = 2.33  # eV, 532 nm laser

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
