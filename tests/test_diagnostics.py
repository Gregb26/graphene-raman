"""
Tests of the model diagnostics (diagnostics): hermiticity_report (F10) on the M0 model and on the
27 x 27 data. Reference values: memoire/EM/EM.md, F10.
"""

import dataclasses

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import hermiticity_report, make_grid_tb

BLOCKS = {'sigma': ([0, 1, 2], [0, 1, 2]),
          'pz':    ([3, 4], [3, 4]),
          'cross': ([0, 1, 2], [3, 4])}

# 27 x 27 data, grid 60^2: block -> (max_R x y, frob_R, max_k x y), Angstrom
HERMITICITY_REFERENCE = {
    'sigma': ((2.5423e-3, 1.6520e-3), 1.1977e-2, (6.9451e-3, 5.8033e-3)),
    'pz':    ((1.2222e-3, 1.3698e-3), 2.5760e-2, (1.2815e-2, 1.4190e-2)),
}


@pytest.fixture(scope="module")
def k60(tb_w90):
    """Cartesian k points of the 60 x 60 grid of the 27 x 27 data."""
    return make_grid_tb(tb_w90, 60).k_cart


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_hermiticity_report_graphene(tb, grid):
    """The M0 r (centres only) is Hermitian: no defect; frob_R undefined (nan), nothing to compare to."""
    with np.errstate(invalid='ignore'):
        rep = hermiticity_report(tb, {'pz': ([0, 1], [0, 1])}, grid.k_cart)['pz']

    assert np.array_equal(rep['max_R'], np.zeros(3)) and np.array_equal(rep['max_k'], np.zeros(3))
    assert np.isnan(rep['frob_R'])


@pytest.mark.parametrize("block", ['sigma', 'pz'])
def test_hermiticity_report_values(tb_w90, k60, block):
    """Defects of the sigma and p_z blocks (EM.md, F10); z at rounding level."""
    rep = hermiticity_report(tb_w90, BLOCKS, k60)[block]
    max_R, frob_R, max_k = HERMITICITY_REFERENCE[block]

    assert np.allclose(rep['max_R'][:2], max_R, rtol=1e-3, atol=0)
    assert np.isclose(rep['frob_R'], frob_R, rtol=1e-3, atol=0)
    assert np.allclose(rep['max_k'][:2], max_k, rtol=1e-3, atol=0)
    assert rep['max_R'][2] < 1e-13 and rep['max_k'][2] < 1e-13


def test_hermiticity_report_cross_block(tb_w90, k60):
    """The sigma-p_z block is zero by the mirror symmetry, hence so is its defect: noise only."""
    rep = hermiticity_report(tb_w90, BLOCKS, k60)['cross']

    assert np.all(rep['max_R'] < 1e-10) and np.all(rep['max_k'] < 1e-9) and rep['frob_R'] < 1e-9


def test_hermiticity_report_origin_independent(tb_w90, k60):
    """Moving the origin (all centres + c) changes nothing: the centres are left out of frob_R."""
    c = np.array([1.3, -0.7, 0.4]) # Angstrom
    i0 = tb_w90.index[(0, 0, 0)]
    nW = tb_w90.r_R.shape[-1]
    r_R = tb_w90.r_R.copy()
    r_R[i0] += c[:, None, None] * np.eye(nW) # tau_n -> tau_n + c
    tb_moved = dataclasses.replace(tb_w90, r_R=r_R)

    rep = hermiticity_report(tb_w90, BLOCKS, k60)
    rep_moved = hermiticity_report(tb_moved, BLOCKS, k60)

    for block in ['sigma', 'pz']:
        for key in ['max_R', 'frob_R', 'max_k']:
            assert np.allclose(rep_moved[block][key], rep[block][key], rtol=1e-10, atol=1e-15)
