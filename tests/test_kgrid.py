"""
Tests of the reciprocal lattice, k grids and paths (kgrid).
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import kpath, reciprocal

PATH = [('G', (0, 0, 0)), ('K', (2/3, 1/3, 0)), ('M', (1/2, 0, 0)), ('G', (0, 0, 0))]

def test_reciprocal(tb, grid):
    """
    a_i . b_j = 2 pi delta_ij, written as lattice @ B.T = 2 pi I.
    """

    assert np.allclose(tb.lattice @ grid.B.T, 2*np.pi*np.eye(3), atol=1e-12)


def test_kpath_lengths(tb):
    """|GK| = 4 pi/3a, |KM| = 2 pi/3a, |MG| = 2 pi/(sqrt(3) a), with a = sqrt(3) a_cc."""
    a = np.sqrt(3) * tb.a_cc
    _, _, ticks, labels = kpath(reciprocal(tb.lattice), PATH, 100)

    assert labels == ['G', 'K', 'M', 'G']
    expected = [4*np.pi/(3*a), 2*np.pi/(3*a), 2*np.pi/(np.sqrt(3)*a)]
    assert np.allclose(np.diff(ticks), expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize("nk", [50, 100, 300])
def test_kpath_steps(tb, nk):
    """
    x is the distance travelled: |k_{j+1} - k_j| = x_{j+1} - x_j at every step, no repeated point,
    near-uniform step, about nk points.
    """
    k, x, ticks, _ = kpath(reciprocal(tb.lattice), PATH, nk)
    dx = np.diff(x)

    assert np.allclose(np.linalg.norm(np.diff(k, axis=0), axis=1), dx, rtol=0, atol=1e-12)
    assert np.all(dx > 0)
    assert dx.max() / dx.min() < 1.1
    assert abs(len(k) - nk) <= len(PATH) - 1 # rounding, at most one point per segment
    assert x[0] == 0 and np.isclose(x[-1], ticks[-1], rtol=0, atol=1e-12)


def test_kpath_vertices(tb):
    """Each vertex is sampled exactly, at x = its tick."""
    B = reciprocal(tb.lattice)
    k, x, ticks, _ = kpath(B, PATH, 100)

    for (_, k_red), tick in zip(PATH, ticks):
        j = np.flatnonzero(np.isclose(x, tick, rtol=0, atol=1e-12))
        assert len(j) == 1
        assert np.allclose(k[j[0]], np.array(k_red) @ B.T, rtol=0, atol=1e-12)
