"""
Tests of the reciprocal lattice and k grids (kgrid).
"""

import numpy as np

def test_reciprocal(tb, grid):
    """
    a_i . b_j = 2 pi delta_ij, written as lattice @ B.T = 2 pi I.
    """

    assert np.allclose(tb.lattice @ grid.B.T, 2*np.pi*np.eye(3), atol=1e-12)
