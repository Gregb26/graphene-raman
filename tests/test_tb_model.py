"""
Tests of the tight-binding model (tb_model): R-space Hermiticity, bonds, relabelling shift_B.
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import make_graphene_tb

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_hermitian(tb, grid):
    """
    R-space Hermiticity X_ij(R) = X_ji(-R)* for H and r (dagger on one side only), R list closed
    under R -> -R, and ndegen(-R) = ndegen(R): together they make H(k) and A(k) Hermitian.
    """

    assert np.allclose(tb.minus[tb.minus], np.arange(len(tb.R_int))), 'R list not closed under R -> -R'
    assert np.allclose(tb.H_R, tb.H_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'TB Hamiltonian not hermitian'
    assert np.allclose(tb.r_R, tb.r_R[tb.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'positian operator not hermitian'
    assert np.allclose(tb.ndegen, tb.ndegen[tb.minus], atol=1e-12), 'R and -R dont have the same weight'


def bond_vectors(tb):
    """
    Bond vectors delta = R_cart + r_BB(0) - r_AA(0) for every R with H_AB(R) != 0, rebuilt from the
    WannierTB alone (a shift L cancels between R and the B centre). Returns (n_hops, 3), Angstrom.
    """
    R_cart = tb.R_int @ tb.lattice.T
    index = {tuple(R): iR for iR, R in enumerate(tb.R_int)}

    iR0 = index[(0,0,0)]
    r_0AA = tb.r_R[iR0,:,0,0].real; r_0BB = tb.r_R[iR0,:,1,1].real # centres are real vectors
    # A -> B hops: rows where the AB element is non-zero
    mask = np.abs(tb.H_R[:,0,1]) > 1e-12

    deltas = R_cart[mask] + r_0BB - r_0AA

    return deltas


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_tb_bond_length(tb):
    """
    Three A -> B hops of length a_cc, in both gauges.
    """
    deltas = bond_vectors(tb)

    assert len(deltas) == 3, 'must have three non-zero hopping for first neighbhour graphene tb'
    assert np.allclose(np.linalg.norm(deltas, axis=1), tb.a_cc, atol=1e-12), 'carbon-carbon bond length not correct'


def sort_by_angle(d):
    """
    Rows of d sorted by polar angle arctan2(y, x), to compare lists of vectors given in any order.
    """

    angles = np.arctan2(d[:,1], d[:,0])
    order = np.argsort(angles)
    return d[order]


def test_tb_shift_same_bonds():
    """
    shift_B = (1,0,0) changes the R list (nR = 5 -> 7) but not the bonds. Builds both models itself.
    """
    tb = make_graphene_tb()
    tb_shifted = make_graphene_tb(shift_B=(1,0,0))

    assert len(tb.R_int) == 5
    assert len(tb_shifted.R_int) == 7

    deltas = bond_vectors(tb)
    deltas_shifted = bond_vectors(tb_shifted)

    assert np.allclose(sort_by_angle(deltas), sort_by_angle(deltas_shifted), atol=1e-12)
