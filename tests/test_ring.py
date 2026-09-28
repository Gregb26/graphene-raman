"""
Tests of the resonant ring (F6) and of the velocity on it: resonance, Dirac limit, symmetries,
node with Berry, isotropy, ring averages, and the effect of dropping the Berry term.
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import ring
from conftest import velocity_from_tb  # temporary: will be replaced by the src chain function

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_on_shell(tb, grid, hw):
    """
    `ring` lies on the resonance: |k - K| = q, and eps_c - eps_v = hw recomputed through
    `velocity_from_tb`, not `_gap` (observed 5e-12 eV).
    """

    k, theta, q, q0 = ring(tb, grid.K, hw)

    assert np.allclose(q, np.linalg.norm(k-grid.K[None, :], axis=1), atol=1e-12)

    eps = velocity_from_tb(tb, k)[1]
    gap = eps[:, 1] - eps[:, 0]

    assert np.allclose(gap, hw, atol=1e-10, rtol=0)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_ring_dirac_limit(tb, grid):
    """
    At hw = 0.1 eV the ring is the Dirac circle q = q0 (observed q/q0 in [0.997, 1.003]).
    """

    k, theta, q, q0 = ring(tb, grid.K, hw=0.1)

    assert np.allclose(q/q0, np.ones(len(q)), rtol=5e-3, atol=0)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_symmetry(tb, grid, hw):
    """
    q(theta) has the C3 symmetry and the mirror symmetry theta -> -theta of the bands around K.
    """

    # test C3 symmetry q(theta + 120) = q(theta)
    k, theta, q, q0 = ring(tb, grid.K, hw)
    ntheta = len(theta)

    assert np.allclose(np.roll(q, ntheta//3), q, atol=1e-12)

    # test mirror symmetry q(-theta) = q(theta): index of -theta_j is (n - j) mod n
    assert np.allclose(np.roll(q[::-1], 1), q, atol=1e-12)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_node(tb, grid, hw):
    """
    With Berry, |hbar v^x_cv|^2 vanishes at theta = 0 and pi (q // x): the node of Gruneis 2003.
    """
    k = ring(tb, grid.K, hw)[0]
    ntheta = k.shape[0]
    hv = velocity_from_tb(tb, k)[-1]
    assert np.allclose(np.abs(hv[0, 0, 1, 0])**2, 0, atol=1e-12) and np.allclose(np.abs(hv[ntheta // 2, 0, 1, 0])**2, 0, atol=1e-12)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_isotropy(tb, grid, hw):
    """
    With Berry, <|v^x_cv|^2> = <|v^y_cv|^2> and <Re(v^x_cv* v^y_cv)> = 0 on the ring (C3), which gives
    sigma_xx = sigma_yy and sigma_xy = 0 in F7. |.|^2 is taken before the mean over theta.
    """
    k = ring(tb, grid.K, hw)[0]
    hv_cv = velocity_from_tb(tb, k)[-1][:, :, 1, 0] # (ntheta, 3): element [c, v] = [1, 0]

    vx2 = np.mean(np.abs(hv_cv[:, 0])**2) # eV^2 Angstrom^2
    vy2 = np.mean(np.abs(hv_cv[:, 1])**2)
    vxy = np.mean(np.real(np.conj(hv_cv[:, 0]) * hv_cv[:, 1]))

    assert np.isclose(vx2, vy2, rtol=1e-10, atol=0), '<|v_x|^2> != <|v_y|^2> on the ring'
    assert abs(vxy) < 1e-10 * vx2, '<Re(v_x* v_y)> != 0 on the ring'


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw, ratio", [(0.1, 1.0), (2.33, 0.98909)])
def test_ring_average(tb, grid, hw, ratio):
    """
    <|hbar v^x_cv|^2> / ((hbar v_F)^2 / 2) on the ring against independent references: 1 in the
    Dirac limit (0.1 eV), 0.989 at 2.33 eV (EM.md, trigonal warping). No reference at 1.0 or 2.54 eV.
    """
    k = ring(tb, grid.K, hw)[0]
    hv = velocity_from_tb(tb, k)[-1] # (ntheta, 3, nW, nW)
    hv_F = 3*tb.t*tb.a_cc/2 # eV Angstrom

    vx2 = np.mean(np.abs(hv[:, 0, 1, 0])**2) # |.|^2 first, then the mean over theta

    assert np.isclose(vx2/(hv_F**2/2), ratio, rtol=1e-4, atol=0), 'ring average of |v_x|^2 off its reference'


def test_ring_without_berry(tb, grid):
    """
    Without Berry, at 2.33 eV with B in cell 0 (EM.md): no node at theta = 0, pi (|v^x_cv|^2 = 2.05);
    two nodes at 193.6 and 343.4 degrees (samples 193.5 and 343.5); <|v^x|^2> without / with 1.125;
    |v_cv| without / with in [0.664, 1.271]. Gauge-dependent numbers, hence the default `tb` only.
    """
    k, theta = ring(tb, grid.K, 2.33)[:2]
    ntheta = len(theta)
    hv = velocity_from_tb(tb, k, berry=True)[-1][:, :2, 1, 0] # (ntheta, 2): x, y of [c, v]
    hv_nb = velocity_from_tb(tb, k, berry=False)[-1][:, :2, 1, 0]

    vx2_nb = np.abs(hv_nb[:, 0])**2 # (ntheta,)

    assert vx2_nb[0] > 1 and vx2_nb[ntheta // 2] > 1, 'node still at theta = 0 or pi without Berry'

    # local minima: lower than both neighbours (np.roll wraps around theta = 2 pi)
    minima = (vx2_nb < np.roll(vx2_nb, 1)) & (vx2_nb < np.roll(vx2_nb, -1))
    nodes = np.degrees(theta[minima])
    assert len(nodes) == 2, f'expected 2 nodes without Berry, found {nodes} degrees'
    assert np.allclose(nodes, [193.5, 343.5], atol=0.5), f'nodes without Berry at {nodes} degrees'

    ratio_x = np.mean(vx2_nb) / np.mean(np.abs(hv[:, 0])**2)
    ratio_norm = np.linalg.norm(hv_nb, axis=1) / np.linalg.norm(hv, axis=1) # |v_cv| in the plane, per theta

    assert np.isclose(ratio_x, 1.125, rtol=1e-3, atol=0), '<|v_x|^2> without / with Berry off 1.125'
    assert np.isclose(ratio_norm.min(), 0.664, atol=1e-3), '|v_cv| without / with Berry: minimum off 0.664'
    assert np.isclose(ratio_norm.max(), 1.271, atol=1e-3), '|v_cv| without / with Berry: maximum off 1.271'
