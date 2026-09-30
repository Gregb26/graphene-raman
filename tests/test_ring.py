"""
Tests of the resonant ring (F6) and of the velocity on it: resonance, Dirac limit, symmetries,
node with Berry, isotropy, ring averages, and the effect of dropping the Berry term; fermi_velocity
(F13), ring_stats (F15), ring_kpoints_crystal (F17) and map_around_K (F19) on the M0 model and on the
27 x 27 data.
"""

import pytest
import numpy as np
from graphene_raman.electron_photon import (centres_only, compute_velocity, fermi_velocity,
                                                         make_graphene_tb, make_grid_tb, map_around_K,
                                                         reciprocal, ring, ring_kpoints_crystal, ring_stats)
from scipy.interpolate import RegularGridInterpolator

@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_on_shell(tb, grid, hw):
    """
    `ring` lies on the resonance: |k - K| = q, and eps_c - eps_v = hw recomputed through
    `compute_velocity`, not `_gap` (observed 5e-12 eV).
    """

    k, theta, q, q0 = ring(tb, grid.K, hw)

    assert np.allclose(q, np.linalg.norm(k-grid.K[None, :], axis=1), atol=1e-12)

    eps = compute_velocity(tb, k)[1]
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
    hv = compute_velocity(tb, k)[-1]
    assert np.allclose(np.abs(hv[0, 0, 1, 0])**2, 0, atol=1e-12) and np.allclose(np.abs(hv[ntheta // 2, 0, 1, 0])**2, 0, atol=1e-12)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
@pytest.mark.parametrize("hw", [0.1, 1.0, 2.33, 2.54])
def test_ring_isotropy(tb, grid, hw):
    """
    With Berry, <|v^x_cv|^2> = <|v^y_cv|^2> and <Re(v^x_cv* v^y_cv)> = 0 on the ring (C3), which gives
    sigma_xx = sigma_yy and sigma_xy = 0 in F7. |.|^2 is taken before the mean over theta.
    """
    k = ring(tb, grid.K, hw)[0]
    hv_cv = compute_velocity(tb, k)[-1][:, :, 1, 0] # (ntheta, 3): element [c, v] = [1, 0]

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
    hv = compute_velocity(tb, k)[-1] # (ntheta, 3, nW, nW)
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
    hv = compute_velocity(tb, k, mode='berry')[-1][:, :2, 1, 0] # (ntheta, 2): x, y of [c, v]
    hv_nb = compute_velocity(tb, k, mode='no_berry')[-1][:, :2, 1, 0]

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


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_fermi_velocity_graphene(tb, grid):
    """M0: both ways give 3 t a_cc / 2, and pi = pi* (no electron-hole asymmetry), in both gauges."""
    res = fermi_velocity(tb, grid.K, 0.0)
    hv_F = 3 * tb.t * tb.a_cc / 2

    for key in ['intra_avg', 'inter_avg', 'pi', 'pi_star']:
        assert np.isclose(res[key], hv_F, rtol=1e-5, atol=0), key


@pytest.fixture(scope="module")
def K_w90(tb_w90):
    """Dirac point of the real lattice, 1/Angstrom."""
    return make_grid_tb(tb_w90, 3).K


def test_fermi_velocity_real(tb_w90, K_w90, w90_ref):
    """Real data: the two ways agree on hbar v_F (EM.md, F13)."""
    res = fermi_velocity(tb_w90, K_w90, w90_ref.E_D)

    assert np.isclose(res['intra_avg'], w90_ref.hv_F, rtol=1e-5, atol=0)
    assert np.isclose(res['inter_avg'], w90_ref.hv_F, rtol=1e-5, atol=0)
    assert np.allclose([res['pi'], res['pi_star']], w90_ref.hv_pi, rtol=1e-5, atol=0)


def test_fermi_velocity_electron_hole(tb_w90, K_w90, w90_ref):
    """
    The pi / pi* asymmetry grows linearly with q (x10 from q = 1e-3 to 1e-2), while their average,
    hbar v_F, stays put (1e-5).
    """
    small = fermi_velocity(tb_w90, K_w90, w90_ref.E_D, q=1e-3)
    large = fermi_velocity(tb_w90, K_w90, w90_ref.E_D, q=1e-2)
    asym = lambda res: res['pi'] - res['pi_star']

    assert np.isclose(asym(large) / asym(small), 10, rtol=0.01, atol=0)
    assert np.isclose(large['intra_avg'], small['intra_avg'], rtol=1e-4, atol=0)


def test_fermi_velocity_modes(tb_w90, K_w90, w90_ref):
    """The intraband way ignores the Berry term (zero on the diagonal); the interband way barely feels it at K."""
    full = fermi_velocity(tb_w90, K_w90, w90_ref.E_D, mode='berry')
    others = [fermi_velocity(tb_w90, K_w90, w90_ref.E_D, mode='no_berry'),
              fermi_velocity(centres_only(tb_w90), K_w90, w90_ref.E_D, mode='berry')]

    for res in others:
        assert np.isclose(res['intra_avg'], full['intra_avg'], rtol=1e-12, atol=0)
        assert np.isclose(res['inter_avg'], full['inter_avg'], rtol=1e-4, atol=0)


def test_fermi_velocity_mu_between_bands(tb_w90, K_w90, w90_ref):
    """
    mu must separate pi and pi* all around the circle: at q = 0.3, pi* spans E_D + 1.2 to 1.9 eV
    (trigonal warping), so mu = E_D + 1.5 eV cuts it and the occupation check fails.
    """
    with pytest.raises(AssertionError):
        fermi_velocity(tb_w90, K_w90, w90_ref.E_D + 1.5, q=0.3)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_ring_stats_graphene(tb, grid):
    """
    M0 at 2.33 eV: the averages of test_ring_average (0.98909, x = y), node of e_x pinned at q // x,
    that of e_y at -4 degrees (trigonal warping); centres only = full (M0 has only centres).
    """
    stats = ring_stats(tb, grid.K, 0.0, 2.33)[2.33]
    full = stats['full']

    assert np.allclose(full['avg'], 0.98909, rtol=1e-4, atol=0)
    assert np.allclose(full['node'], [0.0, -4.0], rtol=0, atol=1e-9) and full['ratio'] == (1.0, 1.0)
    for key in ['avg', 'node', 'ratio']:
        assert np.allclose(stats['centres_only'][key], full[key], rtol=1e-12, atol=0)


def test_ring_stats_gauge():
    """
    Moving the B centre (shift_B) leaves the full velocity unchanged and changes the one without
    Berry; in the default gauge, the no-Berry numbers are those of test_ring_without_berry.
    """
    tb, tb_shifted = make_graphene_tb(), make_graphene_tb(shift_B=(1,0,0))
    K = make_grid_tb(tb, 3).K
    stats, stats_shifted = ring_stats(tb, K, 0.0, 2.33)[2.33], ring_stats(tb_shifted, K, 0.0, 2.33)[2.33]

    for key in ['avg', 'node', 'ratio']:
        assert np.allclose(stats_shifted['full'][key], stats['full'][key], rtol=1e-10, atol=0)
    assert not np.allclose(stats_shifted['no_berry']['avg'], stats['no_berry']['avg'], rtol=1e-2, atol=0)

    no_berry = stats['no_berry']
    assert np.isclose(no_berry['node'][0], -16.5, rtol=0, atol=1e-9) # the node at 343.5 degrees
    assert np.allclose(no_berry['ratio'], (0.664, 1.271), rtol=0, atol=1e-3)
    assert np.isclose(no_berry['avg'][0] / stats['full']['avg'][0], 1.125, rtol=1e-3, atol=0)


@pytest.fixture(scope="module")
def stats_w90(tb_w90, K_w90, w90_ref):
    """ring_stats of the 27 x 27 data at the three laser energies (1.96, 2.33, 2.54 eV)."""
    return ring_stats(tb_w90, K_w90, w90_ref.E_D, list(w90_ref.ring_stats))


def test_ring_stats_values(stats_w90, w90_ref):
    """The table of EM.md (F15): three lasers x three variants."""
    assert set(stats_w90) == set(w90_ref.ring_stats)
    for hw, ref in w90_ref.ring_stats.items():
        for variant, (avg, node, ratio) in ref.items():
            res = stats_w90[hw][variant]
            assert np.allclose(res['avg'], avg, rtol=1e-5, atol=0), (hw, variant)
            assert np.allclose(res['node'], node, rtol=0, atol=1e-9), (hw, variant)
            assert np.allclose(res['ratio'], ratio, rtol=1e-5, atol=0), (hw, variant)


def test_ring_stats_symmetries(stats_w90):
    """
    Node of e_y on the mirror line through K (0) in every variant; x = y averages (C3) with the Berry
    term, full or centres only, not without it; the Berry term moves the node of e_x across q // x.
    Without Berry, <y> = <y> centres only: tau_B - tau_A is along x, so the centres give no y term.
    """
    for hw, stats in stats_w90.items():
        full, centres, no_berry = stats['full'], stats['centres_only'], stats['no_berry']

        assert full['ratio'] == (1.0, 1.0)
        assert all(abs(stats[v]['node'][1]) < 1e-9 for v in stats)
        for res in (full, centres):
            assert np.isclose(res['avg'][0], res['avg'][1], rtol=5e-4, atol=0)
        assert no_berry['avg'][0] / no_berry['avg'][1] > 1.1
        assert np.isclose(no_berry['avg'][1], centres['avg'][1], rtol=1e-12, atol=0)

        assert full['node'][0] > 5 and no_berry['node'][0] < -5
        assert abs(centres['node'][0] - full['node'][0]) < 0.5 + 1e-9 # one grid step


def test_ring_stats_average_direct(tb_w90, K_w90, stats_w90, w90_ref):
    """(avg_x + avg_y) (hbar v_F)^2 / 2 = <|hbar v_cv|^2> computed directly on the 2.33 eV ring."""
    v, c = w90_ref.bands_pi
    k = ring(tb_w90, K_w90, 2.33, mu=w90_ref.E_D)[0]
    hv = compute_velocity(tb_w90, k, mode='berry')[-1]
    direct = np.mean(np.abs(hv[:, 0, c, v])**2 + np.abs(hv[:, 1, c, v])**2) # eV^2 Angstrom^2

    assert np.isclose(np.sum(stats_w90[2.33]['full']['avg']) * w90_ref.hv_F**2 / 2, direct, rtol=1e-4, atol=0)


@pytest.mark.parametrize("variant, mode", [('full', 'berry'), ('no_berry', 'no_berry')])
def test_ring_stats_node_is_minimum(tb_w90, K_w90, stats_w90, w90_ref, variant, mode):
    """
    The node of e_x is the local minimum of |hbar v^x_cv|^2 found independently in the half ring
    facing +x, and a true zero (below 1e-4 of the maximum).
    """
    v, c = w90_ref.bands_pi
    k, theta = ring(tb_w90, K_w90, 2.33, mu=w90_ref.E_D)[:2]
    vx2 = np.abs(compute_velocity(tb_w90, k, mode=mode)[-1][:, 0, c, v])**2 # (ntheta,)

    minima = (vx2 < np.roll(vx2, 1)) & (vx2 < np.roll(vx2, -1))
    angles = np.degrees(np.angle(np.exp(1j * theta[minima]))) # in (-180, 180]
    node = stats_w90[2.33][variant]['node'][0]

    assert np.isclose(angles[np.abs(angles) < 90], node, rtol=0, atol=1e-9).sum() == 1
    i = np.argmin(np.abs(np.degrees(np.angle(np.exp(1j * theta))) - node))
    assert vx2[i] < 1e-4 * vx2.max()


def test_ring_stats_scalar_hw(tb_w90, K_w90, stats_w90, w90_ref):
    """A single energy (float) is accepted and gives the same entry as in the list."""
    stats = ring_stats(tb_w90, K_w90, w90_ref.E_D, 2.33)

    assert list(stats) == [2.33]
    for key in ['avg', 'node', 'ratio']:
        assert np.allclose(stats[2.33]['no_berry'][key], stats_w90[2.33]['no_berry'][key], rtol=1e-12, atol=0)


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_ring_kpoints_crystal_graphene(tb, grid):
    """
    M0: back to Cartesian (k_red @ B.T) gives the points of `ring`, same theta; the ring is centred
    on K = (2/3, 1/3, 0) in crystal coordinates (C3, exact in M0).
    """
    k_red, theta = ring_kpoints_crystal(tb, grid.K, 0.0, {2.33: 12})[2.33]
    k, theta_ring = ring(tb, grid.K, 2.33, ntheta=12)[:2]

    assert k_red.shape == (12, 3) and np.array_equal(theta, theta_ring)
    assert np.allclose(k_red @ reciprocal(tb.lattice).T, k, rtol=0, atol=1e-14)
    assert np.allclose(k_red.mean(axis=0), [2/3, 1/3, 0], rtol=0, atol=1e-12)


@pytest.fixture(scope="module")
def kpoints_w90(tb_w90, K_w90, w90_ref):
    """ring_kpoints_crystal of the 27 x 27 data: the EM2 list, 48 points at 2.33 eV and 12 at 1.96, 2.54."""
    return ring_kpoints_crystal(tb_w90, K_w90, w90_ref.E_D, {2.33: 48, 1.96: 12, 2.54: 12})


def test_ring_kpoints_crystal_real(tb_w90, K_w90, kpoints_w90, w90_ref):
    """
    Structure of the EM2 list; first point at 2.33 eV; the ring points of `ring` converted by an
    independent inverse (solve with B); rings centred on K = (2/3, 1/3, 0) to 2e-7 (C3 of the data).
    """
    B = reciprocal(tb_w90.lattice)

    assert list(kpoints_w90) == [2.33, 1.96, 2.54]
    assert [len(k_red) for k_red, _ in kpoints_w90.values()] == [48, 12, 12]
    assert np.allclose(kpoints_w90[2.33][0][0], [0.73959906, 0.40626573, 0], rtol=0, atol=1e-8)

    for hw, (k_red, theta) in kpoints_w90.items():
        k, theta_ring = ring(tb_w90, K_w90, hw, mu=w90_ref.E_D, ntheta=len(k_red))[:2]
        assert np.array_equal(theta, theta_ring) and np.all(k_red[:, 2] == 0)
        assert np.allclose(k_red, np.linalg.solve(B, k.T).T, rtol=0, atol=1e-14)
        assert np.allclose(k_red.mean(axis=0), [2/3, 1/3, 0], rtol=0, atol=1e-6)


def test_ring_kpoints_crystal_on_shell(tb_w90, kpoints_w90, w90_ref):
    """Each point, back in Cartesian coordinates, absorbs exactly its hw: eps_c - eps_v = hw (5e-12 eV)."""
    v, c = w90_ref.bands_pi
    B = reciprocal(tb_w90.lattice)

    for hw, (k_red, _) in kpoints_w90.items():
        _, eps, _, _ = compute_velocity(tb_w90, k_red @ B.T, mode='no_berry')
        assert np.allclose(eps[:, c] - eps[:, v], hw, rtol=0, atol=1e-10), hw


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_map_around_K_dirac_limit(tb, grid):
    """
    M0 close to K (half width 0.01): a cone, Delta eps = 2 hbar v_F |q| and |hbar v^x_cv|^2 +
    |hbar v^y_cv|^2 = (hbar v_F)^2 (to 4e-3 and 1.7e-2, trigonal corrections), in both gauges.
    """
    deps, P = map_around_K(tb, grid.K, 0.0, 0.01, 20)
    q = np.linspace(-0.01, 0.01, 20)
    qx, qy = np.meshgrid(q, q, indexing='xy')
    hv_F = 3*tb.t*tb.a_cc/2

    assert deps.shape == (20, 20) and P.shape == (20, 20, 2)
    assert np.allclose(deps, 2*hv_F*np.hypot(qx, qy), rtol=1e-2, atol=0)
    assert np.allclose(P.sum(axis=-1), hv_F**2, rtol=3e-2, atol=0)


def test_map_around_K_gauge():
    """Moving the B centre (shift_B) leaves the full map unchanged."""
    maps = []
    for shift in [(0,0,0), (1,0,0)]:
        tb = make_graphene_tb(shift_B=shift)
        maps.append(map_around_K(tb, make_grid_tb(tb, 3).K, 0.0, 0.3, 20))

    assert np.allclose(maps[1][0], maps[0][0], rtol=0, atol=1e-12)
    assert np.allclose(maps[1][1], maps[0][1], rtol=0, atol=1e-9)


@pytest.fixture(scope="module")
def map_w90(tb_w90, K_w90, w90_ref):
    """map_around_K of the 27 x 27 data, half width 0.35 (beyond the 2.54 eV ring), nq = 100."""
    return map_around_K(tb_w90, K_w90, w90_ref.E_D, 0.35, 100)


def test_map_around_K_values(map_w90, w90_ref):
    """Reference values: range of Delta eps and maxima of |hbar v^{x,y}_cv|^2."""
    deps, P = map_w90
    (de_min, de_max), P_max = w90_ref.map_K

    assert np.allclose([deps.min(), deps.max()], [de_min, de_max], rtol=1e-5, atol=0)
    assert np.allclose(P.max(axis=(0, 1)), P_max, rtol=1e-5, atol=0)


def test_map_around_K_orientation(tb_w90, K_w90, map_w90, w90_ref):
    """
    Element [iy, ix] is the point K + (q[ix], q[iy]): checked against compute_velocity at a point off
    the diagonal (q ~ (0.2, 0)), where x and y differ (P = (0.49, 28.2) there, (47.0, 0.0007) transposed).
    """
    deps, P = map_w90
    v, c = w90_ref.bands_pi
    q = np.linspace(-0.35, 0.35, 100)
    ix, iy = np.argmin(np.abs(q - 0.2)), np.argmin(np.abs(q))

    _, eps, _, hv = compute_velocity(tb_w90, K_w90 + np.array([[q[ix], q[iy], 0]]), mode='berry')

    assert np.isclose(deps[iy, ix], eps[0, c] - eps[0, v], rtol=1e-12, atol=0)
    assert np.allclose(P[iy, ix], np.abs(hv[0, :2, c, v])**2, rtol=1e-10, atol=0)
    assert P[iy, ix, 0] < 0.05 * P[iy, ix, 1] # close to the dark line of e_x


def test_map_around_K_rings(tb_w90, K_w90, map_w90, w90_ref):
    """The iso-lines Delta eps = hw of the map are the rings of `ring`: cubic interpolation at the ring points gives hw (5e-6 eV)."""
    deps, _ = map_w90
    q = np.linspace(-0.35, 0.35, 100)
    interp = RegularGridInterpolator((q, q), deps.T, method='cubic') # deps.T is indexed [ix, iy]

    for hw in (1.96, 2.33, 2.54):
        dk = ring(tb_w90, K_w90, hw, mu=w90_ref.E_D)[0] - K_w90
        assert np.allclose(interp(dk[:, :2]), hw, rtol=0, atol=2e-5), hw


def test_map_around_K_modes(tb_w90, K_w90, map_w90, w90_ref):
    """Without Berry: same Delta eps (eigenvalues), different |hbar v_cv|^2; odd nq refused (K on the grid)."""
    deps, P = map_w90
    deps_nb, P_nb = map_around_K(tb_w90, K_w90, w90_ref.E_D, 0.35, 100, mode='no_berry')

    assert np.array_equal(deps_nb, deps)
    assert np.abs(P_nb - P).max() > 10
    with pytest.raises(AssertionError, match='even'): # not the later occupation check, which also fails at K
        map_around_K(tb_w90, K_w90, w90_ref.E_D, 0.35, 101)
