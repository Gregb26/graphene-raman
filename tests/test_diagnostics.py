"""
Tests of the model diagnostics (diagnostics): hermiticity_report (F10) on the M0 model and on the
27 x 27 data, symmetry_report (F11) on the 27 x 27 data, frozen_window_limit (F18) on both (reference
values in w90_ref, conftest.py; see campagnes/EM/EM.md, F10, F11 and F18).
"""

import dataclasses

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import (compute_velocity, frozen_window_limit, hermiticity_report,
                                                         make_graphene_tb, make_grid_tb, pz_block, symmetry_report)


@pytest.fixture(scope="module")
def k60(tb_w90):
    """Cartesian k points of the 60 x 60 grid of the 27 x 27 data."""
    return make_grid_tb(tb_w90, 60).k_cart


@pytest.fixture(scope="module")
def blocks(w90_ref):
    """sigma, p_z and cross blocks of the reference Wannier order."""
    sigma, pz = w90_ref.sigma, w90_ref.pz
    return {'sigma': (sigma, sigma), 'pz': (pz, pz), 'cross': (sigma, pz)}


@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)
def test_hermiticity_report_graphene(tb, grid):
    """The M0 r (centres only) is Hermitian: no defect; frob_R undefined (nan), nothing to compare to."""
    with np.errstate(invalid='ignore'):
        rep = hermiticity_report(tb, {'pz': ([0, 1], [0, 1])}, grid.k_cart)['pz']

    assert np.array_equal(rep['max_R'], np.zeros(3)) and np.array_equal(rep['max_k'], np.zeros(3))
    assert np.isnan(rep['frob_R'])


@pytest.mark.parametrize("block", ['sigma', 'pz'])
def test_hermiticity_report_values(tb_w90, k60, blocks, w90_ref, block):
    """Defects of the sigma and p_z blocks (EM.md, F10); z at rounding level."""
    rep = hermiticity_report(tb_w90, blocks, k60)[block]
    max_R, frob_R, max_k = w90_ref.hermiticity[block]

    assert np.allclose(rep['max_R'][:2], max_R, rtol=1e-3, atol=0)
    assert np.isclose(rep['frob_R'], frob_R, rtol=1e-3, atol=0)
    assert np.allclose(rep['max_k'][:2], max_k, rtol=1e-3, atol=0)
    assert rep['max_R'][2] < 1e-13 and rep['max_k'][2] < 1e-13


def test_hermiticity_report_cross_block(tb_w90, k60, blocks):
    """The sigma-p_z block is zero by the mirror symmetry, hence so is its defect: noise only."""
    rep = hermiticity_report(tb_w90, blocks, k60)['cross']

    assert np.all(rep['max_R'] < 1e-10) and np.all(rep['max_k'] < 1e-9) and rep['frob_R'] < 1e-9


def test_hermiticity_report_origin_independent(tb_w90, k60, blocks):
    """Moving the origin (all centres + c) changes nothing: the centres are left out of frob_R."""
    c = np.array([1.3, -0.7, 0.4]) # Angstrom
    i0 = tb_w90.index[(0, 0, 0)]
    nW = tb_w90.r_R.shape[-1]
    r_R = tb_w90.r_R.copy()
    r_R[i0] += c[:, None, None] * np.eye(nW) # tau_n -> tau_n + c
    tb_moved = dataclasses.replace(tb_w90, r_R=r_R)

    rep = hermiticity_report(tb_w90, blocks, k60)
    rep_moved = hermiticity_report(tb_moved, blocks, k60)

    for block in ['sigma', 'pz']:
        for key in ['max_R', 'frob_R', 'max_k']:
            assert np.allclose(rep_moved[block][key], rep[block][key], rtol=1e-10, atol=1e-15)


def test_symmetry_report_values(tb_w90, w90_ref):
    """Mirror selection rules hold to noise; only r^z between sigma and p_z survives."""
    rep = symmetry_report(tb_w90, w90_ref.sigma, w90_ref.pz)

    assert rep['H_mixed'] < 1e-8 # eV, against |H| up to 15 eV
    assert np.all(rep['r_mixed'][:2] < 1e-10)
    assert np.isclose(rep['r_mixed'][2], w90_ref.rz_sigma_pz, rtol=1e-5, atol=0)
    assert rep['rz_same'] < 1e-10


def perturbed(tb, X, iR, comp, m, n, value):
    """Copy of tb with X_mn(R) += value at row iR (and the Hermitian partner), X = 'H_R' or 'r_R'."""
    arr = getattr(tb, X).copy()
    block = arr[iR] if X == 'H_R' else arr[iR, comp]
    partner = arr[tb.minus[iR]] if X == 'H_R' else arr[tb.minus[iR], comp]
    block[m, n] += value
    partner[n, m] += np.conj(value)
    return dataclasses.replace(tb, **{X: arr})


@pytest.mark.parametrize("X, comp, parity, key", [
    ('H_R', None, 'mixed', 'H_mixed'), # H between sigma and p_z
    ('r_R', 0,    'mixed', 'r_mixed'), # x between sigma and p_z
    ('r_R', 2,    'same',  'rz_same'), # z between the two p_z
], ids=["H_sigma_pz", "x_sigma_pz", "z_pz_pz"])
def test_symmetry_report_detects_breaking(tb_w90, w90_ref, X, comp, parity, key):
    """A forbidden element added by hand (1e-3) shows up in the report."""
    sigma, pz = w90_ref.sigma, w90_ref.pz
    m, n = (sigma[0], pz[0]) if parity == 'mixed' else (pz[0], pz[1])
    iR = tb_w90.index[(1, 0, 0)]

    rep = symmetry_report(perturbed(tb_w90, X, iR, comp, m, n, 1e-3), sigma, pz)
    value = rep[key][comp] if key == 'r_mixed' else rep[key]

    assert np.isclose(value, 1e-3, rtol=1e-6, atol=0)


def test_symmetry_report_sheet_height(tb_w90, w90_ref):
    """Moving the sheet to z0 (every centre + z0 e_z) changes nothing: rz_same excludes the centres."""
    z0 = 7.9 # Angstrom, middle of the cell
    i0 = tb_w90.index[(0, 0, 0)]
    r_R = tb_w90.r_R.copy()
    r_R[i0, 2] += z0 * np.eye(r_R.shape[-1])
    tb_moved = dataclasses.replace(tb_w90, r_R=r_R)

    rep = symmetry_report(tb_w90, w90_ref.sigma, w90_ref.pz)
    rep_moved = symmetry_report(tb_moved, w90_ref.sigma, w90_ref.pz)

    for key in ['H_mixed', 'r_mixed', 'rz_same']:
        assert np.allclose(rep_moved[key], rep[key], rtol=1e-10, atol=1e-14)


@pytest.mark.parametrize("F", [1.0, 2.0])
def test_frozen_window_limit_graphene(F):
    """
    M0 is electron-hole symmetric (eps_c = -eps_v): the limit is 2 F, reached from above as N grows
    (excess 3e-3 eV at N = 300, 6e-4 at N = 600 for F = 1 eV).
    """
    tb = make_graphene_tb()
    coarse, fine = frozen_window_limit(tb, 300, F, 0.0), frozen_window_limit(tb, 600, F, 0.0)

    assert 2*F < fine < coarse < 2*F + 0.01


def test_frozen_window_limit_values(tb_w90, w90_ref):
    """27 x 27 data at froz_max: 4.985 eV (N = 100), 4.979 eV (N = 200), well above the three lasers."""
    for N, hw in w90_ref.hw_froz.items():
        assert np.isclose(frozen_window_limit(tb_w90, N, w90_ref.froz_max, w90_ref.E_D), hw, rtol=1e-6, atol=0)


def test_frozen_window_limit_independent(tb_w90, w90_ref):
    """
    Same value from the p_z model with all k at once (pi, pi* = bands 0, 1, no selection by energy),
    and independent of the block size.
    """
    tb2 = pz_block(tb_w90, w90_ref.pz)
    _, eps, _, _ = compute_velocity(tb2, make_grid_tb(tb2, 100).k_cart, mode='no_berry')
    above = eps[:, 1] > w90_ref.froz_max
    direct = np.min(eps[above, 1] - eps[above, 0])

    assert np.isclose(frozen_window_limit(tb_w90, 100, w90_ref.froz_max, w90_ref.E_D), direct, rtol=1e-12, atol=0)
    assert np.isclose(frozen_window_limit(tb_w90, 100, w90_ref.froz_max, w90_ref.E_D, chunk=100), direct, rtol=1e-12, atol=0)


def test_frozen_window_limit_monotonic(tb_w90, w90_ref):
    """A higher window keeps fewer k, so the limit rises; above every band it is inf."""
    limits = [frozen_window_limit(tb_w90, 100, w90_ref.E_D + d, w90_ref.E_D) for d in (1.5, 2.0, 2.5, 3.0)]

    assert np.all(np.diff(limits) > 0.5)
    assert frozen_window_limit(tb_w90, 30, 100.0, w90_ref.E_D) == np.inf
