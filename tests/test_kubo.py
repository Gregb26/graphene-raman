"""
Tests of the Kubo conductivity (F7) and of the driver `sigma_on_grid`, down to the reference table
of EM.md on the 1800^2 grid (~13 s; skip with -k "not kubo_reference").
"""

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import gaussian_eta, kubo_accumulate, make_graphene_tb, sigma_on_grid

# EM.md table, 1800^2 shifted grid, eta = 0.04 eV: (shift_B, mode) -> {hw: (sigma_xx, sigma_yy, sigma_xy)/sigma_0}
KUBO_REFERENCE = {
    ((0,0,0), 'berry'):    {1.0: (1.0156, 1.0156, 0.0),    2.33: (1.0938, 1.0938, 0.0)},
    ((1,0,0), 'berry'):    {1.0: (1.0156, 1.0156, 0.0),    2.33: (1.0938, 1.0938, 0.0)},
    ((0,0,0), 'no_berry'): {1.0: (1.0389, 1.0234, 0.0134), 2.33: (1.2266, 1.1381, 0.0767)},
    ((1,0,0), 'no_berry'): {1.0: (1.2248, 1.0234, 0.0403), 2.33: (2.2896, 1.1381, 0.2301)},
}

def test_gaussian_eta():
    """
    g_eta is a normalized Gaussian whose standard deviation is eta (not the FWHM = 2.355 eta):
    integral 1, second moment eta^2, peak 1/(eta sqrt(2 pi)).
    """
    eta = 0.04
    x = np.linspace(-1, 1, 200001) # eV, +-25 eta
    dx = x[1] - x[0]
    g = gaussian_eta(x, eta)

    assert np.isclose(np.sum(g)*dx, 1, rtol=0, atol=1e-9), 'gaussian not normalized'
    assert np.isclose(np.sum(x**2 * g)*dx, eta**2, rtol=1e-6, atol=0), 'eta is not the standard deviation'
    assert np.isclose(gaussian_eta(0, eta), 1/(eta*np.sqrt(2*np.pi)), rtol=1e-12, atol=0)


def test_kubo_accumulate_synthetic():
    """
    `kubo_accumulate` on hand-made input, against the closed form S_ab(w) = g(hw - gap) Re(v_a* v_b).

    k0: eps = (-0.6, 0.4) around mu = 0, one allowed transition c = 1 <- v = 0, gap 1 eV, with a
    complex interband velocity (its phase must drop out) and non-zero diagonal elements (which
    must not contribute). k1: both bands below mu, so no transition: occupations by energy.
    """
    eps = np.array([[-0.6, 0.4], [-2.0, -0.5]]) # (2 k, 2 bands), eV
    v_cv = np.array([1.5 + 2.0j, -0.7 + 0.3j, 0.0]) # eV Angstrom
    hv = np.zeros((2, 3, 2, 2), complex)
    hv[0, :, 1, 0] = v_cv
    hv[0, :, 0, 1] = np.conj(v_cv) # Hermitian partner
    hv[0, :, 0, 0] = [3.0, -1.0, 0.0]; hv[0, :, 1, 1] = [-3.0, 1.0, 0.0] # group velocities
    hv[1] = 5.0 # anything: k1 has no allowed transition
    hw = np.array([0.9, 1.0, 1.05])
    eta = 0.04

    S = kubo_accumulate(eps, hv, hw, mu=0, eta=eta)
    expected = gaussian_eta(hw - 1.0, eta)[:, None, None] * np.real(np.conj(v_cv)[:, None] * v_cv[None, :])

    assert S.shape == (3, 3, 3)
    assert np.allclose(S, expected, rtol=1e-12, atol=0), 'S does not match the single-transition closed form'

    # a global phase on the [c, v] element changes nothing
    hv_phase = hv.copy(); hv_phase[0, :, 1, 0] *= np.exp(0.7j); hv_phase[0, :, 0, 1] *= np.exp(-0.7j)
    assert np.allclose(kubo_accumulate(eps, hv_phase, hw, mu=0, eta=eta), S, rtol=1e-12, atol=0)


def test_kubo_accumulate_no_transition():
    """
    A block without any allowed transition (all bands below mu) gives S = 0, not an error.
    """
    eps = np.array([[-2.0, -1.0], [-3.0, -0.5]])
    hv = np.ones((2, 3, 2, 2), complex)

    S = kubo_accumulate(eps, hv, np.array([1.0, 2.0]), mu=0)

    assert S.shape == (2, 3, 3)
    assert np.array_equal(S, np.zeros((2, 3, 3)))


def test_sigma_dirac_limit():
    """
    Universal absorption: sigma/sigma_0 -> 1 for a Dirac cone (pi alpha). At 0.3 eV the trigonal
    correction is 1.4e-3 (it grows like hw^2: 1.0156 at 1 eV), hence rtol = 3e-3. Isotropic, and
    zero out of the plane.
    """
    sigma = sigma_on_grid(make_graphene_tb(), 900, 0.3)[0] # (3, 3)

    assert np.isclose(sigma[0, 0], 1, rtol=3e-3, atol=0), 'no universal absorption at low energy'
    assert np.isclose(sigma[1, 1], sigma[0, 0], rtol=1e-10, atol=0)
    assert abs(sigma[0, 1]) < 1e-10
    assert np.allclose(sigma[2, :], 0, atol=1e-14) and np.allclose(sigma[:, 2], 0, atol=1e-14)


@pytest.mark.parametrize("shift, mode", list(KUBO_REFERENCE), ids=[f"{m}-shift{s[0]}" for s, m in KUBO_REFERENCE])
def test_kubo_reference(shift, mode):
    """
    The 16 values of the EM.md table (two independent implementations, lattice and atomic gauge),
    to 1e-4 (4 decimals; observed 4.8e-5). 1800^2 grid, ~4 s per case.
    """
    ref = KUBO_REFERENCE[(shift, mode)]
    hw = np.array(list(ref))
    sigma = sigma_on_grid(make_graphene_tb(shift_B=shift), 1800, hw, eta=0.04, mode=mode)

    got = np.stack([sigma[:, 0, 0], sigma[:, 1, 1], sigma[:, 0, 1]], axis=1) # (nw, 3)
    assert np.allclose(got, np.array(list(ref.values())), rtol=0, atol=1e-4), f'sigma/sigma_0 = {got}'


def test_kubo_gauge_shift():
    """
    shift_B = (1,0,0) on a small grid: with Berry, sigma is unchanged to round-off (5e-15); without
    Berry, sigma_xx changes (1.23 -> 2.31 at 2.33 eV for N = 300) while sigma_yy does not (L = a1
    along x).
    """
    hw = np.array([1.0, 2.33])
    tb, tb_shifted = make_graphene_tb(), make_graphene_tb(shift_B=(1,0,0))

    assert np.allclose(sigma_on_grid(tb, 300, hw), sigma_on_grid(tb_shifted, 300, hw), rtol=0, atol=1e-12)

    s = sigma_on_grid(tb, 300, hw, mode='no_berry')
    s_shifted = sigma_on_grid(tb_shifted, 300, hw, mode='no_berry')
    assert np.all(np.abs(s_shifted[:, 0, 0] - s[:, 0, 0]) > 0.1), 'shift_B does not reach sigma_xx without Berry'
    assert np.allclose(s_shifted[:, 1, 1], s[:, 1, 1], rtol=0, atol=1e-12)


def test_sigma_on_grid_blocks_and_inputs():
    """
    The result does not depend on the block size (only N_k of the whole grid enters the
    normalization), a scalar hw equals a one-element array, and an unknown mode is refused.
    """
    tb = make_graphene_tb()
    hw = np.array([1.0, 2.33])

    assert np.allclose(sigma_on_grid(tb, 300, hw, chunk=7000), sigma_on_grid(tb, 300, hw, chunk=90000), rtol=0, atol=1e-12)
    assert np.allclose(sigma_on_grid(tb, 300, 2.33), sigma_on_grid(tb, 300, np.array([2.33])), rtol=0, atol=0)
    with pytest.raises((ValueError, TypeError)):
        sigma_on_grid(tb, 30, hw, mode='centres')
