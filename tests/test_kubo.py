"""
Tests of the Kubo conductivity (F7) and of the driver `sigma_on_grid`, down to the reference table
of EM.md on the 1800^2 grid (~13 s; skip with -k "not kubo_reference"), and sigma(omega) of the
27 x 27 data for the three velocity variants (M4; ~7 s); doping and temperature (kT, perspective A);
complex sigma with the complex kernels (perspective B).
"""

import pytest
import numpy as np
from scipy.integrate import quad
from electron_defect_interaction.electron_photon import (centres_only, gaussian_eta, gaussian_complex, kubo_accumulate,
                                                         kubo_normalize, lorentzian_complex,
                                                         kubo_doped_complex_analytical, kubo_doped_finite_T_analytical, make_graphene_tb,
                                                         sigma_on_grid)

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


def test_complex_kernels():
    """
    Complex kernels (perspective B): Re[iK] is pi times the normalized line (the Gaussian exactly, integral pi),
    Im[iK] -> 1/x far away with opposite O(eta^2/x^2) corrections (Lorentzian 1 - eta^2, Gaussian 1 + eta^2 at
    x = 1), K(-x) = -K(x)*, and reference values at x = eta.
    """
    eta = 0.04
    x = np.linspace(-3, 3, 600001)
    dx = x[1] - x[0]
    iKg, iKl = 1j*gaussian_complex(x, eta), 1j*lorentzian_complex(x, eta)

    assert np.allclose(iKg.real, np.pi*gaussian_eta(x, eta), rtol=0, atol=1e-12), 'Re[iK] is not pi g_eta'
    assert np.isclose(np.sum(iKg.real)*dx, np.pi, rtol=0, atol=1e-9)
    assert np.isclose(np.sum(iKl.real)*dx, 2*np.arctan(3/eta), rtol=0, atol=1e-7), 'Lorentzian of half-width eta'
    assert np.isclose((1j*lorentzian_complex(1.0, eta)).imag, 1 - eta**2, rtol=0, atol=1e-5)
    assert np.isclose((1j*gaussian_complex(1.0, eta)).imag, 1 + eta**2, rtol=0, atol=1e-5)
    for K in (gaussian_complex, lorentzian_complex):
        assert np.array_equal(K(-x, eta), -np.conj(K(x, eta)))
    assert np.isclose(1j*lorentzian_complex(eta, eta), 12.5 + 12.5j, rtol=1e-12, atol=0)
    assert np.isclose(1j*gaussian_complex(eta, eta), 19.004336 + 18.119461j, rtol=1e-7, atol=0)


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
    with pytest.raises(ValueError):
        sigma_on_grid(tb, 30, hw, mode='centres')


LASERS = [1.96, 2.33, 2.54] # eV, 633, 532, 488 nm


@pytest.fixture(scope="module")
def sigma_w90(tb_w90, w90_ref):
    """
    sigma(omega)/sigma_0 of the 27 x 27 data, three variants, N = 300 and eta = 0.08 eV (converged to 1e-5
    at the lasers, M4 pilot), at the lasers then 3.80 ... 4.30 eV (van Hove peak).
    """
    hw = np.array(LASERS + list(np.arange(380, 431) / 100))
    variants = {'full': (tb_w90, 'berry'), 'centres_only': (centres_only(tb_w90), 'berry'), 'no_berry': (tb_w90, 'no_berry')}
    return hw, {v: sigma_on_grid(model, 300, hw, mu=w90_ref.E_D, eta=0.08, mode=mode) for v, (model, mode) in variants.items()}


def test_sigma_real_values(sigma_w90, w90_ref):
    """At the lasers, the three variants reproduce the converged pilot (N = 1800, other code path) to 3e-5."""
    _, sigma = sigma_w90
    for v, ref in w90_ref.sigma_eta008.items():
        assert np.allclose(sigma[v][:3, 0, 0], ref, rtol=0, atol=3e-5), v


def test_sigma_real_symmetries(sigma_w90):
    """
    Full and centres only: isotropic (C3; 9e-4 at N = 300, 2e-4 converged) and sigma_xy = 0 (mirror).
    Without Berry: xx != yy, but yy equals centres only to rounding (tau_B - tau_A along x: no y term).
    """
    _, sigma = sigma_w90
    for v in ('full', 'centres_only'):
        assert np.allclose(sigma[v][:, 1, 1], sigma[v][:, 0, 0], rtol=2e-3, atol=0), v
    assert all(np.abs(s[:, 0, 1]).max() < 1e-8 for s in sigma.values())
    assert np.all(sigma['no_berry'][:3, 0, 0] / sigma['no_berry'][:3, 1, 1] > 1.1)
    assert np.allclose(sigma['no_berry'][:, 1, 1], sigma['centres_only'][:, 1, 1], rtol=1e-12, atol=0)


def test_sigma_real_features(sigma_w90):
    """
    van Hove peak at 4.05 eV (transition at M, 4.056 eV) in every variant; centres only / full = 0.951
    at 2.33 eV, the ratio of <|hbar v_cv|^2> on the ring (F15).
    """
    hw, sigma = sigma_w90
    for v, s in sigma.items():
        assert np.isclose(hw[3:][np.argmax(s[3:, 0, 0])], 4.05, rtol=0, atol=1e-9), v
    assert np.isclose(sigma['centres_only'][1, 0, 0] / sigma['full'][1, 0, 0], 0.951, rtol=0, atol=1e-3)


HW_EDGE = np.array([0.2, 0.4, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8, 1.0]) # eV, around the Pauli edge 2 mu = 0.6


def test_dirac_analytical_limits():
    """The reference formula: tanh(mu/kT)/2 at hw = 2 mu, a step as kT -> 0, tanh(hw/4kT) undoped, even in mu."""
    assert np.isclose(kubo_doped_finite_T_analytical(0.6, 0.3, 0.025), np.tanh(0.3 / 0.025) / 2, rtol=0, atol=1e-15)
    assert np.allclose(kubo_doped_finite_T_analytical(np.array([0.5, 0.7]), 0.3, 1e-4), [0, 1], rtol=0, atol=1e-12)
    assert np.allclose(kubo_doped_finite_T_analytical(HW_EDGE, 0.0, 0.025), np.tanh(HW_EDGE / 0.1), rtol=0, atol=1e-15)
    assert np.allclose(kubo_doped_finite_T_analytical(HW_EDGE, 0.3, 0.025), kubo_doped_finite_T_analytical(HW_EDGE, -0.3, 0.025), rtol=0, atol=1e-15)
    with pytest.raises(AssertionError):
        kubo_doped_finite_T_analytical(HW_EDGE, 0.3, 0.0)


def test_kubo_T0_limit_and_inputs():
    """kT -> 0 joins the step (T = 0 unchanged); independent of chunk at kT > 0; negative kT refused."""
    tb = make_graphene_tb()
    step = sigma_on_grid(tb, 300, HW_EDGE, mu=0.3)

    assert np.array_equal(sigma_on_grid(tb, 300, HW_EDGE, mu=0.3, kT=1e-6), step)
    assert np.allclose(sigma_on_grid(tb, 300, HW_EDGE, mu=0.3, kT=0.025, chunk=7000),
                       sigma_on_grid(tb, 300, HW_EDGE, mu=0.3, kT=0.025), rtol=0, atol=1e-12)
    with pytest.raises(AssertionError):
        sigma_on_grid(tb, 30, HW_EDGE, kT=-0.01)


@pytest.mark.parametrize("kT, tol", [(0.025, 4e-3), (0.05, 1.5e-3)])
def test_kubo_pauli_blocking_graphene(kT, tol):
    """
    M0 doped (mu = 0.3 eV) at finite T: sigma / sigma(mu = 0, T = 0) follows the Dirac formula
    (trigonal warping cancels in the ratio); the gap is the eta rounding of the edge (2.6e-3, 8.5e-4).
    Below the edge the thermal tail (weights down to 3e-4) matches in relative terms, to 3 % (eta/kT).
    """
    tb = make_graphene_tb()
    ratio = sigma_on_grid(tb, 900, HW_EDGE, mu=0.3, eta=0.01, kT=kT)[:, 0, 0] / sigma_on_grid(tb, 900, HW_EDGE, eta=0.01)[:, 0, 0]
    ref = kubo_doped_finite_T_analytical(HW_EDGE, 0.3, kT)

    assert np.allclose(ratio, ref, rtol=0, atol=tol)
    tail = HW_EDGE <= 0.4
    assert np.allclose(ratio[tail] / ref[tail], 1, rtol=0, atol=0.05)


def test_kubo_thermal_undoped_graphene():
    """M0 at mu = 0: finite T only thins the lowest transitions, sigma(T)/sigma(0) = tanh(hw / 4kT) (3.5e-4)."""
    tb = make_graphene_tb()
    ratio = sigma_on_grid(tb, 900, HW_EDGE, eta=0.01, kT=0.025)[:, 0, 0] / sigma_on_grid(tb, 900, HW_EDGE, eta=0.01)[:, 0, 0]

    assert np.allclose(ratio, np.tanh(HW_EDGE / 0.1), rtol=0, atol=1e-3)


def test_kubo_electron_hole_graphene():
    """M0 is electron-hole symmetric: n and p doping give the same sigma, to rounding."""
    tb = make_graphene_tb()
    assert np.allclose(sigma_on_grid(tb, 600, HW_EDGE, mu=0.3, eta=0.02, kT=0.025),
                       sigma_on_grid(tb, 600, HW_EDGE, mu=-0.3, eta=0.02, kT=0.025), rtol=0, atol=1e-13)


def test_kubo_pauli_edges_real(tb_w90, w90_ref):
    """
    27 x 27 data, mu = E_D +/- 0.3 eV, kT = 0.025: the half-max edge is 7 meV higher for n than for p doping
    (PBE electron-hole asymmetry; 6.9 meV converged, EM.md); sigma stays positive.
    """
    hw = np.round(np.arange(0.40, 0.801, 0.01), 2)
    base = sigma_on_grid(tb_w90, 300, hw, mu=w90_ref.E_D, eta=0.08)[:, 0, 0]
    edges = []
    for dmu in (0.3, -0.3):
        s = sigma_on_grid(tb_w90, 300, hw, mu=w90_ref.E_D + dmu, eta=0.08, kT=0.025)[:, 0, 0]
        assert np.all(s > 0)
        r = s / base
        i = np.flatnonzero(np.diff(np.sign(r - 0.5)))[0] # r crosses 1/2 between hw[i] and hw[i+1]
        edges.append(hw[i] + (0.5 - r[i]) * (hw[i+1] - hw[i]) / (r[i+1] - r[i]))

    assert np.allclose(edges, w90_ref.pauli_edges, rtol=0, atol=1e-4)
    assert 0.006 < edges[0] - edges[1] < 0.008


# M0, mu = 0.3, T = 0, N = 900, eta = 0.05: sigma_xx/sigma_0 at 0.2, 0.6, 1.0 eV (independent implementation, 2e-16)
COMPLEX_M0 = {gaussian_complex:   (0.000000-0.223532j, 0.501278-1.223066j, 1.015612-0.455188j),
              lorentzian_complex: (0.059750-0.219711j, 0.515714-1.019755j, 0.986309-0.448873j)}


def test_kubo_accumulate_complex_synthetic():
    """
    Complex path on the hand-made input of test_kubo_accumulate_synthetic (one transition, gap 1 eV):
    S_ab = (w/gap) Re(v_a* v_b) [K(hw - gap) + K(hw + gap)] for both kernels; kubo_normalize gives 8i S/(N_k A),
    independent of hw.
    """
    eps = np.array([[-0.6, 0.4], [-2.0, -0.5]])
    v_cv = np.array([1.5 + 2.0j, -0.7 + 0.3j, 0.0])
    hv = np.zeros((2, 3, 2, 2), complex)
    hv[0, :, 1, 0] = v_cv; hv[0, :, 0, 1] = np.conj(v_cv)
    hv[0, :, 0, 0] = [3.0, -1.0, 0.0]; hv[0, :, 1, 1] = [-3.0, 1.0, 0.0]
    hv[1] = 5.0
    hw = np.array([0.9, 1.0, 1.05])
    W = np.real(np.conj(v_cv)[:, None] * v_cv[None, :])
    for K in (gaussian_complex, lorentzian_complex):
        S = kubo_accumulate(eps, hv, hw, mu=0, eta=0.04, kernel=K)
        expected = (K(hw - 1.0, 0.04) + K(hw + 1.0, 0.04))[:, None, None] * W
        assert np.iscomplexobj(S) and S.shape == (3, 3, 3)
        assert np.allclose(S, expected, rtol=1e-13, atol=0), K.__name__
        assert np.allclose(kubo_normalize(S, hw, 10, 2.0, kernel=K), 8j * S / 20, rtol=1e-15, atol=0)
        assert np.array_equal(kubo_normalize(S, hw, 10, 2.0, kernel=K), kubo_normalize(S, 2 * hw, 10, 2.0, kernel=K))


def test_sigma_complex_graphene():
    """
    M0 doped (mu = 0.3), both kernels: values (independent implementation, 1e-6); Gaussian: Re = 0 in the
    Pauli window, Im < 0 below 2 mu; Lorentzian: its tails leave Re = 0.06 at 0.2 eV. Re >= 0 everywhere.
    """
    tb = make_graphene_tb()
    hw = np.array([0.2, 0.6, 1.0])
    s = {K: sigma_on_grid(tb, 900, hw, mu=0.3, eta=0.05, kernel=K)[:, 0, 0] for K in COMPLEX_M0}
    for K, ref in COMPLEX_M0.items():
        assert np.allclose(s[K], ref, rtol=0, atol=1e-6), K.__name__
    assert abs(s[gaussian_complex][0].real) < 1e-6 and s[lorentzian_complex][0].real > 0.05
    s = sigma_on_grid(tb, 300, np.linspace(0.1, 3, 30), mu=0.3, eta=0.05, kernel=gaussian_complex)[:, 0, 0]
    assert np.all(s.real >= 0) and np.all(s.imag[:5] < 0)


def test_sigma_complex_symmetries_and_blocks():
    """
    Exact relations of the complex path on M0: sigma(-w) = sigma(w)* (K(-x) = -K(x)*), n = p doping
    (electron-hole symmetric model), and independence from the block size.
    """
    tb = make_graphene_tb()
    hw = np.array([0.3, 0.7, 2.0])
    for K in (gaussian_complex, lorentzian_complex):
        s = sigma_on_grid(tb, 300, hw, mu=0.3, eta=0.05, kT=0.025, kernel=K)
        assert np.allclose(sigma_on_grid(tb, 300, -hw, mu=0.3, eta=0.05, kT=0.025, kernel=K), np.conj(s), rtol=0, atol=1e-13)
        assert np.allclose(sigma_on_grid(tb, 300, hw, mu=-0.3, eta=0.05, kT=0.025, kernel=K), s, rtol=0, atol=1e-13)
        assert np.allclose(sigma_on_grid(tb, 300, hw, mu=0.3, eta=0.05, kT=0.025, kernel=K, chunk=7000), s, rtol=0, atol=1e-12)


def test_sigma_complex_real(tb_w90, w90_ref, sigma_w90):
    """
    27 x 27 data, gaussian_complex: values at the lasers and at the van Hove peak (kT = 0.025, N = 300,
    eta = 0.08); Re independent of kT and below the absorptive M4 value by the 1/Delta prefactor (-1e-3 at
    eta = 0.08, EM.md decision 9); Im < 0 at the lasers, > 0 above the van Hove peak (Kramers-Kronig).
    """
    hw = np.array(LASERS + [4.05])
    s = sigma_on_grid(tb_w90, 300, hw, mu=w90_ref.E_D, eta=0.08, kT=0.025, kernel=gaussian_complex)[:, 0, 0]
    s0 = sigma_on_grid(tb_w90, 300, hw, mu=w90_ref.E_D, eta=0.08, kernel=gaussian_complex)[:, 0, 0]
    absorptive = sigma_w90[1]['full'][:3, 0, 0]

    assert np.allclose(s, w90_ref.sigma_complex, rtol=0, atol=1e-5)
    assert np.allclose(s.real, s0.real, rtol=0, atol=1e-5)
    assert np.all((-2e-3 < s0.real[:3] - absorptive) & (s0.real[:3] - absorptive < -5e-4))
    assert np.all(s.imag[:3] < 0) and s.imag[3] > 0


def test_complex_analytical_limits():
    """
    Dirac-cone reference: eta -> 0 gives the Pauli step and -(1/pi) ln|(hw + 2mu)/(hw - 2mu)|; the edge value
    1/2 + eta/(4 pi |mu|); exactly 1 at mu = 0; sigma(-w) = sigma(w)*; n = p doping.
    """
    hw = np.array([0.2, 0.5, 0.7, 1.0, 2.0])
    mu = 0.3
    s = kubo_doped_complex_analytical(hw, mu, 1e-9)
    assert np.allclose(s.real, (hw > 2*mu).astype(float), rtol=0, atol=1e-8)
    assert np.allclose(s.imag, -np.log(np.abs((hw + 2*mu) / (hw - 2*mu))) / np.pi, rtol=0, atol=1e-8)
    for eta in (0.01, 0.05):
        assert np.isclose(kubo_doped_complex_analytical(2*mu, mu, eta).real, 0.5 + eta / (4*np.pi*mu), rtol=0, atol=2e-5)
    assert np.allclose(kubo_doped_complex_analytical(hw, 0.0, 0.05), 1, rtol=0, atol=1e-15)
    assert np.allclose(kubo_doped_complex_analytical(-hw, mu, 0.05), np.conj(kubo_doped_complex_analytical(hw, mu, 0.05)), rtol=0, atol=1e-14)
    assert np.array_equal(kubo_doped_complex_analytical(hw, -mu, 0.05), kubo_doped_complex_analytical(hw, mu, 0.05))


def _cone_quad(K, hw, mu, eta):
    """Dirac cone with the kernel K: (i/pi) int_{2|mu|}^inf [K(hw - D) + K(hw + D)] dD, by quadrature."""
    f = lambda D: K(hw - D, eta) + K(hw + D, eta)
    a = 2*abs(mu)
    parts = [(a, hw), (hw, np.inf)] if hw > a else [(a, np.inf)] # split at the resonance
    val = sum(quad(lambda D: f(D).real, *p, limit=400)[0] + 1j*quad(lambda D: f(D).imag, *p, limit=400)[0] for p in parts)
    return 1j * val / np.pi


def test_sigma_complex_dirac_cone_graphene():
    """
    M0 against the Dirac cone: the doping difference sigma(0.3) - sigma(0.2) cancels the lattice background
    (finite bandwidth, van Hove) and matches the cone to 4e-3 (trigonal warping, independent of N) for both
    kernels: closed form for the Lorentzian, quadrature for the Gaussian. The Lorentzian quadrature equals
    the closed form (1e-8).
    """
    tb = make_graphene_tb()
    hw = np.array([0.1, 0.2, 0.3, 0.5, 0.8, 1.2])
    d = {K: (sigma_on_grid(tb, 600, hw, mu=0.3, eta=0.05, kernel=K) - sigma_on_grid(tb, 600, hw, mu=0.2, eta=0.05, kernel=K))[:, 0, 0]
         for K in (lorentzian_complex, gaussian_complex)}
    cone_L = kubo_doped_complex_analytical(hw, 0.3, 0.05) - kubo_doped_complex_analytical(hw, 0.2, 0.05)
    cone_G = np.array([_cone_quad(gaussian_complex, h, 0.3, 0.05) - _cone_quad(gaussian_complex, h, 0.2, 0.05) for h in hw])
    quad_L = np.array([_cone_quad(lorentzian_complex, h, 0.3, 0.05) for h in hw])

    assert np.allclose(quad_L, kubo_doped_complex_analytical(hw, 0.3, 0.05), rtol=0, atol=1e-8) # validates the quadrature
    assert np.abs(cone_L).max() > 0.5
    assert np.allclose(d[lorentzian_complex], cone_L, rtol=0, atol=4e-3)
    assert np.allclose(d[gaussian_complex], cone_G, rtol=0, atol=5e-3)

