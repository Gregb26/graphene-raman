""" Testing of computation of velocity operator """

import pytest
import numpy as np
from electron_defect_interaction.optics.velocity_operator import *
from types import SimpleNamespace

N = 100

@pytest.fixture
def tb():
    return make_graphene_tb()

@pytest.fixture
def grid(tb):
    B = reciprocal(tb.lattice)
    k_red, k_cart, K = k_grid(B,N)
    R_cart = tb.R_int @ tb.lattice.T
    index = {tuple(R): iR for iR, R in enumerate(tb.R_int)}
    minus = np.array([index[tuple(-R)] for R in tb.R_int]) 

    return SimpleNamespace(B=B, k_red=k_red, k_cart=k_cart, K=K, R_cart=R_cart, index=index, minus=minus)


def make_graphene_tb_analytic(tb, grid):
    """
    Inputs:
        -t : float, hopping parameter
        k  : (Nk, 3) ndarray of floats, kpoint grid in cartesian coordinates
        A  : (3, 3)  ndarray of floats, primitive lattice vectors stored in columns A[:,i] = a_i
    Returns:
        TB Hamiltonian : H_AB(k) = -t(1+e^{-ik.a1}+e^{-ik.a2})
    """

    t = tb.t
    A = tb.lattice
    k = grid.k_cart

    return -t*(1+np.exp(-1j*k @ A[:,0].T) + np.exp(-1j*k @ A[:,1].T))

def test_fourier_analytic(tb, grid):
    """ Test that Fourier-transformed TB Hamiltonian reproduces analytical solution """

    H_tb = make_graphene_tb_analytic(tb, grid)
    H_k = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.k_cart)

    assert np.allclose(H_tb, H_k[:,0,1], atol=1e-12), 'tight binding model not correct'

def test_dirac_point(tb, grid):
    """ Test that H_AB(k) = 0 """

    H_K = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.K[None, :])

    assert np.allclose(H_K, 0, atol=1e-12), 'Hamiltonian not zero at Dirac point'

def test_tb_hermitian(tb, grid):
    """  Test that TB Hamiltonian is hermitian """  

    assert np.allclose(grid.minus[grid.minus], np.arange(len(tb.R_int))), 'R list not closed under R -> -R'
    assert np.allclose(tb.H_R, tb.H_R[grid.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'TB Hamiltonian not hermitian'
    assert np.allclose(tb.r_R, tb.r_R[grid.minus].swapaxes(-1, -2).conj(), atol=1e-12), 'positian operator not hermitian'
    assert np.allclose(tb.ndegen, tb.ndegen[grid.minus], atol=1e-12), 'R and -R dont have the same weight'

def test_reciprocal(tb, grid):
    """ Test that a_i.b_j = 2 pi delta_ij """

    assert np.allclose(tb.A @ grid.B.T, 2*np.pi*np.eye(3), atol=1e-12)

def test_fourier_hermitian(tb, grid):
    """ Test that Fourier-transformed objects are hermitian """ 

    H_k = fourier(tb.H_R, grid.R_cart, tb.ndegen, grid.k_cart)
    A_k = fourier(tb.r_R, grid.R_cart, tb.ndegen, grid.k_cart)
    dH_k = fourier(tb.r_R, grid.R_cart, tb.ndegen, grid.k_cart, deriv=True)

    assert np.allclose(H_k, dagger(H_k), atol=1e-12), 'H_k not hermitian'
    assert np.allclose(A_k, dagger(A_k), atol=1e-12), 'A_k not hermitian'
    assert np.allclose(dH_k, dagger(dH_k), atol=1e-12), 'dH_k not hermitian'

def test_hermitize(tb, grid):
    """ Test that hermitize function really hermitizes and is idempotent """

    np.random.default_rng(0)
    A_k = fourier(tb.r_R, grid.R_cart, tb.ndegen, grid.k_cart)

    A = np.random.rand(A_k.shape) + 1j*np.random.rand(A_k.shape)

    assert np.allclose(hermitize(A), dagger(hermitize(A)), atol=1e-12), 'hermitize function does not hermitize the object'
    assert np.allclose(hermitize(A), hermitize(hermitize(A)), atol=1e-12), 'hermitize function not idempotent'
    assert np.allclose(A_k, hermitize(A_k), atol=1e-12), 'hermitize function breaks hermicity of hermitian object'




