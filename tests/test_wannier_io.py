"""
Tests of the Wannier90 `_tb.dat` reader `read_w90_tb` (io/wannier_io.py, F8 of the EM series) on the
27 x 27 wannierisation tracked in results/wannier/27x27/ (fixtures w90_dir, w90_ref and eig_w90 of
conftest.py). Reference values: campagnes/EM/EM.md, section 2.
"""

import re

import pytest
import numpy as np
from graphene_raman.io.wannier_io import read_w90_tb, read_w90_HR


@pytest.fixture(scope="module")
def w90(w90_dir):
    """(HR, R, ndegen, rR, lattice) of the 27 x 27 _tb.dat, read once for the module."""
    return read_w90_tb(w90_dir / "wannier_tb.dat")


@pytest.fixture(scope="module")
def centres(w90_dir):
    """Final Wannier centres of the .wout, (5, 3) Angstrom: independent of the reader."""
    text = (w90_dir / "wannier.wout").read_text()
    final = text[text.rindex("Final State"):]
    rows = re.findall(r"WF centre and spread\s+\d+\s+\(([^)]*)\)", final)
    return np.array([[float(x) for x in row.split(",")] for row in rows])


def minus_index(R):
    """Row of -R for every row of R (KeyError if some -R is missing)."""
    index = {tuple(r): i for i, r in enumerate(R)}
    return np.array([index[tuple(-r)] for r in R])


def test_tb_shapes(w90, w90_ref):
    HR, R, ndegen, rR, lattice = w90
    nR, nW = w90_ref.nR, w90_ref.nW

    assert HR.shape == (nR, nW, nW) and HR.dtype == complex
    assert R.shape == (nR, 3) and np.issubdtype(R.dtype, np.integer)
    assert ndegen.shape == (nR,)
    assert rR.shape == (nR, 3, nW, nW) and rR.dtype == complex
    assert lattice.shape == (3, 3)


def test_tb_lattice(w90, w90_ref):
    """Lattice vectors in the columns of `lattice`, Angstrom (EM.md section 2)."""
    lattice = w90[-1]

    for i, a_i in enumerate(w90_ref.lattice):
        assert np.allclose(lattice[:, i], a_i, rtol=0, atol=1e-6)


def test_tb_ndegen(w90, w90_ref):
    """sum 1/ndegen = number of k points of the grid, and ndegen(-R) = ndegen(R)."""
    _, R, ndegen, _, _ = w90

    assert np.isclose(np.sum(1/ndegen), w90_ref.n_grid**2, rtol=0, atol=1e-10)
    assert {int(d): int(np.sum(ndegen == d)) for d in np.unique(ndegen)} == w90_ref.ndegen_counts
    assert np.array_equal(ndegen[minus_index(R)], ndegen)


def test_tb_H_hermitian_r_not(w90, w90_ref):
    """H_ij(R) = H_ji(-R)*, but r keeps its finite-difference defect: the reader does not hermitize."""
    HR, R, _, rR, _ = w90
    minus = minus_index(R)

    assert np.abs(HR - HR[minus].swapaxes(-1, -2).conj()).max() < 1e-12
    defect = np.abs(rR - rR[minus].swapaxes(-1, -2).conj()).max(axis=(0, 2, 3)) # (3,), x y z
    assert np.allclose(defect[:2], w90_ref.r_defect_R, rtol=0.02, atol=0)


def test_tb_centres(w90, centres):
    """Diagonal of r(0) = .wout centres, real; also checks the Re/Im interleaving of the r lines."""
    _, R, _, rR, _ = w90
    i0 = np.flatnonzero((R == 0).all(axis=1))[0]
    diag = np.diagonal(rR[i0], axis1=-2, axis2=-1).T # (nW, 3)

    assert np.allclose(diag.real, centres, rtol=0, atol=5e-7)
    assert np.abs(diag.imag).max() < 1e-10


def test_tb_row_column_order(w90, centres, w90_ref):
    """
    X_R[iR, m, n] = <0m|X|Rn>, m = first column. Eigenvalues cannot see a transposition (H(k) -> H(k)*),
    hence a geometric test: the three largest |H_{pzA,pzB}(R)| join nearest neighbours, at a_cc.
    """
    HR, R, _, _, lattice = w90
    R_cart = R @ lattice.T
    A, B = w90_ref.pz # p_z of C1 and C2

    top = np.argsort(np.abs(HR[:, A, B]))[::-1][:3]
    dist = np.linalg.norm(centres[B] + R_cart[top] - centres[A], axis=1)

    assert {tuple(int(x) for x in r) for r in R[top]} == w90_ref.nn_R
    assert np.allclose(dist, w90_ref.a_cc, rtol=0, atol=1e-5)
    assert np.allclose(np.abs(HR[top, A, B]), w90_ref.t_nn, rtol=0, atol=1e-4) # eV, ~t


def interpolated_bands(w90, k_red):
    """Eigenvalues of sum_R e^{2 pi i k.R} H(R)/ndegen(R), written out here, eV."""
    HR, R, ndegen, _, _ = w90
    phase = np.exp(2j*np.pi * k_red @ R.T) / ndegen[None, :] # (nk, nR)
    H_k = np.einsum("kr,rmn->kmn", phase, HR)
    return np.linalg.eigvalsh(H_k)


def test_tb_bands_match_eig(w90, eig_w90, w90_ref):
    """
    Interpolated bands = .eig in the frozen window at the coarse k points (residue: use_ws_distance,
    which Wannier90 applies and we do not).
    """
    k_red, E_dft = eig_w90

    eps = interpolated_bands(w90, k_red)
    nfrozen = np.sum(E_dft < w90_ref.froz_max, axis=1)
    assert set(nfrozen) == w90_ref.nfrozen

    err = max(np.abs(eps[k, :n] - E_dft[k, :n]).max() for k, n in enumerate(nfrozen))
    assert err < w90_ref.eig_tol


def test_tb_dirac_point(w90, w90_ref):
    """pi and pi* degenerate at E_D at K = (2/3, 1/3, 0)."""
    eps = interpolated_bands(w90, np.array([[2/3, 1/3, 0]]))[0]
    v, c = w90_ref.bands_pi

    assert abs(eps[c] - eps[v]) < 1e-6
    assert np.allclose(eps[[v, c]], w90_ref.E_D, rtol=0, atol=1e-5)


def test_read_w90_HR_wrapper(w90, w90_dir):
    """The former name returns exactly (HR, R, ndegen) of read_w90_tb."""
    for old, new in zip(read_w90_HR(w90_dir / "wannier_tb.dat"), w90[:3]):
        assert np.array_equal(old, new)
