"""
Tests of the Wannier90 `_tb.dat` reader `read_w90_tb` (io/wannier_io.py, F8 of the EM series) on the
27 x 27 wannierisation tracked in wannier/27x27/ (fixtures w90_dir and eig_w90 of conftest.py).
Reference values: memoire/EM/EM.md, section 2.
"""

import re

import pytest
import numpy as np
from electron_defect_interaction.io.wannier_io import read_w90_tb, read_w90_HR

FROZ_MAX = -1.74 # eV, top of the frozen (inner) window of this wannierisation
E_D = -4.238895  # eV, Dirac point of the .eig at K
PZ_A, PZ_B = 3, 4 # WF 4 = p_z of C1, WF 5 = p_z of C2 (WF 1-3: sigma bonds)


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


def test_tb_shapes(w90):
    HR, R, ndegen, rR, lattice = w90

    assert HR.shape == (741, 5, 5) and HR.dtype == complex
    assert R.shape == (741, 3) and np.issubdtype(R.dtype, np.integer)
    assert ndegen.shape == (741,)
    assert rR.shape == (741, 3, 5, 5) and rR.dtype == complex
    assert lattice.shape == (3, 3)


def test_tb_lattice(w90):
    """Lattice vectors in the columns of `lattice`, Angstrom (EM.md section 2)."""
    lattice = w90[-1]

    assert np.allclose(lattice[:, 0], [2.135490, -1.232926, 0], rtol=0, atol=1e-6)
    assert np.allclose(lattice[:, 1], [2.135490, 1.232926, 0], rtol=0, atol=1e-6)
    assert np.allclose(lattice[:, 2], [0, 0, 15.875316], rtol=0, atol=1e-6)


def test_tb_ndegen(w90):
    """sum 1/ndegen = 27^2 (717 x 1, 24 x 2) and ndegen(-R) = ndegen(R)."""
    _, R, ndegen, _, _ = w90

    assert np.isclose(np.sum(1/ndegen), 27**2, rtol=0, atol=1e-10)
    assert np.sum(ndegen == 1) == 717 and np.sum(ndegen == 2) == 24
    assert np.array_equal(ndegen[minus_index(R)], ndegen)


def test_tb_H_hermitian_r_not(w90):
    """
    H_ij(R) = H_ji(-R)*, but r keeps its finite-difference defect (2.5e-3 Angstrom in x): the reader
    does not hermitize.
    """
    HR, R, _, rR, _ = w90
    minus = minus_index(R)

    assert np.abs(HR - HR[minus].swapaxes(-1, -2).conj()).max() < 1e-12
    defect = np.abs(rR - rR[minus].swapaxes(-1, -2).conj()).max(axis=(0, 2, 3)) # (3,), x y z
    assert np.allclose(defect[:2], [2.54e-3, 1.65e-3], rtol=0.02, atol=0)


def test_tb_centres(w90, centres):
    """Diagonal of r(0) = .wout centres, real; also checks the Re/Im interleaving of the r lines."""
    _, R, _, rR, _ = w90
    i0 = np.flatnonzero((R == 0).all(axis=1))[0]
    diag = np.diagonal(rR[i0], axis1=-2, axis2=-1).T # (nW, 3)

    assert np.allclose(diag.real, centres, rtol=0, atol=5e-7)
    assert np.abs(diag.imag).max() < 1e-10


def test_tb_row_column_order(w90, centres):
    """
    X_R[iR, m, n] = <0m|X|Rn>, m = first column. Eigenvalues cannot see a transposition (H(k) -> H(k)*),
    hence a geometric test: the three largest |H_45(R)| join nearest neighbours (R = 0, -a1, -a2, at a_cc).
    """
    HR, R, _, _, lattice = w90
    R_cart = R @ lattice.T

    top = np.argsort(np.abs(HR[:, PZ_A, PZ_B]))[::-1][:3]
    dist = np.linalg.norm(centres[PZ_B] + R_cart[top] - centres[PZ_A], axis=1)

    assert {tuple(r) for r in R[top]} == {(0, 0, 0), (-1, 0, 0), (0, -1, 0)}
    assert np.allclose(dist, 1.42366, rtol=0, atol=1e-5)
    assert np.allclose(np.abs(HR[top, PZ_A, PZ_B]), 2.9089, rtol=0, atol=1e-4) # eV, ~t


def interpolated_bands(w90, k_red):
    """Eigenvalues of sum_R e^{2 pi i k.R} H(R)/ndegen(R), written out here, eV."""
    HR, R, ndegen, _, _ = w90
    phase = np.exp(2j*np.pi * k_red @ R.T) / ndegen[None, :] # (nk, nR)
    H_k = np.einsum("kr,rmn->kmn", phase, HR)
    return np.linalg.eigvalsh(H_k)


def test_tb_bands_match_eig(w90, eig_w90):
    """
    Interpolated bands = .eig in the frozen window at the 729 coarse points (1.3e-5 eV: use_ws_distance,
    which Wannier90 applies and we do not).
    """
    k_red, E_dft = eig_w90

    eps = interpolated_bands(w90, k_red)
    nfrozen = np.sum(E_dft < FROZ_MAX, axis=1)
    assert set(nfrozen) == {4, 5}

    err = max(np.abs(eps[k, :n] - E_dft[k, :n]).max() for k, n in enumerate(nfrozen))
    assert err < 3e-5


def test_tb_dirac_point(w90):
    """p_z bands degenerate at E_D at K = (2/3, 1/3, 0)."""
    eps = interpolated_bands(w90, np.array([[2/3, 1/3, 0]]))[0]

    assert abs(eps[4] - eps[3]) < 1e-6
    assert np.allclose(eps[3:], E_D, rtol=0, atol=1e-5)


def test_read_w90_HR_wrapper(w90, w90_dir):
    """The former name returns exactly (HR, R, ndegen) of read_w90_tb."""
    for old, new in zip(read_w90_HR(w90_dir / "wannier_tb.dat"), w90[:3]):
        assert np.array_equal(old, new)
