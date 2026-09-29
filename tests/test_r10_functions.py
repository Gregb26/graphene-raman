"""Minimal tests of the R10 functions (2026-09-29): ws_images, ws_phase, Mwr_to_Mwk(ws=), defect_mwr.
Run: PYTHONPATH=src pytest tests/test_r10_functions.py"""
import numpy as np
from electron_defect_interaction.wannier.wannier_interpolation import ws_images


def _graphene_cell(a=2.4659, c=8.0):
    """Graphene unit cell, a1 and a2 at 60 degrees (as the QE cells), vectors in COLUMNS."""
    return np.array([[a * np.sqrt(3) / 2, a * np.sqrt(3) / 2, 0.0], [-a / 2, a / 2, 0.0], [0.0, 0.0, c]])


def _mp_labels(D):
    """R labels of Mwk_to_Mwr for a D x D x 1 grid: arange(D) - D//2 per axis, same (ij) ordering."""
    r = np.arange(D) - D // 2
    return np.array([[i, j, 0] for i in r for j in r])


def test_ws_images_audit_case():
    """T3, audit case (P-c4): D = 20, R_center = 0, R = (9, 9, 0) -> two images at sqrt(103) a with weight 1/2; per-axis norm 9 sqrt(3) a."""
    a = 2.4659; A = _graphene_cell(a); D = 20
    R = _mp_labels(D)
    ws = ws_images(R, np.zeros(3, int), (D, D, 1), A)
    assert ws["R_img"].dtype.kind == "i" and ws["dist"].shape == (len(R),) and ws["n_tie"].shape == (len(R),)
    r = int(np.flatnonzero(np.all(R == [9, 9, 0], axis=1))[0])
    sel = ws["idx"] == r
    assert sorted(map(tuple, ws["R_img"][sel].tolist())) == [(-11, 9, 0), (9, -11, 0)]
    assert np.array_equal(ws["w"][sel], [0.5, 0.5]) and ws["n_tie"][r] == 2
    assert abs(ws["dist"][r] - np.sqrt(103) * a) < 1e-12 * D * a
    assert abs(ws["dist"][r] / (D * a) - 0.507) < 5e-4
    assert abs(np.linalg.norm(A @ R[r]) / (D * a) - 0.779) < 5e-4          # the per-axis norm of the audit bug


def test_ws_images_inside_ws_cell_single_image():
    """T4: a label already inside the Wigner-Seitz cell keeps a single image, weight 1, R_img = R (explicit case, then every |A d| < D a / 2)."""
    a = 2.4659; A = _graphene_cell(a)
    ws = ws_images(np.array([[2, 1, 0]]), np.zeros(3, int), (20, 20, 1), A)
    assert np.array_equal(ws["R_img"], [[2, 1, 0]]) and np.array_equal(ws["idx"], [0]) and np.array_equal(ws["w"], [1.0])
    assert ws["n_tie"][0] == 1
    for D, Rc in ((20, (0, 0, 0)), (27, (4, 4, 0)), (6, (1, -2, 0))):
        R = _mp_labels(D); Rc = np.array(Rc)
        inside = np.flatnonzero(np.linalg.norm((R - Rc) @ A.T, axis=1) < D * a / 2 - 1e-9)
        ws = ws_images(R, Rc, (D, D, 1), A)
        assert np.all(np.bincount(ws["idx"], minlength=len(R))[inside] == 1) and np.all(ws["n_tie"][inside] == 1)
        one = np.isin(ws["idx"], inside)
        assert np.array_equal(ws["R_img"][one], R[ws["idx"][one]]) and np.all(ws["w"][one] == 1.0)
