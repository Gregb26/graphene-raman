"""Minimal tests of the R10 functions (2026-09-29): ws_images, ws_phase, Mwr_to_Mwk(ws=), defect_mwr.
Run: PYTHONPATH=src pytest tests/test_r10_functions.py"""
import numpy as np
import electron_defect_interaction.wannier.wannier_interpolation as wi
import electron_defect_interaction.defects.many_body.local_tmatrix as lt


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
    ws = wi.ws_images(R, np.zeros(3, int), (D, D, 1), A)
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
    ws = wi.ws_images(np.array([[2, 1, 0]]), np.zeros(3, int), (20, 20, 1), A)
    assert np.array_equal(ws["R_img"], [[2, 1, 0]]) and np.array_equal(ws["idx"], [0]) and np.array_equal(ws["w"], [1.0])
    assert ws["n_tie"][0] == 1
    for D, Rc in ((20, (0, 0, 0)), (27, (4, 4, 0)), (6, (1, -2, 0))):
        R = _mp_labels(D); Rc = np.array(Rc)
        inside = np.flatnonzero(np.linalg.norm((R - Rc) @ A.T, axis=1) < D * a / 2 - 1e-9)
        ws = wi.ws_images(R, Rc, (D, D, 1), A)
        assert np.all(np.bincount(ws["idx"], minlength=len(R))[inside] == 1) and np.all(ws["n_tie"][inside] == 1)
        one = np.isin(ws["idx"], inside)
        assert np.array_equal(ws["R_img"][one], R[ws["idx"][one]]) and np.all(ws["w"][one] == 1.0)


def test_ws_images_invariants():
    """T2a: for every label, the weights sum to 1, each image is R_r modulo D, and each image lies at the true distance dist[r] (odd and even D)."""
    a = 2.4659; A = _graphene_cell(a)
    for D, Rc in ((5, (0, 0, 0)), (6, (0, 0, 0)), (6, (2, -1, 0)), (27, (4, 4, 0))):
        R = _mp_labels(D); Rc = np.array(Rc)
        ws = wi.ws_images(R, Rc, (D, D, 1), A)
        assert np.allclose(np.bincount(ws["idx"], weights=ws["w"], minlength=len(R)), 1.0, rtol=0, atol=1e-12)
        assert np.array_equal(np.mod(ws["R_img"] - R[ws["idx"]], [D, D, 1]), np.zeros_like(ws["R_img"]))
        L = np.linalg.norm((ws["R_img"] - Rc) @ A.T, axis=1)
        assert np.allclose(L, ws["dist"][ws["idx"]], rtol=0, atol=1e-12 * D * a)
        assert np.array_equal(np.bincount(ws["idx"], minlength=len(R)), ws["n_tie"])


def test_ws_phase():
    """T2b: ws_phase = plain exp(sign 2 pi i k.R_r) on the MP grid k = m/D; off the grid, the weighted sum over the tied images (audit case)."""
    a = 2.4659; A = _graphene_cell(a)
    for D, Rc in ((5, (0, 0, 0)), (6, (2, -1, 0))):
        R = _mp_labels(D); ws = wi.ws_images(R, np.array(Rc), (D, D, 1), A)
        m = np.arange(D)
        k = np.array([[i / D, j / D, 0.0] for i in m for j in m])
        for sign in (+1, -1):
            P = wi.ws_phase(k, ws, len(R), sign)
            assert P.shape == (len(k), len(R)) and np.iscomplexobj(P)
            assert np.allclose(P, np.exp(sign * 2j * np.pi * (k @ R.T)), rtol=0, atol=1e-12)
    D = 20; R = _mp_labels(D); ws = wi.ws_images(R, np.zeros(3, int), (D, D, 1), A)
    k = np.array([[0.123, 0.456, 0.0], [1 / 3, 2 / 3, 0.0]])
    r = int(np.flatnonzero(np.all(R == [9, 9, 0], axis=1))[0])
    ref = 0.5 * (np.exp(2j * np.pi * (k @ [-11, 9, 0])) + np.exp(2j * np.pi * (k @ [9, -11, 0])))
    Pp, Pm = wi.ws_phase(k, ws, len(R), +1), wi.ws_phase(k, ws, len(R), -1)
    assert np.allclose(Pp[:, r], ref, rtol=0, atol=1e-12) and np.allclose(Pm, Pp.conj(), rtol=0, atol=1e-12)
    assert not np.allclose(Pp[:, r], np.exp(2j * np.pi * (k @ R[r])), atol=1e-3)     # off the grid the label choice matters


def test_Mwr_to_Mwk_ws_option():
    """T1: with ws=, Mwr_to_Mwk and Mwr_to_Mwk_pairs equal the plain FT on the MP grid k = m/D (D odd and even); off the grid, the
    weighted double sum over the images of R and R'."""
    rng = np.random.default_rng(10); A = _graphene_cell(); nw = 2
    for D, Rc in ((5, (1, 2, 0)), (6, (-2, 1, 0))):
        R = _mp_labels(D); nR = len(R)
        Mwr = rng.standard_normal((nw, nR, nw, nR)) + 1j * rng.standard_normal((nw, nR, nw, nR))
        ws = wi.ws_images(R, np.array(Rc), (D, D, 1), A)
        m = np.arange(D); k = np.array([[i / D, j / D, 0.0] for i in m for j in m])
        ref = wi.Mwr_to_Mwk(Mwr, R, k); tol = 1e-12 * np.abs(ref).max()
        assert np.allclose(wi.Mwr_to_Mwk(Mwr, R, k, ws=ws), ref, rtol=0, atol=tol)
        assert np.allclose(wi.Mwr_to_Mwk_pairs(Mwr, R, k[:7], k[7:12], ws=ws), wi.Mwr_to_Mwk_pairs(Mwr, R, k[:7], k[7:12]), rtol=0, atol=tol)
        # off the grid: independent double sum over the images a (bra) and b (ket), M_W(R_img_a, R_img_b) = Mwr[:, idx_a, :, idx_b]
        kb = rng.random((3, 3)); kk = rng.random((4, 3)); kb[:, 2] = kk[:, 2] = 0.0
        M_img = Mwr[:, ws["idx"]][:, :, :, ws["idx"]]
        E = lambda kp, s: ws["w"][None] * np.exp(s * 2j * np.pi * (kp @ ws["R_img"].T))
        off = lambda kbra, kket: np.einsum("ka,waWb,Kb->wkWK", E(kbra, -1), M_img, E(kket, +1), optimize=True)
        assert np.allclose(wi.Mwr_to_Mwk_pairs(Mwr, R, kb, kk, ws=ws), off(kb, kk), rtol=0, atol=tol)
        assert np.allclose(wi.Mwr_to_Mwk(Mwr, R, kb, ws=ws), off(kb, kb), rtol=0, atol=tol)
        assert not np.allclose(wi.Mwr_to_Mwk_pairs(Mwr, R, kb, kk, ws=ws), wi.Mwr_to_Mwk_pairs(Mwr, R, kb, kk), atol=1e-6)


def _random_unitary(nk, n, rng):
    """(nk, n, n) random unitary matrices (QR of complex Gaussians)."""
    Q, _ = np.linalg.qr(rng.standard_normal((nk, n, n)) + 1j * rng.standard_normal((nk, n, n)))
    return Q


def test_defect_mwr_coarse_identity():
    """T5: coarse grid D = N = 4, M = 0, C_N = 1: the box is the whole cell, Mwr = -identity on every (R, R), Mwk = -N_cells delta_kk' identity."""
    rng = np.random.default_rng(20); D = 4; nw = 2; k = lt.mp_grid(D); nk = len(k)
    out = lt.defect_mwr(np.zeros((nw, nk, nw, nk), complex), _random_unitary(nk, nw, rng), None, k, (D, D, 1), D, C_N=1.0)
    nR = len(out["R"])
    assert out["in_box"].dtype == bool and out["in_box"].all() and out["Rn"].shape == out["R"].shape == (nR, 3)
    assert np.allclose(out["Mwr"].reshape(nw * nR, nw * nR), -np.eye(nw * nR), rtol=0, atol=1e-12)
    for R in (out["R"], out["Rn"]):
        assert np.allclose(wi.Mwr_to_Mwk(out["Mwr"], R, k).reshape(nw * nk, nw * nk), -D * D * np.eye(nw * nk), rtol=0, atol=1e-12)


def test_defect_mwr_dense_box():
    """T6: dense grid D = 6, box n_box = N = 3, M = 0, C_N = 1: 9 cells in the box, Mwk[w, k', W, k] = -D_N(k - k') delta_wW,
    D_N(q) = sum_{i, j < N} exp(2 pi i (q1 i + q2 j))."""
    rng = np.random.default_rng(21); D, N = 6, 3; nw = 2; k = lt.mp_grid(D); nk = len(k)
    out = lt.defect_mwr(np.zeros((nw, nk, nw, nk), complex), _random_unitary(nk, nw, rng), None, k, (D, D, 1), N, C_N=1.0)
    assert out["in_box"].sum() == N * N
    box = np.array([[i, j, 0] for i in range(N) for j in range(N)])
    q = k[None, :, :] - k[:, None, :]                                             # q[k', k] = k - k'
    DN = np.exp(2j * np.pi * (q @ box.T)).sum(-1)                                 # (nk', nk)
    ref = -np.einsum("wW,pk->wpWk", np.eye(nw), DN)
    assert np.allclose(wi.Mwr_to_Mwk(out["Mwr"], out["R"], k), ref, rtol=0, atol=1e-12)


def test_defect_mwr_no_shift_is_manual_chain():
    """T7: with C_N = 0, Mwr, R, Rn and R_d are bitwise those of the manual chain Mbk_to_Mwk -> Mwk_to_Mwr -> recenter_mwr."""
    rng = np.random.default_rng(22); nw = 2
    for D, N in ((4, 4), (6, 3)):
        k = lt.mp_grid(D); nk = len(k); MP = (D, D, 1)
        Mbk = rng.standard_normal((nw, nk, nw, nk)) + 1j * rng.standard_normal((nw, nk, nw, nk)); U = _random_unitary(nk, nw, rng)
        out = lt.defect_mwr(Mbk, U, None, k, MP, N)
        Mwr, R = wi.Mwk_to_Mwr(wi.Mbk_to_Mwk(Mbk, U, None), k, MP)
        Rn, R_d = lt.recenter_mwr(Mwr, R, MP)
        for key, ref in (("Mwr", Mwr), ("R", R), ("Rn", Rn), ("R_d", R_d)):
            assert np.array_equal(out[key], ref), key
