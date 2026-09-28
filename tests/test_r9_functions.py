"""Minimal tests of the R9 functions (2026-09-28): P1 atom_sphere_shifts, P2 Mwr_to_Mwk_pairs, P3 cluster_ldos, P4 ldos_from_eigenpairs.
Run: PYTHONPATH=src pytest tests/test_r9_functions.py"""
import numpy as np
from electron_defect_interaction.defects import alignment as al
from electron_defect_interaction.wannier.wannier_interpolation import Mwr_to_Mwk, Mwr_to_Mwk_pairs
from electron_defect_interaction.defects.many_body.local_tmatrix import cluster_ldos
from electron_defect_interaction.wannier.supercell_fold import ldos_from_eigenpairs


def _graphene_supercell(N, a=2.4659, c=8.0):
    """N x N graphene cell (a1, a2 at 60 degrees, as the QE cells), A sites at (i+1/3, j+1/3)/N, B at (i+2/3, j+2/3)/N."""
    A = np.array([[a * N * np.sqrt(3) / 2, a * N * np.sqrt(3) / 2, 0.0], [-a * N / 2, a * N / 2, 0.0], [0.0, 0.0, c]])
    x = [((i + f) / N, (j + f) / N, 0.5) for i in range(N) for j in range(N) for f in (1 / 3, 2 / 3)]
    return A, np.array(x)


def test_atom_sphere_shifts_reproduces_far_atom_alignment():
    rng = np.random.default_rng(0)
    A, x_p = _graphene_supercell(3)
    x_d = np.delete(x_p, 8, axis=0)                                       # vacancy = pristine atom 8
    n = (36, 36, 32)
    g = np.stack(np.meshgrid(*[np.arange(m) / m for m in n], indexing="ij"), -1)
    V_p = np.cos(2 * np.pi * g[..., 0]) + 0.3 * np.sin(2 * np.pi * (g[..., 1] + g[..., 2]))
    V_d = V_p + 0.05 * rng.standard_normal(n)                              # arbitrary difference
    r = al.atom_sphere_shifts(V_d, V_p, x_d, x_p, A, 1.0)
    ref = al.far_atom_alignment(V_d, V_p, x_d, x_p, A, (1.0,))
    i = ref["i_far_d"]
    assert r["i_far_axis"] == i and np.allclose(r["s_vac"], ref["s_vac"])
    assert abs(r["shift"][i] - ref["shifts"][1.0]) < 1e-12 * max(1.0, abs(ref["shifts"][1.0]))
    assert r["npts_d"][i] == ref["means"][1.0][2] and r["npts_p"][i] == ref["means"][1.0][3]
    # a constant difference gives the same shift around every atom; true distances never exceed the per-axis ones
    r2 = al.atom_sphere_shifts(V_p + 0.025, V_p, x_d, x_p, A, 1.0)
    assert np.allclose(r2["shift"], 0.025, atol=1e-14)
    assert np.all(r["dist"] <= r["dist_axis"] + 1e-12)


def test_Mwr_to_Mwk_pairs_matches_square_version():
    rng = np.random.default_rng(1)
    nw, R = 2, np.array([[i, j, 0] for i in range(-1, 2) for j in range(-1, 2)])
    Mwr = rng.standard_normal((nw, len(R), nw, len(R))) + 1j * rng.standard_normal((nw, len(R), nw, len(R)))
    k = rng.random((4, 3)); k[:, 2] = 0.0
    sq = Mwr_to_Mwk(Mwr, R, k)
    assert np.allclose(Mwr_to_Mwk_pairs(Mwr, R, k, k), sq, rtol=0, atol=1e-12)
    rect = Mwr_to_Mwk_pairs(Mwr, R, k[:3], k[3:])                        # (nw, 3, nw, 1) = sub-block of the square array
    assert rect.shape == (nw, 3, nw, 1) and np.allclose(rect, sq[:, :3, :, 3:], rtol=0, atol=1e-12)


def _random_hermitian(n, rng):
    X = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return 0.5 * (X + X.conj().T)


def test_cluster_ldos_is_exact_dyson():
    """Lattice = 12 orbitals, cluster = first 4: G = g0 + g0 T g0 on the cluster must equal the cluster block of (z - H0 - V)^-1."""
    rng = np.random.default_rng(2)
    H0 = _random_hermitian(12, rng); Vc = _random_hermitian(4, rng)
    V = np.zeros((12, 12), complex); V[:4, :4] = Vc
    eta = 0.05; egrid = np.linspace(-3, 3, 7)
    g0 = np.array([np.linalg.inv((e + 1j * eta) * np.eye(12) - H0)[:4, :4] for e in egrid])
    G = np.array([np.linalg.inv((e + 1j * eta) * np.eye(12) - H0 - V)[:4, :4] for e in egrid])
    rho, rho0 = cluster_ldos(g0, Vc, [0, 2, 3])
    assert np.allclose(rho, -np.diagonal(G, axis1=1, axis2=2)[:, [0, 2, 3]].imag / np.pi, atol=1e-12)
    assert np.allclose(rho0, -np.diagonal(g0, axis1=1, axis2=2)[:, [0, 2, 3]].imag / np.pi, atol=1e-12)
    r1, _ = cluster_ldos(g0[3], Vc, [1])                                  # single energy
    assert r1.shape == (1,) and np.isclose(r1[0], -G[3, 1, 1].imag / np.pi)


def test_ldos_from_eigenpairs_equals_resolvent():
    rng = np.random.default_rng(3)
    H = _random_hermitian(10, rng); e, v = np.linalg.eigh(H)
    eta = 0.1; egrid = np.linspace(-4, 4, 9)
    rho = ldos_from_eigenpairs(e, v, [0, 5], egrid, eta)
    ref = np.array([-np.diagonal(np.linalg.inv((E + 1j * eta) * np.eye(10) - H))[[0, 5]].imag / np.pi for E in egrid])
    assert rho.shape == (9, 2) and np.allclose(rho, ref, atol=1e-12)
