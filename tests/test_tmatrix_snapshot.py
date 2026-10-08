"""Non-regression snapshot of chain T before the 2026-09-30 cleanup (cluster_tmatrix, ex local_tmatrix; Mbk_to_Mwk / Mwk_to_Mwr through
defect_mwr, pole_criterion.cluster_t_cache). Synthetic case: make_graphene_tb in the 5-WF layout + a localized Hermitian defect potential put on a
6x6 grid away from the origin, scrambled into a Bloch M with random U and U_dis. The input is built with plain numpy (no
function under test). Reference = regression, code of commit c0fb4d2, stored in tests/snapshots/tmatrix_chain.npz.
Mwr is kept as its on-site blocks and its column towards R0 (what the alignment and mwr_locality touch). Quantities
indexed on the cluster (M_cluster, g0, t) depend on the order of R_cluster; Gamma and the LDOS do not.
Second file, tests/snapshots/wannier_ft.npz: Mwr_to_Mwk (plain, Wigner-Seitz), Mwr_to_Mwk_pairs, Mwk_to_Mbk on random inputs.
Write missing snapshots: .venv/bin/python tests/test_tmatrix_snapshot.py --write (add --force to overwrite, only on purpose)"""
import sys
from pathlib import Path

import numpy as np
import pytest

from graphene_raman.wannier.wannier_hamiltonian import Hwr_to_Hwk
from graphene_raman.wannier import wannier_interpolation as wi
from graphene_raman.defects.many_body import cluster_tmatrix as ct
from graphene_raman.defects.many_body import pole_criterion as pc
from graphene_raman.electron_photon import make_graphene_tb

SNAP = Path(__file__).parent / "snapshots" / "tmatrix_chain.npz"
SNAP_FT = Path(__file__).parent / "snapshots" / "wannier_ft.npz"
D, NW, NB = 6, 5, 7                  # MP grid D x D, Wannier functions, Bloch bands (disentangled)
N_BOX, C_N = 3, -0.025               # supercell box and alignment constant (eV)
R_D = np.array([1, 2, 0])            # defect cell in the raw labels: inside the N_BOX box, off the origin (real recentering)
ETA, R_CUT = 0.1, 1


def _five_wf(e_sigma):
    """make_graphene_tb in the 5-WF layout of the production: sigma 0-2 decoupled at e_sigma, p_z(A), p_z(B) = 3, 4. Also returns the row of R = 0."""
    tb = make_graphene_tb()
    i0 = tb.index[(0, 0, 0)]
    Hwr = np.zeros((len(tb.R_int), NW, NW), complex)
    Hwr[:, 3:5, 3:5] = tb.H_R
    Hwr[i0, [0, 1, 2], [0, 1, 2]] = e_sigma
    return Hwr, tb.R_int, tb.ndegen, i0


def _synthetic_Mbk(rng):
    """Bloch-gauge M (NB, nk, NB, nk) on a permuted D x D grid, with its U (nk, NW, NW) and U_dis (nk, NB, NW)."""
    r = np.arange(D) - D // 2
    R = np.array([[i, j, 0] for i in r for j in r])
    d = ((R - R_D + D // 2) % D) - D // 2
    f = np.exp(-np.linalg.norm(d, axis=1))                            # envelope centred on R_D (periodic distance)
    n = NW * len(R)
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    X = 0.5 * (X + X.conj().T)
    Mwr = X.reshape(NW, len(R), NW, len(R)) * f[None, :, None, None] * f[None, None, None, :]
    i_d = int(np.flatnonzero((R == R_D).all(axis=1))[0])
    Mwr[3, i_d, 3, i_d] += 5.0                                        # strong pz on-site at the defect

    k = np.array([[i / D, j / D, 0.0] for i in range(D) for j in range(D)])[rng.permutation(D * D)]
    bra = np.exp(-2j * np.pi * k @ R.T)
    Mwk = np.einsum("kr,wrWR,KR->wkWK", bra, Mwr, bra.conj())

    U = np.linalg.qr(rng.normal(size=(len(k), NW, NW)) + 1j * rng.normal(size=(len(k), NW, NW)))[0]
    U_dis = np.linalg.qr(rng.normal(size=(len(k), NB, NW)) + 1j * rng.normal(size=(len(k), NB, NW)))[0]
    V = U_dis @ U
    Mbk = np.einsum("kbw,wkWK,KBW->bkBK", V, Mwk, V.conj())
    return Mbk, U, U_dis, k


def _chain():
    rng = np.random.default_rng(20260930)
    Mbk, U, U_dis, k = _synthetic_Mbk(rng)

    d = ct.defect_mwr(Mbk, U, U_dis, k, (D, D, 1), n_box=N_BOX, C_N=C_N)
    Mwr, Rn = d["Mwr"], d["Rn"]
    dist, wt = ct.mwr_locality(Mwr, Rn)
    R_cluster = ct.cluster_cells(Rn, R_CUT)
    M_cluster, herm = ct.cluster_potential(Mwr, Rn, R_cluster)

    Hwr, Rw, nd, _ = _five_wf(e_sigma=-4.0)
    k_int = ct.mp_grid(24)
    Hk_int = Hwr_to_Hwk(Hwr, Rw, k_int, ndegen=nd)[0]
    egrid = np.linspace(-3.0, 3.0, 5)
    g0 = ct.cluster_green_batch(Hk_int, k_int, R_cluster, egrid, ETA, k_chunk=100, e_chunk=2)
    dg0 = ct.cluster_green_batch(Hk_int, k_int, R_cluster, egrid, ETA, k_chunk=100, e_chunk=2, deriv=True)
    g0_single = ct.cluster_green(Hk_int, k_int, R_cluster, 0.7, ETA)
    t = np.array([ct.cluster_t(M_cluster, g) for g in g0])
    i0 = int(np.flatnonzero((R_cluster == 0).all(axis=1))[0])
    rho, rho0 = ct.cluster_ldos(g0, M_cluster, [i0 * NW + 3, i0 * NW + 4])

    k_out = ct.mp_grid(4)
    gamma = ct.scattering_rate(Hwr, Rw, nd, M_cluster, R_cluster, k_out, ETA, k_int=k_int)
    gamma_fast = ct.scattering_rate_fast(Hwr, Rw, nd, M_cluster, R_cluster, k_out, ETA, k_int=k_int, e_window=(-3.0, 3.0), ne_per_eta=4)

    i0_all = int(np.flatnonzero((Rn == 0).all(axis=1))[0])
    return dict(Mwr_onsite=np.einsum("wrWr->wWr", Mwr), Mwr_to_R0=Mwr[:, :, :, i0_all], R=d["R"], Rn=Rn, R_d=d["R_d"], in_box=d["in_box"], dist=dist, wt=wt, R_cluster=R_cluster, M_cluster=M_cluster,
                herm=np.array(herm), g0=g0, dg0=dg0, g0_single=g0_single, t=t, rho=rho, rho0=rho0,
                gamma=gamma, gamma_fast=gamma_fast)


def _wannier_ft():
    """Mwr_to_Mwk (plain and with Wigner-Seitz images), Mwr_to_Mwk_pairs, Mwk_to_Mbk; random Mwr / Mwk, off-grid k."""
    rng = np.random.default_rng(20261001)
    r = np.arange(D) - D // 2
    R = np.array([[i, j, 0] for i in r for j in r])
    Mwr = rng.normal(size=(NW, len(R), NW, len(R))) + 1j * rng.normal(size=(NW, len(R), NW, len(R)))
    k = rng.random((7, 3)); k[:, 2] = 0.0
    k_bra = rng.random((3, 3)); k_bra[:, 2] = 0.0
    a = 2.4659
    A_cols = np.array([[a * np.sqrt(3) / 2, a * np.sqrt(3) / 2, 0.0], [-a / 2, a / 2, 0.0], [0.0, 0.0, 8.0]])
    ws = wi.ws_images(R, R_D, (D, D, 1), A_cols)
    Hwr, Rw, nd, i0 = _five_wf(e_sigma=-4.0)
    Hwr[i0, 0, 0] -= 1.0; Hwr[i0, 1, 1] -= 2.0                       # lift the sigma degeneracy (eigenvectors fixed up to LAPACK's phase)
    Mwk = rng.normal(size=(NW, len(k), NW, len(k))) + 1j * rng.normal(size=(NW, len(k), NW, len(k)))
    return dict(Mwk=wi.Mwr_to_Mwk(Mwr, R, k), Mwk_ws=wi.Mwr_to_Mwk(Mwr, R, k, ws=ws),
                Mwk_pairs=wi.Mwr_to_Mwk_pairs(Mwr, R, k_bra, k), Mwk_pairs_ws=wi.Mwr_to_Mwk_pairs(Mwr, R, k_bra, k, ws=ws),
                Mbk=wi.Mwk_to_Mbk(Mwk, Hwr, Rw, k, ndegen=nd))


@pytest.fixture(scope="module")
def chain():
    return _chain()


@pytest.fixture(scope="module")
def ref():
    with np.load(SNAP) as z:
        return {key: z[key] for key in z.files}


def test_same_keys(chain, ref):
    assert sorted(chain) == sorted(ref)


@pytest.mark.parametrize("key", ["R", "Rn", "R_d", "in_box", "R_cluster"])
def test_labels_exact(chain, ref, key):
    assert np.array_equal(chain[key], ref[key])


@pytest.mark.parametrize("key", ["Mwr_onsite", "Mwr_to_R0", "dist", "wt", "M_cluster", "g0", "dg0", "g0_single", "t", "rho", "rho0",
                                 "gamma", "gamma_fast"])
def test_values(chain, ref, key):
    a, b = np.asarray(chain[key]), ref[key]
    assert a.shape == b.shape
    scale = np.nanmax(np.abs(b))
    assert np.allclose(a, b, rtol=1e-10, atol=1e-12 * scale, equal_nan=True), np.nanmax(np.abs(a - b)) / scale


def test_cluster_t_cache_equals_cluster_t(chain):
    """pole_criterion.cluster_t_cache is the batched cluster_t (two implementations of t, item e of the cleanup)."""
    tc = pc.cluster_t_cache(chain["M_cluster"], chain["g0"])
    assert np.allclose(tc, chain["t"], rtol=1e-10, atol=1e-12 * np.abs(chain["t"]).max())


def test_cluster_cells_matches_inline_truncation():
    """cluster_cells = the inline truncation of the scripts (reduced-label norm <= R_cut + 1e-9, same order); cell counts by enumeration."""
    r = np.arange(27) - 27 // 2
    Rn = np.array([[i, j, 0] for i in r for j in r])
    for rc in (0, 1, 2, 2.5, 3, np.sqrt(8), 5):
        assert np.array_equal(ct.cluster_cells(Rn, rc), Rn[np.linalg.norm(Rn, axis=1) <= rc + 1e-9])
    counts = {rc: sum(i * i + j * j <= rc * rc for i in range(-6, 7) for j in range(-6, 7)) for rc in (1, 2, 3)}
    assert counts == {1: 5, 2: 13, 3: 29}
    assert all(len(ct.cluster_cells(Rn, rc)) == n for rc, n in counts.items())


def test_herm_residual(chain, ref):
    assert abs(float(chain["herm"]) - float(ref["herm"])) <= 1e-12


@pytest.fixture(scope="module")
def ft():
    with np.load(SNAP_FT) as z:
        return _wannier_ft(), {key: z[key] for key in z.files}


@pytest.mark.parametrize("key", ["Mwk", "Mwk_ws", "Mwk_pairs", "Mwk_pairs_ws", "Mbk"])
def test_wannier_ft(ft, key):
    new, ref_ft = ft
    a, b = np.asarray(new[key]), ref_ft[key]
    assert a.shape == b.shape
    scale = np.abs(b).max()
    assert np.allclose(a, b, rtol=1e-10, atol=1e-12 * scale), np.abs(a - b).max() / scale


if __name__ == "__main__" and "--write" in sys.argv:
    SNAP.parent.mkdir(exist_ok=True)
    for path, build in ((SNAP, _chain), (SNAP_FT, _wannier_ft)):
        if path.exists() and "--force" not in sys.argv:
            print(f"kept {path}")
            continue
        np.savez_compressed(path, **build())
        print(f"wrote {path}")
