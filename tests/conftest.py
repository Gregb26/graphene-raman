"""
Shared fixtures of the electron_photon tests (loaded by pytest for every test file).

Most tests run in two gauges, shift_B = (0,0,0) and (1,0,0), through an indirect parametrization of
the `tb` fixture: gauge-independent properties must hold in both. Reference values and the
derivations behind the tolerances: memoire/EM/EM.md (M0, implementation notes). The real data of M1
come from the 27 x 27 wannierisation tracked in wannier/27x27/ (fixtures w90_dir, tb_w90, eig_w90);
everything the tests know about these data is in W90_REF (fixture w90_ref).

Run with:  .venv/bin/python -m pytest tests -v
"""

import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest
import numpy as np
from electron_defect_interaction.electron_photon import make_graphene_tb, make_grid_tb, make_wannier_tb
from electron_defect_interaction.io.wannier_io import read_w90_mat

N = 100    # the k grid has N x N = 1e4 points: fast, yet covers the whole Brillouin zone

HW = 2.33  # eV, 532 nm laser

# Reference data: the 27 x 27 wannierisation of the thesis. The values of W90_REF belong to these
# exact files (sha256 prefixes, checked by w90_dir): a new reference changes both, here only.
W90_DIR = Path(__file__).resolve().parents[1] / "wannier" / "27x27" # repo root from this file
W90_SHA256 = {"wannier_tb.dat": "baa17b88b1b69e51", "wannier.wout": "3db1203ed17c2558",
              "wannier.eig": "8ed92c792c499689", "wannier_u.mat": "64f31f661b4d5c6c"}
W90_REF = SimpleNamespace(
    n_grid=27, nR=741, nW=5,
    ndegen_counts={1: 717, 2: 24},
    lattice=[[2.135490, -1.232926, 0], [2.135490, 1.232926, 0], [0, 0, 15.875316]], # a1, a2, a3, Angstrom
    sigma=[0, 1, 2], pz=[3, 4],       # Wannier order: sigma bonds, then p_z of C1 and C2
    bands_pi=(3, 4),                  # band indices of pi and pi* near K (ascending energies)
    E_D=-4.238895,                    # eV, Dirac point at K (.eig)
    froz_max=-1.74,                   # eV, top of the frozen window
    nfrozen={4, 5},                   # number of frozen bands per k point
    eig_tol=3e-5,                     # eV, interpolated vs .eig (measured 1.3e-5, use_ws_distance)
    nn_R={(0, 0, 0), (-1, 0, 0), (0, -1, 0)}, # R of the three C1 -> C2 nearest-neighbour hops
    a_cc=1.42366, t_nn=2.9089,        # Angstrom, eV: C1-C2 distance, |H_pzA,pzB| of those hops
    r_defect_R=(2.54e-3, 1.65e-3),    # Angstrom, max |r(R) - r(-R)^dagger|, x y (EM1)
    r_defect_K=(2.31e-3, 4.86e-3),    # Angstrom, max |A(K) - A(K)^dagger|, x y
    rz_sigma_pz=0.220394,             # Angstrom, max |r^z| between sigma and p_z (allowed by the mirror)
    ring_q_over_q0=(0.9, 1.35),       # bounds of q/q0 at 2.33 eV with the default t, a_cc
    hv_F=5.46919,                     # eV*Angstrom, fermi_velocity at q = 1e-3 (both ways)
    hv_pi=(5.47029, 5.46809),         # eV*Angstrom, |hbar v_nn| of pi and pi* at q = 1e-3
    hermiticity={                     # hermiticity_report, grid 60^2: (max_R x y, frob_R, max_k x y)
        'sigma': ((2.5423e-3, 1.6520e-3), 1.1977e-2, (6.9451e-3, 5.8033e-3)),
        'pz':    ((1.2222e-3, 1.3698e-3), 2.5760e-2, (1.2815e-2, 1.4190e-2)),
    },
    hw_froz={100: 4.985446, 200: 4.979360}, # eV, frozen_window_limit(N) at froz_max, mu = E_D (4.9595 at N = 800)
    map_K=((0.054565, 6.211405), (59.669842, 47.864443)), # map_around_K, h = 0.35, nq = 100: deps min max (eV), P max x y
    ring_stats={                      # ring_stats at mu = E_D, 720 angles: (avg x y, node x y in degrees, ratio min max)
        1.96: {'full':         ((1.037521, 1.037438), (6.0, 0.0), (1.0, 1.0)),
               'centres_only': ((1.000591, 1.000504), (6.0, 0.0), (0.980848, 0.982751)),
               'no_berry':     ((1.130745, 1.000504), (-8.0, 0.0), (0.661973, 1.231047))},
        2.33: {'full':         ((1.051244, 1.051130), (7.0, 0.0), (1.0, 1.0)),
               'centres_only': ((0.999513, 0.999395), (7.0, 0.0), (0.972914, 0.976287)),
               'no_berry':     ((1.183452, 0.999395), (-9.5, 0.0), (0.578145, 1.283027))},
        2.54: {'full':         ((1.059446, 1.059382), (7.5, 0.0), (1.0, 1.0)),
               'centres_only': ((0.998396, 0.998331), (8.0, 0.0), (0.967812, 0.972348)),
               'no_berry':     ((1.216999, 0.998331), (-10.5, 0.0), (0.527531, 1.314529))},
    },
)

@pytest.fixture
def tb(request):
    """
    Graphene model, shift_B = (0,0,0) by default. Other gauges through
    @pytest.mark.parametrize("tb", [...], indirect=True), which puts each value in request.param.
    """
    shift = getattr(request, "param", (0,0,0))
    return make_graphene_tb(shift_B=shift)


@pytest.fixture
def grid(tb):
    """
    N x N GridTB of `tb`, rebuilt for every gauge of a parametrized test.
    """
    return make_grid_tb(tb, N)


@pytest.fixture(scope="session")
def w90_dir():
    """Directory of the 27 x 27 reference data, after checking the files are those W90_REF describes."""
    for name, sha in W90_SHA256.items():
        digest = hashlib.sha256((W90_DIR / name).read_bytes()).hexdigest()
        if not digest.startswith(sha):
            pytest.fail(f"reference data changed: {name} has sha256 {digest[:16]}, expected {sha}. "
                        "Update W90_SHA256 and W90_REF in tests/conftest.py.", pytrace=False)
    return W90_DIR


@pytest.fixture(scope="session")
def w90_ref():
    """Everything the tests know about the 27 x 27 reference data (W90_REF)."""
    return W90_REF


@pytest.fixture(scope="session")
def tb_w90(w90_dir):
    """
    WannierTB of the real 27 x 27 data (M1), read once per session.
    """
    return make_wannier_tb(w90_dir / "wannier_tb.dat")


@pytest.fixture(scope="session")
def eig_w90(w90_dir):
    """Coarse-grid DFT eigenvalues: k_red (729, 3) in .eig order (from u.mat), E (729, 20) in eV."""
    eig = np.loadtxt(w90_dir / "wannier.eig") # lines: band, k, energy (eV), band fastest
    nb, nk = int(eig[:, 0].max()), int(eig[:, 1].max())
    _, k_red = read_w90_mat(w90_dir / "wannier_u.mat")
    return k_red, eig[:, 2].reshape(nk, nb)
