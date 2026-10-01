"""
M4 production (campagnes/EM/EM.md, section 8) on the 27 x 27 data, with the package functions only:
sigma(omega)/sigma_0 of the three velocity variants at the couple retained by the pilot (N = 1200,
eta = 0.04 eV, see pilote/), the maps around K (map_around_K, F19), the ring statistics (ring_stats,
F15), hbar v_F (fermi_velocity, F13) and hbar omega_froz (frozen_window_limit, F18).

Run from anywhere:  .venv/bin/python campagnes/EM/M4_sigma/m4_prod.py   (outputs next to this script)
"""

import hashlib
import json
import subprocess
import time
from datetime import date
from pathlib import Path

import numpy as np
from graphene_raman.electron_photon import (centres_only, fermi_velocity, frozen_window_limit,
                                                         make_grid_tb, make_wannier_tb, map_around_K,
                                                         ring_stats, sigma_on_grid)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]                   # campagnes/EM/M4_sigma/ -> repo root
E_D = -4.238895                          # eV, Dirac point (.eig at K, W90_REF of tests/conftest.py)
FROZ_MAX = -1.74                         # eV, dis_froz_max of wannier.win
N, ETA = 1200, 0.04                      # retained by the pilot (pilote/sweep.log, EM.md section 8)
HW = np.arange(20, 601) / 100            # eV, 0.20 ... 6.00, exact decimals (lasers on the grid)
LASERS = [1.96, 2.33, 2.54]              # eV, 633, 532, 488 nm
H_MAP, NQ_MAP = 0.35, 300                # 1/Angstrom, points per axis (beyond the 2.54 eV ring)
N_FROZ = 800                             # grid of frozen_window_limit (converged to 1e-3 eV)
HBAR = 6.582119569e-16                   # eV s

t_start = time.time()
tb = make_wannier_tb(REPO / "results" / "wannier" / "27x27" / "wannier_tb.dat")
K = make_grid_tb(tb, 3).K
variants = {"full": (tb, "berry"), "centres_only": (centres_only(tb), "berry"), "no_berry": (tb, "no_berry")}

# provenance: commit, and sha256 of the modules actually used (the tree may hold uncommitted changes)
src = REPO / "src" / "graphene_raman"
modules = sorted((src / "electron_photon").glob("*.py")) + [src / "io" / "wannier_io.py"]
provenance = {
    "date": date.today().isoformat(),
    "commit": subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip(),
    "src_dirty": bool(subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "src"], capture_output=True, text=True).stdout.strip()),
    "sha256": {str(p.relative_to(REPO)): hashlib.sha256(p.read_bytes()).hexdigest()[:16] for p in modules},
    "data_sha256": hashlib.sha256((REPO / "results" / "wannier" / "27x27" / "wannier_tb.dat").read_bytes()).hexdigest()[:16],
}

# sigma(omega) / sigma_0, one file per variant
iL = [int(np.flatnonzero(HW == e)[0]) for e in LASERS]
sigma_lasers = {}
for name, (model, mode) in variants.items():
    t0 = time.time()
    sigma = sigma_on_grid(model, N, HW, mu=E_D, eta=ETA, mode=mode, chunk=int(5e4)) # (nw, 3, 3)
    np.savez(HERE / f"em_sigma_{name}_N{N}_eta{ETA}.npz", hw=HW, sigma=sigma, N=N, eta=ETA, mu=E_D, variant=name)
    sigma_lasers[name] = {str(e): {"xx": sigma[i, 0, 0], "yy": sigma[i, 1, 1], "xy": sigma[i, 0, 1]} for e, i in zip(LASERS, iL)}
    print(f"sigma {name:13s} {time.time() - t0:6.1f} s | xx at lasers {np.round(sigma[iL, 0, 0], 5)}", flush=True)

# maps around K (panel (a)): Delta eps is the same for the three variants
maps = {}
for name, (model, mode) in variants.items():
    deps, P = map_around_K(model, K, E_D, H_MAP, NQ_MAP, mode=mode)
    maps[f"P_{name}"] = P
np.savez_compressed(HERE / "em_map_K.npz", q=np.linspace(-H_MAP, H_MAP, NQ_MAP), deps=deps, **maps,
                    K=K, half_width=H_MAP, nq=NQ_MAP, mu=E_D)
print(f"maps done, deps {deps.min():.5f}-{deps.max():.3f} eV", flush=True)

# ring statistics (F15)
stats = ring_stats(tb, K, E_D, LASERS)
stats_json = {str(float(hw)): {v: {key: np.asarray(val, dtype=float).tolist() for key, val in res.items()}
                               for v, res in d.items()} for hw, d in stats.items()}
(HERE / "em_ring_stats.json").write_text(json.dumps(stats_json, indent=1))

# scalars: hbar v_F, hbar omega_froz, sigma at the lasers
fv = fermi_velocity(tb, K, E_D)
hw_froz = float(frozen_window_limit(tb, N_FROZ, FROZ_MAX, E_D))
scalars = {
    "hv_F_eV_A": {"inter": float(fv["inter_avg"]), "intra": float(fv["intra_avg"])},
    "v_F_m_s": float(fv["inter_avg"]) * 1e-10 / HBAR,
    "hw_froz_eV": hw_froz, "N_froz": N_FROZ, "froz_max_eV": FROZ_MAX,
    "E_D_eV": E_D, "N": N, "eta_eV": ETA, "sigma_over_sigma0_at_lasers": sigma_lasers,
    "provenance": provenance,
}
(HERE / "em_scalars.json").write_text(json.dumps(scalars, indent=1, default=float))
print(f"hbar v_F = {fv['inter_avg']:.5f} eV A, hbar omega_froz = {hw_froz:.4f} eV", flush=True)

# consistency with the pilot (different code path: three eta per diagonalization)
pilot = np.load(HERE / "pilote" / f"em_sigma_full_N{N}.npz")
i_eta = int(np.flatnonzero(np.isclose(pilot["eta"], ETA))[0])
full = np.load(HERE / f"em_sigma_full_N{N}_eta{ETA}.npz")["sigma"]
print(f"production vs pilot (full, N = {N}, eta = {ETA}): max |diff| = {np.abs(full - pilot['sigma'][i_eta]).max():.1e}")
print(f"DONE in {time.time() - t_start:.0f} s", flush=True)
