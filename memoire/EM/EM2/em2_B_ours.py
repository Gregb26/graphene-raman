"""
EM2, step B (decomposition): our sigma(omega)/sigma_0 (full velocity, package functions only) recomputed with the
conventions of postw90 3.1.0 switched on one at a time, so that postw90 - M4 splits into its causes.

Differences between kubo.py (M4) and postw90 berry_get_kubo_k (read in the 3.1.0 source):
  1. prefactor: M4 weights each transition by 1/hw (delta identity at resonance), postw90 by 1/(eps_c - eps_v)
     (|A_cv|^2 (eps_c - eps_v) = |hbar v_cv|^2 / (eps_c - eps_v)); same eta -> 0 limit, O(eta^2) apart;
  2. grid: M4 = shifted N x N (k_grid shift 0.5), postw90 = Gamma-centred (i/N);
  3. use_ws_distance (postw90 default T) and r(R) rebuilt by postw90 from the .mmn (not read from _tb.dat): not
     reproducible here; what remains after 1 and 2.
Spin g_s = 2 and the Gaussian of std eta are those of M4 (postw90 is run with width sqrt(2) eta and no spin factor).

Runs (all full velocity, eta = 0.04 eV, mu = E_D, hw = 0.20 ... 6.00 by 0.01):
  m4_check   : N = 1200, shift 0.5, 1/hw   -> must equal memoire/EM/M4_sigma/em_sigma_full_N1200_eta0.04.npz
  shift_delta: N = 1200, shift 0.5, 1/Delta
  gam1201    : N = 1201, shift 0,   1/Delta (grid and formula of postw90 berry_kmesh = 1201)
  gam301     : N = 301,  shift 0,   1/Delta (same for berry_kmesh = 301)
Writes em2_B_ours.npz and prints the pairwise differences. Run as a job (~1.5 M k points per N = 1200 run).
"""

import os
import time
from pathlib import Path

import numpy as np
from electron_defect_interaction.electron_photon import compute_velocity, make_wannier_tb
from electron_defect_interaction.electron_photon.kgrid import k_grid, reciprocal
from electron_defect_interaction.electron_photon.kubo import gaussian_eta, kubo_normalize

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
E_D, ETA = -4.238895, 0.04
HW = np.arange(20, 601) / 100
CHUNK = int(5e4)


def sigma_two_prefactors(tb, N, shift):
    """sigma/sigma_0 (nw, 3, 3) with the 1/hw prefactor (as kubo.py) and with the 1/Delta prefactor (as postw90)."""
    B = reciprocal(tb.lattice)
    _, k_cart, _ = k_grid(B, N, shift=shift)
    nk = k_cart.shape[0]
    S_w = np.zeros((len(HW), 9)); S_d = np.zeros((len(HW), 9))
    for b in range(0, nk, CHUNK):
        _, eps, _, hv = compute_velocity(tb, k_cart[b:b + CHUNK], "berry")
        occ = eps < E_D                                           # same selection as kubo_accumulate
        ik, ic, iv = np.nonzero((~occ)[:, :, None] & occ[:, None, :])
        gaps = eps[ik, ic] - eps[ik, iv]
        hv_cv = hv[ik, :, ic, iv]
        W = np.real(np.conj(hv_cv[:, :, None]) * hv_cv[:, None, :]).reshape(-1, 9)
        G = gaussian_eta(HW[:, None] - gaps[None, :], ETA)        # (nw, nt)
        S_w += G @ W
        S_d += (G / gaps[None, :]) @ W
    A_cell = np.linalg.norm(np.cross(tb.lattice[:, 0], tb.lattice[:, 1]))
    sig_w = kubo_normalize(S_w.reshape(-1, 3, 3), HW, nk, A_cell)
    sig_d = kubo_normalize(S_d.reshape(-1, 3, 3) * HW[:, None, None], HW, nk, A_cell)   # hw cancels the 1/hw
    return sig_w, sig_d


tb = make_wannier_tb(REPO / "wannier" / "27x27" / "wannier_tb.dat")
res = {}
t0 = time.time()
res["m4_check"], res["shift_delta"] = sigma_two_prefactors(tb, 1200, 0.5)
print(f"N = 1200 shifted: {time.time() - t0:.0f} s", flush=True)
t0 = time.time()
_, res["gam1201"] = sigma_two_prefactors(tb, 1201, 0.0)
print(f"N = 1201 Gamma: {time.time() - t0:.0f} s", flush=True)
_, res["gam301"] = sigma_two_prefactors(tb, 301, 0.0)
np.savez(HERE / "em2_B_ours.npz", hw=HW, eta=ETA, mu=E_D, **res)

m4 = np.load(REPO / "memoire" / "EM" / "M4_sigma" / "em_sigma_full_N1200_eta0.04.npz")["sigma"]
iL = [int(np.flatnonzero(HW == e)[0]) for e in (1.96, 2.33, 2.54)]
def report(name, a, b):
    d = a - b
    k = int(np.argmax(np.abs(d[:, 0, 0])))
    print(f"{name:32s} xx max |d| {np.abs(d[:, 0, 0]).max():.2e} at {HW[k]:.2f} eV; yy {np.abs(d[:, 1, 1]).max():.2e}; "
          f"xy {np.abs(d[:, 0, 1]).max():.1e}; at lasers xx " + ", ".join(f"{d[i, 0, 0]:+.2e}" for i in iL)
          + f"; at 0.2 eV {d[0, 0, 0]:+.1e}")
report("m4_check - M4 npz", res["m4_check"], m4)
report("shift_delta - M4 (prefactor)", res["shift_delta"], m4)
report("gam1201 - shift_delta (grid)", res["gam1201"], res["shift_delta"])
report("gam301 - gam1201 (N)", res["gam301"], res["gam1201"])
print("DONE")
