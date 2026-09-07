#!/usr/bin/env python
"""
analyze_Ved.py -- checks of V_ed^L = V_d - V_p (local defect potential), existing pp.x data only:
  1. 2D map in the graphene plane (z = z_C exactly, grid plane), 5x5 and 9x9
  2. profile along a1 WITH and WITHOUT the G=0 (3D mean) subtraction; 3D mean and in-plane mean
  3. azimuthal average in the plane vs distance to the vacancy, four N
Boundary values (max |V| on the half-box planes along a1/a2, all z and in-plane only), raw and mean-subtracted.
Output: results/M/ved_analysis.npz
"""
import numpy as np
from scipy.ndimage import map_coordinates
from electron_defect_interaction.io import qe_io
from electron_defect_interaction.config import HA2EV
BOHR = 0.529177210903; DATA = "data/graphene"; out = {}
import sys
SIZES = sys.argv[1].split(",") if len(sys.argv) > 1 else ["5x5", "7x7", "8x8", "9x9"]
if len(sys.argv) > 2 and sys.argv[2] == "--merge": out.update(dict(np.load("results/M/ved_analysis.npz")))
for S in SIZES:
    scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; scd = f"{DATA}/supercell/qe/defect_{S}_d.save"
    A, _ = qe_io.get_A_volume(scd); xp = np.mod(qe_io.get_x_red(scp), 1.0); xd = np.mod(qe_io.get_x_red(scd), 1.0)
    dmin = np.array([np.min(np.linalg.norm(np.mod(xd - p + 0.5, 1) - 0.5, axis=1)) for p in xp]); s_vac = xp[int(np.argmax(dmin))]
    Vp, _ = qe_io.get_pot(f"{scp}/Vks_{S}_p", subtract_mean=False, to_hartree=True); Vd, _ = qe_io.get_pot(f"{scd}/Vks_{S}_d", subtract_mean=False, to_hartree=True)
    dV = (Vd - Vp).transpose(2, 1, 0) * HA2EV; nr = np.array(dV.shape); del Vp, Vd
    mean3d = float(dV.mean()); iz = int(np.round(s_vac[2] * nr[2])) % nr[2]; plane = dV[:, :, iz]; mean_plane = float(plane.mean())
    zC = s_vac[2] * A[2, 2] * BOHR; zgrid = iz / nr[2] * A[2, 2] * BOHR
    # in-plane Cartesian coordinates (A) relative to the vacancy, minimum image
    i = np.arange(nr[0]); j = np.arange(nr[1]); s1 = (i / nr[0] - s_vac[0] + 0.5) % 1 - 0.5; s2 = (j / nr[1] - s_vac[1] + 0.5) % 1 - 0.5
    S1, S2 = np.meshgrid(s1, s2, indexing="ij"); X = (S1 * A[0, 0] + S2 * A[0, 1]) * BOHR; Y = (S1 * A[1, 0] + S2 * A[1, 1]) * BOHR
    R = np.sqrt(X ** 2 + Y ** 2)
    # 1. map (roll so that the vacancy is centred)
    i0 = int(np.round(s_vac[0] * nr[0])); j0 = int(np.round(s_vac[1] * nr[1])); sh = (nr[0] // 2 - i0, nr[1] // 2 - j0)
    # continuous coordinates of the rolled (vacancy-centred) grid: no wrap-around cells in the mesh
    ic = (np.arange(nr[0]) - nr[0] // 2) / nr[0]; jc = (np.arange(nr[1]) - nr[1] // 2) / nr[1]; IC, JC = np.meshgrid(ic, jc, indexing="ij")
    out[f"{S}_map"] = np.roll(plane, sh, axis=(0, 1)); out[f"{S}_map_X"] = (IC * A[0, 0] + JC * A[0, 1]) * BOHR; out[f"{S}_map_Y"] = (IC * A[1, 0] + JC * A[1, 1]) * BOHR
    # 2. profile along a1 (in-plane, z = z_C plane), raw; subtracted = raw - mean3d (what get_pot(subtract_mean=True) does)
    t = np.linspace(0, 0.5, 401); s = s_vac[None, :] + t[:, None] * np.array([1.0, 0, 0])[None, :]
    line = map_coordinates(dV, (s * nr[None, :]).T, order=1, mode="wrap"); a1_len = np.linalg.norm(A[:, 0]) * BOHR
    # 3. azimuthal average in the plane
    rmax = 0.5 * a1_len; edges = np.arange(0, rmax + 0.05, 0.05); m = R.ravel() < rmax
    cnt, _ = np.histogram(R.ravel()[m], edges); sm, _ = np.histogram(R.ravel()[m], edges, weights=plane.ravel()[m])
    rad = np.where(cnt > 0, sm / np.maximum(cnt, 1), np.nan); rc = 0.5 * (edges[1:] + edges[:-1])
    # boundary planes (half box along a1 and a2), all z and in-plane only
    i1 = int(np.round((s_vac[0] + 0.5) * nr[0])) % nr[0]; i2 = int(np.round((s_vac[1] + 0.5) * nr[1])) % nr[1]
    b_all = max(np.abs(dV[i1]).max(), np.abs(dV[:, i2]).max()); b_pl = max(np.abs(plane[i1]).max(), np.abs(plane[:, i2]).max())
    b_all_s = max(np.abs(dV[i1] - mean3d).max(), np.abs(dV[:, i2] - mean3d).max()); b_pl_s = max(np.abs(plane[i1] - mean3d).max(), np.abs(plane[:, i2] - mean3d).max())
    site = float(map_coordinates(dV, (s_vac * nr)[:, None], order=1, mode="wrap")[0])
    print(f"[{S}] z_C = {zC:.3f} A (grid plane {zgrid:.3f} A, iz={iz}); <V_ed>_3D = {mean3d*1e3:+.2f} meV, <V_ed>_plane = {mean_plane*1e3:+.2f} meV; site {site:+.2f} eV; "
          f"a1 end (x=1): raw {line[-1]*1e3:+.1f} meV, subtracted {(line[-1]-mean3d)*1e3:+.1f} meV; boundary max|V| (all z): raw {b_all*1e3:.1f} / sub {b_all_s*1e3:.1f} meV; "
          f"(plane): raw {b_pl*1e3:.1f} / sub {b_pl_s*1e3:.1f} meV; plane min {plane.min():+.3f} eV at r={R.ravel()[np.argmin(plane.ravel())]:.2f} A; "
          f"radial: r=1.42 A -> {np.interp(1.42, rc, np.nan_to_num(rad))*1e3:+.1f} meV, r=3a={3*a1_len/int(S.split('x')[0]):.2f} A -> {np.interp(3*a1_len/int(S.split('x')[0]), rc, np.nan_to_num(rad))*1e3:+.1f} meV, r=half box -> {np.nanmean(rad[-3:])*1e3:+.1f} meV", flush=True)
    out.update(**{f"{S}_x": t / 0.5, f"{S}_line_raw": line, f"{S}_line_sub": line - mean3d, f"{S}_mean3d": mean3d, f"{S}_mean_plane": mean_plane, f"{S}_site": site,
                  f"{S}_rc": rc, f"{S}_rad": rad, f"{S}_a": a1_len / int(S.split('x')[0]), f"{S}_halfbox": rmax, f"{S}_b_all": b_all, f"{S}_b_all_sub": b_all_s, f"{S}_b_pl": b_pl, f"{S}_b_pl_sub": b_pl_s, f"{S}_zC": zC, f"{S}_zgrid": zgrid})
    del dV
np.savez("results/M/ved_analysis.npz", **out); print("saved results/M/ved_analysis.npz")

# ---------------- anomaly at r = 1.42 A: sublattice, grid alignment, neighbour displacements, core-masked radial profile
print("\n=== anomaly checks ===", flush=True)
out2 = dict(np.load("results/M/ved_analysis.npz"))
for S in SIZES:
    N = int(S.split('x')[0]); scp = f"{DATA}/supercell/qe/defect_{S}_p.save"; scd = f"{DATA}/supercell/qe/defect_{S}_d.save"
    A, _ = qe_io.get_A_volume(scd); xp = np.mod(qe_io.get_x_red(scp), 1.0); xd = np.mod(qe_io.get_x_red(scd), 1.0)
    dmin = np.array([np.min(np.linalg.norm(np.mod(xd - p + 0.5, 1) - 0.5, axis=1)) for p in xp]); iv = int(np.argmax(dmin)); s_vac = xp[iv]
    s_uc = np.mod(N * s_vac, 1.0); sub = "A" if np.allclose(s_uc[:2], 1/3, atol=1e-3) else ("B" if np.allclose(s_uc[:2], 2/3, atol=1e-3) else "?")
    Vd, _ = qe_io.get_pot(f"{scd}/Vks_{S}_d", subtract_mean=False, to_hartree=True); nr = np.array(Vd.shape[::-1]); del Vd
    g = s_vac * nr; gn = np.round(g); dgrid = (gn - g) / nr; dcart = np.linalg.norm((A @ dgrid) * BOHR)
    # neighbours: pristine atoms within 1.6 A of the site (minimum image), their positions in p and d and displacement
    def cart(s): return (A @ s) * BOHR
    dp = np.mod(xp - s_vac + 0.5, 1) - 0.5; rp = np.linalg.norm((A @ dp.T).T * BOHR, axis=1); nn = [i for i in np.argsort(rp) if 0.5 < rp[i] < 1.6]
    disp = []
    for i in nn:
        dd = np.mod(xd - xp[i] + 0.5, 1) - 0.5; j = int(np.argmin(np.linalg.norm(dd, axis=1))); u = (A @ dd[j]) * BOHR
        rd = np.linalg.norm((A @ (np.mod(xd[j] - s_vac + 0.5, 1) - 0.5)) * BOHR); disp.append((float(rp[i]), float(rd), float(np.linalg.norm(u))))
    print(f"[{S}] removed atom: s_red = {np.round(s_vac, 5).tolist()} -> unit-cell reduced {np.round(s_uc, 4).tolist()} = sublattice {sub}; "
          f"site vs nearest FFT point: grid coords {np.round(g, 3).tolist()} -> offset {dcart*1e3:.2f} mA; neighbours (r_p, r_d, |u|) A: {[tuple(round(x, 4) for x in t) for t in disp]}", flush=True)
    # core-masked radial profile (mask r < 0.5 A around every atom present in the DEFECTIVE cell), plane z = z_C
    Vp, _ = qe_io.get_pot(f"{scp}/Vks_{S}_p", subtract_mean=False, to_hartree=True); Vd, _ = qe_io.get_pot(f"{scd}/Vks_{S}_d", subtract_mean=False, to_hartree=True)
    dV = (Vd - Vp).transpose(2, 1, 0) * HA2EV; del Vp, Vd; iz = int(np.round(s_vac[2] * nr[2])) % nr[2]; plane = dV[:, :, iz]; del dV
    i = np.arange(nr[0]); j = np.arange(nr[1]); S1, S2 = np.meshgrid((i / nr[0] - s_vac[0] + 0.5) % 1 - 0.5, (j / nr[1] - s_vac[1] + 0.5) % 1 - 0.5, indexing="ij")
    X = (S1 * A[0, 0] + S2 * A[0, 1]) * BOHR; Y = (S1 * A[1, 0] + S2 * A[1, 1]) * BOHR; R = np.sqrt(X ** 2 + Y ** 2)
    mask = np.ones(plane.shape, bool)
    for s_at in xd:
        d1 = (i / nr[0] - s_at[0] + 0.5) % 1 - 0.5; d2 = (j / nr[1] - s_at[1] + 0.5) % 1 - 0.5; D1, D2 = np.meshgrid(d1, d2, indexing="ij")
        ra = np.sqrt(((D1 * A[0, 0] + D2 * A[0, 1]) * BOHR) ** 2 + ((D1 * A[1, 0] + D2 * A[1, 1]) * BOHR) ** 2); mask &= ra >= 0.5
    rmax = 0.5 * np.linalg.norm(A[:, 0]) * BOHR; edges = np.arange(0, rmax + 0.05, 0.05); m = (R < rmax) & mask
    cnt, _ = np.histogram(R[m], edges); sm, _ = np.histogram(R[m], edges, weights=plane[m]); rad = np.where(cnt > 0, sm / np.maximum(cnt, 1), np.nan); rc = 0.5 * (edges[1:] + edges[:-1])
    vals = [float(np.interp(r0, rc[np.isfinite(rad)], rad[np.isfinite(rad)])) for r0 in (1.0, 2.0, 3.0)]
    print(f"[{S}] core-masked azimuthal average: r = 1.0 / 2.0 / 3.0 A -> {vals[0]*1e3:+.1f} / {vals[1]*1e3:+.1f} / {vals[2]*1e3:+.1f} meV; masked fraction of the plane {1 - mask.mean():.3f}", flush=True)
    out2.update(**{f"{S}_rad_masked": rad, f"{S}_rc_masked": rc, f"{S}_sublattice": sub, f"{S}_s_vac": s_vac, f"{S}_grid_offset_A": dcart, f"{S}_nn": np.array(disp),
                   f"{S}_atoms_xy": np.array([[(A @ (np.mod(s - s_vac + 0.5, 1) - 0.5))[0] * BOHR, (A @ (np.mod(s - s_vac + 0.5, 1) - 0.5))[1] * BOHR] for s in xd])})
np.savez("results/M/ved_analysis.npz", **out2); print("saved results/M/ved_analysis.npz (+anomaly)")
