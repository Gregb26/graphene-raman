#!/usr/bin/env python
"""
epw_ring_check.py -- P11 point 2: decides the convention of EPW's phonselfen gamma___ (Im Pi = HWHM, or 2 Im Pi = FWHM) from a
30-second prtgkk run with k on rings around K (24k-24q/epw_g_ring: filkf = ring.kpt, 12 angles on two rings, filqf = GK.kpt = Gamma, K).
On each ring the interband |g_{45}|^2 (band 4 = pi below E_D, 5 = pi* above) is read directly, so the Dirac-cone Im Pi needs NO
assumption on the vertex angular structure:
  E2g (q = Gamma, modes 5+6, 2 valleys): Im Pi per mode = A_c hw <sum_nu |g_45^nu|^2>_ring / (4 v^2)
  A1' (q = K, mode 3, 1 valley):          Im Pi          = A_c hw_K <|g_45|^2>_ring / (4 v^2)
(sum_k wkf = 2 for spin; int d^2k delta(2vk - hw) = pi hw/(2v^2); A_c = sqrt3 a^2/2; v = Fermi velocity of the EPW bands.)
Prints the ring-averaged block sums, interband fractions, predicted Im Pi and the ratio gamma___(EPW)/Im Pi.
Usage: python scripts/epw_ring_check.py [--dir <epw_g_ring>] [--v 5.464] [--gamma-G 1.3346 --gamma-K 2.6448]
"""
import argparse, os, re, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--dir", default="/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/epw/24k-24q/epw_g_ring")
ap.add_argument("--v", type=float, default=5.464); ap.add_argument("--gamma-G", type=float, default=1.3346); ap.add_argument("--gamma-K", type=float, default=2.6448); ap.add_argument("--E_D", type=float, default=-4.2389); ap.add_argument("--tag", default="24k24q")
a = ap.parse_args(); A_c = np.sqrt(3) / 2 * 2.4659 ** 2; v = a.v; ED = a.E_D
data = {}; iq = ik = None
for l in open(os.path.join(a.dir, "epw.out")):
    m = re.match(r"\s*iq =\s*(\d+) coord", l)
    if m: iq = int(m.group(1)); continue
    m = re.match(r"\s*ik =\s*(\d+) coord", l)
    if m: ik = int(m.group(1)); data[(iq, ik)] = []; continue
    p = l.split()
    if iq and ik and len(p) == 10 and p[0].isdigit() and p[1].isdigit() and p[2].isdigit(): data[(iq, ik)].append([float(x) for x in p])
data = {k: np.array(vv) for k, vv in data.items() if len(vv)}
print(f"[ring] {len(data)} (iq, ik) blocks in {a.dir}/epw.out ; A_c = {A_c:.4f} A^2, v = {v} eV.A")
out = {}
def pick_modes(iq, ik, kind):
    """kind 'E2g': the two highest-frequency modes (checked degenerate); 'A1p': the mode with the largest pi/pi* block sum (character). Never by index."""
    r = data[(iq, ik)]; w = np.array([r[(r[:, 0] == 4) & (r[:, 1] == 4) & (r[:, 2] == nu)][0][5] for nu in range(1, 7)])
    if kind == "E2g":
        o = np.argsort(w); m = (int(o[-2]) + 1, int(o[-1]) + 1); assert abs(w[m[0] - 1] - w[m[1] - 1]) < 0.5, f"E2g not degenerate: {w}"; return m
    S = [sum(r[(r[:, 0] == ib) & (r[:, 1] == jb) & (r[:, 2] == nu)][0][7] ** 2 for ib in (4, 5) for jb in (4, 5)) for nu in range(1, 7)]
    return (int(np.argmax(S)) + 1,)
mG, mK = pick_modes(1, 1, "E2g"), pick_modes(2, 13, "A1p")
print(f"[modes] E2g(q=Gamma) = modes {mG} (two highest), A1'(q=K) = mode {mK[0]} (largest pi/pi* block sum; branch 3 at degauss 0.002 Ry, 6 at 0.02)")
for iq, name, modes, ring, nval, gepw in ((1, "E2g, q = Gamma", mG, range(1, 13), 2, a.gamma_G), (2, "A1', q = K", mK, range(13, 25), 1, a.gamma_K)):
    inter, intra, blk, e, om = [], [], [], [], []
    for ik in ring:
        r = data[(iq, ik)]; sel = lambda ib, jb, nu: r[(r[:, 0] == ib) & (r[:, 1] == jb) & (r[:, 2] == nu)][0]
        g2 = lambda ib, jb: sum(sel(ib, jb, nu)[7] ** 2 for nu in modes) * 1e-6
        om.append(np.mean([sel(4, 4, nu)[5] for nu in modes]) * 1e-3); inter.append(g2(4, 5)); intra.append(g2(4, 4) + g2(5, 5)); blk.append(g2(4, 4) + g2(5, 5) + g2(4, 5) + g2(5, 4))
        e.append((sel(4, 4, 1)[3] - ED, sel(5, 5, 1)[3] - ED))
    inter, blk, e, w = np.array(inter), np.array(blk), np.array(e), float(np.mean(om))
    impi = nval * A_c * w * inter.mean() / (8 * v ** 2)          # per mode: pi*2(spin)*nval*A_c/(2pi)^2 * pi hw/(2v^2) * <|g_45|^2 per mode>; sum over modes /len(modes)
    impi = impi * 2 / len(modes) if len(modes) == 2 else impi * 2  # -> both cases reduce to A_c hw <inter>/(4 v^2) for E2g (2 valleys, 2 modes) and A1' (1 valley, 1 mode)
    print(f"[{name}] hw = {w*1e3:.2f} meV ; ring energies E_k(4,5) - E_D = {e[:,0].mean()*1e3:+.1f} / {e[:,1].mean()*1e3:+.1f} meV")
    print(f"   block sum S = {blk.mean():.5f} eV^2 (spread {blk.min():.5f}..{blk.max():.5f}) ; interband <|g_45|^2> = {inter.mean():.5f} eV^2 = S/{blk.mean()/inter.mean():.2f} ; per-angle fractions {np.round(inter/blk, 3).tolist()}")
    print(f"   Dirac-cone Im Pi (HWHM) = {impi*1e3:.4f} meV ; EPW gamma___ = {gepw:.4f} meV ; gamma___/Im Pi = {gepw*1e-3/impi:.3f}")
    out[name] = (w, blk.mean(), inter.mean(), impi, gepw)
r = [out[k][4] * 1e-3 / out[k][3] for k in out]
print(f"[verdict] gamma___/Im Pi = {r[0]:.3f} (Gamma) and {r[1]:.3f} (K): EPW's gamma___ is 2 Im Pi = the FULL width (FWHM)" if all(abs(x - 2) < 0.05 for x in r) else f"[verdict] ratios {r}: no clean factor")
np.savez(f"results/epw/ring_check_{a.tag}.npz", **{k.replace(", ", "_").replace(" = ", "").replace("'", "p"): np.array(vv) for k, vv in out.items()}, v=v, A_c=A_c, units="hw eV; S, inter eV^2; ImPi eV (HWHM); gamma_epw meV")
print(f"saved results/epw/ring_check_{a.tag}.npz")
