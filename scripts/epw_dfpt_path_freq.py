#!/usr/bin/env python
"""
epw_dfpt_path_freq.py -- phonon frequencies of a ph.x run with qplot/q_in_band_form (e.g. the DFPT 'prt' run dfpt_g_K, q on
Gamma-K-M-Gamma in CARTESIAN 2pi/a), read from the 'Diagonalizing the dynamical matrix' blocks of ph.out. No new run.
Output: results/epw/dfpt_path_freq_<tag>.npz : q_cart (nq,3) [2pi/a], s (nq,) path coordinate normalised to 1 on
Gamma-K-M-Gamma (same convention as validation_<tag>.npz: cumulative Cartesian distance), freq (nq, nmodes) cm^-1 (ph.x order, ascending).
Usage: python scripts/epw_dfpt_path_freq.py --phout <ph.out> --tag 24k24q
"""
import argparse, re, os, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("--phout", required=True); ap.add_argument("--tag", required=True); a = ap.parse_args()
q, F, cur = [], [], None
for l in open(a.phout):
    if "Diagonalizing the dynamical matrix" in l: cur = "diag"; continue
    if cur == "diag":
        m = re.match(r"\s*q = \(\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\)", l)
        if m: q.append([float(m.group(i)) for i in (1, 2, 3)]); F.append([]); cur = "freq"; continue
    if cur == "freq":
        m = re.match(r"\s*freq \(\s*\d+\)\s*=\s*([-\d.]+)\s*\[THz\]\s*=\s*([-\d.]+)\s*\[cm-1\]", l)
        if m: F[-1].append(float(m.group(2)))
        elif F[-1]: cur = None
q = np.array(q); F = np.array(F); s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]); s /= s[-1]
out = f"results/epw/dfpt_path_freq_{a.tag}.npz"; os.makedirs("results/epw", exist_ok=True)
np.savez(out, q_cart=q, s=s, freq=F, units="q cartesian 2pi/a; s normalised Gamma-K-M-Gamma; freq cm^-1", source=os.path.abspath(a.phout))
iK = int(np.argmin(np.abs(s - (2 / 3) / (2 / 3 + 1 / 3 + 1 / np.sqrt(3))))); iM = int(np.argmin(np.abs(s - 1 + (1 / np.sqrt(3)) / (2 / 3 + 1 / 3 + 1 / np.sqrt(3)))))
print(f"[dfpt {a.tag}] {len(q)} q, {F.shape[1]} modes ; Gamma: {F[0]} ; K (s={s[iK]:.4f}, q={q[iK]}): {F[iK]} ; M (s={s[iM]:.4f}, q={q[iM]}): {F[iM]} cm^-1 ; saved {out}")
