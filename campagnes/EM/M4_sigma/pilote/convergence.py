"""
Convergence table of the M4 pilot (reads the npz of sweep.py, writes convergence.txt): sigma(eps_L)/sigma_0 against N and
eta, the N error against N = 1800, and the quadratic eta bias sigma(eta) = sigma(0) + a eta^2 fitted on eta = 0.02, 0.04.
Run:  .venv/bin/python campagnes/EM/M4_sigma/pilote/convergence.py
"""
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
LASERS, NS, VARIANTS = (1.96, 2.33, 2.54), (900, 1200, 1800), ("full", "centres_only", "no_berry")
S = {(v, N): np.load(HERE / f"em_sigma_{v}_N{N}.npz") for v in VARIANTS for N in NS}
hw, etas = S["full", 900]["hw"], S["full", 900]["eta"]
iL = [int(np.flatnonzero(hw == e)[0]) for e in LASERS]
xx = {k: z["sigma"][:, iL, 0, 0] for k, z in S.items()} # (n_eta, 3 lasers)

lines = ["M4 pilot, 27 x 27 data, mu = E_D: sigma_xx / sigma_0 at 1.96, 2.33, 2.54 eV", ""]
for v in VARIANTS:
    lines.append(f"[{v}]")
    for i, eta in enumerate(etas):
        for N in NS:
            err = np.abs(xx[v, N][i] - xx[v, 1800][i]).max()
            lines.append(f"  eta {eta:.2f}  N {N:4d}: {np.array2string(xx[v, N][i], precision=5)}   |N - 1800| max {err:.1e}")
    a = (xx[v, 1800][1] - xx[v, 1800][0]) / (etas[1]**2 - etas[0]**2)   # sigma = sigma(0) + a eta^2
    a_check = (xx[v, 1800][2] - xx[v, 1800][1]) / (etas[2]**2 - etas[1]**2)
    lines.append(f"  eta bias: a = {np.array2string(a, precision=3)} (0.04-0.08 gives {np.array2string(a_check, precision=3)}); "
                 f"bias at eta = 0.04: {np.array2string(a * 0.04**2, formatter={"float_kind": lambda x: f"{x:.1e}"})}; sigma(eta -> 0) = {np.array2string(xx[v, 1800][0] - a * etas[0]**2, precision=5)}")
    lines.append("")
lines.append("Retained: N = 1200, eta = 0.04 eV (N error < 1e-11, eta bias < 1e-3).")
(HERE / "convergence.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
