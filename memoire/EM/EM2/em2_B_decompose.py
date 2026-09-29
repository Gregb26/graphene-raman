"""
EM2, step B: postw90 (default, berry_kmesh 1201) - M4 split into its causes, telescoping so that the terms add up
exactly to the total:
  prefactor : ours 1/Delta, shifted 1200        - M4 (ours 1/hw, shifted 1200)
  grid      : ours 1/Delta, Gamma 1201          - ours 1/Delta, shifted 1200
  residual  : postw90 transl_inv, no ws         - ours 1/Delta, Gamma 1201   (r(R) from .mmn vs _tb.dat, H(R), ...)
  ws        : postw90 transl_inv                - postw90 transl_inv, no ws  (use_ws_distance)
  diagonal  : postw90 default                   - postw90 transl_inv         (diagonal of r: i sum w b M_nn vs -sum w b Im ln M_nn)
Table at 0.2 eV, the three lasers and the van Hove peak, and max |term| over the curve, for xx and yy; also |yy - xx|.
Writes em2_B_decomposition.txt.
"""

import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get("GRAPHENE_RAMAN", Path(os.environ["PROJECTS"]) / "graphene-raman"))
sys.path.insert(0, str(HERE))
from em2_B_compare import read_postw90, at   # noqa: E402

m4 = np.load(REPO / "memoire" / "EM" / "M4_sigma" / "em_sigma_full_N1200_eta0.04.npz")["sigma"]
o = np.load(HERE / "em2_B_ours.npz")
hw, p_def = read_postw90("k1201")
p_ti = read_postw90("k1201_ti")[1]
p_ti_nows = read_postw90("k1201_ti_nows")[1]

terms = [("prefactor 1/Delta", o["shift_delta"] - m4),
         ("grid Gamma 1201", o["gam1201"] - o["shift_delta"]),
         ("residual", p_ti_nows - o["gam1201"]),
         ("use_ws_distance", p_ti - p_ti_nows),
         ("diagonal of r", p_def - p_ti)]
total = p_def - m4
assert np.abs(sum(t for _, t in terms) - total).max() < 1e-12

pts = [0.2, 1.96, 2.33, 2.54, 4.05]
lines = ["EM2 B — postw90 (default, 1201) - M4 = sum of the terms below (sigma/sigma_0, g_s = 2 applied to postw90)", ""]
for comp, (i, j) in (("xx", (0, 0)), ("yy", (1, 1))):
    lines.append(f"{comp}: term               | " + " | ".join(f"{e:5.2f} eV " for e in pts) + " | max |term| (at eV)")
    for name, t in terms + [("TOTAL postw90 - M4", total)]:
        k = int(np.argmax(np.abs(t[:, i, j])))
        lines.append(f"    {name:18s} | " + " | ".join(f"{t[at(hw, e), i, j]:+.1e}" for e in pts)
                     + f" | {np.abs(t[:, i, j]).max():.2e} ({hw[k]:.2f})")
    lines.append("")
lines.append("|sigma_yy - sigma_xx| max over the curve (at the lasers):")
for name, s in (("M4", m4), ("ours 1/Delta Gamma 1201", o["gam1201"]), ("postw90 default", p_def),
                ("postw90 transl_inv", p_ti), ("postw90 transl_inv, no ws", p_ti_nows)):
    d = np.abs(s[:, 1, 1] - s[:, 0, 0])
    lines.append(f"    {name:26s} {d.max():.1e}  ({', '.join(f'{d[at(hw, e)]:.1e}' for e in pts[1:4])})")
text = "\n".join(lines) + "\n"
(HERE / "em2_B_decomposition.txt").write_text(text)
print(text)
