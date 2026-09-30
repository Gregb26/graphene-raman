#!/usr/bin/env python
"""
epw_extract_gkk.py -- split the prtgkk blocks of an EPW output (one block per (iq, ik)) into per-(ibnd, jbnd, imode)
files  <prefix>_g{ib}{jb}_{nu}.dat  (one line per q, same columns as epw.out:
  ibnd jbnd imode enk[eV] enk+q[eV] omega(q)[meV] |g_sym|[meV] |g|[meV] Re(g)[meV] Im(g)[meV]),
the layout expected by epw_validate.py (--g-prefix). Same convention as the 16k-8q files epw_16k8q_g44_1.dat.
Usage: epw_extract_gkk.py <dir> [--prefix epw_g] [--out epw.out]
"""
import argparse, os, re, collections
ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("--prefix", default="epw_g"); ap.add_argument("--out", default="epw.out"); a = ap.parse_args()
rows = collections.defaultdict(list); nq = 0
for line in open(os.path.join(a.dir, a.out)):
    if re.match(r"\s*iq\s*=", line): nq += 1; continue
    p = line.split()
    if len(p) == 10 and p[0].isdigit() and p[1].isdigit() and p[2].isdigit():
        rows[(int(p[0]), int(p[1]), int(p[2]))].append(line.rstrip("\n"))
for (ib, jb, nu), ls in sorted(rows.items()):
    with open(os.path.join(a.dir, f"{a.prefix}_g{ib}{jb}_{nu}.dat"), "w") as f: f.write("\n".join(ls) + "\n")
print(f"{a.dir}: {nq} q blocks, {len(rows)} (ib,jb,nu) files, {min(len(v) for v in rows.values())}-{max(len(v) for v in rows.values())} lines each")
