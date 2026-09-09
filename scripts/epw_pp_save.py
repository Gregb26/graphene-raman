#!/usr/bin/env python
"""Equivalent of EPW's pp.py (not shipped with the quantumespresso/7.5 module): gather the ph.x outputs into save/.
save/<prefix>.dvscf_q<i>  <- _ph0/<prefix>.dvscf1 (i=1) or _ph0/<prefix>.q_<i>/<prefix>.dvscf1 (i>1)
save/<prefix>.dyn_q<i>.xml <- <prefix>.dyn<i>.xml
save/<prefix>.phsave/     <- _ph0/<prefix>.phsave/
Verified against 16k-8q/phonons/save produced by the original pp.py (same names, same sizes).
Usage: python epw_pp_save.py <phonons_dir> <prefix> <nq_irr> [--check]"""
import os, shutil, sys, filecmp
d, prefix, nq = sys.argv[1], sys.argv[2], int(sys.argv[3]); check = "--check" in sys.argv
save = os.path.join(d, "save"); os.makedirs(save, exist_ok=True); n = 0
for i in range(1, nq + 1):
    src = os.path.join(d, "_ph0", f"{prefix}.dvscf1") if i == 1 else os.path.join(d, "_ph0", f"{prefix}.q_{i}", f"{prefix}.dvscf1")
    dyn = os.path.join(d, f"{prefix}.dyn{i}.xml")
    for s, t in ((src, os.path.join(save, f"{prefix}.dvscf_q{i}")), (dyn, os.path.join(save, f"{prefix}.dyn_q{i}.xml"))):
        if not os.path.exists(s): raise SystemExit(f"missing {s}")
        if check:
            ok = os.path.exists(t) and filecmp.cmp(s, t, shallow=False); print(f"{t}: {'OK' if ok else 'DIFF/MISSING'}"); continue
        shutil.copyfile(s, t); n += 1
ps, pt = os.path.join(d, "_ph0", f"{prefix}.phsave"), os.path.join(save, f"{prefix}.phsave")
if check: print(f"{pt}: {'OK' if os.path.isdir(pt) and not filecmp.dircmp(ps, pt).diff_files else 'DIFF/MISSING'}")
else:
    if os.path.isdir(pt): shutil.rmtree(pt)
    shutil.copytree(ps, pt); print(f"copied {n} files + {pt} ({len(os.listdir(pt))} entries) into {save}")
