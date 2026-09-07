"""Band-resolved version of the coincident-k check: eigenvalue agreement coarse vs dense nscf at shared k,
and the singular-value mismatch of the M^L blocks restricted to the lowest nb_sub bands (excludes the
16-band truncation edge)."""
import sys, numpy as np
from electron_defect_interaction.io import qe_io
size = sys.argv[1]; N = int(size.split("x")[0]); D = {"5x5": 25, "6x6": 24, "7x7": 28, "8x8": 32, "9x9": 27, "12x12": 24}[size]
uc_c = f"data/graphene/unit_cell/qe/defect_{size}.save"
uc_d = f"/home/gregb26/links/scratch/qe_tmp/defect_uc_dense_{D}/defect_uc_dense_{D}.save"
kc, ec = qe_io.get_k_eigenvalues(uc_c, False); kd, ed = qe_io.get_k_eigenvalues(uc_d, False)
key = lambda k: tuple(np.round(np.mod(k + 1e-9, 1.0), 6))
idx_d = {key(k): i for i, k in enumerate(kd)}
pairs = [(ic, idx_d[key(k)]) for ic, k in enumerate(kc) if key(k) in idx_d]
ec = np.asarray(ec); ed = np.asarray(ed)
if ec.shape[0] != len(kc): ec = ec.T
if ed.shape[0] != len(kd): ed = ed.T
nbc = min(ec.shape[1], ed.shape[1]); ec = ec[:, :nbc]; ed = ed[:, :nbc]   # coarse nscf may have fewer bands (16) than the dense one (20)
de = np.array([ec[ic] - ed[id_] for ic, id_ in pairs]) * 27.211386
print(f"coincident k: {len(pairs)}; eps shapes {ec.shape} {ed.shape}; raw sample: {ec[0][:3]}")
print("max |eps_coarse - eps_dense| per band (raw units x 27.2114):", np.round(np.max(np.abs(de), axis=0), 5))
gap16 = None
print("min gap band16-band17 (dense, raw units):", None if gap16 is None else gap16.min())
Mc = np.load(sys.argv[2]); Md = np.load(sys.argv[3])
for nb_sub in (4, 8, 12, 15, 16):
    worst = 0.0
    for ic, id_ in pairs:
        for jc, jd in pairs:
            sa = np.linalg.svd(Mc[:nb_sub, ic, :nb_sub, jc], compute_uv=False); sb = np.linalg.svd(Md[:nb_sub, id_, :nb_sub, jd], compute_uv=False)
            worst = max(worst, np.max(np.abs(sa - sb)) / max(1e-30, sa[0]))
    print(f"bands 0..{nb_sub-1}: max rel singular-value mismatch = {worst:.3e}")
