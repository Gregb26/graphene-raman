import numpy as np, sys
from electron_defect_interaction.config import load_production, results_dir
RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)
S=sys.argv[1]
import os
if not os.path.exists(f"{RES}/M_dense_{S}_nb16.npy"): raise SystemExit(f"[non rejouable] {RES}/M_dense_{S}_nb16.npy absent (artefacts nbnd 16 supprimés au ménage de septembre 2026)")
from electron_defect_interaction.io import matrix_io; matrix_io.check_manifest(f"{RES}/M_dense_{S}.npy", require_normalization=matrix_io.M_NORM_V2)
A=np.load(f"{RES}/M_dense_{S}_nb16.npy", mmap_mode="r"); B=np.load(f"{RES}/M_dense_{S}.npy", mmap_mode="r")
nk=A.shape[1]; rng=np.random.default_rng(0); ks=rng.choice(nk, 40, replace=False)
w15=w16=0.0; wd=0.0
for i in ks:
    for j in ks:
        a=np.array(A[:,i,:,j]); b=np.array(B[:16,i,:16,j])
        s15=np.linalg.svd(a[:15,:15],compute_uv=False); t15=np.linalg.svd(b[:15,:15],compute_uv=False)
        s16=np.linalg.svd(a,compute_uv=False); t16=np.linalg.svd(b,compute_uv=False)
        w15=max(w15,np.abs(s15-t15).max()/s15[0]); w16=max(w16,np.abs(s16-t16).max()/s16[0]); wd=max(wd,np.abs(a-b).max()/np.abs(a).max())
print(f"[{S}] nb20[:16] vs nb16, 40x40 random k pairs: rel SV mismatch bands0..14={w15:.2e}  bands0..15={w16:.2e}  raw max|a-b|/max|a|={wd:.2e} (gauge-dependent)")
