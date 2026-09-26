"""eta sensitivity of median Gamma*N_cells (guardrail 3c). Usage: python _eta_scan.py 7x7"""
import sys
import numpy as np
from electron_defect_interaction.io import qe_io, matrix_io
from electron_defect_interaction.defects.many_body.single_defect import compute_T
from electron_defect_interaction.config import load_production, results_dir
raise SystemExit("[obsolète, R6 2026-09-26] ce script lit M_ed_*_norm.npy / gamma_*.npz (supprimés) et appartient à la chaîne d'avant le 2026-09-05 ; chaîne courante : compute_spectral_wannier.py, resonance_metrics.py, resonance_criteria.py")
RES = results_dir(load_production(verbose=False))          # R6 : results/M2 (results/M gelé)
HA2EV = 27.211386245988
N = sys.argv[1]
etas = [float(x) for x in sys.argv[2].split(",")] if len(sys.argv) > 2 else [0.05,0.10,0.15,0.20,0.30,0.40]
uc = f"data/graphene/unit_cell/qe/defect_{N}.save"
M = matrix_io.load_M_checked(f"{RES}/M_ed_{N}_norm.npy", units=matrix_io.EV)
eigs = qe_io.aligned_eigenvalues(uc, nk_expected=M.shape[1], shift_Fermi=True) * HA2EV
nb, nk = eigs.shape; NN = nb*nk; Ncells = nk; eflat = eigs.reshape(NN)
lvl = 6.0 / max(1, int((np.abs(eigs) <= 3.0).sum()))
print(f"{N}: nk={nk}, 2x level spacing = {2*lvl:.3f} eV (eta below this is undersampled)", flush=True)
for eta in etas:
    g = np.zeros(NN)
    for s in range(0, NN, 16):
        e = eflat[s:s+16]; T = compute_T(e, eigs, M, eta).reshape(len(e), NN, NN)
        for j in range(len(e)): g[s+j] = -2.0*T[j, s+j, s+j].imag
    print(f"{N}  eta={eta:.3f}  median G*Ncells = {np.median(np.abs(g))*Ncells*1e3:.2f} meV", flush=True)
