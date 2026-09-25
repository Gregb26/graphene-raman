#!/bin/bash
#SBATCH --job-name=nkint
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# Convergence test of the INTERNAL k-grid N_k^int used for g0 (never varied before 2026-09-18: every run used 300^2).
# Same chain as production (rcut_resigma.py: frozen R_cut/grid/eta/window/ne_per_eta, hard gauge gate), only nk_int swept.
# Usage: sbatch scripts/submit_nkint_check.sh 9x9 "150 300 450 600"
# Compare afterwards (median on-shell Gamma*N_cells, Re Sigma, E_res per nk_int; states inside +-e_window only):
#   .venv/bin/python - <<'PY'
#   import numpy as np, glob
#   for f in sorted(glob.glob("results/M/resigma_9x9_rc3_nk*.npz"), key=lambda s: int(s.split("nk")[-1][:-4])):
#       z = np.load(f); S = z["Sigma_rc3"]; E = z["E_out"]; ED = float(z["E_D"]); m = np.isfinite(S.real)
#       G = -2 * S.imag; mm = m & (np.abs(E - ED) <= 1.5); e_res = float(E[mm][np.argmax(G[mm])] - ED)
#       print(f"nk_int={int(z['nk_int']):4d}  med Gamma = {np.median(G[m])*1e3:8.2f} meV  med ReSigma = {np.median(S.real[m])*1e3:7.2f} meV  E_res-E_D = {e_res:+.3f} eV")
#   PY
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:?size}; NKS=${2:-"150 300 450 600"}
for NK in $NKS; do
  echo "[$(date)] nk_int = $NK"
  "$PROJ/.venv/bin/python" -u scripts/rcut_resigma.py --size "$SIZE" --rcut 3 --nk-int "$NK" \
      --out "results/M/resigma_${SIZE}_rc3_nk${NK}.npz" || exit 1
done
echo "[$(date)] done"
