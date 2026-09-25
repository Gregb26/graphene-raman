#!/bin/bash
#SBATCH --job-name=golden_dense
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# BLOCKING golden test on the dense 5x5 (25x25) M: local Wannier t-matrix vs dense compute_T in the
# same 5-WF subspace. Then the dense-vs-coarse coincident-k check on M^L (9x9) / M (others).
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load mpi4py/4.0.3 scipy-stack
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:-5x5}
echo "[$(date)] golden dense $SIZE"
.venv/bin/python -u scripts/test_local_tmatrix_real.py "$SIZE" --dense; rc=$?
echo "[$(date)] golden exit code $rc"
echo "[$(date)] coincident-k check (full M dense vs M_ed coarse)"
.venv/bin/python -u scripts/check_M_dense_vs_coarse.py "$SIZE" "results/M/M_ed_$SIZE.npy" "results/M/M_dense_$SIZE.npy"
exit $rc
