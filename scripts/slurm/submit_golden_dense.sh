#!/bin/bash
#SBATCH --job-name=golden_dense
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=results/M2/logs/%x_%A.out
#SBATCH --error=results/M2/logs/%x_%A.err
# BLOCKING golden test on the dense 5x5 (25x25) M: local Wannier t-matrix vs dense compute_T in the
# same 5-WF subspace. Exit code = result of the test (0 PASS, 1 FAIL). The dense-vs-coarse coincident-k check that
# used to follow (check_M_dense_vs_coarse.py, obsolete) is done by analyze_M.py.
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load mpi4py/4.0.3 scipy-stack
export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:-5x5}
echo "[$(date)] golden dense $SIZE"
.venv/bin/python -u scripts/validation/test_local_tmatrix_real.py "$SIZE" --dense; rc=$?
echo "[$(date)] golden exit code $rc"
exit $rc
