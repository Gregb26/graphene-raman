#!/bin/bash
#SBATCH --job-name=lnlfrob
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=results/M2/logs/%x_%A.out
#SBATCH --error=results/M2/logs/%x_%A.err
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
RES=$("$PROJ/.venv/bin/python" -c 'from electron_defect_interaction.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')   # results/M2 (results/M gelé, R6)
mkdir -p "$RES/logs"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 FLEXIBLAS_NUM_THREADS=4
"$PROJ/.venv/bin/python" -u scripts/t/lnl_frobenius_all.py "${1:-5x5,6x6,7x7,8x8,9x9,12x12}"
