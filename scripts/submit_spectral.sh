#!/bin/bash
#SBATCH --job-name=spectral
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
"$PROJ/.venv/bin/python" -u scripts/compute_spectral.py --size "${1:-5x5}" --plot
