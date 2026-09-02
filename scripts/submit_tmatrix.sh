#!/bin/bash
#SBATCH --job-name=tmatrix
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=results/M/logs/%x_%a_%A.out
#SBATCH --error=results/M/logs/%x_%a_%A.err
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
SIZE=${1:-5x5}
"$PROJ/.venv/bin/python" -u scripts/compute_tmatrix.py --size "$SIZE" --plot
