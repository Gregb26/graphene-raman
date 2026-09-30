#!/bin/bash
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=em2_Bours
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00-00:30
#SBATCH --output=em2_B_ours_%j.out
#SBATCH --error=em2_B_ours_%j.err
# EM2 B (decomposition): our sigma with the postw90 prefactor and grid, em2_B_ours.py (package functions, no src change).
source /etc/profile
module restore qe
module load mpi4py/4.0.3 scipy-stack
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
export PYTHONPATH=$PROJ/src:$PYTHONPATH
cd "$SLURM_SUBMIT_DIR" || exit 1
$PROJ/.venv/bin/python em2_B_ours.py
