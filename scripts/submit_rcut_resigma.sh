#!/bin/bash
#SBATCH --job-name=resigma
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# On-shell Re/Im Sigma per R_cut on the dense reference M (complement to m_rcut_convergence): rcut list as arg 2.
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:?size}; RC=${2:?rcut list}; TAG=$(echo "$RC" | tr -d ,)
"$PROJ/.venv/bin/python" -u scripts/rcut_resigma.py --size "$SIZE" --rcut "$RC" --out "results/M/resigma_${SIZE}_rc${TAG}.npz"
