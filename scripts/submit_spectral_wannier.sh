#!/bin/bash
#SBATCH --job-name=specw
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:-7x7}
"$PROJ/.venv/bin/python" -u scripts/compute_spectral_wannier.py --size "$SIZE" \
   --manifest "$PROJ/wannier/$SIZE/wannier_manifest.json" \
   --grids 60,120,240 --etas 0.05,0.02,0.01 --rcut 0,1,2 --nk-int 300 \
   --out "results/M/specw_$SIZE.npz"
