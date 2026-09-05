#!/bin/bash
#SBATCH --job-name=specwd
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# Level-1 map on the DENSE (zero-padded) M: R_cut x output grid x eta, one dense wannierization per size.
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:?size}; TAG=${2:+_$2}
declare -A DD=([5x5]=25 [7x7]=28 [8x8]=32 [9x9]=27); D=${DD[$SIZE]}
"$PROJ/.venv/bin/python" -u scripts/compute_spectral_wannier.py --size "$SIZE" --dense \
   --manifest "$PROJ/wannier/${D}x${D}/wannier_manifest.json" \
   --grids 60,120,240 --etas 0.05,0.02,0.01 --rcut 0,1,2,3 \
   --out "results/M/specwd_${SIZE}${TAG}.npz"
