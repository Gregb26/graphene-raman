#!/bin/bash
#SBATCH --job-name=Mdense
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# Dense M on the (p*N)^2 grid: ml (MPI, zero-padded M^L) -> nl (serial, fresh process) -> combine.
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load mpi4py/4.0.3 scipy-stack
SIZE=${1:?size e.g. 5x5}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
echo "[$(date)] $SIZE stage ml (64 ranks)"
srun --cpu-bind=cores "$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage ml --size "$SIZE" --out "results/M/M_L_dense_$SIZE.npy" || exit 2
echo "[$(date)] $SIZE stage nl (serial, fresh process)"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
"$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage nl --size "$SIZE" --out "results/M/M_NL_dense_$SIZE.npy" || exit 3
echo "[$(date)] $SIZE stage combine"
"$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage combine --size "$SIZE" \
   --ml "results/M/M_L_dense_$SIZE.npy" --nl "results/M/M_NL_dense_$SIZE.npy" \
   --out "results/M/M_dense_$SIZE.npy" --out-norm "results/M/M_dense_${SIZE}_norm.npy" || exit 4
echo "[$(date)] DONE $SIZE"
