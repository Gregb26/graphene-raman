#!/bin/bash
#SBATCH --job-name=Mdense
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=12:00:00
#SBATCH --output=results/M/logs/%x_%A.out
#SBATCH --error=results/M/logs/%x_%A.err
# Dense M on the (p*N)^2 grid: ml (MPI real-space, node-shared u_nk) -> nl (serial, fresh process) -> combine.
# Usage: sbatch scripts/submit_M_dense.sh <size> [extra ml args, e.g. --coarse]
PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1
module restore qe; module load mpi4py/4.0.3 scipy-stack
SIZE=${1:?size e.g. 5x5}; shift
EXTRA="$@"
TAG=$SIZE; [[ "$EXTRA" == *coarse* ]] && TAG="${SIZE}_coarsecheck"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "[$(date)] $SIZE stage ml ($SLURM_NTASKS ranks x $SLURM_CPUS_PER_TASK threads) $EXTRA"
srun --cpu-bind=cores "$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage ml --size "$SIZE" --block-size 2000 $EXTRA --out "results/M/M_L_dense_$TAG.npy" || exit 2
[[ "$EXTRA" == *coarse* ]] && { echo "[$(date)] coarse check done"; exit 0; }
echo "[$(date)] $SIZE stage nl (serial, fresh process)"
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 FLEXIBLAS_NUM_THREADS=32
"$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage nl --size "$SIZE" --out "results/M/M_NL_dense_$SIZE.npy" || exit 3
echo "[$(date)] $SIZE stage combine"
"$PROJ/.venv/bin/python" -u scripts/compute_M_dense_stages.py --stage combine --size "$SIZE" \
   --ml "results/M/M_L_dense_$SIZE.npy" --nl "results/M/M_NL_dense_$SIZE.npy" \
   --out "results/M/M_dense_$SIZE.npy" --out-norm "results/M/M_dense_${SIZE}_norm.npy" || exit 4
echo "[$(date)] DONE $SIZE"
