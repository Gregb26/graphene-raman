#!/bin/bash
#SBATCH --job-name=M_ed
#SBATCH --account=rrg-cotemich-ac
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=06:00:00
#SBATCH --output=results/M/logs/M_ed_%a_%A.out
#SBATCH --error=results/M/logs/M_ed_%a_%A.err

# Compute the full electron-defect scattering matrix M = M^L + M^NL for one supercell per
# array task, V_ed = V_d - V_p, all bands. Works for ANY size via the two-stage split:
#   stage 1 (MPI, srun)   -> M^L              (compute_M.py --stage ml)
#   stage 2 (serial, fresh process) -> M^NL + combine -> M  (compute_M.py --stage nl)
# Stage 2 MUST be a fresh process: compute_ML_R_mpi leaves a latent heap corruption for
# nk >= 81 that would abort compute_M_NL if run in the same process.
#
# Edit SIZES to whatever you need and set --array=0-$((len-1)) to match.

PROJ=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects
cd "$PROJ" || exit 1

SIZES=(5x5 6x6 7x7 8x8 9x9 10x10 11x11 12x12)
N=${SIZES[$SLURM_ARRAY_TASK_ID]}

module restore qe
module load mpi4py/4.0.3 scipy-stack

UC=data/graphene/unit_cell/qe/defect_${N}.save
SCP=data/graphene/supercell/qe/defect_${N}_p.save
SCD=data/graphene/supercell/qe/defect_${N}_d.save
PY="$PROJ/.venv/bin/python"

ML="results/M/M_L_${N}.npy"
NL="results/M/M_NL_${N}.npy"
MED="results/M/M_ed_${N}.npy"

echo "[$(date)] size=$N  stage 1/3: M^L (MPI, $SLURM_NTASKS ranks)"
# Pin all threading layers to 1: parallelism is over MPI ranks; oversubscribing BLAS
# corrupts the heap ("double free") in M^NL for nk >= 81.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
srun --cpu-bind=cores "$PY" -u scripts/compute_M.py --stage ml \
    --uc "$UC" --sc-p "$SCP" \
    --pot-p "$SCP/Vks_${N}_p" --pot-d "$SCD/Vks_${N}_d" \
    --bands all --block-size 50000 \
    --out "$ML"

echo "[$(date)] size=$N  stage 2/3: M^NL (serial, fresh process)"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 FLEXIBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMEXPR_NUM_THREADS=8
"$PY" -u scripts/compute_M.py --stage nl \
    --uc "$UC" --sc-p "$SCP" --sc-d "$SCD" --upf "$UC/C.upf" \
    --bands all \
    --out "$NL"

echo "[$(date)] size=$N  stage 3/3: combine M = M^L + M^NL"
"$PY" -u scripts/compute_M.py --stage combine --ml "$ML" --nl "$NL" --out "$MED"

echo "[$(date)] done size=$N"
