#!/bin/bash
# R5 J4 -- M grossier à 128 bandes : M^L (srun MPI, noyau R partagé) -> M^NL (processus frais) -> somme + vérifications.
# Ressources calquées sur scripts/submit_M_dense.sh (32 rangs x 6 fils, nœud exclusif).
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=r5m128
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=slurm-r5-%x-%j.out
#SBATCH --error=slurm-r5-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"
cd "$WORK" || exit 1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK FLEXIBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "[$(date)] b2ml ($SLURM_NTASKS rangs x $SLURM_CPUS_PER_TASK fils)"
srun --cpu-bind=cores "$PROJ/.venv/bin/python" -u "$WORK/r5_driver.py" b2ml || exit 2
export OMP_NUM_THREADS=32 OPENBLAS_NUM_THREADS=32 FLEXIBLAS_NUM_THREADS=32 MKL_NUM_THREADS=32
echo "[$(date)] b2nl (processus frais)"
"$PROJ/.venv/bin/python" -u "$WORK/r5_driver.py" b2nl || exit 3
echo "[$(date)] b2sum"
"$PROJ/.venv/bin/python" -u "$WORK/r5_driver.py" b2sum || exit 4
echo "[$(date)] done"
