#!/bin/bash
# R7 -- une soumission par sous-commande du pilote r7_driver.py. Usage : sbatch [-J nom] [--dependency=...] submit_r7.sh d1 [--sizes 6,9,12] [--out d1_reg]
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-r7-%x-%j.out
#SBATCH --error=slurm-r7-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}   # racine du depot (jamais codee en dur)
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack     # mpi4py : importe au chargement par defects/non_local.py (via r4_driver)
NT=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT FLEXIBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
export GRAPHENE_RAMAN="$PROJ"
cd "$WORK" || exit 1
echo "[$(date)] r7_driver $*"
"$PROJ/.venv/bin/python" -u "$WORK/r7_driver.py" "$@"; rc=$?
echo "[$(date)] done (exit $rc)"; exit $rc
