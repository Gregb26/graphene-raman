#!/bin/bash
# R5 -- une soumission par sous-commande du pilote r5_driver.py (J1 a, J2 c, J5 b1 b3). Usage : sbatch [--dependency=...] submit_r5.sh <cmd> [<cmd> ...]
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-r5-%x-%j.out
#SBATCH --error=slurm-r5-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}   # racine du dépôt (jamais codée en dur)
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack     # mpi4py : importé au chargement par defects/non_local.py
NT=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT FLEXIBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
export GRAPHENE_RAMAN="$PROJ"
cd "$WORK" || exit 1
for CMD in "$@"; do
  echo "[$(date)] r5_driver $CMD"
  "$PROJ/.venv/bin/python" -u "$WORK/r5_driver.py" "$CMD"; rc=$?
  if [ $rc -ne 0 ]; then echo "[$(date)] $CMD FAILED (exit $rc)"; exit 2; fi
done
echo "[$(date)] done: $*"
