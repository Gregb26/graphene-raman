#!/bin/bash
# R4 -- une soumission par sous-commande du pilote r4_driver.py (J1 prep, J2 d6, J3 d1+d5+r, J4 d4, J5 d3).
# Usage : sbatch [--dependency=afterok:...] submit_r4.sh <cmd> [<cmd> ...]
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-r4-%x-%j.out
#SBATCH --error=slurm-r4-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}   # racine du dépôt (jamais codée en dur)
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack     # mpi4py : importé au chargement par defects/non_local.py (compute_M_NL)
NT=${R4_OMP:-$SLURM_CPUS_PER_TASK}                        # R4_OMP=1 pour les balayages de diagonalisations par lots (d3) : 1 fil BLAS + fils Python
export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT FLEXIBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT R4_EIG_WORKERS=$SLURM_CPUS_PER_TASK
export GRAPHENE_RAMAN="$PROJ"
cd "$WORK" || exit 1
for CMD in "$@"; do
  echo "[$(date)] r4_driver $CMD"
  "$PROJ/.venv/bin/python" -u "$WORK/r4_driver.py" "$CMD"; rc=$?
  if [ $rc -ne 0 ]; then echo "[$(date)] $CMD FAILED (exit $rc)"; exit 2; fi
done
echo "[$(date)] done: $*"
