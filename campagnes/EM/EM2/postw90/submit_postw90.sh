#!/bin/bash
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=em2_pw90
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00-01:00
#SBATCH --output=em2_pw90_%j.out
#SBATCH --error=em2_pw90_%j.err
# EM2 B: postw90.x 3.1.0 (module qe) berry_task = kubo in each directory given as argument (k301, k301_nows, k1201).
# The module postw90.x is a SERIAL executable ('Running in serial (with serial executable)'): one task only
# (job 22035375 with 8 tasks ran 8 identical copies writing the same files; cancelled).
module restore qe
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR" || exit 1
set -e
for d in "$@"; do
  cd "$SLURM_SUBMIT_DIR/$d"
  start=$(date +%s)
  srun -n 1 postw90.x wannier
  echo "$d: postw90 $(( $(date +%s) - start )) s"
done
echo "DONE em2 postw90 $*"
