#!/bin/bash
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=em2_bands
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00-00:30
#SBATCH --output=em2_bands_%j.out
#SBATCH --error=em2_bands_%j.err
# EM2 A.2-A.3: pw.x 'bands' on the 84 k of em2_kpoints_crystal.txt (potential of defect_uc_dense_27, copied
# density), then bands.x with lp = .true. (<c|p|v>, commutator [V_NL, r] included, PP/src/compute_ppsi.f90).
module restore qe
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR" || exit 1
set -e
srun pw.x -nk 4 < bands.in > bands.out
srun -n 8 bands.x < bands_pp.in > bands_pp.out
echo "DONE em2 bands"
