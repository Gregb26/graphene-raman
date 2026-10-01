#!/bin/bash
# R5 (B.1) -- nscf pw.x 128 bandes, maille unitaire, grille 9x9 (81 k). NE PAS LANCER SANS LE GO « nscf » (un seul calcul QE autorisé).
# Le .save de production (defect_unit_cell_9x9) n'est jamais écrit : charge-density.hdf5 + data-file-schema.xml sont copiés
# dans un outdir NEUF (R5_uc9x9_nb128) avant le nscf. Ressources calquées sur unit_cell/27x27/submit_dense_nb20.sh.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --job-name=r5nscf
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-r5-%x-%j.out
#SBATCH --error=slurm-r5-%x-%j.err
set -e
module restore qe
export OMP_NUM_THREADS=1
SCRATCH=/home/gregb26/links/scratch/qe_tmp
SRC=$SCRATCH/defect_unit_cell_9x9/defect_unit_cell_9x9.save          # production, lecture seule
OUT=$SCRATCH/R5_uc9x9_nb128/defect_unit_cell_9x9.save                 # neuf (même prefix, outdir différent)
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
cd "$WORK"
[ -e "$OUT/wfc1.hdf5" ] && { echo "outdir déjà rempli ($OUT) : refus d'écraser"; exit 3; }
mkdir -p "$OUT"
cp -p "$SRC/charge-density.hdf5" "$SRC/data-file-schema.xml" "$OUT/"
md5sum "$SRC/charge-density.hdf5" "$OUT/charge-density.hdf5" "$SRC/data-file-schema.xml" "$OUT/data-file-schema.xml"
echo "[$(date)] nscf 128 bandes -> $OUT"
srun pw.x < nscf_nb128.in > nscf_nb128.out
grep -E "PWSCF.*WALL|number of k points|not converged" nscf_nb128.out | head -20
ls "$OUT" | wc -l; du -sh "$OUT"
echo "[$(date)] done"
