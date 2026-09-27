#!/bin/bash
# R6 étape 3, fin automatique (sans session interactive) : ligne porte A.2 dans M_tests_summary.csv, table v1 -> v2 (3.7), copie versionnée,
# miroir md5 du .save à 128 bandes de R5 (manifeste du prompt). À soumettre avec --dependency=afterok:<job figures>.
#SBATCH --job-name=r6final
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm-r6-%x-%j.out
#SBATCH --error=slurm-r6-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:$PYTHONPATH"
PY="$PROJ/.venv/bin/python"
cd "$PROJ" || exit 1
echo "[$(date)] 3.6 : ligne porte A.2"; "$PY" -u "$WORK/r6_tests_gate_row.py" || echo "gate row FAILED"
echo "[$(date)] 3.7 : table v1 -> v2"; "$PY" -u "$WORK/r6_compare_v1_v2.py" || echo "compare FAILED"
echo "[$(date)] copie versionnée article/"
A=article/R6_production_corrigee/etape3; mkdir -p "$A/figures_v2" "$A/csv_v2"
cp -p "$WORK"/etape3/*.md "$WORK"/etape3/*.sh "$WORK"/etape3/*.py "$WORK"/etape3/*.txt "$A"/ 2>/dev/null
cp -p "$WORK"/r6_level1_gate.py "$WORK"/r6_compare_v1_v2.py "$WORK"/r6_tests_gate_row.py "$WORK"/submit_r6_final.sh "$A"/
cp -p results/M2/*.csv "$A/csv_v2/"; cp -p figures/*.pdf figures/*.png "$A/figures_v2/" 2>/dev/null
cp -p results/M2/logs/post_*.out results/M2/logs/specwd_218572*.out results/M2/logs/nkint_21857278.out results/M2/logs/resigma_218572*.out results/M2/logs/r6golden_21852238.out "$WORK/etape3/" 2>/dev/null
chgrp -R rrg-cotemich-ac "$A" "$WORK" 2>/dev/null
echo "[$(date)] manifeste : miroir md5 du .save 128 bandes de R5 (1,5 Go) -> qe_tmp_backup/R5_uc9x9_nb128/"
SRC=/home/gregb26/links/scratch/qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.save
DST=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/qe_tmp_backup/R5_uc9x9_nb128
mkdir -p "$DST" && rsync -a -r --no-o --no-g --open-noatime --chmod=Dg+s "$SRC" "$DST/" && ( cd "$SRC/.." && find defect_unit_cell_9x9.save -type f -exec md5sum {} + > "$DST/MD5SUMS_2026-09-26.txt" ) && ( cd "$DST" && md5sum -c MD5SUMS_2026-09-26.txt | grep -c -v ": OK$" | sed 's/^/écarts md5 : /' ) && chgrp -R rrg-cotemich-ac "$DST"
echo "[$(date)] md5 de results/M2 (relecture) : $(cd results/M2 && md5sum -c MD5SUMS_2026-09-25.txt 2>/dev/null | grep -c ': OK$') OK sur $(wc -l < results/M2/MD5SUMS_2026-09-25.txt)"
echo "[$(date)] done r6final"
