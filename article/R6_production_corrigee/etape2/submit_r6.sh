#!/bin/bash
# R6 -- tâches série (16 cœurs) : assemble (J2, étape 1.2), gate (J3, étape 1.3), states (J4, identité des 30 états 9x9).
# Usage : sbatch [--dependency=afterok:ID] submit_r6.sh <assemble|gate|states|d4|b3|d3>
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-r6-%x-%j.out
#SBATCH --error=slurm-r6-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
R5=$(dirname "$WORK")/R5_base_vs_M
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:$PYTHONPATH"
NT=$SLURM_CPUS_PER_TASK; export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT FLEXIBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
PY="$PROJ/.venv/bin/python"
cd "$PROJ" || exit 1
TASK=${1:?tâche : assemble | gate | states}
echo "[$(date)] R6 $TASK (HEAD $(git rev-parse --short HEAD))"
case "$TASK" in
  assemble)
    M=results/M; M2=results/M2; ML=$WORK/ml; SUM=$WORK/assemble_summary.jsonl; mkdir -p "$M2/logs"
    A="$PY -u scripts/assemble_M2.py --summary-file $SUM"
    fail=0
    # grossiers : parties sur disque (7) ; M^NL = M_ed - M_L de juin (9, 10, 12) ; M^L recalculé (5, 6, 8 ; 11 avec M_NL = M_ed - M_L juin)
    $A --ml $M/M_L_7x7.npy   --nl $M/M_NL_7x7.npy           --n-cells 49  --out-l $M2/M_L_7x7.npy   --out-nl $M2/M_NL_7x7.npy   --out-m $M2/M_ed_7x7.npy   --label "coarse 7x7 : L (sept., 217) et NL sur disque" || fail=1
    for S in 9x9:81 10x10:100 12x12:144; do N=${S%%:*}; C=${S##*:}
      $A --ml $M/M_L_$N.npy --nl-from-diff $M/M_ed_$N.npy --n-cells $C --out-l $M2/M_L_$N.npy --out-nl $M2/M_NL_$N.npy --out-m $M2/M_ed_$N.npy --label "coarse $N : M_NL = M_ed(juin) - M_L(juin)" || fail=1
    done
    for S in 5x5:25 6x6:36 8x8:64; do N=${S%%:*}; C=${S##*:}
      $A --ml-v2 $ML/M_L_${N}_v2.npy --nl-from-diff $M/M_ed_$N.npy --n-cells $C --out-l $M2/M_L_$N.npy --out-nl $M2/M_NL_$N.npy --out-m $M2/M_ed_$N.npy --label "coarse $N : M_L recalculé (ee051ec), M_NL = M_ed(juin) - M_L_v2/N_cells" || fail=1
    done
    echo "[assemble] 11x11 exclu : M^L de juin sur grille non commensurable, .save de maille écrasé (XML 181 k) -> pas de M2 11x11"
    # denses : N_cells = N^2 ; M_NL en lien symbolique vers results/M (inchangé)
    for S in 5x5:25 6x6:36 7x7:49 8x8:64 9x9:81 12x12:144; do N=${S%%:*}; C=${S##*:}
      $A --ml $M/M_L_dense_$N.npy --nl $M/M_NL_dense_$N.npy --nl-link --n-cells $C --out-l $M2/M_L_dense_$N.npy --out-nl $M2/M_NL_dense_$N.npy --out-m $M2/M_dense_$N.npy --label "dense $N" || fail=1
    done
    $A --ml $M/M_L_dense_9x9_coarsecheck.npy --n-cells 81 --out-l $M2/M_L_dense_9x9_coarsecheck.npy --label "coarsecheck 9x9 (L seul)" || fail=1
    # R5, 128 bandes x 81 k (étape 2, c-all)
    $A --ml $R5/b/M_L_9x9_nb128.npy --nl $R5/b/M_NL_9x9_nb128.npy --n-cells 81 --out-l $M2/M_L_9x9_nb128.npy --out-nl $M2/M_NL_9x9_nb128.npy --out-m $M2/M_ed_9x9_nb128.npy --label "coarse 9x9 nb128 (R5 B.2)" || fail=1
    ls -la $M2 | tee "$WORK/assemble_ls.txt"; du -sh $M2
    echo "[$(date)] assemble done fail=$fail"; exit $fail ;;
  gate)
    G="$PY -u scripts/gate_M_normalization.py"; mkdir -p "$WORK/gate"; fail=0
    for N in 5x5 6x6 7x7 8x8 9x9 10x10 12x12; do
      $G --size $N --level coarse --out "$WORK/gate/gate_${N}_coarse.json" || { echo "[gate] $N coarse rc=$?"; fail=1; }
    done
    for N in 5x5 6x6 7x7 8x8 9x9 12x12; do
      $G --size $N --level dense --out "$WORK/gate/gate_${N}_dense.json" || { echo "[gate] $N dense rc=$?"; fail=1; }
    done
    $G --size 9x9 --level coarse --uc /home/gregb26/links/scratch/qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.save \
       --ml results/M2/M_L_9x9_nb128.npy --nl results/M2/M_NL_9x9_nb128.npy --out "$WORK/gate/gate_9x9_coarse-nb128.json" || { echo "[gate] nb128 rc=$?"; fail=1; }
    # documentation : les fichiers v1 gelés (attendu : rapport L = N_cells, REFUSÉ)
    $G --size 9x9 --level dense --ml results/M/M_L_dense_9x9.npy --nl results/M/M_NL_dense_9x9.npy --allow-v1 --out "$WORK/gate/v1_gate_9x9_dense.json"; echo "[gate] v1 9x9 dense rc=$? (attendu 3)"
    $G --summary "$WORK/gate" | tee "$WORK/gate/gate_table.md"
    echo "[$(date)] gate done fail=$fail"; exit $fail ;;
  states)
    "$PY" -u "$WORK/r6_states_identity.py"; rc=$?; echo "[$(date)] states rc=$rc"; exit $rc ;;
  d4)
    "$PY" -u "$WORK/r6_d4.py" prep d4 d3; rc=$?; echo "[$(date)] d4 rc=$rc"; exit $rc ;;
  b3)
    "$PY" -u "$WORK/r6_d4.py" b3; rc=$?; echo "[$(date)] b3 rc=$rc"; exit $rc ;;
  d3)
    # diagonalisations par lots (pole_criterion._batched) : 1 fil BLAS + 16 fils Python (R4 : 4,7 ms contre 61 ms par matrice avec 16 fils BLAS)
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 R4_EIG_WORKERS=$SLURM_CPUS_PER_TASK
    "$PY" -u "$WORK/r6_d4.py" d3; rc=$?; echo "[$(date)] d3 rc=$rc"; exit $rc ;;
  *) echo "tâche inconnue $TASK"; exit 1 ;;
esac
