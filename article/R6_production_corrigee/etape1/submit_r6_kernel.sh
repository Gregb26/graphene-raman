#!/bin/bash
# R6 J1 (étape 1.1) : M^L grossier avec le noyau corrigé (ee051ec) pour 9x9 (contrôle), 7x7 (contrôle, grille rééchantillonnée),
# 5x5, 6x6, 8x8, 11x11 (parties manquantes) ; test V_p ; noyau série et partagé ; reconstruction KS rejouée. Sorties dans ml/, kernel/, ksrec/.
#SBATCH --job-name=r6kernel
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --output=slurm-r6-%x-%j.out
#SBATCH --error=slurm-r6-%x-%j.err
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
WORK=${SLURM_SUBMIT_DIR:-$(dirname "$(readlink -f "$0")")}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:$PYTHONPATH"     # pas d'install éditable ; PYTHONPATH préfixé (h5py de scipy-stack conservé)
PY="$PROJ/.venv/bin/python"
cd "$PROJ" || exit 1
mkdir -p "$WORK/ml" "$WORK/kernel"
echo "[$(date)] HEAD $(git rev-parse --short HEAD) ; noyau : $(grep -n 'inv_sqrtO =' src/electron_defect_interaction/defects/local_R.py | tr '\n' ' ')"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
for N in 9x9 7x7 5x5 6x6 8x8; do      # 11x11 exclu : .save de maille écrasé par un run bands (XML 181 k, wfc 121 k)
  [ -f "$WORK/ml/M_L_${N}_v2.npy" ] && { echo "[$(date)] M^L $N déjà calculé, sauté"; continue; }
  UC=data/graphene/unit_cell/qe/defect_${N}.save; SCP=data/graphene/supercell/qe/defect_${N}_p.save; SCD=data/graphene/supercell/qe/defect_${N}_d.save
  echo "[$(date)] M^L $N grossier, noyau corrigé, 32 rangs"
  srun --cpu-bind=cores "$PY" -u scripts/compute_M.py --stage ml --uc "$UC" --sc-p "$SCP" --pot-p "$SCP/Vks_${N}_p" --pot-d "$SCD/Vks_${N}_d" \
       --bands all --block-size 50000 --out "$WORK/ml/M_L_${N}_v2.npy" || { echo "[$(date)] M^L $N FAILED"; exit 2; }
done
echo "[$(date)] test V_p (super-cellule 9x9, 4 bandes, MPI)"
srun --cpu-bind=cores "$PY" -u "$WORK/r6_kernel_check.py" vp || echo "[$(date)] vp FAILED rc=$?"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16
"$PY" -u "$WORK/r6_kernel_check.py" compare --sizes 9x9,7x7,5x5,6x6,8x8 || echo "compare FAILED rc=$?"
echo "[$(date)] porte A.2 sur le M^L 9x9 recalculé (NL = _test_mnl/M_NL_serial.npy de juin, sans sidecar)"
"$PY" -u scripts/gate_M_normalization.py --size 9x9 --level coarse --ml "$WORK/ml/M_L_9x9_v2.npy" --nl results/M/_test_mnl/M_NL_serial.npy --allow-v1 \
      --out "$WORK/kernel/gate_9x9_coarse_v2kernel.json"; echo "[$(date)] gate rc=$?"
echo "[$(date)] noyau série compute_ML_R (5x5, 2 bandes)"
"$PY" -u "$WORK/r6_kernel_check.py" serial || echo "serial FAILED rc=$?"
echo "[$(date)] reconstruction KS rejouée (9x9 grossier + dense), sortie hors results/"
mkdir -p "$WORK/ksrec/results/M"; ln -sfn "$PROJ/data" "$WORK/ksrec/data"; ln -sfn "$PROJ/scripts" "$WORK/ksrec/scripts"
( cd "$WORK/ksrec" && "$PY" -u "$PROJ/scripts/ks_reconstruction_all.py" --sizes 9x9 ) || echo "ksrec FAILED rc=$?"
"$PY" -u "$WORK/r6_kernel_check.py" ksrec || echo "ksrec compare FAILED rc=$?"
echo "[$(date)] noyau partagé compute_ML_R_mpi_shared (--coarse 5x5, 2 bandes ; attendu : NameError 'Ndiag' hors rang 0 du nœud)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FLEXIBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
timeout 20m srun --cpu-bind=cores "$PY" -u scripts/compute_M_dense_stages.py --stage ml --size 5x5 --coarse --bands 0,1 --block-size 2000 \
      --out "$WORK/ml/M_L_5x5_shared_coarsecheck_v2.npy"; echo "[$(date)] shared rc=$?"
export OMP_NUM_THREADS=16
"$PY" -u "$WORK/r6_kernel_check.py" shared
echo "[$(date)] done J1"
