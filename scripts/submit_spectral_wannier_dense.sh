#!/bin/bash
#SBATCH --job-name=specwd
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=results/M2/logs/%x_%A.out
#SBATCH --error=results/M2/logs/%x_%A.err
# Level-1 map on the DENSE (zero-padded) M: R_cut x output grid x eta, one dense wannierization per size.
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
RES=$("$PROJ/.venv/bin/python" -c 'from electron_defect_interaction.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')   # results/M2 (results/M gelé, R6)
mkdir -p "$RES/logs"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:?size}; TAG=${2:+_$2}
declare -A DD=([5x5]=25 [6x6]=24 [7x7]=28 [8x8]=32 [9x9]=27 [12x12]=24); D=${DD[$SIZE]}
GRIDS=${GRIDS:-60,120,240}; ETAS=${ETAS:-0.05,0.02,0.01}; RCUTS=${RCUTS:-0,1,2,3}
"$PROJ/.venv/bin/python" -u scripts/compute_spectral_wannier.py --size "$SIZE" --dense \
   --manifest "$PROJ/wannier/${D}x${D}/wannier_manifest.json" \
   --grids $GRIDS --etas $ETAS --rcut $RCUTS \
   --out "$RES/specwd_${SIZE}${TAG}.npz"
