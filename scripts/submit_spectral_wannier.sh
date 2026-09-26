#!/bin/bash
#SBATCH --job-name=specw
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=results/M2/logs/%x_%A.out
#SBATCH --error=results/M2/logs/%x_%A.err
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load scipy-stack
export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
RES=$("$PROJ/.venv/bin/python" -c 'from electron_defect_interaction.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')   # results/M2 (results/M gelé, R6)
mkdir -p "$RES/logs"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 FLEXIBLAS_NUM_THREADS=16
SIZE=${1:-7x7}
"$PROJ/.venv/bin/python" -u scripts/compute_spectral_wannier.py --size "$SIZE" \
   --manifest "$PROJ/wannier/$SIZE/wannier_manifest.json" \
   --grids 60,120,240 --etas 0.05,0.02,0.01 --rcut 0,1,2 --nk-int 300 \
   --out "$RES/specw_$SIZE.npz"
