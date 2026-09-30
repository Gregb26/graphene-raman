#!/bin/bash
# Post-traitements de production (chaîne M2, R6) : une tâche par soumission, sorties dans results_dir(config) via les scripts eux-mêmes.
# Usage : sbatch [--time=…] scripts/slurm/submit_post.sh <locality|analyze|ksrec|resonance SIZE|criteria SIZE|c14|figures>
#SBATCH --job-name=post
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=08:00:00
#SBATCH --output=results/M2/logs/%x_%A.out
#SBATCH --error=results/M2/logs/%x_%A.err
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
cd "$PROJ" || exit 1
module restore qe; module load mpi4py/4.0.3 scipy-stack
export PYTHONPATH="$PROJ/src:$PYTHONPATH"          # pas d'installation éditable dans .venv (R6) ; préfixe : h5py de scipy-stack conservé
RES=$("$PROJ/.venv/bin/python" -c 'from graphene_raman.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')   # results/M2 (results/M gelé, R6)
mkdir -p "$RES/logs"
NT=${SLURM_CPUS_PER_TASK:-16}; export OMP_NUM_THREADS=$NT OPENBLAS_NUM_THREADS=$NT FLEXIBLAS_NUM_THREADS=$NT MKL_NUM_THREADS=$NT
PY="$PROJ/.venv/bin/python"; TASK=${1:?tâche}; shift
echo "[$(date)] post $TASK $* (HEAD $(git rev-parse --short HEAD), RES=$RES)"
run() { echo "[$(date)] >>> $*"; "$PY" -u "$@" || { echo "[$(date)] FAILED ($?): $*"; exit 2; }; }
case "$TASK" in
  locality)   run scripts/t/mwr_locality_coarse_vs_dense.py
              run scripts/t/m_rcut_convergence.py --size 9x9 --rcuts 0,1,2,3,4,5,6
              run scripts/t/m_rcut_convergence.py --size 12x12 --rcuts 0,1,2,3,4,5,6 ;;
  analyze)    run scripts/t/analyze_M.py
              run scripts/t/lnl_frobenius_all.py 5x5,6x6,7x7,8x8,9x9,12x12 ;;
  ksrec)      run scripts/t/ks_reconstruction_all.py --sizes 5x5,6x6,7x7,8x8,9x9,12x12
              # ved_analysis.npz ne dépend pas de M (potentiels V_ed seulement) : copie du produit v1 gelé, tracée
              [ -f "$RES/ved_analysis.npz" ] || { cp -p results/M/ved_analysis.npz "$RES/ved_analysis.npz" && chmod u+w "$RES/ved_analysis.npz" && echo "[$(date)] ved_analysis.npz copié de results/M (indépendant de M ; V_ed inchangés)"; }
              run scripts/t/sampling_table.py ;;
  resonance)  run scripts/t/resonance_metrics.py --size "${1:?size}"
              run scripts/t/resonance_criteria.py --size "$1" ;;
  criteria)   run scripts/t/resonance_criteria.py --size "${1:?size}" ;;
  c14)        run scripts/t/resonance_metrics.py --size 9x9 --shift-L-meV 25,-25 --out "$RES/resonance_9x9_shiftL.npz" ;;
  figures)    run scripts/fig/make_figures.py
              run scripts/t/level2_families.py
              run scripts/t/nkint_check_post.py --size 9x9 --nk 150,300,450,600
              run scripts/fig/make_figures_memoire.py
              # chapitre 5 : options de production de NOTES_EPW (chaîne mv0.02). Sans elles, les défauts (chaîne mv0.002, --prod-tag prod
              # absent) réécrivent les figures génériques avec l'ancienne chaîne et ne produisent pas fig_epw_vs_ed (R6, job 21872955)
              run scripts/t/epw_ed_vs_ep.py --selfen results/epw/selfen_240_dg0.02_mv0.02_T300.npz --tag 24k24q_mv0.02
              run scripts/fig/make_figures_epw.py --prod-tag 240_dg0.02_mv0.02 --phself-tag path_1200_dg0.02_mv0.02 --val-tag 24k24q_mv0.02 \
                  --sel-suffix _mv0.02 --control --kohn-val-tags 24k24q,24k24q_mv0.02 --dfpt-tag 24k24q ;;
  *) echo "tâche inconnue : $TASK"; exit 1 ;;
esac
echo "[$(date)] done $TASK $*"
