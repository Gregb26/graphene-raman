#!/bin/bash
#SBATCH --job-name=epw_p1_post
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=results/epw/logs/%x_%A.out
#SBATCH --error=results/epw/logs/%x_%A.err
# P1 post-processing (24k-24q): split prtgkk outputs into per-(ib,jb,nu) files, then run the validations
# (bands, phonons, decay ratios, |g| DFPT vs EPW by degenerate subspaces). Run AFTER the three EPW interpolation jobs.
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
EPW=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/epw/${1:-24k-24q}; TAG=${2:-24k24q}
cd "$PROJ" || exit 1; mkdir -p results/epw/logs
module restore qe; module load scipy-stack
for d in epw_g_G epw_g_K; do "$PROJ/.venv/bin/python" scripts/epw_extract_gkk.py "$EPW/$d" --prefix epw_g || exit 1; done
"$PROJ/.venv/bin/python" -u scripts/epw_validate.py --root "$EPW" --tag "$TAG"
