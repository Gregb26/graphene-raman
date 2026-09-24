#!/bin/bash
#SBATCH --job-name=epw_p2_post
#SBATCH --account=rrg-cotemich-ac
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=results/epw/logs/%x_%A.out
#SBATCH --error=results/epw/logs/%x_%A.err
# P18 post-processing of a re-done EPW chain (default: 24k-24q_mv0.02, tag suffix _mv0.02): selfen (3 runs), phself (path + zoom,
# K convergence table), <D^2>, ring check (gamma___ / Im Pi with THIS chain's gamma), Gamma^ed / Gamma^ep ratio. Run AFTER the P2/P6
# EPW jobs and the P1 post (validation_<tag>.npz gives E_D). Usage: sbatch scripts/submit_epw_p2_post_mv.sh [chain dir] [tag] [E_D]
PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}   # racine du dépôt : variable d'environnement, sinon dépôt git du répertoire de soumission
CH=${1:-24k-24q_mv0.02}; TAG=${2:-24k24q_mv0.02}; SUF=${3:-_mv0.02}
EPW=/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/graphene/qe/epw/$CH
cd "$PROJ" || exit 1; module restore qe; module load scipy-stack; PY="$PROJ/.venv/bin/python"
ED=${4:-$($PY -c "import numpy as np; print(round(float(np.load('results/epw/validation_$TAG.npz')['bands_ED_epw']), 4))")}
echo "[p2 post] chain $CH tag $TAG E_D = $ED eV"
for r in 240_T300 240_T10 120_T300; do n=${r%_T*}; T=${r#*_T}; d=$EPW/epw2_selfen_${n}_T${T}_dg0.02
  [ -f $d/epw.out ] && $PY -u scripts/epw_selfen_post.py --dir $d --tag ${n}_dg0.02$SUF --E_D $ED || echo "[skip] $d"; done
$PY -u scripts/epw_phself_post.py --dir $EPW/epw6_phself_path_1200_dg0.02 $EPW/epw6_phself_zoom_1200_dg0.02 --tag path_1200_dg0.02$SUF
$PY -u scripts/epw_phself_post.py --table $EPW/epw6_phself_K_*_dg* --ref 1200_dg0.02
$PY -u scripts/epw_d2_extract.py --root $EPW --tag $TAG --phself-tag path_1200_dg0.02$SUF --E_D $ED
read GG GK V <<<$($PY -c "
import numpy as np; from electron_defect_interaction.electron_phonon import phself
P=np.load('results/epw/phself_path_1200_dg0.02$SUF.npz', allow_pickle=True); sp=phself.special_points(P['q']); i10=int(np.argmin(np.abs(P['T']-10)))
Rp=dict(T=P['T'], omega=P['omega'], gamma_epw=P['gamma_epw']); mE=list(phself.modes_E2g(Rp, sp['G'][0])); mA=phself.mode_A1p(Rp, sp['K'][0])
D=np.load('results/epw/d2_extract_$TAG.npz', allow_pickle=True); print(P['gamma_epw'][i10, sp['G'][0], mE].mean(), P['gamma_epw'][i10, sp['K'][0], mA], float(D['v_F_epw']))")
echo "[ring] gamma___(E2g, 10 K) = $GG meV, gamma___(A1', 10 K) = $GK meV, v_F EPW = $V eV.A"
$PY -u scripts/epw_ring_check.py --dir $EPW/epw_g_ring --gamma-G $GG --gamma-K $GK --v $V --E_D $ED --tag $TAG
$PY -u scripts/epw_ed_vs_ep.py --selfen results/epw/selfen_240_dg0.02${SUF}_T300.npz --tag $TAG
echo "[p2 post] done"
