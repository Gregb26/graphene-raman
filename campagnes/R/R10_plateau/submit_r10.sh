#!/bin/bash
# R10 -- lanceur unique de la campagne (entrée de campagne, CLAUDE.md « Campagnes de calcul ») ; la tâche est choisie par argument.
#
# Depuis le nœud de connexion :   bash submit_r10.sh <tâche> [arguments]      (DEP=afterok:ID[:ID…] pour une dépendance ; DRY=1 : diff seul)
#   1. affiche le diff de r10_driver.py et de ce fichier depuis la dernière soumission (première soumission : md5 et nombre de lignes) ;
#   2. sbatch avec les ressources de la tâche ;
#   3. garde une copie des deux fichiers et du diff dans submitted/<jobid>/ ; une ligne par job dans JOBID.
# Dans le job (SLURM_JOB_ID défini) : exécute la tâche.
#
# Tâches (GO 1) :
#   a           J1 : a0 (porte A.0 ; échec : STOP sur A) puis a1, a2, a3. 1 tâche, 8 cœurs, 64 Go, 1 h.
#   b [args]    J2 : B.0 (porte), B.1, B.2 (args passés à « r10_driver.py b », ex. --parts b1,b2). 16 cœurs, 128 Go, 3 h.
# Tâches (GO 2, partie C ; produits dans results_dir de la config = results/M2_plateau, journaux de C1 et C1f dans <results_dir>/logs) :
#   c0          portes C.0 (r10_driver.py c0 ; échec : code 3). 16 cœurs, 128 Go, 2 h.
#   c1          macro : soumet la production C1 avec les lanceurs de scripts/ tels quels (mêmes réglages que R6, D13 ; --output/--error vers
#               <results_dir>/logs) : specwd ×6 (grilles 60/120/240, η 0,05/0,02/0,01, R_cut 0…4, tag prod), nkint 9x9, resigma 9x9 (0,1,2,3 et 4),
#               post locality, post analyze (GOLDEN_RESULT de R6), post resonance 9x9 / 6x6 / 12x12 ; plus p_ved et p_c14 (ce lanceur).
#   p_ved       analyze_Ved.py (8 tailles, P-c2), copie md5 de ks_reconstruction.npz depuis results/M2, sampling_table.py. 16 cœurs, 64 Go, 1 h.
#   p_c14       resonance_metrics.py 9x9 --shift-L-meV 9.05,-9.05 (C14, D10) -> resonance_9x9_shiftL.npz. 16 cœurs, 96 Go, 2 h.
#   c1f         r10_driver.py c1post (m_rcut_resigma.csv, lignes de tab:tests_M) puis make_figures.py, level2_families.py, nkint_check_post.py,
#               make_figures_memoire.py, epw_ed_vs_ep.py, make_figures_epw.py (options NOTES_EPW), figures dans fig/ de ce répertoire. 16 cœurs, 96 Go, 1 h.
#   c2          r10_driver.py c2 (tab:rcut_M à trois colonnes, Kaasbjerg). 16 cœurs, 128 Go, 2 h.
#   c3          r10_driver.py c3, c4, c5 (contrôles de l'audit, planches, table de correspondance, README/MD5 de results_dir). 8 cœurs, 32 Go, 1 h.
#   c6          r10_driver.py c6 (anomalie ⟨ΔV⟩_3D : SCF et profils ⟨ΔV⟩(z)). 8 cœurs, 32 Go, 1 h.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --output=slurm-r10-%x-%j.out
#SBATCH --error=slurm-r10-%x-%j.err

if [ -z "${SLURM_JOB_ID:-}" ]; then
  # ------------------------------------------------------------------ côté nœud de connexion : diff puis soumission
  WORK=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
  TASK=${1:?"tâche : a | b [args] | c0 | c1 | p_ved | p_c14 | c1f | c2 | c3 | c6"}; shift
  PROJ=${GRAPHENE_RAMAN:-$(cd "$WORK/../../../../graphene-raman" 2>/dev/null && pwd)}
  LOGRES=0
  case "$TASK" in
    a)     RES=(--nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G --time=01:00:00); NAME="r10a" ;;
    b)     RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=128G --time=03:00:00); NAME="r10b" ;;
    c0)    RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=128G --time=02:00:00); NAME="r10c0" ;;
    p_ved) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=01:00:00); NAME="r10ved"; LOGRES=1 ;;
    p_c14) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=02:00:00); NAME="r10c14"; LOGRES=1 ;;
    c1f)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=01:00:00); NAME="r10c1f"; LOGRES=1 ;;
    c2)    RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=128G --time=02:00:00); NAME="r10c2" ;;
    c3)    RES=(--nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G --time=01:00:00); NAME="r10c3" ;;
    c6)    RES=(--nodes=1 --ntasks=1 --cpus-per-task=8 --mem=32G --time=01:00:00); NAME="r10c6" ;;
    c1)    NAME="r10c1" ;;
    *) echo "tâche inconnue : $TASK"; exit 1 ;;
  esac
  mkdir -p "$WORK/submitted"; LAST="$WORK/submitted/last"; DIFF="$WORK/submitted/pending_diff.txt"
  {
    echo "=== R10 : soumission de « $TASK $* » le $(date '+%F %T') ; HEAD $(git -C "$PROJ" rev-parse --short HEAD 2>/dev/null)${DEP:+ ; dépendance $DEP}"
    for f in r10_driver.py submit_r10.sh; do
      if [ -f "$LAST/$f" ]; then
        if diff -q "$LAST/$f" "$WORK/$f" >/dev/null; then echo "=== $f : inchangé depuis la soumission $(cat "$LAST/JOBID_OF_SNAPSHOT" 2>/dev/null)"
        else diff -u "$LAST/$f" "$WORK/$f"; fi
      else
        echo "=== $f : première soumission ($(wc -l < "$WORK/$f") lignes, md5 $(md5sum < "$WORK/$f" | cut -c1-12)) — fichier entier dans submitted/<jobid>/"
      fi
    done
  } | tee "$DIFF"
  [ "${DRY:-0}" = 1 ] && { echo "(DRY=1 : rien soumis)"; exit 0; }
  RESD=$(cd "$PROJ" && PYTHONPATH="$PROJ/src:${PYTHONPATH:-}" "$PROJ/.venv/bin/python" -c 'from graphene_raman.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))') \
    || { echo "results_dir illisible"; exit 1; }
  archive() {   # $1 = jobid de référence ; $2… = lignes « id nom » des jobs soumis
    local J=$1; shift
    mkdir -p "$WORK/submitted/$J"; cp -p "$WORK/r10_driver.py" "$WORK/submit_r10.sh" "$WORK/submitted/$J/"; mv "$DIFF" "$WORK/submitted/$J/diff.txt"
    printf '%s\n' "$@" > "$WORK/submitted/$J/jobs.txt"
    echo "$J" > "$WORK/submitted/$J/JOBID_OF_SNAPSHOT"; rm -rf "$LAST"; cp -rp "$WORK/submitted/$J" "$LAST"
  }
  if [ "$TASK" != c1 ]; then
    OUT=(); [ "$LOGRES" = 1 ] && { mkdir -p "$RESD/logs"; OUT=(--output="$RESD/logs/%x_%j.out" --error="$RESD/logs/%x_%j.err"); }
    JID=$(sbatch --parsable --chdir="$WORK" --job-name="$NAME" ${DEP:+--dependency=$DEP} "${RES[@]}" "${OUT[@]}" --export=ALL,GRAPHENE_RAMAN="$PROJ" "$WORK/submit_r10.sh" "$TASK" "$@") \
      || { echo "sbatch a échoué"; exit 1; }
    archive "$JID" "$JID $NAME"
    echo "$JID $NAME $TASK $* ${DEP:+dep=$DEP} $(date '+%F %T')" >> "$WORK/JOBID"
    echo "soumis : $JID ($NAME)"; exit 0
  fi
  # ---- macro c1 : production dans results_dir, lanceurs de scripts/ tels quels (sorties redirigées vers <results_dir>/logs)
  mkdir -p "$RESD/logs"
  C=(--parsable --chdir="$PROJ" --output="$RESD/logs/%x_%j.out" --error="$RESD/logs/%x_%j.err" ${DEP:+--dependency=$DEP} --export=ALL,GRAPHENE_RAMAN="$PROJ")
  L=(); IDS=()
  sub() { local name=$1; shift; local j; j=$("$@") || { echo "sbatch a échoué pour $name (déjà soumis : ${L[*]})"; exit 1; }; L+=("$j $name"); IDS+=("$j"); echo "soumis : $j ($name)"; }
  for S in 5x5 6x6 7x7 8x8 9x9 12x12; do
    sub "specwd_$S" env GRIDS=60,120,240 ETAS=0.05,0.02,0.01 RCUTS=0,1,2,3,4 sbatch "${C[@]}" --job-name="specwd_$S" "$PROJ/scripts/slurm/submit_spectral_wannier_dense.sh" "$S" prod
  done
  sub nkint     sbatch "${C[@]}" --job-name=nkint "$PROJ/scripts/slurm/submit_nkint_check.sh" 9x9
  sub resigma0123 sbatch "${C[@]}" --job-name=resigma "$PROJ/scripts/slurm/submit_rcut_resigma.sh" 9x9 0,1,2,3
  sub resigma4  sbatch "${C[@]}" --job-name=resigma "$PROJ/scripts/slurm/submit_rcut_resigma.sh" 9x9 4
  sub post_loc  sbatch "${C[@]}" --job-name=post_loc "$PROJ/scripts/slurm/submit_post.sh" locality
  sub post_ana  env GOLDEN_RESULT="1.80e-13,r6golden_21852238" sbatch "${C[@]}" --job-name=post_ana "$PROJ/scripts/slurm/submit_post.sh" analyze
  for S in 9x9 6x6 12x12; do
    sub "post_res${S%%x*}" sbatch "${C[@]}" --job-name="post_res${S%%x*}" "$PROJ/scripts/slurm/submit_post.sh" resonance "$S"
  done
  CW=(--parsable --chdir="$WORK" --output="$RESD/logs/%x_%j.out" --error="$RESD/logs/%x_%j.err" ${DEP:+--dependency=$DEP} --export=ALL,GRAPHENE_RAMAN="$PROJ")
  sub r10ved sbatch "${CW[@]}" --job-name=r10ved --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=01:00:00 "$WORK/submit_r10.sh" p_ved
  sub r10c14 sbatch "${CW[@]}" --job-name=r10c14 --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=02:00:00 "$WORK/submit_r10.sh" p_c14
  archive "${IDS[0]}" "${L[@]}"
  for x in "${L[@]}"; do echo "$x c1 ${DEP:+dep=$DEP} $(date '+%F %T')" >> "$WORK/JOBID"; done
  echo "C1_ALL=$(IFS=:; echo "${IDS[*]}")"
  echo "C1_LOC=$(printf '%s\n' "${L[@]}" | awk '$2=="post_loc"{print $1}')"
  exit 0
fi

# ------------------------------------------------------------------ dans le job
WORK=$(pwd)                                                                        # = --chdir de la soumission (répertoire de campagne)
PROJ=${GRAPHENE_RAMAN:?GRAPHENE_RAMAN non défini}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:${PYTHONPATH:-}"      # pas d'install éditable ; h5py de scipy-stack conservé
PY="$PROJ/.venv/bin/python"; DRV="$WORK/r10_driver.py"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 FLEXIBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
TASK=${1:?tâche}; shift
echo "[$(date)] R10 $TASK $* ; job $SLURM_JOB_ID ; pilote md5 $(md5sum < "$DRV" | cut -c1-12) ; HEAD $(git -C "$PROJ" rev-parse --short HEAD)"
threads "$SLURM_CPUS_PER_TASK"
RESD=$(cd "$PROJ" && "$PY" -c 'from graphene_raman.config import load_production, results_dir; print(results_dir(load_production(verbose=False)))')
run() { echo "[$(date)] >>> $*"; (cd "$PROJ" && "$PY" -u "$@") || { rc=$?; echo "[$(date)] ÉCHEC ($rc) : $*"; exit 2; }; }
case "$TASK" in
  a)
    "$PY" -u "$DRV" a0 || { echo "[$(date)] porte A.0 : ÉCHEC, STOP sur A"; exit 3; }
    "$PY" -u "$DRV" a1 || { echo "[$(date)] a1 ÉCHEC"; exit 2; }
    "$PY" -u "$DRV" a2 || { echo "[$(date)] a2 ÉCHEC"; exit 2; }
    "$PY" -u "$DRV" a3; rc=$?; echo "[$(date)] a3 rc=$rc"; exit $rc ;;
  b)
    "$PY" -u "$DRV" b "$@"; rc=$?; echo "[$(date)] b rc=$rc"; exit $rc ;;
  c0|c2|c6)
    "$PY" -u "$DRV" "$TASK"; rc=$?; echo "[$(date)] $TASK rc=$rc"; exit $rc ;;
  c3)
    "$PY" -u "$DRV" c3; r1=$?; echo "[$(date)] c3 rc=$r1"
    "$PY" -u "$DRV" c4; r2=$?; echo "[$(date)] c4 rc=$r2"
    "$PY" -u "$DRV" c5; r3=$?; echo "[$(date)] c5 rc=$r3"
    rc=$r1; [ $r2 -gt $rc ] && rc=$r2; [ $r3 -gt $rc ] && rc=$r3; exit $rc ;;
  p_ved)
    echo "[$(date)] results_dir = $RESD"
    run scripts/t/analyze_Ved.py 5x5,6x6,7x7,8x8,9x9,10x10,11x11,12x12
    cp -p "$PROJ/results/M2/ks_reconstruction.npz" "$RESD/ks_reconstruction.npz" && chmod u+w "$RESD/ks_reconstruction.npz"
    m1=$(md5sum < "$PROJ/results/M2/ks_reconstruction.npz" | cut -d' ' -f1); m2=$(md5sum < "$RESD/ks_reconstruction.npz" | cut -d' ' -f1)
    echo "[$(date)] ks_reconstruction.npz copié de results/M2 : md5 $m1 / $m2"; [ "$m1" = "$m2" ] || { echo "md5 différent : STOP"; exit 2; }
    run scripts/t/sampling_table.py
    echo "[$(date)] p_ved terminé" ;;
  p_c14)
    run scripts/t/resonance_metrics.py --size 9x9 --shift-L-meV 9.05,-9.05 --out "$RESD/resonance_9x9_shiftL.npz"
    echo "[$(date)] p_c14 terminé" ;;
  c1f)
    "$PY" -u "$DRV" c1post || { echo "[$(date)] c1post ÉCHEC (porte de m_rcut_resigma.csv ou lignes de tab:tests_M)"; exit 3; }
    run scripts/fig/make_figures.py --outdir "$WORK/fig"
    run scripts/t/level2_families.py
    run scripts/t/nkint_check_post.py --size 9x9 --nk 150,300,450,600
    run scripts/fig/make_figures_memoire.py --outdir "$WORK/fig"
    # chapitre 5 : options de production de NOTES_EPW (chaîne mv0.02), comme submit_post.sh figures
    run scripts/t/epw_ed_vs_ep.py --selfen results/epw/selfen_240_dg0.02_mv0.02_T300.npz --tag 24k24q_mv0.02
    run scripts/fig/make_figures_epw.py --outdir "$WORK/fig" --prod-tag 240_dg0.02_mv0.02 --phself-tag path_1200_dg0.02_mv0.02 --val-tag 24k24q_mv0.02 \
        --sel-suffix _mv0.02 --control --kohn-val-tags 24k24q,24k24q_mv0.02 --dfpt-tag 24k24q
    echo "[$(date)] c1f terminé" ;;
  *) echo "tâche inconnue $TASK"; exit 1 ;;
esac
