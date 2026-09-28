#!/bin/bash
# R9 -- lanceur unique de la campagne (entrée de campagne, CLAUDE.md « Campagnes de calcul ») ; la tâche est choisie par argument.
#
# Depuis le nœud de connexion :   bash submit_r9.sh <tâche> [arguments]      (DEP=afterok:ID pour une dépendance ; DRY=1 : diff seul)
#   1. affiche le diff de r9_driver.py et de ce fichier depuis la dernière soumission (première soumission : fichiers entiers) ;
#   2. sbatch avec les ressources de la tâche (tableau RES ci-dessous), répertoire de travail = ce répertoire ;
#   3. garde une copie datée des deux fichiers et du diff dans submitted/<jobid>/ ; ligne dans JOBID.
# Dans le job (SLURM_JOB_ID défini) : exécute la tâche.
#
# Tâches :
#   box <S>     A.2 : M^L[1_boîte] dense de la taille S (9x9 ou 5x5) par le noyau partagé de production (ΔV = 1 via le module io substitut),
#               puis contrôles et F_W (a2dpost) ; pour 9x9, d'abord a2c (M^L d'une constante, 9x9 grossière, 4 bandes, contre 81·𝕀).
#               Nœud exclusif 32 rangs × 6 fils (comme scripts/submit_M_dense.sh).
#   b [args]    B : résonance contre N_k^int (args passés à « r9_driver.py b », ex. --variants brut ; --variants aligne,exact). 16 cœurs, 96 Go.
#   a           A.0 (porte) puis A.1 (C_N par taille). 16 cœurs, 128 Go.
#   a3 [args]   A.3 (variantes alignées ; args de « r9_driver.py a3 ») puis a3pole à 1 fil BLAS + 16 fils Python. 16 cœurs, 160 Go.
#   c [args]    C (porte C.0, C.1, C.2 bloc π jusqu'à N = 81, C.3). 16 cœurs, 160 Go.
#   d           D (porte D.0, D.2, D.1). 16 cœurs, 96 Go.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --output=slurm-r9-%x-%j.out
#SBATCH --error=slurm-r9-%x-%j.err

if [ -z "${SLURM_JOB_ID:-}" ]; then
  # ------------------------------------------------------------------ côté nœud de connexion : diff puis soumission
  WORK=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
  TASK=${1:?"tâche : box <S> | b [args] | a | a3 [args] | c [args] | d"}; shift
  case "$TASK" in
    box) S=${1:?"taille : 9x9 | 5x5"}
         case "$S" in 9x9) T=03:00:00 ;; 5x5) T=02:00:00 ;; *) echo "box : taille $S non prévue (9x9, 5x5)"; exit 1 ;; esac
         RES=(--nodes=1 --ntasks=32 --cpus-per-task=6 --exclusive --mem=0 --time=$T); NAME="r9box$S" ;;
    b)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=04:00:00); NAME="r9b" ;;
    a)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=128G --time=03:00:00); NAME="r9a" ;;
    a3)  RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=160G --time=08:00:00); NAME="r9a3" ;;
    c)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=160G --time=08:00:00); NAME="r9c" ;;
    d)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=02:00:00); NAME="r9d" ;;
    *)   echo "tâche inconnue : $TASK"; exit 1 ;;
  esac
  mkdir -p "$WORK/submitted"; LAST="$WORK/submitted/last"; DIFF="$WORK/submitted/pending_diff.txt"
  {
    echo "=== R9 : soumission de « $TASK $* » le $(date '+%F %T') ; HEAD $(git -C "${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}" rev-parse --short HEAD 2>/dev/null)"
    for f in r9_driver.py submit_r9.sh; do
      if [ -f "$LAST/$f" ]; then
        if diff -q "$LAST/$f" "$WORK/$f" >/dev/null; then echo "=== $f : inchangé depuis la soumission $(cat "$LAST/JOBID_OF_SNAPSHOT" 2>/dev/null)"
        else diff -u "$LAST/$f" "$WORK/$f"; fi
      else
        echo "=== $f : première soumission ($(wc -l < "$WORK/$f") lignes, md5 $(md5sum < "$WORK/$f" | cut -c1-12)) — fichier entier dans submitted/<jobid>/"
      fi
    done
  } | tee "$DIFF"
  [ "${DRY:-0}" = 1 ] && { echo "(DRY=1 : rien soumis)"; exit 0; }
  JID=$(sbatch --parsable --chdir="$WORK" --job-name="$NAME" ${DEP:+--dependency=$DEP} "${RES[@]}" "$WORK/submit_r9.sh" "$TASK" "$@") || { echo "sbatch a échoué"; exit 1; }
  mkdir -p "$WORK/submitted/$JID"; cp -p "$WORK/r9_driver.py" "$WORK/submit_r9.sh" "$WORK/submitted/$JID/"; mv "$DIFF" "$WORK/submitted/$JID/diff.txt"
  echo "$JID" > "$WORK/submitted/$JID/JOBID_OF_SNAPSHOT"; rm -rf "$LAST"; cp -rp "$WORK/submitted/$JID" "$LAST"
  echo "$JID $NAME $TASK $* ${DEP:+dep=$DEP} $(date '+%F %T')" >> "$WORK/JOBID"
  echo "soumis : $JID ($NAME)"; exit 0
fi

# ------------------------------------------------------------------ dans le job
WORK=$(pwd)                                                                        # = --chdir de la soumission
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:${PYTHONPATH:-}"      # pas d'install éditable ; h5py de scipy-stack conservé
PY="$PROJ/.venv/bin/python"; DRV="$WORK/r9_driver.py"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 FLEXIBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
TASK=${1:?tâche}; shift
echo "[$(date)] R9 $TASK $* ; job $SLURM_JOB_ID ; pilote md5 $(md5sum < "$DRV" | cut -c1-12) ; HEAD $(git -C "$PROJ" rev-parse --short HEAD) ; noyau : $(grep -n 'inv_sqrtO =' "$PROJ/src/electron_defect_interaction/defects/local_R.py" | tr '\n' ' ')"
case "$TASK" in
  box)
    S=$1
    if [ "$S" = 9x9 ]; then
      threads 1                                                                   # comme le J1 de R6 (compute_ML_R_mpi, 32 rangs)
      srun --cpu-bind=cores "$PY" -u "$DRV" a2c || echo "[$(date)] a2c ÉCHEC rc=$?"
    fi
    threads "$SLURM_CPUS_PER_TASK"                                                # comme scripts/submit_M_dense.sh (32 rangs × 6 fils)
    srun --cpu-bind=cores "$PY" -u "$DRV" a2d --size "$S" --block-size 2000 || { echo "[$(date)] a2d ÉCHEC"; exit 2; }
    threads 32
    "$PY" -u "$DRV" a2dpost --size "$S" || { echo "[$(date)] a2dpost ÉCHEC"; exit 3; }
    echo "[$(date)] box $S terminé" ;;
  b)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" b "$@"; rc=$?; echo "[$(date)] b rc=$rc"; exit $rc ;;
  a)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" a0 || { echo "[$(date)] porte A.0 : ÉCHEC, STOP sur A"; exit 3; }
    "$PY" -u "$DRV" a1; rc=$?; echo "[$(date)] a1 rc=$rc"; exit $rc ;;
  a3)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" a3 "$@" || { echo "[$(date)] a3 ÉCHEC"; exit 2; }
    threads 1; export R4_EIG_WORKERS=$SLURM_CPUS_PER_TASK                          # eig par lots : 1 fil BLAS + fils Python (leçon R4/R6)
    "$PY" -u "$DRV" a3pole; rc=$?; echo "[$(date)] a3pole rc=$rc"; exit $rc ;;
  c)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" c "$@"; rc=$?; echo "[$(date)] c rc=$rc"; exit $rc ;;
  d)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" d; rc=$?; echo "[$(date)] d rc=$rc"; exit $rc ;;
  *) echo "tâche inconnue $TASK"; exit 1 ;;
esac
