#!/bin/bash
# R10 -- lanceur unique de la campagne (entrée de campagne, CLAUDE.md « Campagnes de calcul ») ; la tâche est choisie par argument.
#
# Depuis le nœud de connexion :   bash submit_r10.sh <tâche> [arguments]      (DEP=afterany:ID pour une dépendance ; DRY=1 : diff seul)
#   1. affiche le diff de r10_driver.py et de ce fichier depuis la dernière soumission (première soumission : md5 et nombre de lignes) ;
#   2. sbatch avec les ressources de la tâche, répertoire de travail = ce répertoire ;
#   3. garde une copie des deux fichiers et du diff dans submitted/<jobid>/ ; ligne dans JOBID.
# Dans le job (SLURM_JOB_ID défini) : exécute la tâche.
#
# Tâches (GO 1) :
#   a           J1 : a0 (porte A.0 ; échec : STOP sur A) puis a1, a2, a3. 1 tâche, 8 cœurs, 64 Go, 1 h.
#   b [args]    J2 : B.0 (porte), B.1, B.2 (args passés à « r10_driver.py b », ex. --parts b1,b2). 16 cœurs, 128 Go, 3 h.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --output=slurm-r10-%x-%j.out
#SBATCH --error=slurm-r10-%x-%j.err

if [ -z "${SLURM_JOB_ID:-}" ]; then
  # ------------------------------------------------------------------ côté nœud de connexion : diff puis soumission
  WORK=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
  TASK=${1:?"tâche : a | b [args]"}; shift
  case "$TASK" in
    a) RES=(--nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G --time=01:00:00); NAME="r10a" ;;
    b) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=128G --time=03:00:00); NAME="r10b" ;;
    *) echo "tâche inconnue : $TASK"; exit 1 ;;
  esac
  mkdir -p "$WORK/submitted"; LAST="$WORK/submitted/last"; DIFF="$WORK/submitted/pending_diff.txt"
  PROJ=${GRAPHENE_RAMAN:-$(cd "$WORK/../../../../graphene-raman" 2>/dev/null && pwd)}
  {
    echo "=== R10 : soumission de « $TASK $* » le $(date '+%F %T') ; HEAD $(git -C "$PROJ" rev-parse --short HEAD 2>/dev/null)"
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
  JID=$(sbatch --parsable --chdir="$WORK" --job-name="$NAME" ${DEP:+--dependency=$DEP} "${RES[@]}" --export=ALL,GRAPHENE_RAMAN="$PROJ" "$WORK/submit_r10.sh" "$TASK" "$@") \
    || { echo "sbatch a échoué"; exit 1; }
  mkdir -p "$WORK/submitted/$JID"; cp -p "$WORK/r10_driver.py" "$WORK/submit_r10.sh" "$WORK/submitted/$JID/"; mv "$DIFF" "$WORK/submitted/$JID/diff.txt"
  echo "$JID" > "$WORK/submitted/$JID/JOBID_OF_SNAPSHOT"; rm -rf "$LAST"; cp -rp "$WORK/submitted/$JID" "$LAST"
  echo "$JID $NAME $TASK $* ${DEP:+dep=$DEP} $(date '+%F %T')" >> "$WORK/JOBID"
  echo "soumis : $JID ($NAME)"; exit 0
fi

# ------------------------------------------------------------------ dans le job
WORK=$(pwd)                                                                        # = --chdir de la soumission
PROJ=${GRAPHENE_RAMAN:?GRAPHENE_RAMAN non défini}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:${PYTHONPATH:-}"      # pas d'install éditable ; h5py de scipy-stack conservé
PY="$PROJ/.venv/bin/python"; DRV="$WORK/r10_driver.py"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 FLEXIBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
TASK=${1:?tâche}; shift
echo "[$(date)] R10 $TASK $* ; job $SLURM_JOB_ID ; pilote md5 $(md5sum < "$DRV" | cut -c1-12) ; HEAD $(git -C "$PROJ" rev-parse --short HEAD)"
threads "$SLURM_CPUS_PER_TASK"
case "$TASK" in
  a)
    "$PY" -u "$DRV" a0 || { echo "[$(date)] porte A.0 : ÉCHEC, STOP sur A"; exit 3; }
    "$PY" -u "$DRV" a1 || { echo "[$(date)] a1 ÉCHEC"; exit 2; }
    "$PY" -u "$DRV" a2 || { echo "[$(date)] a2 ÉCHEC"; exit 2; }
    "$PY" -u "$DRV" a3; rc=$?; echo "[$(date)] a3 rc=$rc"; exit $rc ;;
  b)
    "$PY" -u "$DRV" b "$@"; rc=$?; echo "[$(date)] b rc=$rc"; exit $rc ;;
  *) echo "tâche inconnue $TASK"; exit 1 ;;
esac
