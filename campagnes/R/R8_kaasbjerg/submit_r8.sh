#!/bin/bash
# R8 -- lanceur unique de la campagne (CLAUDE.md « Campagnes de calcul ») ; la tâche est choisie par argument.
#
# Depuis le nœud de connexion :   bash submit_r8.sh <tâche>      (DEP=afterok:ID pour une dépendance ; DRY=1 : diff seul)
#   1. affiche le diff de r8_driver.py et de ce fichier depuis la dernière soumission (première soumission : fichiers entiers) ;
#   2. sbatch avec les ressources de la tâche, répertoire de travail = ce répertoire ;
#   3. garde une copie des deux fichiers et du diff dans submitted/<jobid>/ ; ligne dans JOBID.
# Dans le job (SLURM_JOB_ID défini) : exécute la tâche.
#
# Tâches :
#   g0     prep (M2 9x9 et 12x12 -> V_loc) puis tous les g0 de la campagne (local_green_batch, 16 fils BLAS).
#   gate   étape 1 : portes P1-P4 et normes d'information (1 fil BLAS + 16 fils Python) ; code de sortie 3 si une porte échoue.
#   dos    étapes 2 et 2 bis.      spec   étape 3.      sens   étape 5.      fig    7b, figures, tables.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --output=slurm-r8-%x-%j.out
#SBATCH --error=slurm-r8-%x-%j.err

if [ -z "${SLURM_JOB_ID:-}" ]; then
  # ------------------------------------------------------------------ côté nœud de connexion : diff puis soumission
  WORK=$(cd "$(dirname "$(readlink -f "$0")")" && pwd)
  TASK=${1:?"tâche : g0 | gate | dos | spec | sens | fig"}; shift
  case "$TASK" in
    g0)   RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=160G --time=05:00:00); NAME="r8g0" ;;
    gate) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=03:00:00); NAME="r8gate" ;;
    dos)  RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=02:00:00); NAME="r8dos" ;;
    spec) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=01:00:00); NAME="r8spec" ;;
    sens) RES=(--nodes=1 --ntasks=1 --cpus-per-task=16 --mem=96G --time=03:00:00); NAME="r8sens" ;;
    fig)  RES=(--nodes=1 --ntasks=1 --cpus-per-task=4 --mem=32G --time=00:30:00); NAME="r8fig" ;;
    *)    echo "tâche inconnue : $TASK"; exit 1 ;;
  esac
  mkdir -p "$WORK/submitted"; LAST="$WORK/submitted/last"; DIFF="$WORK/submitted/pending_diff.txt"
  {
    echo "=== R8 : soumission de « $TASK » le $(date '+%F %T') ; HEAD $(git -C "${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}" rev-parse --short HEAD 2>/dev/null)"
    for f in r8_driver.py submit_r8.sh; do
      if [ -f "$LAST/$f" ]; then
        if diff -q "$LAST/$f" "$WORK/$f" >/dev/null; then echo "=== $f : inchangé depuis la soumission $(cat "$LAST/JOBID_OF_SNAPSHOT" 2>/dev/null)"
        else diff -u "$LAST/$f" "$WORK/$f"; fi
      else
        echo "=== $f : première soumission ($(wc -l < "$WORK/$f") lignes, md5 $(md5sum < "$WORK/$f" | cut -c1-12)) — fichier entier dans submitted/<jobid>/"
      fi
    done
  } | tee "$DIFF"
  [ "${DRY:-0}" = 1 ] && { echo "(DRY=1 : rien soumis)"; exit 0; }
  JID=$(sbatch --parsable --chdir="$WORK" --job-name="$NAME" ${DEP:+--dependency=$DEP} "${RES[@]}" "$WORK/submit_r8.sh" "$TASK" "$@") || { echo "sbatch a échoué"; exit 1; }
  mkdir -p "$WORK/submitted/$JID"; cp -p "$WORK/r8_driver.py" "$WORK/submit_r8.sh" "$WORK/submitted/$JID/"; mv "$DIFF" "$WORK/submitted/$JID/diff.txt"
  echo "$JID" > "$WORK/submitted/$JID/JOBID_OF_SNAPSHOT"; rm -rf "$LAST"; cp -rp "$WORK/submitted/$JID" "$LAST"
  echo "$JID $NAME $TASK ${DEP:+dep=$DEP} $(date '+%F %T')" >> "$WORK/JOBID"
  echo "soumis : $JID ($NAME)"; exit 0
fi

# ------------------------------------------------------------------ dans le job
WORK=$(pwd)                                                                        # = --chdir de la soumission
PROJ=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}
module restore qe; module load mpi4py/4.0.3 scipy-stack
export GRAPHENE_RAMAN="$PROJ"; export PYTHONPATH="$PROJ/src:${PYTHONPATH:-}"      # pas d'install éditable ; h5py de scipy-stack conservé
PY="$PROJ/.venv/bin/python"; DRV="$WORK/r8_driver.py"
threads() { export OMP_NUM_THREADS=$1 OPENBLAS_NUM_THREADS=$1 FLEXIBLAS_NUM_THREADS=$1 MKL_NUM_THREADS=$1; }
pool() { threads 1; export R8_WORKERS=$SLURM_CPUS_PER_TASK R4_EIG_WORKERS=$SLURM_CPUS_PER_TASK; }   # boucles d'énergie : leçon R4/R6
TASK=${1:?tâche}; shift
echo "[$(date)] R8 $TASK ; job $SLURM_JOB_ID ; pilote md5 $(md5sum < "$DRV" | cut -c1-12) ; HEAD $(git -C "$PROJ" rev-parse --short HEAD)"
case "$TASK" in
  g0)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" prep || { echo "[$(date)] prep ÉCHEC"; exit 2; }
    "$PY" -u "$DRV" g0; rc=$?; echo "[$(date)] g0 rc=$rc"; exit $rc ;;
  gate|dos|spec|sens)
    pool
    "$PY" -u "$DRV" "$TASK"; rc=$?; echo "[$(date)] $TASK rc=$rc"; exit $rc ;;
  fig)
    threads "$SLURM_CPUS_PER_TASK"
    "$PY" -u "$DRV" fig; rc=$?; echo "[$(date)] fig rc=$rc"; exit $rc ;;
  *) echo "tâche inconnue $TASK"; exit 1 ;;
esac
