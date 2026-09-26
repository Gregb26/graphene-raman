#!/bin/bash
# R7 -- miroir des wfc1.hdf5 15x15/18x18 vers qe_tmp_backup (decision de Greg 2026-09-25), md5 source et miroir, groupe du projet.
#SBATCH --account=rrg-cotemich-ac
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=slurm-r7-%x-%j.out
#SBATCH --error=slurm-r7-%x-%j.err
S=$HOME/links/scratch/qe_tmp; B=${PROJECTS}/graphene/qe/qe_tmp_backup; W=${SLURM_SUBMIT_DIR}
cd "$B" || exit 1
for c in defect_15x15_p defect_15x15_d defect_18x18_p defect_18x18_d; do
  echo "[$(date)] rsync $c/wfc1.hdf5"
  rsync -a -r --no-o --no-g --open-noatime --chmod=Dg+s "$S/$c/$c.save/wfc1.hdf5" "$B/$c/$c.save/" || { echo "RSYNC ECHEC $c"; exit 2; }
  chgrp rrg-cotemich-ac "$B/$c/$c.save/wfc1.hdf5"
done
echo "[$(date)] md5 source"
for c in defect_15x15_p defect_15x15_d defect_18x18_p defect_18x18_d; do (cd "$S/$c" && md5sum "$c.save/wfc1.hdf5" | sed "s|^\([0-9a-f]*\)  |\1  $c/|"); done > "$W/MD5SUMS_R7_wfc_source.txt"
echo "[$(date)] md5 miroir"
for c in defect_15x15_p defect_15x15_d defect_18x18_p defect_18x18_d; do md5sum "$c/$c.save/wfc1.hdf5"; done > "$B/MD5SUMS_R7_wfc_2026-09-25.txt"
diff <(sort "$W/MD5SUMS_R7_wfc_source.txt") <(sort "$B/MD5SUMS_R7_wfc_2026-09-25.txt") && echo "[$(date)] MD5 SOURCE = MIROIR : OK (4/4)" || { echo "[$(date)] MD5 DIFFERENT"; exit 3; }
cat "$B/MD5SUMS_R7_wfc_2026-09-25.txt" >> "$B/MD5SUMS_R7_2026-09-25.txt"
du -sh "$B"/defect_15x15_? "$B"/defect_18x18_?
echo "[$(date)] done"
