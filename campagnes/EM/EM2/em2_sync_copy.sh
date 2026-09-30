#!/bin/bash
# EM2: (re)synchronise the versioned copy campagnes/EM/EM2/ from this working directory (rule 5 of CLAUDE.md), cp -p,
# explicit list: inputs, scripts, reports, analyses, small outputs. Not copied: .chk/.eig (= reference), .mmn link,
# slurm outputs (em2_*_<jobid>.out/.err), JOBID_*, bands_pp.out and em2_bands.dat* (bands.x side outputs),
# postw90 jdos / kubo_A files and the cancelled run's .wpout.
set -e
cd "$(dirname "$0")"
DST=${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}/campagnes/EM/EM2
mkdir -p "$DST/postw90"
cp -p README.md EM2_rapport.md MD5SUMS_EM2_2026-09-29.txt em2_sync_copy.sh \
      em2_kpoints.py em2_kpoints.log em2_kpoints_crystal.txt em2_kpoints_table.txt \
      bands.in bands.in.diff bands_pp.in submit_bands.sh bands.out em2_p_avg.dat \
      em2_A_compare.py em2_A_compare.log em2_A_summary.txt em2_A_table.txt em2_A.npz em2_A_figure.py em2_A_anneaux.pdf em2_A_anneaux.png \
      em2_B_ours.py submit_B_ours.sh em2_B_ours.log em2_B_ours.npz em2_B_compare.py em2_B_compare.log em2_B_compare_*.txt \
      em2_B_decompose.py em2_B_decomposition.txt em2_B_figure.py em2_B_sigma.pdf em2_B_sigma.png em2_postw90_sigma.npz "$DST/"
cp -p postw90/mkwin.sh postw90/submit_postw90.sh postw90/wannier.win.ref "$DST/postw90/"
for d in k301 k301_nows k301_ti k301_ti_nows k1201 k1201_ti k1201_ti_nows; do
  mkdir -p "$DST/postw90/$d"
  cp -p postw90/$d/wannier.win postw90/$d/wannier.win.diff postw90/$d/wannier.wpout postw90/$d/wannier-kubo_S_*.dat "$DST/postw90/$d/"
done
find "$DST" -type f -size +5M -print | sed 's/^/TOO LARGE: /'
echo "synced to $DST: $(find "$DST" -type f | wc -l) files, $(du -sh --apparent-size "$DST" | cut -f1)"
