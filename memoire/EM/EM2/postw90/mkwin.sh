#!/bin/bash
# EM2 B: build DIR/wannier.win = wannier.win.ref (reference .win, read only copy) + postw90 kubo block after 'write_hr',
# with berry_kmesh = N N 1 and optional extra lines; links .chk/.eig/.mmn; writes DIR/wannier.win.diff.
# Usage: ./mkwin.sh DIR N ["extra line" ...]
set -e
d=$1; N=$2; shift 2
[[ $N =~ ^[0-9]+$ ]] || { echo "N must be an integer, got '$N'"; exit 1; }
mkdir -p "$d"
{
  sed -n '1,/^ *write_hr = .true./p' wannier.win.ref
  cat <<BLOCK
! --- EM2 B: postw90 kubo; Gaussian of std 0.04 eV = W90 width sqrt(2)*0.04 (utility_w0gauss: exp(-x^2)/sqrt(pi)) ---
 berry = true
 berry_task = kubo
 berry_kmesh = $N $N 1
 fermi_energy = -4.238895
 kubo_freq_min = 0.2
 kubo_freq_max = 6.0
 kubo_freq_step = 0.01
 kubo_adpt_smr = false
 kubo_smr_type = gauss
 kubo_smr_fixed_en_width = 0.0565685425
 kubo_eigval_max = 1000.0
BLOCK
  for l in "$@"; do echo " $l"; done
  sed '1,/^ *write_hr = .true./d' wannier.win.ref
} > "$d/wannier.win"
(cd "$d" && ln -sf ../wannier.chk ../wannier.eig ../wannier.mmn .)
diff wannier.win.ref "$d/wannier.win" > "$d/wannier.win.diff" || true
