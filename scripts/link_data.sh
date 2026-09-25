#!/bin/bash
# Assemble the data/ tree expected by the scripts (run.py, test_ks_reconstruction.py)
# entirely out of symlinks, so the multi-GB .save directories are never copied.
#
# Layout produced (relative to the project root):
#   data/graphene/unit_cell/qe/defect_<N>.save/   -> wavefunctions (scratch) + Vks_uc + C.upf
#   data/graphene/supercell/qe/defect_<N>_p.save/ -> wavefunctions (scratch) + Vks_<N>_p + C.upf
#   data/graphene/supercell/qe/defect_<N>_d.save/ -> wavefunctions (scratch) + Vks_<N>_d + C.upf
#
# Usage:  bash scripts/link_data.sh 5x5      (repeat for 6x6, 7x7, ... 12x12)
set -euo pipefail

N=${1:?usage: link_data.sh NxN  (e.g. 5x5)}

ROOT=${PROJECTS:-/home/gregb26/links/projects/rrg-cotemich-ac/gregb26}
PROJ=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)   # racine du dépôt, déduite de l'emplacement du script
SCRATCH=/home/gregb26/links/scratch/qe_tmp
GRAPHENE=$ROOT/graphene/qe/defects

# Find the .save directory inside a scratch run dir (naming is not 100% consistent).
find_save () {
    local rundir="$1"
    local save
    save=$(find "$rundir" -maxdepth 2 -name '*.save' -type d 2>/dev/null | head -1)
    [ -n "$save" ] && echo "$save"
}

# $1 = scratch run dir, $2 = target data .save dir, $3 = Vks source path, $4 = Vks link name
link_save () {
    local rundir="$1" dst="$2" vks_src="$3" vks_name="$4"
    local src
    src=$(find_save "$rundir")
    if [ -z "$src" ]; then echo "SKIP  $dst : pas de .save dans $rundir"; return; fi
    mkdir -p "$dst"
    # Link every file of the scratch .save (wfc*.hdf5, charge-density, data-file-schema.xml, C.upf)
    local f
    for f in "$src"/*; do ln -sfn "$f" "$dst/$(basename "$f")"; done
    # Link the local KS potential (pp.x dump); warn if it has not been generated yet
    if [ -e "$vks_src" ]; then
        ln -sfn "$vks_src" "$dst/$vks_name"
    else
        echo "WARN  $vks_name absent ($vks_src) -- genere-le via pp.x"
    fi
    echo "OK    $dst"
}

UC=$PROJ/data/graphene/unit_cell/qe
SC=$PROJ/data/graphene/supercell/qe

link_save "$SCRATCH/defect_unit_cell_${N}" "$UC/defect_${N}.save" \
          "$GRAPHENE/unit_cell/${N}/Vks_uc_${N}" "Vks_uc"
link_save "$SCRATCH/defect_${N}_p" "$SC/defect_${N}_p.save" \
          "$GRAPHENE/super_cell/${N}/pristine/Vks_${N}_p" "Vks_${N}_p"
link_save "$SCRATCH/defect_${N}_d" "$SC/defect_${N}_d.save" \
          "$GRAPHENE/super_cell/${N}/defective/Vks_${N}_d" "Vks_${N}_d"

echo "Termine pour ${N}."
