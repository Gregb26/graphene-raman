#!/bin/bash
# R6 étape 3 — installation des figures corrigées dans figures/ du dépôt. À lancer par Greg (Code n'écrase pas figures/).
# Les versions de etape3/figures_fix/ viennent des scripts du dépôt après correction des limites d'axe (2026-09-27) :
#   fig_spectral_final   (a) : bornes tirées des données (4e2–8e4 meV, calées sur v1, coupaient Born v2) ; légende au-dessus des courbes
#   fig_M_scaling_final      : 0 → 1,25 max (5–9 eV, calé sur v1 : figure vide avec M2)
#   fig_epw_vs_ed (ch. 5)    : Γ^ed de M2 (options de production NOTES_EPW) ; place au-dessus des courbes pour la légende
# Les quatre figures EPW génériques réécrites par l'ancienne chaîne (job 21872955) ont déjà été rétablies (git checkout, 2026-09-27 14:39).
# Usage : bash install_figures_R6.sh [--pdf-dates]
#   --pdf-dates : rétablit aussi depuis HEAD les PDF dont seule /CreationDate diffère (PNG identiques au bit : fig_Ved*, fig_ks_reconstruction,
#                 fig_epw_kohn_degauss), pour que le commit ne contienne que les figures qui changent.
set -euo pipefail
FX="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)/figures_fix"
cd "${GRAPHENE_RAMAN:?GRAPHENE_RAMAN non défini}"
for b in fig_spectral_final fig_M_scaling_final fig_epw_vs_ed; do cp -p "$FX/$b.pdf" "$FX/$b.png" figures/; echo "installé : figures/$b.{pdf,png}"; done
if [ "${1:-}" = "--pdf-dates" ]; then
  git checkout -- figures/fig_epw_kohn_degauss.pdf figures/fig_ks_reconstruction.pdf figures/fig_Ved.pdf \
                  figures/fig_Ved_{boundary,map,profile_mean,radial,radial_masked,zoom}.pdf
  echo "PDF à date seule rétablis depuis HEAD"
fi
A=article/R6_production_corrigee/etape3/figures_v2
mkdir -p "$A" && cp -p figures/*.pdf figures/*.png "$A/" && echo "copie versionnée : $A"
git status --short figures/
