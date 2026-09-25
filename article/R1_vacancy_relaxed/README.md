# R1 / R1b — relaxation de la supercellule 9x9 avec lacune (préparation article, hors mémoire)

Copie versionnée des inputs, scripts et rapports de `graphene/qe/defects/super_cell_relaxed/9x9/` (déplacé le 2026-09-23 depuis `graphene/qe/vacancy_relaxed/`) (hors dépôt, à côté de celui-ci, même
convention que `graphene/qe/epw/` pour le chapitre 5). Les `relax.out`/`scf.out` de pw.x (< 5 Mo) sont versionnés depuis le 2026-09-23 ; les autres sorties (`projwfc.out`, `slurm-*`, `pdos_*`,
`proj_*`) et les `.save` ne sont pas versionnés : `.save` de R1 sur scratch `qe_tmp/vacancy_relaxed/nspin{1,2}` et miroir
`graphene/qe/qe_tmp_backup/vacancy_relaxed/` (md5) ; `.save` de R1b (`k3x3/`) sur scratch seulement (tests reproductibles).
Les scripts `make_inputs.py`, `analyze_relax.py`, `k3x3/make_k3x3.py`, `k3x3/analyze_k3x3.py` portent des chemins absolus
vers `graphene/qe/` et s'exécutent depuis le répertoire d'origine.
- 2026-09-25 (R4, J0) : potentiels locaux pp.x (plot_num 1) des géométries relaxées, calculés pour la campagne R4 (`graphene/qe/defects/R4_quasi_lie/`) : `nspin1/Vks_R1_nspin1` (spin_component 0) et `nspin2/Vks_R1_nspin2_{up,dw}` (spin_component 1, 2), entrées `pp.in`, `pp_up.in`, `pp_dw.in`, `submit.pp`, jobs 21796633 / 21796634 ; fichiers Vks (241 Mo) et `pp*.out` (184 Mo) non versionnés.
