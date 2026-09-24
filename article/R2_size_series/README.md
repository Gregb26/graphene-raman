# R2 — série en taille 5×5 → 12×12 de la lacune relaxée (préparation article, hors mémoire)

Copie versionnée (règle 5 de CLAUDE.md, § « Campagnes de calcul ») du répertoire de travail
`graphene/qe/defects/super_cell_relaxed/series/` (hors dépôt, à côté de celui-ci). Statut : PRODUCTION.
Versionné ici : `make_inputs_series.py`, `analyze_relax_series.py`, `R2_phase0_table.txt`, `R2_jobs.md`, `R2_rapport.md`,
`analyse_R2_2026-09-24.txt`, `NxN/perturbation.log`, `NxN/nspin{1,2}/{relax.in,relax.out,submit.relax}` (relax.out de pw.x ≤ 2.3 Mo),
`NxN/nspin2/projwfc/{projwfc.in,submit.projwfc}`.
Non versionné : `JOBID`, `slurm-*`, `projwfc.out`, `pdos_*`, `proj_*`, et les `.save` (scratch `qe_tmp/vacancy_relaxed/series/NxN/nspinN/`,
miroir `graphene/qe/qe_tmp_backup/vacancy_relaxed/series/` avec `MD5SUMS_series_2026-09-24.txt`).
Le 9×9 de la série est celui de R1 (`article/R1_vacancy_relaxed/`, répertoire de travail `super_cell_relaxed/9x9/`).
Les scripts portent des chemins absolus vers `graphene/qe/` et s'exécutent depuis le répertoire de travail.
