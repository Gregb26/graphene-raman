# R4 — Diagnostic de l'état quasi-lié de la lacune (préparation article, hors mémoire)

- But : savoir où la DFT place les états de la lacune (π et σ) de la 9×9 et quelle étape de la chaîne
  M → Wannier → matrice T s'en écarte. Le facteur α de D3 est un diagnostic, jamais un paramètre.
- Prompt d'origine : R4 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) D6 → D1 → D4 → D2 → D3 → D5 → R → rapport → STOP.
- Statut : **TEST** (post-traitement seulement ; aucun calcul pw.x ; données de production en lecture seule).
- Dates : phase 0 le 2026-09-25 (matin) ; GO reçu le 2026-09-25 ; jobs J0–J6 exécutés et rapport terminé le 2026-09-25 (STOP).
- Répertoire de travail (celui-ci, hors dépôt, règle 5 de CLAUDE.md) : `graphene/qe/defects/R4_quasi_lie/`, à côté de
  `super_cell/9x9/` (données du ch. 4) et de `super_cell_relaxed/9x9/` (R1) qu'il prolonge. Copie versionnée prévue :
  `graphene-raman/article/R4_quasi_lie/` (rapport `R4_rapport.md`, pilote, tables, figures ; jamais les npz > 5 Mo ni les `slurm-*`).
- Contenu : `R4_rapport.md` (rapport complet : phase 0, méthodes, D6, D1, D4, D2, D3, D5, R), pilote unique `r4_driver.py`
  (sous-commandes prep, d6, d1, d5, r, d4, g0, d3, tables, figs) + `submit_r4.sh`, résultats `prep/ d6/ d1/ d4/ d3/ d5/ r/` (json, npz),
  `sections/` (rédaction + `tables.md` généré), `fig/`, caches `cache/` (3.4G, non versionnés), `slurm-*`, `JOBID`, `failed_runs/`.
- Jobs : J0 21796633/4 (pp.x, dans R1), J1 21796852, J2 21796853, J3 21796854, J4 21797745 → 21797985 (J6, avec d1), J5a 21797930, J5b 21797931.
- Aucune suppression ; tout nettoyage passe par un manifeste et un GO séparé. Git en lecture seulement.
