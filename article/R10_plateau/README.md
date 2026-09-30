# R10 — Base cohérente pour les chiffres du chapitre 4 (alignement par la moyenne du plateau, argument E_res, familles, rejeu de production aligné)

- But : mettre tous les chiffres du ch. 4 sur une seule base — ΔV aligné par la moyenne du plateau (i) de R9 pour toutes les tailles (A), argument E_res et familles au niveau 1 (B), rejeu de la production avec M_W(R,R) − C_N·𝕀₅ (approximation (i)) et les corrections de l'audit de l'image minimale (C).
- Prompt d'origine : « R10 — Base cohérente pour les chiffres du chapitre 4 » (Greg, 2026-09-29).
- Statut : **TEST** (post-traitement seul, aucun calcul QE, aucun `.save`) ; sorties de C = production du ch. 4 dans `results/M2_plateau` (produits seulement ; matrices M2 brutes dans `results/M2`).
- Dates : phase 0, étape G et GO 1 (A, B) 2026-09-29 ; (a) `c61c118`, (b) `ea91d61`, (c) `23ee3bb` ; GO 2 (C) soumis le 2026-09-29, terminé le 2026-09-30 ; figures installées et **campagne close** le 2026-09-30.
- Jobs (`JOBID`, diffs dans `submitted/<jobid>/`) : 22041340 `r10a` (A), 22041341 `r10b` (B) ; 22058838 `r10c0` (C.0), 22058839–54 (C.1, lanceurs de `scripts/`), 22058855 `r10c1f`, 22058856 `r10c2`, 22058857 `r10c3` (C.3–C.5 ; porte des anneaux refusée à 1e-13, acceptée au seuil de 1e-10 fixé par Greg, sans recalcul), 22058858 `r10c6` (C.6).
- Pilote et lanceur : `r10_driver.py`, `submit_r10.sh` ; fonctions de 0.5 (a) écrites par Greg (D1) ; changements (b) de `scripts/` écrits par Code (section « (b) » du rapport).
- Sorties : `R10_rapport.md` (phase 0, étape G, A et B, (b), (c), GO 2, rapport de C) ; `a/` (C_N des 13 tailles : `a1_results.json`), `b/`, `c/` (portes, tab:rcut_M à trois colonnes, Kaasbjerg, contrôles de l'audit, `table_v2_plateau.md`, C.6), `fig/` (figures de A, B, C, figures de production régénérées, `avant_apres/`) ; `manifeste_R10.md` (cache/).
- Copie versionnée : `article/R10_plateau/` (rapport, README, pilote, lanceur, tables, json, npz < 5 Mo, figures ; jamais slurm ni `cache/`).
- Lecture :
