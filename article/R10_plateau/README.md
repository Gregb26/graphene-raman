# R10 — Base cohérente pour les chiffres du chapitre 4 (alignement par la moyenne du plateau, argument E_res, familles, rejeu de production aligné)

- But : mettre tous les chiffres du ch. 4 sur une seule base — ΔV aligné par la moyenne du plateau (i) de R9 pour toutes les tailles (A), argument E_res et familles au niveau 1 (B), rejeu de la production avec M_W(R,R) − C_N·𝕀₅ (approximation (i)) et les corrections de l'audit de l'image minimale (C).
- Prompt d'origine : « R10 — Base cohérente pour les chiffres du chapitre 4 » (Greg, 2026-09-29).
- Statut : **TEST** (post-traitement seul, aucun calcul QE, aucun `.save`) ; les sorties de C deviendront la production du ch. 4 dans un nouveau répertoire de résultats (proposé : `results/M2_plateau`, à confirmer au GO 2).
- Dates : phase 0, étape G et GO 1 (A, B) 2026-09-29.
- Jobs (`JOBID`, diffs dans `submitted/<jobid>/`) : 22041340 `r10a` (A.0–A.3), 22041341 `r10b` (B.0–B.2).
- Pilote et lanceur : `r10_driver.py`, `submit_r10.sh` ; fonctions de 0.5 (a) réservées à Greg (signatures acceptées, D1).
- Sorties : `R10_rapport.md` (phase 0, étape G, rapport intermédiaire A et B) ; `a/` (C_N des 13 tailles : `a1_results.json`), `b/`, `fig/` (offset_profiles_13, levels_vs_invN, eres_vs_grid).
- Copie versionnée : `article/R10_plateau/` (rapport, README, pilote, lanceur, tables, json, npz < 5 Mo, figures ; jamais slurm ni `cache/`).
- Lecture :
