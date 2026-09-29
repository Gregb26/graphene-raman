# R9 — Contrôles avant le chapitre 4 (alignement de ΔV, résonance contre N_k^int, chaîne repliée contre R7, grandeur de Kaasbjerg)

- But : (A) décalage uniforme de ΔV^L dans M2 et variantes alignées ; (B) position de la résonance quand N_k^int converge ; (C) chaîne repliée N×N contre R7 et limite diluée ; (D) grandeur comparée aux ~70 eV Å² de Kaasbjerg. Clôture : rejeu avec C_N = moyenne du plateau (i), test intensif hors k = k′. Audit : image minimale dans la cellule à 60° (écarts 1–3).
- Prompts d'origine : « R9 — Contrôles avant le chapitre 4 » (Greg, 2026-09-27) ; « R9 — Clôture » et « Audit — image minimale » (Greg, 2026-09-28).
- Statut : **TEST** (post-traitement seul, aucun calcul QE ; résultats consignés dans `R9_rapport.md` ; rien promu en production).
- Dates : phase 0 2026-09-27 ; étape G et GO 2026-09-28 ; A, B, C, D terminés 2026-09-28 ; clôture (phase 0, GO, rejeu) et audit 2026-09-28 ; `cache/` supprimé 2026-09-29 (GO de Greg) ; R9 clos 2026-09-29.
- Jobs (`JOBID`, diffs dans `submitted/<jobid>/`) : 21954787, 21954788 (M^L[1_boîte] 9×9, 5×5), 21954789 (B), 21955512 (A, échec d'exécution), 21955644 (A), 21955645 (A.3), 21955646 (B aligné), 21955647 (C), 21955648 (D) ; clôture 21976369 (R.0), 21976370 (A.3 plateau), 21976371 (B plateau).
- Pilote et lanceur : `r9_driver.py`, `submit_r9.sh` ; fonctions P1–P4 dans `src/` et `tests/test_r9_functions.py` (dépôt, commités dans cb7241d) ; sous-commande `audit` (nœud de connexion, sans job).
- Sorties principales : `R9_rapport.md` ; `a/` (A.1, A.3, porte R.0), `a2/`, `b/`, `c/`, `d/` (json, tables md) ; `cloture/synthese.md` ; `audit_image_minimale.md`, `audit/audit_results.json` ; `fig/` (offset_profiles, resonance_vs_nkint, rcut_aligned, folded_vs_R7, kaasbjerg_fig3_map) ; `manifeste_R9.md`.
- Copie versionnée : `article/R9_controles/` (rapport, README, manifeste, pilote, lanceur, tables, json, figures ; ni npz, ni slurm ; `cache/` supprimé, reconstructible par le manifeste).
- Lecture :
