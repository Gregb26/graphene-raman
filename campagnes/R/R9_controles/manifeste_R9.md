# Manifeste de fin — R9

**Mise à jour 2026-09-29 : `cache/` supprimé sur GO de Greg ; R9 clos (voir la dernière section).**

Campagne R9 (contrôles avant le chapitre 4, puis clôture « rejeu C_N = plateau »), statut **TEST**. Établi le 2026-09-28 après la clôture.
**Aucun fichier n'est promu en production** : rien n'est écrit dans `results/`, `config/`, `scripts/` ; les fonctions P1–P4 de `src/` sont des fonctions de bibliothèque
(commit cb7241d), pas des résultats. La suppression de `cache/` fait l'objet d'un GO séparé sur ce manifeste.

| catégorie | contenu | taille |
|---|---|---|
| versionné | `article/R9_controles/` (41 fichiers), P1–P4 (4 fichiers de `src/`), `tests/test_r9_functions.py` | 5,30 Mo |
| gardé dans le répertoire de travail | npz de a/, b/, c/, d/, npy de a2/, `slurm-r9-*`, `submitted/`, `JOBID`, `r9_log.txt`, `cloture/submit_r0.log`, `__pycache__/` | 9,69 Mo |
| supprimable (TEST consigné, CLAUDE.md règle 2) | `cache/` (72 fichiers) | 50,04 Go (taille apparente ; `du` : 12 Go sur le disque) |

## Versionné (dépôt `graphene-raman`)

| chemin | taille | état git |
|---|---|---|
| `article/R9_controles/R9_rapport.md` | 76,03 ko | modifié depuis cb7241d |
| `article/R9_controles/README.md` | 1,71 ko | modifié depuis cb7241d |
| `article/R9_controles/r9_driver.py` | 124,74 ko | modifié depuis cb7241d |
| `article/R9_controles/submit_r9.sh` | 7,02 ko | modifié depuis cb7241d |
| `article/R9_controles/a2/a2c_9x9.json` | 405 o | commité (cb7241d), inchangé |
| `article/R9_controles/a2/a2d_5x5.json` | 29,78 ko | commité (cb7241d), inchangé |
| `article/R9_controles/a2/a2d_9x9.json` | 29,91 ko | commité (cb7241d), inchangé |
| `article/R9_controles/a/A1_tables.md` | 2,14 ko | commité (cb7241d), inchangé |
| `article/R9_controles/a/A3_tables.md` | 8,72 ko | commité (cb7241d), inchangé (restauré, voir le rapport) |
| `article/R9_controles/a/A3_tables_plateau.md` | 6,58 ko | non suivi |
| `article/R9_controles/a/a0_gate.json` | 446 o | commité (cb7241d), inchangé |
| `article/R9_controles/a/a1_results.json` | 9,58 ko | commité (cb7241d), inchangé |
| `article/R9_controles/a/a3_results.json` | 36,46 ko | commité (cb7241d), inchangé |
| `article/R9_controles/a/a3_results_plateau.json` | 182,21 ko | non suivi |
| `article/R9_controles/a/r0_gate.json` | 544 o | non suivi |
| `article/R9_controles/b/B_tables.md` | 7,29 ko | commité (cb7241d), inchangé |
| `article/R9_controles/b/B_tables_plateau.md` | 5,46 ko | non suivi |
| `article/R9_controles/b/b_results.json` | 21,72 ko | commité (cb7241d), inchangé |
| `article/R9_controles/b/b_results_plateau.json` | 16,03 ko | non suivi |
| `article/R9_controles/fig/folded_vs_R7.pdf` | 304,36 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/folded_vs_R7.png` | 61,79 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/kaasbjerg_fig3_map.pdf` | 851,80 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/kaasbjerg_fig3_map.png` | 383,97 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/offset_profiles.pdf` | 423,86 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/offset_profiles.png` | 71,81 ko | commité (cb7241d), inchangé |
| `article/R9_controles/fig/rcut_aligned.pdf` | 481,63 ko | non suivi |
| `article/R9_controles/fig/rcut_aligned.png` | 83,30 ko | non suivi |
| `article/R9_controles/fig/resonance_vs_nkint.pdf` | 522,41 ko | modifié depuis cb7241d |
| `article/R9_controles/fig/resonance_vs_nkint.png` | 72,74 ko | modifié depuis cb7241d |
| `article/R9_controles/fig/resonance_vs_nkint_R9.pdf` | 521,33 ko | non suivi |
| `article/R9_controles/fig/resonance_vs_nkint_R9.png` | 61,86 ko | non suivi |
| `article/R9_controles/fig/resonance_vs_nkint_plateau.pdf` | 519,19 ko | non suivi |
| `article/R9_controles/fig/resonance_vs_nkint_plateau.png` | 58,81 ko | non suivi |
| `article/R9_controles/cloture/phase0_r9_driver.diff` | 37,07 ko | non suivi |
| `article/R9_controles/cloture/phase0_submit_r9.diff` | 2,80 ko | non suivi |
| `article/R9_controles/cloture/r9_driver_propose.py` | 123,35 ko | non suivi |
| `article/R9_controles/cloture/submit_r9_propose.sh` | 7,02 ko | non suivi |
| `article/R9_controles/cloture/synthese.md` | 15,66 ko | non suivi |
| `article/R9_controles/c/C_tables.md` | 6,80 ko | commité (cb7241d), inchangé |
| `article/R9_controles/c/c_results.json` | 75,92 ko | commité (cb7241d), inchangé |
| `article/R9_controles/d/D_tables.md` | 1,67 ko | commité (cb7241d), inchangé |
| `article/R9_controles/d/d_results.json` | 3,66 ko | commité (cb7241d), inchangé |
| `article/R9_controles/manifeste_R9.md` | ce fichier | non suivi |
| `src/electron_defect_interaction/defects/alignment.py` | 10,77 ko | commité (cb7241d), inchangé |
| `src/electron_defect_interaction/wannier/wannier_interpolation.py` | 12,24 ko | commité (cb7241d), inchangé |
| `src/electron_defect_interaction/defects/many_body/local_tmatrix.py` | 15,02 ko | commité (cb7241d), inchangé |
| `src/electron_defect_interaction/wannier/supercell_fold.py` | 8,35 ko | commité (cb7241d), inchangé |
| `tests/test_r9_functions.py` | 4,40 ko | commité (cb7241d), inchangé |

## Gardé dans le répertoire de travail (`graphene/qe/defects/R9_controles/`)

| chemin | taille |
|---|---|
| `a/a1_profiles_12x12.npz` | 33,28 ko |
| `a/a1_profiles_5x5.npz` | 7,58 ko |
| `a/a1_profiles_6x6.npz` | 9,95 ko |
| `a/a1_profiles_7x7.npz` | 12,76 ko |
| `a/a1_profiles_8x8.npz` | 16,00 ko |
| `a/a1_profiles_9x9.npz` | 19,67 ko |
| `a2/ML_const_9x9_coarse_4b.npy` | 1,68 Mo |
| `b/b_curves_12x12.npz` | 488,07 ko |
| `b/b_curves_12x12_plateau.npz` | 254,10 ko |
| `b/b_curves_9x9.npz` | 722,18 ko |
| `b/b_curves_9x9_plateau.npz` | 488,34 ko |
| `c/c1_ldos_N12_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N12_brut.npz` | 9,47 ko |
| `c/c1_ldos_N15_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N15_brut.npz` | 9,47 ko |
| `c/c1_ldos_N18_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N18_brut.npz` | 9,47 ko |
| `c/c1_ldos_N21_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N21_brut.npz` | 9,47 ko |
| `c/c1_ldos_N24_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N24_brut.npz` | 9,47 ko |
| `c/c1_ldos_N27_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N27_brut.npz` | 9,47 ko |
| `c/c1_ldos_N6_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N6_brut.npz` | 9,47 ko |
| `c/c1_ldos_N9_aligne.npz` | 9,47 ko |
| `c/c1_ldos_N9_brut.npz` | 9,47 ko |
| `c/c2_ldos_N27_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N27_brut.npz` | 9,47 ko |
| `c/c2_ldos_N36_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N36_brut.npz` | 9,47 ko |
| `c/c2_ldos_N45_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N45_brut.npz` | 9,47 ko |
| `c/c2_ldos_N54_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N54_brut.npz` | 9,47 ko |
| `c/c2_ldos_N63_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N63_brut.npz` | 9,47 ko |
| `c/c2_ldos_N72_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N72_brut.npz` | 9,47 ko |
| `c/c2_ldos_N81_aligne.npz` | 9,47 ko |
| `c/c2_ldos_N81_brut.npz` | 9,47 ko |
| `c/c3_ldos_nk600_aligne.npz` | 136,09 ko |
| `c/c3_ldos_nk600_brut.npz` | 136,09 ko |
| `c/c3_ldos_nk900_aligne.npz` | 136,09 ko |
| `c/c3_ldos_nk900_brut.npz` | 136,09 ko |
| `d/d1_maps.npz` | 2,77 Mo |
| `slurm-r9-r9a-21955512.err` | 752 o |
| `slurm-r9-r9a-21955512.out` | 1,31 ko |
| `slurm-r9-r9a-21955644.err` | 33 o |
| `slurm-r9-r9a-21955644.out` | 3,71 ko |
| `slurm-r9-r9a3-21955645.err` | 2,54 ko |
| `slurm-r9-r9a3-21955645.out` | 14,50 ko |
| `slurm-r9-r9a3-21976370.err` | 451 o |
| `slurm-r9-r9a3-21976370.out` | 11,42 ko |
| `slurm-r9-r9b-21954789.err` | 33 o |
| `slurm-r9-r9b-21954789.out` | 5,55 ko |
| `slurm-r9-r9b-21955646.err` | 33 o |
| `slurm-r9-r9b-21955646.out` | 7,83 ko |
| `slurm-r9-r9b-21976371.err` | 2,28 ko |
| `slurm-r9-r9b-21976371.out` | 7,79 ko |
| `slurm-r9-r9box5x5-21954788.err` | 33 o |
| `slurm-r9-r9box5x5-21954788.out` | 5,24 ko |
| `slurm-r9-r9box9x9-21954787.err` | 33 o |
| `slurm-r9-r9box9x9-21954787.out` | 5,85 ko |
| `slurm-r9-r9c-21955647.err` | 33 o |
| `slurm-r9-r9c-21955647.out` | 8,43 ko |
| `slurm-r9-r9d-21955648.err` | 33 o |
| `slurm-r9-r9d-21955648.out` | 1,76 ko |
| `slurm-r9-r9r0-21976369.err` | 33 o |
| `slurm-r9-r9r0-21976369.out` | 1,96 ko |
| `JOBID` | 1,05 ko |
| `r9_log.txt` | 67,91 ko |
| `cloture/submit_r0.log` | 50,05 ko |
| `submitted/` (copies du pilote, du lanceur et diff par soumission, 16 soumissions + `last/`) | 1,88 Mo |
| `__pycache__/` | 278,47 ko |

Les fichiers du répertoire de travail qui ont une copie dans `article/R9_controles/` (rapport, README, manifeste, pilote, lanceur, json, tables, figures) ne sont pas répétés.

## Supprimable : `cache/` (proposition du 2026-09-28 ; supprimé le 2026-09-29, voir la dernière section)

Tous ces fichiers se reconstruisent à l'identique depuis `results/M2/`, les wannierisations et les `.save` (lecture seule) par le pilote de R9 ; les g₀ et M_W se
recréent d'eux-mêmes à la première sous-commande qui en a besoin.

| fichier | taille | rôle | reconstruction (commande, durée mesurée) |
|---|---|---|---|
| `cache/Fw_5x5.npz` | 156,29 Mo | F_W = V† M^L[1_boîte] V, 5x5 (variantes exactes) | `r9_driver.py a2dpost --size 5x5` (≈ 1 min ; exige ML_box_5x5) |
| `cache/Fw_9x9.npz` | 212,62 Mo | F_W = V† M^L[1_boîte] V, 9x9 (variantes exactes) | `r9_driver.py a2dpost --size 9x9` (≈ 1 min ; exige ML_box_9x9) |
| `cache/ML_box_5x5.json` | 839 o | M^L[1_boîte] 5x5 dense (A.2 exact ; sans dimension) — métadonnées | `bash submit_r9.sh box 5x5` (a2d, nœud exclusif 32 × 6) : 13,4 min (job 14 min) |
| `cache/ML_box_5x5.npy` | 2,50 Go | M^L[1_boîte] 5x5 dense (A.2 exact ; sans dimension) | `bash submit_r9.sh box 5x5` (a2d, nœud exclusif 32 × 6) : 13,4 min (job 14 min) |
| `cache/ML_box_9x9.json` | 840 o | M^L[1_boîte] 9x9 dense (A.2 exact ; sans dimension) — métadonnées | `bash submit_r9.sh box 9x9` (a2d, nœud exclusif 32 × 6) : 53,1 min (job 54 min) |
| `cache/ML_box_9x9.npy` | 3,40 Go | M^L[1_boîte] 9x9 dense (A.2 exact ; sans dimension) | `bash submit_r9.sh box 9x9` (a2d, nœud exclusif 32 × 6) : 53,1 min (job 54 min) |
| `cache/Mwr_12x12.npz` | 132,74 Mo | M_W du M2 dense 12x12 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 5 s |
| `cache/Mwr_5x5.npz` | 156,28 Mo | M_W du M2 dense 5x5 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 5 s |
| `cache/Mwr_6x6.npz` | 132,74 Mo | M_W du M2 dense 6x6 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 5 s |
| `cache/Mwr_7x7.npz` | 245,90 Mo | M_W du M2 dense 7x7 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 10 s |
| `cache/Mwr_8x8.npz` | 419,48 Mo | M_W du M2 dense 8x8 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 19 s |
| `cache/Mwr_9x9.npz` | 212,61 Mo | M_W du M2 dense 9x9 (chaîne de production) | automatique (`mwr_cached`, toute sous-commande) : rotation 8 s |
| `cache/g0_12x12_nk300.json` | 986 o | g₀ sur l'amas, 12x12, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 1,0 min |
| `cache/g0_12x12_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 12x12, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 1,0 min |
| `cache/g0_12x12_nk300_rc3_sig.json` | 985 o | g₀ sur l'amas, 12x12, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 3,4 min |
| `cache/g0_12x12_nk300_rc3_sig.npy` | 812,41 Mo | g₀ sur l'amas, 12x12, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 3,4 min |
| `cache/g0_12x12_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 12x12, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_12x12_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 12x12, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_12x12_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 12x12, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_12x12_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 12x12, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_12x12_nk450.json` | 985 o | g₀ sur l'amas, 12x12, 450², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 2,2 min |
| `cache/g0_12x12_nk450.npy` | 813,08 Mo | g₀ sur l'amas, 12x12, 450², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 2,2 min |
| `cache/g0_12x12_nk600.json` | 986 o | g₀ sur l'amas, 12x12, 600², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 4,1 min |
| `cache/g0_12x12_nk600.npy` | 813,08 Mo | g₀ sur l'amas, 12x12, 600², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 4,1 min |
| `cache/g0_12x12_nk900.json` | 985 o | g₀ sur l'amas, 12x12, 900², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 8,2 min |
| `cache/g0_12x12_nk900.npy` | 813,08 Mo | g₀ sur l'amas, 12x12, 900², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 8,2 min |
| `cache/g0_5x5_nk300.json` | 984 o | g₀ sur l'amas, 5x5, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 0,8 min |
| `cache/g0_5x5_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 5x5, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande a3) : 0,8 min |
| `cache/g0_5x5_nk300_rc3_sig.json` | 983 o | g₀ sur l'amas, 5x5, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 0,8 min |
| `cache/g0_5x5_nk300_rc3_sig.npy` | 812,41 Mo | g₀ sur l'amas, 5x5, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 0,8 min |
| `cache/g0_5x5_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 5x5, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 6,7 min |
| `cache/g0_5x5_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 5x5, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 6,7 min |
| `cache/g0_5x5_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 5x5, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 6,8 min |
| `cache/g0_5x5_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 5x5, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 6,8 min |
| `cache/g0_6x6_nk300.json` | 984 o | g₀ sur l'amas, 6x6, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 2,0 min |
| `cache/g0_6x6_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 6x6, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande a3) : 2,0 min |
| `cache/g0_6x6_nk300_rc3_sig.json` | 983 o | g₀ sur l'amas, 6x6, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 3,8 min |
| `cache/g0_6x6_nk300_rc3_sig.npy` | 812,41 Mo | g₀ sur l'amas, 6x6, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 3,8 min |
| `cache/g0_6x6_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 6x6, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 2,5 min |
| `cache/g0_6x6_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 6x6, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 2,5 min |
| `cache/g0_6x6_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 6x6, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,4 min |
| `cache/g0_6x6_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 6x6, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 5,4 min |
| `cache/g0_7x7_nk300.json` | 983 o | g₀ sur l'amas, 7x7, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 3,3 min |
| `cache/g0_7x7_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 7x7, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande a3) : 3,3 min |
| `cache/g0_7x7_nk300_rc3_sig.json` | 982 o | g₀ sur l'amas, 7x7, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 2,2 min |
| `cache/g0_7x7_nk300_rc3_sig.npy` | 812,74 Mo | g₀ sur l'amas, 7x7, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 2,2 min |
| `cache/g0_7x7_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 7x7, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_7x7_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 7x7, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_7x7_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 7x7, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,1 min |
| `cache/g0_7x7_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 7x7, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 5,1 min |
| `cache/g0_8x8_nk300.json` | 983 o | g₀ sur l'amas, 8x8, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 3,3 min |
| `cache/g0_8x8_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 8x8, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande a3) : 3,3 min |
| `cache/g0_8x8_nk300_rc3_sig.json` | 982 o | g₀ sur l'amas, 8x8, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 2,9 min |
| `cache/g0_8x8_nk300_rc3_sig.npy` | 811,73 Mo | g₀ sur l'amas, 8x8, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 2,9 min |
| `cache/g0_8x8_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 8x8, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 4,0 min |
| `cache/g0_8x8_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 8x8, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 4,0 min |
| `cache/g0_8x8_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 8x8, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 1,9 min |
| `cache/g0_8x8_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 8x8, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 1,9 min |
| `cache/g0_9x9_nk300.json` | 982 o | g₀ sur l'amas, 9x9, 300², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 1,0 min |
| `cache/g0_9x9_nk300.npy` | 813,08 Mo | g₀ sur l'amas, 9x9, 300², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 1,0 min |
| `cache/g0_9x9_nk300_rc3_sig.json` | 982 o | g₀ sur l'amas, 9x9, 300², R_cut 3, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a0 / a3) : 1,1 min |
| `cache/g0_9x9_nk300_rc3_sig.npy` | 812,41 Mo | g₀ sur l'amas, 9x9, 300², R_cut 3, grille « sig » | automatique (`g0_cached`, sous-commande a0 / a3) : 1,1 min |
| `cache/g0_9x9_nk300_rc4_res.json` | 1,52 ko | g₀ sur l'amas, 9x9, 300², R_cut 4, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_9x9_nk300_rc4_res.npy` | 2,32 Go | g₀ sur l'amas, 9x9, 300², R_cut 4, grille « res » | automatique (`g0_cached`, sous-commande a3) : 5,2 min |
| `cache/g0_9x9_nk300_rc4_sig.json` | 1,52 ko | g₀ sur l'amas, 9x9, 300², R_cut 4, grille « sig » — métadonnées | automatique (`g0_cached`, sous-commande a3) : 4,1 min |
| `cache/g0_9x9_nk300_rc4_sig.npy` | 2,32 Go | g₀ sur l'amas, 9x9, 300², R_cut 4, grille « sig » | automatique (`g0_cached`, sous-commande a3) : 4,1 min |
| `cache/g0_9x9_nk450.json` | 982 o | g₀ sur l'amas, 9x9, 450², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 2,3 min |
| `cache/g0_9x9_nk450.npy` | 813,08 Mo | g₀ sur l'amas, 9x9, 450², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 2,3 min |
| `cache/g0_9x9_nk600.json` | 983 o | g₀ sur l'amas, 9x9, 600², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 3,9 min |
| `cache/g0_9x9_nk600.npy` | 813,08 Mo | g₀ sur l'amas, 9x9, 600², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 3,9 min |
| `cache/g0_9x9_nk900.json` | 983 o | g₀ sur l'amas, 9x9, 900², R_cut 3, grille « res » — métadonnées | automatique (`g0_cached`, sous-commande b) : 8,5 min |
| `cache/g0_9x9_nk900.npy` | 813,08 Mo | g₀ sur l'amas, 9x9, 900², R_cut 3, grille « res » | automatique (`g0_cached`, sous-commande b) : 8,5 min |

## git status (lecture, 2026-09-28 18:00) — branche main à jour avec origin/main (dernier commit cb7241d « R9 checkpoint »)

```
 M article/R9_controles/R9_rapport.md
 M article/R9_controles/README.md
 M article/R9_controles/fig/resonance_vs_nkint.pdf
 M article/R9_controles/fig/resonance_vs_nkint.png
 M article/R9_controles/r9_driver.py
 M article/R9_controles/submit_r9.sh
?? article/R9_controles/a/A3_tables_plateau.md
?? article/R9_controles/a/a3_results_plateau.json
?? article/R9_controles/a/r0_gate.json
?? article/R9_controles/b/B_tables_plateau.md
?? article/R9_controles/b/b_results_plateau.json
?? article/R9_controles/cloture/
?? article/R9_controles/fig/rcut_aligned.pdf
?? article/R9_controles/fig/rcut_aligned.png
?? article/R9_controles/fig/resonance_vs_nkint_R9.pdf
?? article/R9_controles/fig/resonance_vs_nkint_R9.png
?? article/R9_controles/fig/resonance_vs_nkint_plateau.pdf
?? article/R9_controles/fig/resonance_vs_nkint_plateau.png
?? article/R9_controles/manifeste_R9.md
```

## À commiter par Greg

- `article/R9_controles/` : fichiers modifiés depuis cb7241d (`R9_rapport.md`, `README.md`, `r9_driver.py`, `submit_r9.sh`, `fig/resonance_vs_nkint.{pdf,png}`) et nouveaux
  (`manifeste_R9.md`, `a/a3_results_plateau.json`, `a/A3_tables_plateau.md`, `a/r0_gate.json`, `b/b_results_plateau.json`, `b/B_tables_plateau.md`, `cloture/` (5 fichiers), `fig/rcut_aligned.{pdf,png}`,
  `fig/resonance_vs_nkint_R9.{pdf,png}`, `fig/resonance_vs_nkint_plateau.{pdf,png}`).
- `src/` (P1–P4) et `tests/test_r9_functions.py` : déjà commités dans cb7241d, inchangés depuis ; rien à commiter.
- Les autres lignes de `git status` (s'il y en a hors `article/R9_controles/`) ne viennent pas de R9.

## Clôture définitive (2026-09-29)

- **`cache/` supprimé** le 2026-09-29 sur GO de Greg. Avant suppression, le contenu (72 fichiers) était identique à la liste de la section
  « Supprimable » ; aucun autre script ne le lisait, aucun job R9 n'était en file. Place libérée : 50,04 Go en taille apparente, 12 Go sur le
  disque. Les fichiers se reconstruisent par les commandes de cette liste. Le répertoire de travail fait maintenant 14 Mo.
- **Audit « image minimale »** (2026-09-28, après ce manifeste) : `audit_image_minimale.md`, `audit/audit_results.json` et la sous-commande
  `audit` de `r9_driver.py`. Les trois sont versionnés dans `article/R9_controles/`. Le json fait 196,6 ko. Il n'y a pas de cache, rien de
  supprimable.
- **État git** : tout `article/R9_controles/` est commité, y compris l'audit, dans cff65cd « R9 done » (numéros après la réécriture de
  l'historique du 2026-09-29 ; cb7241d = ancien 5a4bc94, 1883861 = ancien ab2d884). Cette mise à jour du README et du manifeste reste à
  commiter.
- **Corrections proposées par l'audit** (§7 d'`audit_image_minimale.md`) : gardées par Greg pour plus tard, non appliquées.
- **Ligne « Lecture » du README** : laissée vide (texte à fournir par Greg).
- Aucun fichier n'est promu en production.
