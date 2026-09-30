# R10 — Base cohérente pour les chiffres du chapitre 4 (alignement par la moyenne du plateau, argument E_res, familles, rejeu de production aligné)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 12 figures hors mémoire (`fig_locality`, `fig_M_map`, `fig_M_scaling`, `fig_spectral`, `fig_Ved_{map,radial,radial_masked,boundary,profile_mean}`, `dV_z_profiles`, `eres_vs_grid`, `M_map_brut_aligne`) et les 28 planches de `fig/avant_apres/`. `r10_driver.py` et `submit_r10.sh` sont gardés (base vivante du ch. 4). Gardés : rapport, tables et json/npz de `a/ b/ c/`, `fig/` (22 figures du mémoire), `manifeste_R10.md`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R10_plateau/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R10 — Base unique du chapitre 4 : alignement par la moyenne du plateau (Kumagai–Oba)
- **Question** : mettre tous les chiffres du ch. 4 sur une seule base :
  - (A) ΔV aligné par la moyenne du plateau (i) pour 13 tailles ;
  - (B) E_res et familles ;
  - (C) rejeu de la production avec M_W(R, R) − C_N·𝕀₅ et les corrections de l'audit, dans `results/M2_plateau`.
- **Pourquoi** : c'est le but du README (une base cohérente pour le ch. 4). R9 avait laissé des chiffres tel quel, avec Lu, et au plateau pour six tailles seulement. L'audit a invalidé Lu, et le critère « sans plateau » est abandonné (D14).
- **Méthode** :
  - A : moyenne des décalages de sphère de 1,0 Å sur les atomes à distance vraie ≥ 0,75 r_max (5…12, 15, 18, 21, 24, 27) ; niveaux QE réalignés.
  - B.1 : E_res en fonction de la grille de sortie, de 120² à 960² (9×9, N_k^int 900). B.2 : familles 5…12.
  - C : portes (C.0), production (C.1), tab:rcut_M et Kaasbjerg avec Wigner-Seitz (C.2), contrôles de l'audit (C.3), avant/après (C.4), table de correspondance (C.5), ⟨ΔV⟩_3D (C.6).
  - Jobs : 22041340 `r10a`, 22041341 `r10b` ; 22058838 `r10c0`, 22058839–54 (C.1), 22058855 `r10c1f`, 22058856 `r10c2`, 22058857 `r10c3`, 22058858 `r10c6`.
- **Résultat** :
  - C_N plateau (i) : −57,56 (5×5), −25,14 (9×9), −18,69 (12×12), −13,01 meV (27×27). Le rms du plateau passe de 94,14 à 0,53 meV ; C_N − Lu va de +64,77 (5×5) à −0,24 meV (27×27).
  - E_res 9×9 (120² / 240² / 480² / 960²) : tel quel −0,22689 / −0,20184 / −0,19993 / −0,19050 eV ; plateau −0,13441 / −0,17485 / −0,17485 / −0,17210 eV. Le pic de la courbe Γ_T au plateau reste à −0,1800 eV de 240² à 960².
  - Niveau 1 (240², R_cut 3), médiane Γ·N_cells tel quel → plateau : 3 288,90 → 2 997,14 (5×5), 3 132,60 → 3 189,01 (9×9), 3 162,78 → 3 264,51 meV (12×12). E_res 6×6 passe de −0,227 à −0,175 eV ; en 9×9 et 12×12 elle reste à −0,175 eV.
  - tab:rcut_M (R_cut 3) : 9×9 6,6789e-2 → 2,4704e-2 (Wigner-Seitz 2,5879e-2) ; 12×12 9,0206e-2 → 3,0966e-2. Kaasbjerg : valence en K 76,956 → 81,605 eV Å², ½ Tr (K, K) 75,87 → 86,59 eV Å².
  - tab:tests_M par famille : 3m 2,257e-2, non-3m 2,057e-2. C14 (±9,05 meV) : écart relatif max de Γ_T 6,36e-2. 12 figures changent au pixel, 16 sont identiques.
  - C.6 : les 26 SCF ont convergé. N²⟨ΔV⟩_3D vaut 713,1 à 749,7 meV pour 5…12, mais −36 137,7 meV (24×24) et −72 836,6 meV (27×27).
- **Interprétation / décision** :
  - C_N = moyenne du plateau pour toutes les tailles, sans critère d'arrêt (D14).
  - `results/M2_plateau` devient la production du ch. 4 ; les M2 brutes restent dans `results/M2`.
  - La porte des anneaux a été refusée au seuil de 1e-13 (1,21e-13 et 1,01e-13 eV), puis acceptée au seuil de 1e-10 fixé par Greg, sans recalcul.
  - 12 figures installées dans `figures/` ; commit f2dae43.
  - Interprétation de l'anomalie C.6 : non consignée. R10 clos le 2026-09-30.
- **Où sont les données** (relatif à `campagnes/R/R10_plateau/`) :
  - Rapport : `R10_rapport.md`.
  - A : `a/A1_tables.md`, `a/a1_results.json` (`C_N_eV`), `a/A2_tables.md`, `a/a2_results.json`, `a/A3_tables.md`, `a/a3_results.json`, `a/profiles_<N>x<N>.npz`.
  - B : `b/B_tables.md`, `b/b_results.json`, `b/b1_curves.npz`.
  - C : `c/C0_tables.md`, `c/C2_tables.md`, `c/C3_tables.md`, `c/C6_tables.md`, `c/table_v2_plateau.md`, `c/c0_results.json`, `c/c1post_results.json`, `c/c2_results.json`, `c/c3_results.json`, `c/c4_results.json`, `c/c6_results.json`, `c/c2_kaasbjerg_maps.npz`, `c/c6_profile_<N>x<N>.npz`.
  - Figures du mémoire : `fig/` (22 figures : les 19 de production, identiques à `figures/`, plus `offset_profiles_13`, `kaasbjerg_plateau_ws`, `levels_vs_invN`).
  - Production : `results/M2_plateau/` (33 fichiers, `MD5SUMS_2026-09-30.txt`).

---

- But : mettre tous les chiffres du ch. 4 sur une seule base — ΔV aligné par la moyenne du plateau (i) de R9 pour toutes les tailles (A), argument E_res et familles au niveau 1 (B), rejeu de la production avec M_W(R,R) − C_N·𝕀₅ (approximation (i)) et les corrections de l'audit de l'image minimale (C).
- Prompt d'origine : « R10 — Base cohérente pour les chiffres du chapitre 4 » (Greg, 2026-09-29).
- Statut : **TEST** (post-traitement seul, aucun calcul QE, aucun `.save`) ; sorties de C = production du ch. 4 dans `results/M2_plateau` (produits seulement ; matrices M2 brutes dans `results/M2`).
- Dates : phase 0, étape G et GO 1 (A, B) 2026-09-29 ; (a) `c61c118`, (b) `ea91d61`, (c) `23ee3bb` ; GO 2 (C) soumis le 2026-09-29, terminé le 2026-09-30 ; figures installées et **campagne close** le 2026-09-30.
- Jobs (`JOBID`, diffs dans `submitted/<jobid>/`) : 22041340 `r10a` (A), 22041341 `r10b` (B) ; 22058838 `r10c0` (C.0), 22058839–54 (C.1, lanceurs de `scripts/`), 22058855 `r10c1f`, 22058856 `r10c2`, 22058857 `r10c3` (C.3–C.5 ; porte des anneaux refusée à 1e-13, acceptée au seuil de 1e-10 fixé par Greg, sans recalcul), 22058858 `r10c6` (C.6).
- Pilote et lanceur : `r10_driver.py`, `submit_r10.sh` ; fonctions de 0.5 (a) écrites par Greg (D1) ; changements (b) de `scripts/` écrits par Code (section « (b) » du rapport).
- Sorties : `R10_rapport.md` (phase 0, étape G, A et B, (b), (c), GO 2, rapport de C) ; `a/` (C_N des 13 tailles : `a1_results.json`), `b/`, `c/` (portes, tab:rcut_M à trois colonnes, Kaasbjerg, contrôles de l'audit, `table_v2_plateau.md`, C.6), `fig/` (figures de A, B, C, figures de production régénérées, `avant_apres/`) ; `manifeste_R10.md` (cache/).
- Copie versionnée : `campagnes/R/R10_plateau/` (rapport, README, pilote, lanceur, tables, json, npz < 5 Mo, figures ; jamais slurm ni `cache/`).
- Lecture :
