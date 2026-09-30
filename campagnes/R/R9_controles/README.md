# R9 — Contrôles avant le chapitre 4 (alignement de ΔV, résonance contre N_k^int, chaîne repliée contre R7, grandeur de Kaasbjerg)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 3 figures hors mémoire (`offset_profiles`, `kaasbjerg_fig3_map`, `resonance_vs_nkint_R9`), `submit_r9.sh`, `cloture/r9_driver_propose.py` et `cloture/submit_r9_propose.sh` (proposition remplacée). `r9_driver.py` est gardé : `synth` et `cfig` redessinent les quatre figures du mémoire depuis les json. Gardés : rapports, tables et json de `a/ a2/ b/ c/ d/ audit/`, `cloture/synthese.md` et diffs, `fig/` (4 figures). Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R9_controles/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R9 — Contrôles avant le chapitre 4 (alignement de ΔV, résonance, chaîne repliée, grandeur de Kaasbjerg)
- **Question** : (A) quel décalage uniforme C_N de ΔV^L retirer de M2, et que deviennent les chiffres alignés ? (B) la résonance bouge-t-elle quand N_k^int converge ? (C) la chaîne repliée N×N redonne-t-elle R7 et la limite diluée ? (D) notre grandeur se compare-t-elle aux ~70 eV Å² de Kaasbjerg ? S'y ajoutent la clôture (rejeu avec C_N = moyenne du plateau) et l'audit de l'image minimale.
- **Pourquoi** :
  - M2 est produit sans aucune soustraction. Or un ΔV constant C ajoute exactement N_cells·C·𝕀 à M : N_cells·|C_N| vaut 1,47 à 4,21 eV, contre max|M| de 23 à 27 eV.
  - Les E_res de R6 variaient avec N_k^int (−0,175 ; −0,202 ; −0,238 eV).
  - Deux valeurs « Kaasbjerg » déjà citées (107,22 et 76 eV Å²) ont des définitions différentes.
  - Lu reposait sur un seul atome.
- **Méthode** :
  - A.1 : C_N pour 5, 6, 7, 8, 9 et 12 = moyenne des décalages de sphère de 1,0 Å sur les atomes à ≥ 0,75 r_max. A.2 : identité M^L[1_boîte] (9×9, 5×5). A.3 : variantes alignées sans recalcul de M.
  - B : pic de la courbe Γ_T, pic de −Im T̄(K) et E_res à N_k^int 300 / 450 / 600 / 900.
  - C : V_loc 9×9 replié en N×N (6…27, puis 36…81 en bloc π), comparé à R7.
  - D : Ṽ = A_cell·|M| près de K.
  - Clôture : rejeu de A.3 et B avec C_N plateau. Audit géométrique sans job.
  - Jobs : 21954787, 21954788, 21954789, 21955512 (échec), 21955644 à 21955648 ; clôture 21976369, 21976370, 21976371.
- **Résultat** :
  - A : aucune taille n'a max|écart| ≤ 5 meV sur son plateau (de 24,10 meV en 9×9 à 233,70 meV en 5×5), donc R9 aligne avec Lu. C_N plateau : −57,56 / −50,51 / −26,91 / −14,41 / −25,14 / −18,69 meV. L'écart plateau − Lu va de +64,77 à −0,49 meV. M^L[1] 9×9 = 81·𝕀 à 1,56e-13.
  - Clôture (plateau, R_cut 3) : tab:rcut_M 9×9 6,679e-2 → 2,470e-2 ; 12×12 9,021e-2 → 3,097e-2. Médiane Γ·N_cells 9×9 : 3 132,60 → 3 189,01 meV. Dispersion des familles : 3m 0,96 % → 3,60 %, non-3m 8,60 % → 1,83 %.
  - B (9×9) : pic de la courbe Γ_T −0,1800 / −0,1825 / −0,1825 / −0,1825 eV (tel quel). Tel quel, E_res saute de −0,1748 (x = 7) à −0,2018 eV (x = 9) dès N_k^int 450 ; aligné, elle reste à −0,1748 aux quatre N_k^int.
  - C : π quasi-lié, modèle / QE : −0,6774 / −0,7373 eV (N = 9) et −0,2812 / −0,2785 eV (N = 27) ; ε∞ en 1/N −0,0486 / −0,0523 eV. Maximum de LDOS dans la limite diluée : −0,1925 eV (tel quel), −0,1775 eV (aligné).
  - D : ½ Tr du bloc (K, K) 75,87 eV Å² (tel quel), 86,38 (aligné exact) ; disques autour de K 76,96 (π) / 78,16 (π*) eV Å². La valeur 107,22 est la norme d'une ligne.
  - Audit : Lu était choisi par une image minimale axe par axe (9×9 : 18,51 Å affiché, 11,12 Å vrai). Lu est donc faux de +298,47 (5×5), +67,75 (8×8) et −21,23 meV (9×9). Corrections proposées, non appliquées.
- **Interprétation / décision** :
  - C_N = moyenne du plateau (i) est retenu à la clôture, choisi avant tout résultat aligné ; R8 et R10 le reprennent.
  - Observables de la résonance : pic de la courbe Γ_T et pic de −Im T̄(K). E_res est rapportée avec sa couronne.
  - C et D ne sont pas rejoués au plateau (écart estimé ≤ 0,33 meV et +0,09 eV Å²). Les corrections de l'audit sont reprises par R10.
  - Statut TEST, rien n'est promu ; R9 clos le 2026-09-29.
- **Où sont les données** (relatif à `campagnes/R/R9_controles/`) :
  - Rapports : `R9_rapport.md`, `cloture/synthese.md`, `audit_image_minimale.md`.
  - A : `a/A1_tables.md`, `a/a1_results.json` (`C_i_eV`, `C_retenu_eV`), `a/A3_tables.md`, `a/a3_results.json`, `a/A3_tables_plateau.md`, `a/a3_results_plateau.json`, `a/a0_gate.json`, `a/r0_gate.json`, `a2/a2c_9x9.json`, `a2/a2d_{5x5,9x9}.json`.
  - B : `b/B_tables.md`, `b/b_results.json`, `b/B_tables_plateau.md`, `b/b_results_plateau.json`.
  - C et D : `c/C_tables.md`, `c/c_results.json`, `d/D_tables.md`, `d/d_results.json`.
  - Audit : `audit/audit_results.json`.
  - Figures du mémoire : `fig/` (`resonance_vs_nkint`, `resonance_vs_nkint_plateau`, `rcut_aligned` par `r9_driver.py synth` ; `folded_vs_R7` par `r9_driver.py cfig` ; depuis les json, sans calcul). Aucun npz dans la copie (`cache/` supprimé, voir `manifeste_R9.md`).

---

- But : (A) décalage uniforme de ΔV^L dans M2 et variantes alignées ; (B) position de la résonance quand N_k^int converge ; (C) chaîne repliée N×N contre R7 et limite diluée ; (D) grandeur comparée aux ~70 eV Å² de Kaasbjerg. Clôture : rejeu avec C_N = moyenne du plateau (i), test intensif hors k = k′. Audit : image minimale dans la cellule à 60° (écarts 1–3).
- Prompts d'origine : « R9 — Contrôles avant le chapitre 4 » (Greg, 2026-09-27) ; « R9 — Clôture » et « Audit — image minimale » (Greg, 2026-09-28).
- Statut : **TEST** (post-traitement seul, aucun calcul QE ; résultats consignés dans `R9_rapport.md` ; rien promu en production).
- Dates : phase 0 2026-09-27 ; étape G et GO 2026-09-28 ; A, B, C, D terminés 2026-09-28 ; clôture (phase 0, GO, rejeu) et audit 2026-09-28 ; `cache/` supprimé 2026-09-29 (GO de Greg) ; R9 clos 2026-09-29.
- Jobs (`JOBID`, diffs dans `submitted/<jobid>/`) : 21954787, 21954788 (M^L[1_boîte] 9×9, 5×5), 21954789 (B), 21955512 (A, échec d'exécution), 21955644 (A), 21955645 (A.3), 21955646 (B aligné), 21955647 (C), 21955648 (D) ; clôture 21976369 (R.0), 21976370 (A.3 plateau), 21976371 (B plateau).
- Pilote et lanceur : `r9_driver.py`, `submit_r9.sh` ; fonctions P1–P4 dans `src/` et `tests/test_r9_functions.py` (dépôt, commités dans cb7241d) ; sous-commande `audit` (nœud de connexion, sans job).
- Sorties principales : `R9_rapport.md` ; `a/` (A.1, A.3, porte R.0), `a2/`, `b/`, `c/`, `d/` (json, tables md) ; `cloture/synthese.md` ; `audit_image_minimale.md`, `audit/audit_results.json` ; `fig/` (offset_profiles, resonance_vs_nkint, rcut_aligned, folded_vs_R7, kaasbjerg_fig3_map) ; `manifeste_R9.md`.
- Copie versionnée : `campagnes/R/R9_controles/` (rapport, README, manifeste, pilote, lanceur, tables, json, figures ; ni npz, ni slurm ; `cache/` supprimé, reconstructible par le manifeste).
- Lecture :
