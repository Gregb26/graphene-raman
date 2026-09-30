# R6 — Production corrigée : normalisation de la partie locale de M (chapitre 4)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : toutes les figures (`etape2/fig/`, `etape3/figures_v2/`, `etape3/figures_fix/` : versions R6 des figures du ch. 4, remplacées par R10 dans `figures/`), les pilotes et lanceurs (`etape1/r6_kernel_check.py`, `r6_states_identity.py`, `etape2/r6_d4.py`, `etape3/r6_compare_v1_v2.py`, `r6_level1_gate.py`, `r6_m_rcut_resigma.py`, `r6_tests_gate_row.py`, `refactor_3_0.py`, `runbook_3.sh`, `install_figures_R6.sh`, `submit_r6*.sh`). L'outil de reconstruction des M2, `scripts/assemble_M2.py`, reste dans `scripts/`. Gardés : rapport, tables, json/jsonl, `csv_v2/`, `MD5SUMS`, diffs de la phase 0 et `diff_3.0.patch`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R6_production_corrigee/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R6 — Correction de la normalisation de M^L (M2) et régénération du chapitre 4
- **Question** : dans tous les M de production, M^L était calculé avec des états normés sur la super-cellule (ψ = u e^{ik·r}/√Ω_sc) et M^NL avec des états normés sur la maille (4π/√Ω_uc). M^L est-il donc N_cells fois trop petit ? Il s'agissait de le corriger (M2 = N_cells·M^L + M^NL), de vérifier chaque M2 (porte A.2, ≤ 1e-6 eV), de revalider la chaîne contre la super-cellule 9×9 (escalier D4), puis de régénérer tout le chapitre 4.
- **Pourquoi** : R5 (A.2) avait mesuré un rapport de 81,0000 entre l'application directe de ΔV et M^L sur la 9×9. La reconstruction KS ne pouvait pas voir ce facteur : elle est calculée sur la maille, et pour un V_p périodique le 1/Ω_sc est exactement compensé par les N_cells copies. Toute la chaîne du chapitre 4, ainsi que Γ^ed du chapitre 5, lisait les M v1.
- **Méthode** : phase 0 (inventaire, diffs proposés) → étape G (Greg corrige `local_R.py`, commit 99d64da ; les deux écarts des autres noyaux sont clos par a223687) → étape 1 : M^L grossiers 5×5, 6×6 et 8×8 recalculés, 15 M2 réassemblés dans `results/M2/`, porte A.2 par fichier (J1b 21820489, J2b 21820490, J3b 21820491, J4b 21820492) → étape 2 : escalier D4, (a1) à n bandes et critère de pôle sur la 9×9 (J5 21833649, J6 21833650, J7 21850495) → étape 3 : lecteurs de M passés au `results_dir` avec porte v2, test d'or 5×5 (21852238), niveaux 1 et 2, résonance 6/9/12, C14, tests, figures (21857272–21857284, 21862371–21862375, 21872955). Aucun calcul QE.
- **Résultat** :
  - Noyau corrigé : M^L 9×9 v2 contre 81 × v1 à 8,97e-16 relatif (7×7 : 8,72e-16). Les 15 M2 sont exacts au bit. Porte A.2 OK pour 8 grossiers et 6 denses (L ≤ 3,3e-14 / 1,2e-8 eV, NL ≤ 8,1e-15 / 4,0e-9 eV, seuil 1e-6 eV). 11×11 exclu : son `.save` de maille a été écrasé par un run `bands`.
  - D4 9×9, état π : (a1) −0,727 eV (w₂ 0,236), (a3) −0,734, (c-3) −0,677, contre QE −0,737 (v1 : (a1) −1,350). Décalage rigide résiduel +24,3 meV. Paire σ à −0,812 eV à 5 WF contre QE +0,101 ×2 ; elle est retrouvée à +0,116 eV avec 128 bandes.
  - Pôle 9×9 avec M2 : minimum de det et de λ à −0,812 eV (λ = 0,0019 i), porté par le bloc σ. Le bloc π n'a pas de zéro (min |λ| 0,40 vers −0,17). Pic de −Im T̄(K) à −0,177 eV (v1 : −1,292).
  - 9×9, v1 → v2 : médiane Γ·N_cells 2 473,55 → 3 132,60 meV ; E_res −1,238 → −0,175 eV ; p_z–p_z sur la lacune 6,617 → 31,521 eV ; ⟨‖M^NL‖_F⟩/⟨‖M^L‖_F⟩ 20,44 → 0,252 (0,25 pour les six tailles) ; Friedel ∫δρ −0,057 → −1,0007 ; Born/T médian 3,30 → 46,5.
  - Les paramètres gelés du 2026-09-05 restent valides pour les six tailles (9×9 : C11 grille 0,49 %, C11 η 0,10 %, C10 0,47 %). Test d'or 5×5 : 1,80e-13, PASS. Chapitre 5 : Γ^ed/Γ^ep médian 0,361 → 0,537.
- **Interprétation / décision** : M^L était bien N_cells fois trop petit. M2 devient la convention de production : `production.json` porte `M_normalization = v2` (et, jusqu'à R10, `results_dir = results/M2`), et une porte v2 refuse les v1 dans tous les lecteurs. `results/M/` est gelé, jamais modifié. D4 reproduit à l'identique la variante « M^L × 81 » de R5. Le pic de −Im T̄ suit maintenant la partie locale. Le test « convention intensive » (7,7e-2 pour un seuil de 5e-2) est marqué « À VOIR » et rapporté sans jugement. Pas de miroir des M2 denses (Greg, 2026-09-27) : ils se reconstruisent en ≈ 11 min avec `scripts/assemble_M2.py`. Suite : R10 a remplacé les produits de l'étape 3 (`results/M2_plateau`) et le ch. 4 est réécrit depuis `campagnes/M/ch4/table_v1_final.md`.
- **Où sont les données** :
  - Campagne : `R6_rapport.md` ; `etape1/gate/gate_table.md` et `gate_*.json`, `etape1/kernel/*.json`, `etape1/ml/M_L_*_v2.json`, `etape1/assemble_summary.jsonl`, `etape1/MD5SUMS_2026-09-25.txt`, `etape1/states/states_tables.md` ; `etape2/d4/d4_tables.md` (+ `d4_results.json`, `prep_results.json`), `etape2/b3/b3_tables_M2.md`, `etape2/d3/d3_tables_M2.md` (+ json) ; `etape3/table_v1_v2.md`, `etape3/level1_gate.md`, `etape3/csv_v2/*.csv`, `etape3/diff_3.0.patch` ; `phase0/*.diff`.
  - Dépôt : `results/M2/` (npz et csv de l'étape 3, `README.md` avec la recette de reconstruction, `MD5SUMS_2026-09-25.txt`). Les 43 fichiers M `.npy` sont sur `/project` et ne sont pas dans le dépôt. `results/M/` = v1 gelé. Pas de `.in`/`.out` (aucun calcul QE).

---

- But : R5 (A.2) a montré que dans tous les fichiers M de production, M^L est calculé avec des états normés sur la
  super-cellule (noyau `compute_ML_R*`, ψ = u e^{ik·r}/√Ω_sc) et M^NL avec des états normés sur la maille (`compute_M_NL`,
  4π/√Ω_uc) : M^L est N_cells fois trop petit. Corriger le noyau local (Greg), réassembler les M existants avec le facteur
  N_cells sur M^L (M2 = N_cells·M^L + M^NL, `results/M2/`), vérifier chaque M2 par application directe de ΔV (porte A.2,
  ≤ 1e-6 eV), revalider la chaîne contre la super-cellule 9×9 (escalier D4), puis régénérer tous les résultats du chapitre 4.
- Prompt d'origine : R6 (Greg, 2026-09-25). Ordre : phase 0 → STOP → correction du noyau par Greg (étape G) → (GO 1)
  assemblage et portes → STOP → (GO 2) validation super-cellule → STOP → (GO 3) production du chapitre 4 → rapport → STOP.
- Statut : **PRODUCTION**. Aucun calcul QE. Anciens résultats `results/M/` gelés (README), jamais modifiés ni supprimés ;
  nouveaux résultats dans `results/M2/` (mêmes noms). `config/production.json` reçoit une clé de version (diff proposé,
  appliqué par Greg). Git en lecture. Aucune suppression sans manifeste et GO séparé.
- Dates : phase 0 rédigée le 2026-09-25 (STOP, attente de l'étape G puis du GO 1).
- Répertoire de travail (celui-ci, hors dépôt, règle 5 de CLAUDE.md) : `graphene/qe/defects/R6_production_corrigee/`, à côté
  de `R5_base_vs_M/` qu'il prolonge. Copie versionnée : `graphene-raman/campagnes/R/R6_production_corrigee/`.
- Contenu (phase 0) : `R6_rapport.md` (inventaire 0.1–0.5, plan, coûts), `phase0/*.diff` (diffs proposés, NON appliqués :
  `local_R.py`, `local_G.py`, `config.py`, `matrix_io.py`, `compute_M.py`, `production.json`),
  `phase0/README_results_M_gele.md` (texte du README à déposer dans `results/M/` au GO 1).
- État du dépôt au moment de la phase 0 : HEAD a5afdd4 (EM1), `campagnes/R/R5_base_vs_M/` non suivi, `campagnes/R/R4_quasi_lie/R4_rapport.md`
  modifié (erratum R5), rien commité par Code.
