# R6 — Production corrigée : normalisation de la partie locale de M (chapitre 4)

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
  de `R5_base_vs_M/` qu'il prolonge. Copie versionnée : `graphene-raman/article/R6_production_corrigee/`.
- Contenu (phase 0) : `R6_rapport.md` (inventaire 0.1–0.5, plan, coûts), `phase0/*.diff` (diffs proposés, NON appliqués :
  `local_R.py`, `local_G.py`, `config.py`, `matrix_io.py`, `compute_M.py`, `production.json`),
  `phase0/README_results_M_gele.md` (texte du README à déposer dans `results/M/` au GO 1).
- État du dépôt au moment de la phase 0 : HEAD f4b7ec3 (EM1), `article/R5_base_vs_M/` non suivi, `article/R4_quasi_lie/R4_rapport.md`
  modifié (erratum R5), rien commité par Code.
