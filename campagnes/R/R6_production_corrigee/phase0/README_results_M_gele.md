# results/M/ — GELÉ le 2026-09-25 (R6)

Résultats obtenus avec M^L non normalisé (facteur N_cells manquant), remplacés par R6.

- Tous les fichiers M de ce répertoire (`M_ed_*`, `M_L_*`, `M_NL_*`, `M_L_dense_*`, `M_NL_dense_*`, `M_dense_*`, `_test_mnl/`,
  `*_coarsecheck.npy`) et tout ce qui en dérive (`specwd_*_prod.npz`, `resonance_*.npz`, `resonance_criteria_*.npz`,
  `resigma_*.npz`, `mwr_locality.npz`, `M_analysis.npz`, `lnl_frobenius.csv`, `level1_summary.csv`, `level2_*.csv`,
  `m_rcut_*.csv`, `nkint_check_9x9.csv`, `M_tests_summary.csv`, `gamma_*`, `dos_*`, `convergence.*`, `logs/`) ont été
  calculés avec la partie locale M^L des noyaux `compute_ML_R*` en norme super-cellule (ψ = u e^{ik·r}/√Ω_sc) alors que
  M^NL est en norme maille (4π/√Ω_uc) : M^L y est N_cells = N² fois trop petit par rapport à M^NL (constat R5-A.2,
  `graphene/qe/defects/R5_base_vs_M/R5_rapport.md` §A.2, 2026-09-25). Les sidecars `bloch_norm = unit_cell` de ces
  fichiers ne valent que pour M^NL.
- Ne pas modifier, ne pas supprimer (manifeste et GO séparé pour tout nettoyage). Les scripts de production lisent
  désormais `results/M2/` (clé `results_dir` de `config/production.json`, version `M_normalization = v2`).
- Fichiers NON concernés par le facteur (ils ne lisent pas M) : `ks_reconstruction.npz` (reconstruction KS sur la maille,
  norme Ω_uc), `ved_analysis.npz`, `sampling_table.csv`.
- Remplacés par : `results/M2/` (R6, rapport `graphene/qe/defects/R6_production_corrigee/R6_rapport.md`, copie
  `article/R6_production_corrigee/`).
