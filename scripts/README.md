# scripts/ — scripts de production, d'analyse, de figures et de validation (rangés le 2026-09-30, ménage étape A)

Tous se lancent depuis la racine du dépôt (`.venv/bin/python scripts/<dossier>/<script>.py`, `sbatch scripts/slurm/<lanceur>.sh`) ;
les lanceurs font `cd $PROJ`. Les trois helpers partagés restent ici : `_palette.py` (palette et colormaps des figures),
`_bands.py` (chemins de bandes, points de symétrie), `_paths.py` (chemins des `.save` locaux, `EDI_DATA`) ; les pilotes de
`campagnes/` les importent par `sys.path.insert(…, "scripts")`.

| Dossier | Contenu | Lancé par |
|---|---|---|
| `m/` | chaîne M : `compute_M.py` (grossier, MPI), `compute_M_dense_stages.py` (dense, zero-padding), `assemble_M2.py` (M2 = N_cells·M^L + M^NL), `gate_M_normalization.py` (porte A.2), `tag_vacancy_sublattice.py`, `finalize_wannier.py`, `link_data.sh` | `slurm/submit_M.sh`, `submit_M_dense.sh` ; manuel |
| `t/` | chaîne T et analyses du ch. 4 : `compute_spectral_wannier.py`, `rcut_resigma.py`, `resonance_metrics.py`, `resonance_criteria.py`, `nkint_check_post.py`, `m_rcut_convergence.py`, `mwr_locality_coarse_vs_dense.py`, `level2_families.py`, `analyze_M.py`, `analyze_Ved.py`, `ks_reconstruction_all.py`, `sampling_table.py`, `lnl_frobenius_all.py`, `epw_ed_vs_ep.py` (pont ch. 5) | `slurm/submit_*.sh`, `submit_post.sh`, `campagnes/R/R10_plateau/submit_r10.sh` |
| `epw/` | post-traitement EPW (ch. 5) : `epw_pp_save`, `epw_extract_gkk`, `epw_validate`, `epw_selfen_post`, `epw_phself_post`, `epw_d2_extract`, `epw_ring_check`, `epw_dfpt_path_freq`, `epw_phdos_extract` | `slurm/submit_epw_p1_post.sh`, `submit_epw_p2_post_mv.sh` |
| `fig/` | figures du mémoire : `make_figures.py` (ch. 4, travail), `make_figures_memoire.py` (ch. 4, finales), `make_figures_epw.py` (ch. 5 ; options de production dans `slurm/submit_post.sh figures`), `make_figures_em.py` (§2.5) | `slurm/submit_post.sh figures`, `submit_r10.sh` ; local |
| `validation/` | scripts autonomes PASS/FAIL (code 0/1), **tous enveloppés par pytest** (`tests/test_scripts_{M,tmatrix}.py`, marqueurs `slow`, `needs_data`, `cluster`) : `test_ks_reconstruction`, `test_wannier`, `test_zero_pad_dense`, `test_pad_vs_full_supercell` (cluster), `test_local_tmatrix`, `test_local_rcut`, `test_local_green_batch`, `test_local_tmatrix_real` (cluster) ; `validate_wannier_bands`, `compare_bands_qe`, `compare_bands_w90_qe`, `run_test_A_batch` | local ; `slurm/submit_golden_dense.sh`, `submit_test_A.sh` |
| `slurm/` | les 11 lanceurs SLURM (rorqual, MPICH `srun`) | `sbatch scripts/slurm/…` |

Les tests pytest sont dans `tests/` (`python -m pytest tests`). Scripts obsolètes retirés le 2026-09-30 (dans l'historique git,
dernier commit qui les contient : `d9a3588`) : `compute_spectral`, `compute_tmatrix`, `compute_convergence`, `_eta_scan`, `_normtest`,
`_mcheck`, `migrate_M_norm`, `_old_vs_new_7x7`, `check_M_dense_nb20_vs_nb16`, `check_M_dense_vs_coarse{,_bands}`, `check_ML_coarse_kernel`,
`check_onsite_and_NL`, `validate_ML_grid_7x7`, `summarize_level1_maps`, `_bz_ratio_LNL`, `_diag_mnl_mpi`, `_test_mnl_mpi.sh`,
`submit_spectral.sh`, `submit_tmatrix.sh`, `submit_spectral_wannier.sh`.
