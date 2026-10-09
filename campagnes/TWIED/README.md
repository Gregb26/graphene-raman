# TWIED — référence 5×5 de la chaîne complète, pour les tests de twied

- But : référence produite par le code de graphene-raman, de bout en bout (fichiers QE et Wannier90 → M → M_W → amas → g₀ → t → Γ),
  qui stocke l'entrée et la sortie de chaque étape, pour tester les fonctions de `twied` une à une (`../twied`, ROADMAP Phase 1 :
  jeu de test 5×5 commis, exécutable sans le cluster).
- Prompt d'origine : Greg (2026-10-08), session locale ; pilote écrit par Greg, script de vérification et fiche par Code.
- Statut : PRODUCTION (2026-10-09), calcul local de quelques minutes (M^L domine), aucun `.save` produit ; seul emplacement.
- Entrées : `data/export_chaine_5x5/` (git-ignoré, rapatrié de rorqual ; son README et son MD5SUMS sont copiés ici,
  `export_chaine_5x5_{README.md,MD5SUMS.txt}`) : maille `defect_unit_cell_5x5_wann.save` (le NSCF qui a nourri la wannierisation),
  paire `defect_5x5_{p,d}.save` avec `Vks_5x5_{p,d}` et `C.upf` ; wannierisation `results/wannier/5x5` (`config.wannier_dir`).
  **Ne pas** prendre `data/graphene/unit_cell/qe/defect_unit_cell_5x5.save` : c'est un SCF (bande 16 non convergée) dont la jauge
  n'est pas celle de `wannier/5x5` (étalement recalculé 107 Å² contre 3,48 Å²).
- `make_ref_5x5.py` → `config.twied_dir` (`results/twied/`) : `ref_5x5_inputs.h5` (140 Mo, reste local) et `ref_5x5_chain.h5`
  (22 Mo, à copier dans `twied/tests/data/`) ; même provenance dans les deux (`run_id`, commit, md5 des 44 fichiers sources),
  mode strict (`git_clean`). Chaque dataset porte `units` et `axes` ; M en Ha jusqu'à `M_coarse`, eV ensuite (`provenance/HA2EV`).
- Paramètres : deux variantes d'alignement, `unaligned` (C_N = 0) et `kumagai_oba` (C_N = −0,05756 eV, config) ; R_cut = 3 (config ;
  toute la boîte 5×5, 25 mailles) et R_cut = 1 (5 mailles, test de troncature, étape 5 seulement) ; nk_int = 300, η = 0,02 eV ;
  g₀ et t sur 21 énergies E_D ± 1,5 eV ; Γ sur 60 × 60, fenêtre E_D ± 1 eV, ne_per_eta = 8 (production : 240², ± 3 eV, réduits
  pour que les tests de `twied` restent rapides).
- Arbre :
  ```
  ref_5x5_inputs.h5   provenance/{source_files/}      input/unit_cell/   C_nk, nG, G_red, k_red, eps, A_cols, Omega, x_red, ecut
                      input/supercell_{p,d}/  A_cols, Omega, x_red, V      input/pseudo/   ekb_li, fr_li, rgrid, V_L
                      input/wannier/  U, U_dis, k_w90, perm, H_R, R, ndegen, r_R, lattice      M_coarse/prep/  Ved, Omega_sc
  ref_5x5_chain.h5    provenance/{source_files/}      M_coarse/  M_L, M_NL, M
                      M_W/         U, U_dis, k_coarse (ordre QE), Mwk, Mwr_raw, R          (commun aux variantes)
                      G0/          H_R, R_w, ndegen, E_D, gap, R_cluster, egrid, g0         (commun : g₀ ne dépend pas de M)
                      Gamma_grid/  E_out, e_window                                          (commun)
                      variants/{unaligned,kumagai_oba}/  C_N, Mwr, in_box, Rn, R_d, dist, wt, R_cluster, M_cluster,
                                                         R_cut_1/{R_cluster, M_cluster}, t, Gamma
  ```
  Indice plat de M_cluster, g₀ et t : L·nw + w, dans l'ordre de `R_cluster` (non trié). Γ est rangé (bande, k), E_out (k, bande).
- `check_ref_5x5.py` : 59 contrôles en moins d'une seconde (structure, unités, appariement des provenances, chaque étape recalculée
  depuis les entrées stockées) ; `--sources --full` : 68 contrôles en 37 s (md5 et lecteurs, g₀ et Γ recalculés). Tous PASS le
  2026-10-09, recalculs exacts au bit.
- Relancer : `.venv/bin/python campagnes/TWIED/make_ref_5x5.py && .venv/bin/python campagnes/TWIED/check_ref_5x5.py --sources --full`
  (arbre de travail propre exigé pour `src`, `config` et `campagnes/TWIED`).
- Constat en route : sur la paire cohérente, le M_W de la grille grossière est localisé (41,5 eV sur site, voisins 1,4 à 2,1 eV,
  comme la grille dense). Les courbes « grille grossière » de `fig_locality_final` (`défauts.tex`, l. 419) sont donc un artefact de
  jauge (SCF + U du `_wann`), pas de repliement ; établi pour la 5×5, la 8×8 est à refaire sur rorqual avec son `_wann`.
- Suite : copie de `ref_5x5_chain.h5` dans `twied/tests/data/` (Greg) ; tests de `twied` sur ce fichier.
