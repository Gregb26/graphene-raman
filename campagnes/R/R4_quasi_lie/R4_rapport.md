# R4 — Diagnostic de l'état quasi-lié de la lacune : rapport de campagne

Statut TEST. Prompt d'origine : R4 (Greg, 2026-09-25). Répertoire de travail : `graphene/qe/defects/R4_quasi_lie/` (hors dépôt).
Copie versionnée : `graphene-raman/article/R4_quasi_lie/` (au rapport final). Chemins relatifs à `/project/6004866/gregb26/`
(`$PROJECTS`) ; le dépôt est `graphene-raman/`. Unités : énergies en eV sauf mention ; les XML de QE sont en Hartree.

## Phase 0 — Localisation et interfaces (2026-09-25)

Rien de la campagne n'a été calculé. Seule évaluation faite, sur le nœud de connexion : g⁽⁰⁾ à une énergie (§0.3), grille
interne 120², 2 s. Le miroir `qe_tmp_backup/` et le scratch `qe_tmp/` portent les mêmes `.save` (md5 2026-09-17 et 2026-09-22).

### 0.1 Chemins

| Donnée | Chemin | Constat |
|---|---|---|
| 9×9 non relaxée, défaut (ch. 4) | `graphene/qe/defects/super_cell/9x9/defective/` : `scf.in`, `scf.out`, `pp.in`, `submit.pp`, `Vks_9x9_d` (241 Mo) | 161 atomes, lacune A à (13/27, 13/27), `K_POINTS gamma`, nbnd 352, nelec 644, E_F = −0,1551537 Ha = −4,2219 eV |
| — `.save` | scratch `/scratch/gregb26/qe_tmp/defect_9x9_d/defect_9x9_d.save/` ; miroir `graphene/qe/qe_tmp_backup/defect_9x9_d/` ; lien `graphene-raman/data/graphene/supercell/qe/defect_9x9_d.save` | `wfc1.hdf5` 2,18 G : **gamma_only = .TRUE.**, 352 bandes, 385 710 ondes planes (demi-sphère), `charge-density.hdf5`, `data-file-schema.xml`, `C.upf` |
| 9×9 non relaxée, parfaite | `super_cell/9x9/pristine/` (`Vks_9x9_p`) ; `.save` scratch/miroir `defect_9x9_p` | 162 atomes, nbnd 354, nelec 648, E_F = −0,1550643 Ha = −4,2195 eV ; gamma_only |
| Maille grossière 9×9 | scratch `qe_tmp/defect_unit_cell_9x9/defect_unit_cell_9x9.save` ; lien `data/graphene/unit_cell/qe/defect_9x9.save` ; miroir `qe_tmp_backup/defect_unit_cell_9x9/` | 16 bandes, 81 k (grille 9×9 complète, nosym + noinv), gamma_only = .FALSE., E_F = −0,1550643 Ha (= parfaite 9×9) |
| Maille dense 27×27 (p = 3) | scratch `qe_tmp/defect_uc_dense_27/defect_uc_dense_27.save` ; miroir `qe_tmp_backup/defect_uc_dense_27/` ; inputs `graphene/qe/defects/unit_cell/27x27/` | 20 bandes, 729 k, `diago_full_acc` ; la grille 9×9 est le sous-ensemble d'indices {0, 3, 6, …} de la 27×27 (81 points, vérifié) |
| M grossier 9×9 | `graphene-raman/results/M/M_ed_9x9.npy` (+ `.json` : unit_cell, hartree, « raw, pre-migration »), `M_L_9x9.npy` (2026-06-18, **sans sidecar**), `M_L_dense_9x9_coarsecheck.npy` (+ `.json`, part `M_L_coarse_check`, noyau R, 2026-09-04) | forme (16, 81, 16, 81) ; M^NL grossier n'existe pas en fichier (= M_ed − M_L) |
| M dense 9×9 | `results/M/M_dense_9x9.npy`, `M_L_dense_9x9.npy`, `M_NL_dense_9x9.npy` (+ `.json` : unit_cell, hartree, p 3, D 27) | 3,4 G chacun, forme (20, 729, 20, 729) |
| Wannier 27×27 | `graphene-raman/wannier/27x27/` : `wannier_tb.dat` (H(R), 741 R, ndegen ∈ {1, 2}), `wannier_u.mat`, `wannier_u_dis.mat`, `wannier.wout`, `wannier.eig`, `wannier_manifest.json` (run_id `baa17b88b1b69e51`) | `_hr.dat` seulement dans `graphene/qe/defects/unit_cell/27x27/wannier_hr.dat` (hors manifeste) ; centres : voir §0.3 |
| Sorties matrice T de production | `results/M/resonance_9x9.npz`, `resonance_criteria_9x9.npz`, `specwd_9x9_prod.npz`, `resigma_9x9_rc{0123,3_nk150..600,3_npe2,3_npe32,4}.npz`, `mwr_locality.npz`, `M_analysis.npz`, `ved_analysis.npz`, `ks_reconstruction.npz`, `lnl_frobenius.csv` ; logs `results/M/logs/resonance_criteria_9x9_20414406.out`, `resonance_9x9_20421370.out`, `specwd_20294202.out`, `three_nb20_20226399.out`, `nkint_21337627.out` | tous présents |
| R1 (relaxée 9×9) | `graphene/qe/defects/super_cell_relaxed/9x9/` : `nspin1/`, `nspin2/` (`relax.in`, `relax.out`), `R1_rapport.md`, `R1b_rapport.md`, `k3x3/`, `projwfc/` ; copie `graphene-raman/article/R1_vacancy_relaxed/` | `.save` scratch `qe_tmp/vacancy_relaxed/nspin1/vac_9x9_relax_nspin1.save` (`wfc1.hdf5`, gamma_only, 352 bandes, E_F −0,1552269 Ha = −4,2239 eV) et `nspin2/vac_9x9_relax_nspin2.save` (`wfcup1.hdf5` ispin 1 + `wfcdw1.hdf5` ispin 2, lsda, nbnd_up = nbnd_dw = 352, E_F −0,1552444 Ha = −4,2244 eV, `atomic_proj.xml`) ; miroirs `qe_tmp_backup/vacancy_relaxed/nspin{1,2}/` (md5 2026-09-22) |
| R1b projwfc | `super_cell_relaxed/9x9/projwfc/gamma_nspin2/` (R1 nspin2 à Γ) et `projwfc/k3x3_relax2_nspin2/` | par répertoire : `projwfc.out` (décomposition par état « ==== e(n) ==== psi = … », \|psi\|²), `proj_*.projwfc_{up,down}` (filproj), 322 fichiers `pdos_atm#N(C\|C1)_wfc#{1(s),2(p)}` avec colonnes m-résolues (E, ldos↑↓, pdos↑↓ × 3 m), `pdos_tot`. Atomes 64, 80, 81 = espèce C1. Les `.save` de `k3x3/` ont été supprimés le 2026-09-23 (xml seuls). **Aucun projwfc pour la 9×9 non relaxée ni pour R1 nspin1.** |
| C.upf | `abinit_processing/pseudo/C.upf` (md5 `34a24e64c0a39f27c6c36b90a16ac686`) | copies identiques (md5) dans `defect_9x9_d.save`, `defect_uc_dense_27.save`, `vac_9x9_relax_nspin1.save` |
| Config figée | `graphene-raman/config/production.json` | R_cut 3, grid 240, η 0,02 eV, nk_int 300, e_window 3, ne_per_eta 8, K_red (2/3, 1/3), dense 9x9 : p 3, D 27 |
| Potentiels pp.x (D5) | `super_cell/NxN/{defective,pristine}/Vks_NxN_{d,p}` pour N = 5, 6, 7, 8, 9, 10, 11, 12 | **aucun potentiel pp.x pour R1** (nspin1, nspin2) |

Valeurs de référence (lecture seule) :

| Valeur attendue | Retrouvée | Source |
|---|---|---|
| M_W p_z–p_z site lacunaire = 6.617 eV | 6.616636913763318 eV | `mwr_locality.npz` clé `9x9_dense_onsite_pzvac` (= `9x9_dense_onsite_pzA`, `vac_sublattice` A) |
| ‖M_W(0,0)‖ = 9.365 eV | 9.36469608639188 eV | `mwr_locality.npz` clé `9x9_dense_onsite_norm` |
| min \|det\|/max = 2.09e-4 à −2.530 eV | 2.091e-04 à −2.530 eV | log `resonance_criteria_9x9_20414406.out` l. 6 |
| min_i \|λ_i\| = 0.0108, λ = +0.0003 + 0.0108i au même point | 0.0108, +0.0003+0.0108i à −2.530 eV | idem l. 13–14 |
| minima secondaires entre −2.19 et −1.97 eV | −2.185, −2.132, −2.082, −2.027, −1.972 eV (\|det\|/max 2.6e-2 … 4.1e-2) | idem l. 7–11 |
| pic de Γ_T à −1.24 eV | −1.2400 eV | `resonance_9x9.npz` clé `peak_GT` |
| pic de −Im T̄ à −1.29 eV | −1.2925 eV | clé `peak_ImTbar` |
| Re T̄(E_D) = 2.521 eV | 2.521257922195672 eV | clé `ReTbar_at_ED` |
| C14 : décalage G = 0 de 67.0 meV → écart relatif max 6.4e-4 | 67.00 meV ; max rel 6.433e-04 (médian 1.817e-04) | log `resonance_9x9_20421370.out` l. 3 et 25 ; npz clé `ML_diag_mean` = 0.06700 eV |

Toutes les valeurs de référence sont retrouvées. E_D (Wannier, 90×90) = −4,238896 eV (`resonance_9x9.npz` clé `E_D`).
Remarque : deux logs plus anciens (`resonance_9x9_20294415.out`, `_20294682.out`) utilisent K = (1/3, 1/3) (index 19280) et
donnent un pic de −Im T̄ à −1,185 eV ; la production (job 20421370, npz) utilise K = (2/3, 1/3) (index 38480).

### 0.2 Interface de la matrice T

Chaîne de production (scripts au niveau module, sans fonction `main` : `resonance_criteria.py`, `resonance_metrics.py`,
`rcut_resigma.py`, `compute_spectral_wannier.py`) :

1. `config.load_production()` → `dense_paths(cfg, "9x9")` (chemins du `.save` dense, de `M_dense_9x9.npy`, de `wannier/27x27/`, du manifeste).
2. Porte de jauge `wannier_provenance.load_wannier_checked(manifest)` (sha256 de tb / u / u_dis ; refuse sinon).
3. `matrix_io.load_M_checked(mfile, require_bloch_norm="unit_cell", units="eV")` : exige le sidecar `.json` (`bloch_norm = unit_cell`, `units = hartree`), multiplie par HA2EV une fois.
4. `qe_io.get_k_red(uc_dense)` → grille k de M ; `read_w90_mat(u)`, `read_w90_mat(u_dis)` réordonnés par `_match_kpoint_order` ; `read_w90_HR(tb)` → `(Hwr (741,5,5), Rw (741,3), ndegen)`.
5. `Mbk_to_Mwk(M, U, U_dis)` (V†MV) → `Mwk_to_Mwr(Mwk, k, MP)` (double TF, boîte 27×27) → `lt.recenter_mwr` (R_d = [4, 4, 0]) → `lt.mwr_locality` (garde-fou a) → `Rloc = R[‖R‖ ≤ R_cut]` → `lt.extract_V_loc(Mwr, R, Rloc)` (garde-fou b) : V_loc (145 × 145), index L·nw + w.
6. `lt.Hwr_to_Hwk(Hwr, Rw, k_int, ndegen)` ; `lt.local_green_batch(Hwk_int, k_int, Rloc, egrid, eta)` → g⁽⁰⁾ (nE, 145, 145) ; `lt.local_t(V, g0)` ; critère de pôle : lignes 39–54 de `resonance_criteria.py` (`slogdet`, `eigvals` de 𝕀 − V g⁽⁰⁾, minima locaux) ; T̄(K) : lignes 86–90 de `resonance_metrics.py` (φ_K = U_out[K] e^{2πiK·R}, bandes 3–4).

Entrées des routines du module `defects/many_body/local_tmatrix.py` : **tableaux seulement**. `scattering_rate(Hwr, Rw, ndegen,
V_loc, R_local, k_out, eta, k_int)`, `scattering_rate_fast(…)`, `local_green(Hwk, k_int, R_local, eps, eta)`,
`local_green_batch(Hwk, k_int, R_local, egrid, eta, deriv)`, `local_t(V_loc, g0)`, `extract_V_loc(Mwr, R_mwr, R_local)`,
`Hwr_to_Hwk(Hwr, Rw, k, ndegen)`. Aucun fichier, aucun manifeste, aucune config à ce niveau ; nw est déduit de `Hwr.shape[1]`
et V_loc doit être `(nL·nw, nL·nw)` hermitienne.

Réponse à 0.2 : **un modèle de liaisons fortes arbitraire est accepté sans modification** (H(R) + R + ndegen et V_loc en
tableaux) ; `scripts/test_local_tmatrix.py` le fait déjà avec un H aléatoire. Les seuls éléments non appelables tels quels
sont les 16 lignes du critère (det, λ, minima) et les 5 lignes de T̄(K), inline dans les scripts ; elles seront recopiées
dans une fonction nouvelle (§0.6) sans toucher aux scripts. Le chemin « fichiers » (manifeste, `.npy` + sidecar, `.save`)
n'est requis que pour obtenir V_loc de production ; il sera exécuté une fois (job J1) et le résultat mis en cache dans ce
répertoire.

### 0.3 Blocs miroir

Ordre des 5 fonctions de Wannier (`wannier.win` : `C1:sp2;pz`, `C2:pz` ; C1 = A à (1/3, 1/3), C2 = B à (2/3, 2/3) ; `wannier.wout` « Final State ») :

| index (0-based) | fonction | centre (Å) | étalement (Å²) | parité z → −z |
|---|---|---|---|---|
| 0, 1, 2 | sp² de A (centres sur les trois liaisons, à 0,71 Å de A) | (1,0677, ±0,6165, 0), (2,1355, 0, 0) | 0,6109 ×3 | paire |
| 3 | p_z de A (atome de la lacune de production) | (1,4237, 0, 0) = A | 0,9670 | impaire |
| 4 | p_z de B | (2,8473, 0, 0) = B | 0,9670 | impaire |

Tous les centres sont à z = 0 (plan des atomes) ; le plan miroir est z = 0. H(R = 0) (eV, `wannier_tb.dat`) : sp² sur site
−15,085, sp²–sp² −2,154, p_z sur site −3,9125, p_z(A)–p_z(B) −2,909 (à R = 0, (−1, 0, 0), (0, −1, 0) : les trois premiers
voisins de A), p_z–p_z même sous-réseau 0,222 (R = ±a₁, ±a₂), p_z(A)–p_z(B) 0,021 aux autres R.

Éléments hors bloc σ–π (max de la valeur absolue, indices w ∈ {0,1,2} contre w ∈ {3,4}) :

| objet | max \|hors bloc\| | max \|élément\| | comment |
|---|---|---|---|
| H(R), 741 R | 5,5e-10 eV | 15,08 eV | lecture de `wannier_tb.dat` |
| g⁽⁰⁾(ε = E_D − 1 eV), η 0,02, grille interne 120², amas R_cut 3 (29 mailles, dim 145) | 1,1e-11 eV⁻¹ | 0,147 eV⁻¹ | `local_green_batch`, 2 s (nœud de connexion). g⁽⁰⁾ p_zA–p_zA sur site = 0,1274 − 0,0737i eV⁻¹ |
| M_W(0, 0) (bloc sur site 5 × 5, 9×9 dense) | < 5e-4 eV (imprimé 0.000 à trois décimales) | 6,617 eV | log `three_nb20_20226399.out` l. 44–49 (`scripts/check_onsite_and_NL.py`) ; σ diag 0,631, σ–σ 2,653, p_zA–p_zB 0,491, p_zB 0,044 |
| M_loc complet (29 mailles) | **non calculé en phase 0** (rotation de M dense : 8–10 min, 16 cœurs) | — | première sortie du job J1 |

Sous-blocs : `local_green_batch` renvoie (nE, nL·nw, nL·nw) avec l'index L·nw + w ; le bloc π (29 × 2 = 58) et le bloc σ
(29 × 3 = 87) s'obtiennent par sélection d'indices (`np.ix_`) sans modification ; idem pour V_loc. Le critère de pôle étant
inline (§0.2), la fonction nouvelle prendra l'ensemble d'indices en argument.

### 0.4 Parties locale et non locale de M

| niveau | stockage séparé | fichiers | constat |
|---|---|---|---|
| grossier 9×9 (16 bandes, 81 k) | partiel | `M_ed_9x9.npy` (total, sidecar), `M_L_9x9.npy` (sans sidecar : `load_M_checked` refuse, `np.load` direct), `M_L_dense_9x9_coarsecheck.npy` (M^L recalculé le 2026-09-04, noyau R, sidecar) | M^NL = M_ed − M_L ; les deux M^L grossiers seront comparés (D4) |
| dense 27×27 (20 bandes, 729 k) | oui | `M_L_dense_9x9.npy`, `M_NL_dense_9x9.npy`, somme `M_dense_9x9.npy` (`compute_M_dense_stages.py`, étapes ml / nl / combine) | 3 × 3,4 G |
| Wannier (M_W(R, R′), V_loc) | non | rien : chaque script de production refait la rotation V†MV et la double TF (8–10 min, MaxRSS 10–15 G) | — |

Coût d'un recalcul séparé écrit dans le répertoire de travail : trois rotations (M, M^L, M^NL) ≈ 30 min sur 16 cœurs, sorties
`Mwr_{tot,L,NL}_9x9.npz` de 213 Mo chacune ((5, 729, 5, 729) complexe) + les blocs V_loc (145 × 145) ; c'est le job J1.

### 0.5 C.upf (tel quel, unités du fichier)

`PP_HEADER` : ONCVPSP 3.3.0 (D. R. Hamann, 2017-10-31), `pseudo_type NC`, scalaire-relativiste, PBE, Z_val 4,00,
`core_correction T`, `l_max 1`, `l_local −1`, `mesh_size 1248` (pas 0,01 bohr, r_max 12,47), `number_of_wfc 2`,
`number_of_proj 4`. Projecteurs : β1 ℓ = 0, β2 ℓ = 0, β3 ℓ = 1, β4 ℓ = 1, `cutoff_radius_index 132` (r = 1,31 bohr) pour les quatre.

`PP_DIJ` (4 × 4, Ry, diagonale ; hors diagonale 0) :

| i | ℓ | D_ii (Ry) |
|---|---|---|
| 1 | 0 | +12,963096312 |
| 2 | 0 | +0,77100728970 |
| 3 | 1 | −8,3999228241 |
| 4 | 1 | −1,7550333227 |

`pseudo_io.read_upf` lit ces valeurs, les groupe par ℓ (2 canaux par ℓ) et divise par 2 (Ry → Ha) ; `compute_M_NL`
utilise le préfacteur 4π/√Ω_uc.

### 0.6 Plan d'exécution

Constats qui conditionnent le plan :

- Les wfc à Γ (9×9 défaut/parfaite, R1 nspin1/2) sont en **gamma_only** (demi-sphère G, ψ(−G) = ψ*(G)) ; `qe_io._read_all_wfc`
  ne gère ni ce cas ni la paire `wfcup1/wfcdw1` (glob `wfc*.hdf5`, tri par attribut `ik` = 1 pour les deux). D1 exige un
  lecteur nouveau (parité ⟨σ_h⟩ = Σ_G C*(g₁, g₂, −g₃) C(g₁, g₂, g₃), complétée par conjugaison). Fenêtre
  [E_D − 3, E_D + 1] eV ≈ quelques dizaines de bandes autour de 322 : lecture par tranche de bandes, ≤ 6 G.
- Aucun potentiel pp.x pour R1 : l'alignement de D1 sur les géométries relaxées et la décroissance de ΔV de R demandent
  **pp.x** (exécutable QE de post-traitement, 44 s sur 1 cœur pour la 9×9 du ch. 4, job 14406255). Listé en J0, GO explicite.
- Pas de projwfc pour la 9×9 non relaxée ni R1 nspin1 : D1 préparera `projwfc.in` sans les lancer ; les poids par état
  viendront de `projwfc.out` (R1 nspin2 Γ) et, pour les autres cas, de la parité + du lecteur d'ondes planes.
- D4 (a) « nbnd 16 puis 20 » : le M grossier à 20 bandes est la restriction de `M_dense_9x9` aux 81 k du sous-réseau 9×9
  de la 27×27 (indices 0, 3, 6, …) ; c'est aussi le M de (b) (matrices U aux mêmes k).

Coûts estimés (16 cœurs, 64 G, comme `scripts/submit_rcut_resigma.sh` ; références : rotation 8–10 min, g⁽⁰⁾ 300² sur 2 417
énergies ≈ 3 min, `resonance_criteria` 1 h 34 pour 2 401 + 6 790 énergies) :

| étape | contenu | coût | job |
|---|---|---|---|
| prep | 3 rotations (M, M^L, M^NL) → `Mwr_*.npz`, V_loc (R_cut 3) total/L/NL, max hors bloc de M_loc, D2 base de Wannier, D2 base de Bloch à K et K′ (lignes mmap de M dense) | 40 min, ~15 G | J1 |
| D6.1 | modèle synthétique (vérification bandes, van Hove, g₀ sur site vs somme directe), 10 valeurs de U + site retiré, g⁽⁰⁾ 300² (2 min) et 1200² (~30 min, W par bloc ≈ 1,8 G) | 40 min | J2 |
| D6.2 | même protocole sur le vrai H(R) (g⁽⁰⁾ 300²/1200² recalculés pour ce H), + bloc π de M_loc | 40 min | J2 (après J1) |
| D1 | valeurs propres XML (Ha → eV), parité par ondes planes (4 fichiers de 2,2 G, par tranches), alignement Lu 2019 (`get_pot` sur `Vks_9x9_{d,p}`), PDOS/projwfc.out, tableau + figure | 20 min, ≤ 8 G | J3 |
| D4 | (a) diag(ε_uc) + M/81 en base repliée (dim 1296 et 1620, variantes L/NL) ; (b) idem projeté sur 5 WF (U aux 81 k) ; (c) M_W dense replié (R_cut 3 puis toutes mailles), H_SC 405 × 405 à Γ | 15 min | J4 (après J1) |
| D2 | (fait dans J1) comparaison à t, U_c, Kaasbjerg ; ⟨ΔV^L⟩ = 67,0 meV, A_sc = 81 × 5,24 Å² | — | J1 |
| D3 | g⁽⁰⁾ 300² (cache de J2) et 600² (~12 min) sur [E_D − 3, E_D + 1] (1 601 énergies) ; 8 α × 3 blocs × eig 145 (~10 min) ; T̄(K) ; variantes L/NL | 40 min | J5 (après J1, J2) |
| D5 | `ved_analysis.npz` existe déjà (moyenne 3D, profil a₁, bord, 8 tailles) ; ajout de la valeur au bord après alignement D1 (9×9) | 5 min | J3 |
| R | lecture de code (faite en partie, §ci-dessous) + décroissance de ΔV pour R1 nspin1 si J0 | 5 min | J3 (après J0) |

Lecture de code pour R (à confirmer dans le rapport final) : `compute_M_NL` lit les positions de **tous** les atomes de
`sc_p` et de `sc_d` séparément (`tau_s_p`, `tau_s_d`) et fait M_d − M_p → une géométrie relaxée est traitée par construction ;
aucun indice de spin nulle part (`compute_M_NL`, `compute_ML_R`, `get_pot`, `_read_all_wfc` : `npol = 1` seulement) → nspin = 2
exigerait deux potentiels pp.x (`spin_component` 1 et 2) et deux passages.

sbatch envisagés (tous `--account=rrg-cotemich-ac`, 1 nœud, `module restore qe; module load scipy-stack`, `.venv/bin/python`,
sorties dans ce répertoire) :

| job | contenu | ressources | dépendance | écart à `config/production.json` |
|---|---|---|---|---|
| J0 (optionnel, exécutable QE) | `pp.x` plot_num 1 sur `vac_9x9_relax_nspin1` (spin_component 0) et `vac_9x9_relax_nspin2` (spin_component 1, 2) → `Vks_R1_nspin1`, `Vks_R1_nspin2_{up,dw}` (3 × 241 Mo, ici) | 1 cœur, 32 G, 15 min | — | hors config (données nouvelles, lecture seule des `.save`) |
| J1 `r4_prep` | rotations, V_loc, D2 | 16 cœurs, 64 G, 1 h 30 | — | aucun (R_cut 3, η 0,02) |
| J2 `r4_d6` | D6.1 + D6.2 | 16 cœurs, 64 G, 2 h | afterok J1 | nk_int 1200 en plus de 300 ; fenêtre [−3, +1] |
| J3 `r4_d1` | D1 + D5 + R | 16 cœurs, 32 G, 1 h | afterok J0 si lancé | — |
| J4 `r4_d4` | D4 | 16 cœurs, 64 G, 1 h | afterok J1 | — |
| J5 `r4_d3` | D3 | 16 cœurs, 64 G, 2 h | afterok J1, J2 | nk_int 600 en plus de 300 ; fenêtre [−3, +1] ; grille de sortie 240² non utilisée (T̄ à K seulement) ; α ≠ 1 = diagnostic |

Paramètres inchangés : η = 0,02 eV, R_cut = 3 (29 mailles, norme des indices réduits), ne_per_eta = 8 (pas 2,5 meV),
E_D Wannier (min du gap sur 90×90), K = (2/3, 1/3), porte de jauge C1, sidecars C2, recentrage R_d = [4, 4, 0].
Total ≈ 6 h de jobs à 16 cœurs, séquencés par dépendances ; aucun job n'est lancé en phase 0.

Fonctions à ajouter (modules **nouveaux** dans `src/electron_defect_interaction/` ; aucun fichier existant modifié) :

| module | signature | rôle |
|---|---|---|
| `io/qe_gamma_io.py` | `read_wfc_gamma(save_dir, spin=None, bands=None) -> (C, mill)` | lit `wfc1`/`wfcup1`/`wfcdw1` gamma_only, tranche de bandes, demi-sphère |
| idem | `mirror_parity_z(C, mill) -> parity (nb,)` | ⟨σ_h⟩ par état, g₃ → −g₃, complétion par conjugaison |
| idem | `get_eigenvalues_spin(save_dir) -> (nspin, nbnd) Ha` | découpe du bloc 2·nbnd des XML lsda |
| `io/projwfc_io.py` | `read_projwfc_states(projwfc_out) -> list[(e, {state: w}, psi2)]` | blocs « ==== e(n) ==== » |
| idem | `read_pdos_m(path) -> (E, ldos, pdos_m)` | colonnes m-résolues des `pdos_atm#…` |
| `defects/alignment.py` | `far_atom_alignment(V_d, V_p, x_d, x_p, A, radius) -> shift, atom` | moyenne sphérique de V_d − V_p autour de l'atome le plus loin de la lacune (Lu 2019) |
| idem | `rigid_shift_fit(eps_ref, eps_model, exclude_window) -> shift` | décalage rigide ajusté hors fenêtre (D4) |
| `defects/many_body/pole_criterion.py` | `block_indices(nL, nw, wfs) -> idx` | indices L·nw + w d'un bloc de parité |
| idem | `det_eig_criterion(V, g0, idx=None) -> (logdet_rel, minlam, lam_min, vec_min)` | réplique de `resonance_criteria.py` l. 41–44 sur un sous-bloc, avec vecteur propre |
| idem | `local_minima(y, x, n)` ; `sign_changes(y, x)` | copie de l. 45–47 ; racines de Re λ |
| idem | `tbar_pi_K(t_cache, phi_K) -> (Tbar (nE,2,2), tr)` | réplique de `resonance_metrics.py` l. 86–90 |
| `defects/many_body/tb_models.py` | `graphene_pz_tb(t=2.7, e_sigma=-1e3) -> (Hwr, Rw, ndegen)` | modèle π premiers voisins dans le bloc p_z, σ découplé, ordre de production |
| idem | `check_pz_model(Hwr, Rw, ndegen, t) -> dict` | bandes Γ-K-M, van Hove ±t, E_D = 0, g₀ sur site vs somme directe |
| idem | `onsite_vloc(R_local, nw, wf, U)` ; `removed_site_vloc(Hwr, Rw, R_local, nw, wf, U)` | V_loc des cas « U » et « site retiré » |
| `wannier/supercell_fold.py` | `bloch_folded_hamiltonian(eps_uc, M, N_cells) -> (eigs, vecs)` | D4 (a)/(b) : diag(ε) + M/N_cells |
| idem | `fold_hwr_to_supercell(Hwr, Rw, ndegen, N)` ; `fold_mwr_to_supercell(Mwr, R, N, R_cut=None)` | D4 (c) : H_SC et M_SC 405 × 405 à Γ |
| idem | `weights_near_vacancy(vec, R_local, nw, shells)` | poids p_z / sp² sur la lacune et ses voisins |

Pilote unique de campagne : `R4_quasi_lie/r4_driver.py` (sous-commandes `prep`, `d6`, `d1`, `d4`, `d2`, `d3`, `d5`, `r`,
chargement de `config/production.json`, porte de jauge, caches npz dans ce répertoire) et `submit_r4.sh` (une soumission par
sous-commande, dépendances ci-dessus). Aucun script jetable dans le dépôt.

**STOP — phase 0 terminée le 2026-09-25 ; attente du GO (couvre J1–J5 ; J0 = pp.x demande une décision explicite).**

## Phase D — Méthodes (GO du 2026-09-25 ; J0 à J5)

Précisions du GO appliquées telles quelles (sorties pp.x dans R1 ; vérifications J1 ; valeurs propres non triviales en D6 ;
E_D de super-cellule := quadruplet de la parfaite ; plan miroir z₀ ; deux rayons d'alignement ; échelle D4 (a1)→(a2)→(a3)→(b)→(c) ;
porte de régression D3 ; suivi de branche ; ΔV^NL indépendant du spin).

### Jobs (rrg-cotemich-ac, 1 nœud, `module restore qe; module load mpi4py/4.0.3 scipy-stack`, `.venv/bin/python`)

| job | sous-commande du pilote | ressources | dépendance | identifiant |
|---|---|---|---|---|
| J0 | `pp.x` plot_num 1 sur R1 : nspin1 (spin_component 0), nspin2 (1 puis 2) ; sorties dans `super_cell_relaxed/9x9/nspin{1,2}/` | 1 cœur, 32 G, 20/40 min | — | 21796633 (49 s), 21796634 (1 min 36) ; COMPLETED |
| J1 | `prep` : vérifications, rotations M / M^L / M^NL, V_loc, D2 | 16 cœurs, 64 G, 2 h | — | 21796852 (1re soumission 21796819 : `mpi4py` absent du job, corrigé) |
| J2 | `d6` : D6.1 puis D6.2 | 16 cœurs, 64 G, 3 h 30 | afterok J1 | 21796853 |
| J3 | `d1 d5 r` | 16 cœurs, 64 G, 2 h | afterok J2 | 21796854 |
| J4 | `d4` | 16 cœurs, 64 G, 2 h | afterok J3 | 21796855 |
| J5 | `d3` | 16 cœurs, 64 G, 4 h | afterok J2 | 21796856 |

Un échec de D6 (STOP demandé) annule automatiquement J3, J4, J5 (afterok). Aucun job n'est relancé automatiquement.

### Code

Modules nouveaux (`graphene-raman/src/electron_defect_interaction/`, aucun fichier existant modifié, non commités) :
`io/qe_gamma_io.py` (wfc gamma_only : lecture par tranche de bandes, parité ⟨σ_h⟩, densité en plan, disque), `io/projwfc_io.py`
(états atomiques et décomposition par bande de `projwfc.out`, pdos m-résolus), `defects/alignment.py` (site de la lacune, atome le
plus éloigné, moyenne sphérique, alignement Lu 2019, décalage rigide ajusté), `defects/many_body/pole_criterion.py` (critère det /
λ sur un sous-bloc = copie des lignes 41–47 de `resonance_criteria.py` avec vecteur propre ; valeurs propres non triviales sur le
support de V ; changements de signe ; recouvrements de branche ; T̄ de la paire π = lignes 87–89 de `resonance_metrics.py`),
`defects/many_body/tb_models.py` (modèle π, vérifications, V_loc « U » et « site retiré » ; matériel de test),
`wannier/supercell_fold.py` (H et M repliés dans la 9×9 à Γ, bases |nk⟩, |wk⟩, |wr⟩, poids en espace réel des vecteurs propres).
Tests unitaires du 2026-09-25 (nœud de connexion, 2 min) : modèle π exact (écart de bandes 0, E(M) = ∓2,7, g₀ batché = somme
directe à 1e-16) ; λ non triviale = 1 − U g₀ à 0 ; porte 1 du repliement 4,9e-14 eV ; repliement de M_W = F†(M_k/N²)F à 1e-15 ;
parité ±1 exacte sur des wfc gamma_only.
Pilote : `r4_driver.py` (sous-commandes `prep d6 d1 d5 r d4 d3`) ; `submit_r4.sh` ; journal `r4_log.txt` ; résultats
`prep/`, `d6/`, `d1/`, `d4/`, `d3/`, `d5/`, `r/` (json + npz), caches `cache/` (M_W, V_loc, g⁽⁰⁾ par grille), figures `fig/`.

### Définitions

- Fenêtre : ε − E_D ∈ [−3, +1] eV. E_D de super-cellule = moyenne du quadruplet dégénéré à Γ de la parfaite 9×9 (K et K′ repliés,
  états 323–326) ; E_D Wannier = milieu du gap minimal π/π* sur 90×90 (−4,238896 eV, production) ; chaque modèle de D4 est rapporté
  à son propre E_D (valeurs à K de la maille grossière 16 bandes, de la maille dense 20 bandes, de V†εV, de H(R)).
- Parité : ⟨σ_h⟩ = Σ_G C*(G) C(g₁, g₂, −g₃) e^{−4πi g₃ z₀} (z₀ = z moyen des atomes, réduit) ; demi-sphère gamma_only complétée par
  C(−G) = C*(G) ; « paire » = σ (⟨σ_h⟩ = +1), « impaire » = π (−1). Dans la base de Wannier : bloc σ = fonctions 0–2, bloc π = 3–4.
- Poids près de la lacune : w₂ = fraction de |ψ|² (sommée sur z) dans le disque en plan de rayon 2,0 Å (image minimale) centré sur
  le site de la lacune (contient les trois premiers voisins à 1,42 Å, pas les seconds à 2,47 Å) ; w₁ idem avec 1,0 Å. Fraction d'aire
  du disque de 2 Å : 2,96 %. État « localisé » : w₂ > 3 × ⟨w₂⟩ des états de la fenêtre de la parfaite (seuil explicite dans le rapport).
  Pour les modèles de D4, |Ψ|² est reconstruit en ondes planes de la super-cellule (Ψ = Σ d_nk ψ_nk, g_sc = 9(k + G)).
- Alignement (Lu 2019) : décalage = ⟨V_d⟩ − ⟨V_p⟩, moyenne du potentiel local pp.x (plot_num 1) dans une sphère (0,5 et 1,0 Å) autour
  de l'atome le plus éloigné de la lacune ; valeur de travail : 1,0 Å ; pour nspin2, potentiel moyen (↑ + ↓)/2. Énergies rapportées
  en ε − E_D (parfaite alignée) et relatives à E_F du calcul défectueux.
- D4 : (a) H = diag(ε_nk) + M/81 dans la base |nk⟩ repliée (M en norme unit_cell, eV) ; (b) H = V†εV ⊕ M_W(k,k′)/81 dans |wk⟩ ;
  (c) H_SC = Σ_{R≡r′−r} H(R)/ndegen + Σ_{R≡r,R′≡r′} M_W(R,R′) dans |wr⟩ (repliement modulo 9 des étiquettes recentrées) ; (c-3) : M_W
  restreint à |R|,|R′| ≤ R_cut = 3 avant repliement. Porte 1 : valeurs propres de H_SC = valeurs de Wannier aux 81 k (≤ 1e-8 eV).
  Porte 2 : valeurs propres de F†(V†εV)F + fold(M_W) = celles de (b) (≤ 1e-9 eV) — même partie à un corps des deux côtés, pour
  isoler le repliement de M ; l'écart littéral (b) [V†εV] contre (c-all) [H(R)] est rapporté à part (il contient le résidu
  d'interpolation H(R) contre V†εV sur la grille grossière).
- D3 : A(ε) = 𝕀 − αV g⁽⁰⁾(ε) sur le bloc π (58), σ (87) ou complet (145) ; |det|/max normalisé sur [−3, +1] ; λ_min = valeur propre de
  plus petit module, vecteur propre droit normé ; racines = changements de signe de Re λ_min ; saut de branche = recouvrement
  |⟨v_j|v_{j+1}⟩| < 0,9 entre énergies consécutives ; poids du vecteur propre sur p_z(lacune), p_z(3 voisins B), sp²(A, R = 0), autres.
  T̄(K) = trace/2 de ⟨nK| t(ε) |n′K⟩ sur la paire π (bandes 3, 4), t = V[𝕀 − g⁽⁰⁾V]⁻¹ (bloc ou complet).
- D6 : V à support restreint ⇒ valeurs propres non triviales de 𝕀 − V g⁽⁰⁾ = celles de 𝕀_S − V_SS g⁽⁰⁾_SS ; |det| non normalisé =
  exp(log|det|) (slogdet) ; racines = changements de signe de Re de la valeur propre non triviale de plus petit module.
- Figures : style `figures/memoire.mplstyle` et palette `scripts/_palette.py`, mais `text.usetex = False` (pas de LaTeX garanti sur
  les nœuds de calcul) ; textes en français.

## Phase D — Résultats (2026-09-25, jobs du matin ; rapport écrit le même jour)

Ordre du prompt : D6 → D1 → D4 → D2 → D3 → D5 → R. Tous les jobs ont abouti ; aucun calcul QE (pw.x) ; pp.x seulement (J0).
Chiffres bruts, aucune interprétation physique ; les constats demandés sont marqués en gras.

Problèmes rencontrés et traitement (règle « deux sortes de problèmes ») :

| n° | problème | sorte | traitement |
|---|---|---|---|
| 1 | J1 (21796819) : `mpi4py` absent de l'environnement du job (`defects/non_local.py` l'importe au chargement) | exécution | `module load mpi4py/4.0.3` ajouté à `submit_r4.sh` ; chaîne resoumise (21796852–56) |
| 2 | D4, porte 2 littérale : 3,7e-7 eV contre seuil 1e-9 | exécution / précision des données | cause : k du `.save` dense écrits à 1,3e-7 près (XML) ; porte réévaluée avec M_W recalculé à k = m/27 exacts : 2,4e-12 eV ; effet sur V_loc de production 1,6e-7 eV ; rapporté, production inchangée |
| 3 | J5 (21796856) : diagonalisation par lots 13 × trop lente (16 fils BLAS sur 145 × 145) ; annulé après 20 min | exécution | `_batched` (1 fil BLAS + fils Python) dans `pole_criterion.py`, g⁽⁰⁾ 600² calculé à part (J5a) ; résultats identiques à 1e-14 |
| 4 | D4 (21797745) : poids w₂ des variantes (c) faux (étiquettes recentrées non retraduites en mailles réelles) | exécution (bug du pilote) | phase e^{−2πi k·R_d} ajoutée ; (b) et (c-all) identiques après correction (J6, 21797985) |
| 5 | D1 : E_D de la parfaite ≠ E_F (−18,96 meV) ; alignement Lu à 0,5 Å instable pour R1 (déplacement des atomes de 1–3e-3 Å) | écart | rapporté tel quel |
| 6 | D4 : M grossier (juin) et M dense restreint diffèrent hors des valeurs singulières (16ᵉ bande dégénérée avec la 17ᵉ) | écart | rapporté (§D4) |

## D6 — Banc d'essai liaisons fortes (J2, job 21796853, 23 min 47, MaxRSS 36,8 G)

Mêmes routines que la production (`local_green_batch`, `extract_V_loc`-compatible V_loc 145 × 145, amas R_cut = 3 de 29 mailles,
η = 20 meV, pas 2,5 meV, fenêtre [−3, +1] eV), g⁽⁰⁾ calculé une fois par grille et par Hamiltonien (300² : 38–52 s ; 1200² : 604–616 s).
Pour un V à un seul élément, det[𝕀 − V g⁽⁰⁾] = 1 − U g₀,vv et la seule valeur propre non triviale est λ = 1 − U g₀,vv ; c'est
elle qui est rapportée (les 144 autres valent 1 exactement ; site retiré : 4 valeurs propres non triviales sur le support
{(R = 0, p_z A), (R = 0, p_z B), (R = (−1,0,0), p_z B), (R = (0,−1,0), p_z B)}).

### D6.1 — modèle synthétique (t = 2,7 eV, E_D = 0, σ découplé à −100 eV)

Vérification du modèle (erreur d'exécution si échec) : bandes sur Γ-K-M-Γ contre ±t|1 + e^{−2πik₁} + e^{−2πik₂}| : écart
1,8e-15 eV ; E(K) = ∓1,2e-15 eV ; E(M) = −2,7000 / +2,7000 eV (van Hove à ±t) ; gap minimal π/π* sur 300² : 1,8e-15 eV, E_D = 0 ;
g₀ p_z(A) sur site à ε = −1 eV, amas {0}, `local_green_batch` contre somme directe Σ_k [(ε + iη − H(k))⁻¹]₃₃/N :
0,105448 − 0,084741 i pour les deux (écart 2e-16). Le modèle a passé sa vérification.

Grille interne et g₀ : à 120² (§0.3) g₀ = 0,1164 − 0,0664 i (jouet) contre 0,1274 − 0,0737 i (H(R)) ; à 300² 0,10545 − 0,08474 i
contre 0,11446 − 0,09130 i ; à 1200² 0,10717 − 0,08439 i contre 0,11616 − 0,09308 i eV⁻¹. Maximum de Re g₀ (jouet) : 0,1787 eV⁻¹
à −2,525 eV (300²), 0,1683 à −2,4725 eV (1200²), soit U_c = 1/max Re g₀ = 5,60 (300²) et 5,94 eV (1200² ; attendu 5,94 à −2,47).

Verdict du critère d'échec (racine à 1200² à plus de 0,02 eV de l'attendu pour U ∈ {8, 12, 20, 27, 40}) : **aucun écart**,
écart maximal 0,0039 eV (U = 12). Écarts sur argmin |det| ≤ 0,005 eV, sur |det| au minimum ≤ 0,005 (U = 1e3 : 6,073 contre 6,1).
Sous U_c (U = 3, 5) : aucune racine, |det| à −2,5 eV = 1,13 / 1,70 (≈ 2 attendu pour un Im g₀ ≈ −0,35 : ici Im g₀ vaut −0,084 à −1 eV
et |det|(−2,5) dépend de U). Pour U ≥ 6,617, seconde racine entre −2,644 et −2,711 eV (van Hove ; non signalée comme écart).
Site retiré : minimum de |det| exactement à E_D (argmin −6e-14 eV), |det| = 6,05 (1200²) / 6,47 (300²) ; aucune valeur propre ne
s'annule : min sur les énergies de |λ_i| = 0,926, 0,927, 1,000 et 6,22 (1200²) ; 0,924 (300²).

À 300², racines multiples rapprochées pour U ≤ 12 (bruit de grille) : U = 6,617 : douze changements de signe entre −2,648 et −1,789 eV
(dont neuf entre −1,990 et −1,789 ; attendu « entre −1,99 et −1,79 ») ; U = 8 : six (trois entre −1,334 et −1,289 ; attendu −1,33 à −1,29) ;
U = 12 : quatre (trois entre −0,682 et −0,645). À 1200², deux racines par U (van Hove + près de E_D).

#### D6.1 — jouet, 1200² (valeurs attendues du prompt entre parenthèses ; écarts = calculé − attendu)

| U (eV) | racines de Re λ (eV) | racine près de E_D (attendu) | écart | argmin \|det\| (attendu) | écart | \|det\| min (attendu) | écart | \|det\| à −2,5 eV | λ à E_D |
|---|---|---|---|---|---|---|---|---|---|
| 3 | aucune | — (aucune) | — | -1.2850 (-1.280) | -0.0050 | 0.713 (0.71) | +0.0026 | 1.129 | +1.0000 +0.0182 i |
| 5 | aucune | — (aucune) | — | -0.9050 (-0.900) | -0.0050 | 0.624 (0.62) | +0.0037 | 1.697 | +1.0000 +0.0303 i |
| 6.617 | -2.644, -1.879 | -1.879 (-1.880) | +0.0011 | -0.7100 (-0.710) | -0.0000 | 0.577 (0.58) | -0.0027 | 2.239 | +1.0000 +0.0401 i |
| 8 | -2.674, -1.310 | -1.310 (-1.310) | +0.0004 | -0.5925 (-0.590) | -0.0025 | 0.549 (0.55) | -0.0014 | 2.726 | +1.0000 +0.0484 i |
| 12 | -2.690, -0.664 | -0.664 (-0.660) | -0.0039 | -0.3900 (-0.390) | -0.0000 | 0.498 (0.5) | -0.0017 | 4.182 | +1.0000 +0.0726 i |
| 20 | -2.699, -0.313 | -0.313 (-0.310) | -0.0034 | -0.2200 (-0.220) | -0.0000 | 0.463 (0.46) | +0.0030 | 7.161 | +1.0000 +0.1211 i |
| 27 | -2.702, -0.209 | -0.209 (-0.210) | +0.0011 | -0.1550 (-0.160) | +0.0050 | 0.461 (0.46) | +0.0009 | 9.789 | +1.0000 +0.1635 i |
| 40 | -2.704, -0.126 | -0.126 (-0.130) | +0.0037 | -0.0975 (-0.100) | +0.0025 | 0.486 (0.49) | -0.0045 | 14.683 | +1.0000 +0.2422 i |
| 80 | -2.707, -0.055 | -0.055 (-0.055) | -0.0003 | -0.0450 (-0.044) | -0.0010 | 0.643 (0.64) | +0.0029 | 29.770 | +1.0000 +0.4843 i |
| 1000 | -2.711, -0.004 | -0.004 (-0.004) | +0.0000 | -0.0025 (-0.004) | +0.0015 | 6.073 (6.1) | -0.0268 | 377.080 | +1.0000 +6.0542 i |

Site retiré (U = 1e3 + −H sur les trois liaisons ; support [44, 69, 73, 74]) à 1200² : argmin \|det\| = -6.4e-14 eV, \|det\| = 6.054 ; min sur les énergies de \|λ_i\| pour les quatre valeurs propres non triviales : 6.2248, 0.9264, 0.9267, 1.0000 ; aucun changement de signe de Re λ_min. À 300² : \|det\| min 6.472 à -6.4e-14 eV, min \|λ\| 0.9241.


#### D6.1 — jouet, 300² (bruit de grille : toutes les racines rapportées)

| U (eV) | racines de Re λ (eV) | argmin \|det\| | \|det\| min | min \|λ\| (position) |
|---|---|---|---|---|
| 3 | aucune | -1.3450 | 0.703 | 0.703 (-1.345) |
| 5 | aucune | -0.8600 | 0.612 | 0.612 (-0.860) |
| 6.617 | -2.648, -2.614, -2.601, -1.990, -1.983, -1.940, -1.921, -1.889, -1.861, -1.839, -1.799, -1.789 | -0.6950 | 0.566 | 0.566 (-0.695) |
| 8 | -2.685, -2.680, -2.657, -1.334, -1.310, -1.289 | -0.5975 | 0.535 | 0.535 (-0.598) |
| 12 | -2.695, -0.682, -0.664, -0.645 | -0.3825 | 0.490 | 0.490 (-0.383) |
| 20 | -2.699, -0.320 | -0.2300 | 0.445 | 0.445 (-0.230) |
| 27 | -2.701, -0.216 | -0.1350 | 0.444 | 0.444 (-0.135) |
| 40 | -2.702, -0.127 | -0.0850 | 0.468 | 0.468 (-0.085) |
| 80 | -2.704, -0.056 | -0.0375 | 0.624 | 0.624 (-0.038) |
| 1000 | -2.705, -0.005 | -0.0075 | 6.402 | 6.402 (-0.008) |



### D6.2 — vrai H(R) de production (U sur la p_z du site lacunaire ; fenêtre ε − E_D ∈ [−3, +1], E_D Wannier −4,2389)

Pas de valeur attendue. Maximum de Re g₀ (vrai H(R)) : 0,1979 eV⁻¹ à −2,0925 eV (300²) et 0,1765 à −2,1875 eV (1200²), soit
U_c = 5,05 (300²) et 5,67 eV (1200²) contre 5,94 pour le jouet. Racine près de E_D à 1200² : jouet −1,879 / vrai −1,563 (U = 6,617),
−1,310 / −1,128 (8), −0,664 / −0,599 (12), −0,313 / −0,297 (20), −0,209 / −0,205 (27), −0,126 / −0,132 (40), −0,055 / −0,067 (80),
−0,004 / −0,019 (1e3) ; seconde racine (van Hove) entre −2,330 et −2,387 eV au lieu de −2,64 à −2,71.
Bloc π du M_loc complet de production (58 × 58, α = 1) : |det|/max ne descend pas sous 0,244 (300², −0,905 eV) / 0,248 (1200², −0,833 eV) ;
min |λ| = 0,775 / 0,790, λ = 0,581 + 0,513 i / 0,627 + 0,481 i ; aucun changement de signe de Re λ_min sur la fenêtre.
Le minimum de production (2,09e-4 à −2,530 eV) n'est pas dans le bloc π (voir D3 : vecteur propre entièrement σ).

#### D6.2 — vrai H(R), 300² (U sur la p_z du site lacunaire ; énergies ε − E_D, E_D Wannier)

g⁽⁰⁾ p_z–p_z sur site à E_D − 1 eV : vrai H(R) +0.11446 -0.09130 i eV⁻¹ ; jouet à −1 eV, même grille : +0.10545 -0.08474 i eV⁻¹. Maximum de Re g⁽⁰⁾ (vrai) : 0.1979 eV⁻¹ à -2.0925 eV (U_c = 5.05 eV) ; jouet : 0.1787 à -2.5250 eV (U_c = 5.60 eV).

| U (eV) | racines de Re λ (eV) | argmin \|det\| | \|det\| min | \|det\| à −2,5 eV | λ à E_D | jouet même grille : racine près de E_D / argmin \|det\| |
|---|---|---|---|---|---|---|
| 3 | aucune | -1.1850 | 0.691 | 1.609 | +1.0126 +0.0216 i | — / -1.3450 |
| 5 | aucune | -0.9050 | 0.603 | 2.259 | +1.0209 +0.0361 i | — / -0.8600 |
| 6.617 | -2.331, -1.669, -1.662, -1.619, -1.598, -1.568, -1.536, -1.516 | -0.6375 | 0.560 | 2.826 | +1.0277 +0.0477 i | -1.789 / -0.6950 |
| 8 | -2.344, -1.161, -1.146, -1.118 | -0.5250 | 0.539 | 3.325 | +1.0335 +0.0577 i | -1.289 / -0.5975 |
| 12 | -2.362, -0.619, -0.612, -0.587 | -0.3775 | 0.499 | 4.803 | +1.0502 +0.0866 i | -0.645 / -0.3825 |
| 20 | -2.372, -0.301 | -0.2200 | 0.478 | 7.813 | +1.0837 +0.1443 i | -0.320 / -0.2300 |
| 27 | -2.375, -0.209 | -0.1700 | 0.489 | 10.467 | +1.1130 +0.1948 i | -0.216 / -0.1350 |
| 40 | -2.379, -0.128 | -0.1200 | 0.541 | 15.409 | +1.1674 +0.2886 i | -0.127 / -0.0850 |
| 80 | -2.384, -0.069 | -0.0675 | 0.785 | 30.647 | +1.3348 +0.5773 i | -0.056 / -0.0375 |
| 1000 | -2.388, -0.020 | -0.0200 | 7.036 | 381.477 | +5.1855 +7.2157 i | -0.005 / -0.0075 |

Bloc π de M_loc complet (58 × 58), 300² : min \|det\|/max = 0.2442 à -0.9050 eV (\|det\| non normalisé 0.7787) ; minima locaux de \|det\|/max : -1.015 (0.2487), -0.955 (0.2473), -0.905 (0.2442), -0.850 (0.2456), -0.797 (0.2443), -0.747 (0.2471) ; min \|λ\| = 0.7753 à -0.9050 eV, λ = +0.5811 +0.5132 i ; minima locaux de \|λ_min\| : -1.015 (0.7897), -0.955 (0.7852), -0.905 (0.7753), -0.850 (0.7798), -0.797 (0.7755), -0.747 (0.7843) ; changements de signe de Re λ_min : aucun.



#### D6.2 — vrai H(R), 1200² (U sur la p_z du site lacunaire ; énergies ε − E_D, E_D Wannier)

g⁽⁰⁾ p_z–p_z sur site à E_D − 1 eV : vrai H(R) +0.11616 -0.09308 i eV⁻¹ ; jouet à −1 eV, même grille : +0.10717 -0.08439 i eV⁻¹. Maximum de Re g⁽⁰⁾ (vrai) : 0.1765 eV⁻¹ à -2.1875 eV (U_c = 5.67 eV) ; jouet : 0.1683 à -2.4725 eV (U_c = 5.94 eV).

| U (eV) | racines de Re λ (eV) | argmin \|det\| | \|det\| min | \|det\| à −2,5 eV | λ à E_D | jouet même grille : racine près de E_D / argmin \|det\| |
|---|---|---|---|---|---|---|
| 3 | aucune | -1.2025 | 0.701 | 1.606 | +1.0126 +0.0205 i | — / -1.2850 |
| 5 | aucune | -0.8450 | 0.615 | 2.247 | +1.0209 +0.0341 i | — / -0.9050 |
| 6.617 | -2.330, -1.563 | -0.6650 | 0.573 | 2.805 | +1.0277 +0.0452 i | -1.879 / -0.7100 |
| 8 | -2.348, -1.128 | -0.5550 | 0.548 | 3.297 | +1.0335 +0.0546 i | -1.310 / -0.5925 |
| 12 | -2.361, -0.599 | -0.3675 | 0.508 | 4.754 | +1.0502 +0.0820 i | -0.664 / -0.3900 |
| 20 | -2.370, -0.297 | -0.2125 | 0.493 | 7.723 | +1.0837 +0.1366 i | -0.313 / -0.2200 |
| 27 | -2.373, -0.205 | -0.1550 | 0.506 | 10.341 | +1.1130 +0.1844 i | -0.209 / -0.1550 |
| 40 | -2.377, -0.132 | -0.1025 | 0.555 | 15.216 | +1.1674 +0.2732 i | -0.126 / -0.0975 |
| 80 | -2.381, -0.067 | -0.0550 | 0.777 | 30.249 | +1.3348 +0.5463 i | -0.055 / -0.0450 |
| 1000 | -2.387, -0.019 | -0.0150 | 7.176 | 376.362 | +5.1854 +6.8293 i | -0.004 / -0.0025 |

Bloc π de M_loc complet (58 × 58), 1200² : min \|det\|/max = 0.2482 à -0.8325 eV (\|det\| non normalisé 0.7937) ; minima locaux de \|det\|/max : -0.832 (0.2482) ; min \|λ\| = 0.7901 à -0.8325 eV, λ = +0.6268 +0.4811 i ; minima locaux de \|λ_min\| : -2.482 (0.9996), -2.295 (0.9997), -1.617 (0.9997), -1.605 (0.9997), -0.832 (0.7901), -0.127 (0.9900) ; changements de signe de Re λ_min : aucun.




Figures : `fig/d6_toy_1200`, `fig/d6_real_1200` (|det| et Re λ non triviale contre ε − E_D, dix U, 1200²), `fig/d6_toy_grid`
(Re λ à 300² et 1200², U = 6,617, 12, 27).

## D1 — Vérité DFT (J3, job 21796854, 1 min 25 ; J6, 21797985 : figure)

Géométrie (lecture des `.save`) : lacune au site s = (13/27, 13/27, 0) ; premiers voisins (1-based, numérotation de `scf.in`) 64, 80, 81
à 1,4237 Å ; seconds voisins 63, 65, 79, 82, 96, 98 à 2,4659 Å ; grille FFT 270 × 270 × 192 ; tous les atomes à z = 0 (9×9) et
|z| ≤ 5e-7 Å (R1). Les wfc à Γ sont en gamma_only (demi-sphère de 385 710 ondes planes).

**E_D de la super-cellule** : quadruplet des états 323–326 de la parfaite 9×9, dégénéré à 1,4e-7 eV : **E_D = −4,23847 eV**.
E_F de la parfaite = −4,21952 eV : E_D − E_F = **−18,96 meV** (attendu « E_D = E_F = −4,2195 » : non retrouvé ; les quatre
états sont à −4,23847). E_D Wannier = −4,23890 : E_D(SC) − E_D(Wannier) = +0,43 meV. La maille grossière 16 bandes donne à K
−4,238470 eV (= quadruplet à 7e-7 eV) ; la maille dense 27×27 donne −4,238895 (= Wannier à 1e-6) : l'écart de 0,43 meV est entre les
deux calculs DFT (grossier / dense), pas dans l'interpolation. Parfaite 9×9 (354 états) contre maille repliée (16 × 81) : écart
1,3e-6 eV dans la fenêtre, 45 meV au maximum sur les 354 états les plus bas (états du haut de la liste).

**Alignement** (Lu 2019 ; atome 161, le plus loin de la lacune, 18,51 Å) : 9×9 non relaxée −24,65 meV (1,0 Å) et −21,59 (0,5 Å) ;
géométries R1 : +20,22 / +92,65 (nspin1), +25,73 / +330,25 (nspin2, potentiel moyen), +15,89 / +317,07 (↑), +35,57 / +343,43 (↓) meV.
Pour R1 l'atome 161 est déplacé de 1,3e-3 (nspin1) et 2,6e-3 Å (nspin2) par rapport à la parfaite et la sphère de 0,5 Å ne contient
pas le même ensemble de points (1045 / 1057 contre 1039) ; la valeur de travail est celle à 1,0 Å. Moyenne 3D ⟨V_d⟩ − ⟨V_p⟩ :
+9,09 (9×9), +9,07 (nspin1), −17,29 (nspin2 moyen ; ↑ −18,61, ↓ −15,97) meV.

**Parité** : max |1 − |⟨σ_h⟩|| ≤ 6e-10 sur tous les états lus (structures planes, z₀ = 0). **Comptage dans la fenêtre**
(ε − E_D ∈ [−3, +1] après alignement) : parfaite 28 états, tous impairs (π), bandes 299–326 ; 9×9 lacune 30 états = 28 impairs + 2 pairs,
bandes 297–326 ; R1 nspin1 29 = 28 + 1 ; R1 nspin2 ↑ 29 = 28 + 1, ↓ 28 = 28 + 0. Aucun état pair (σ) de la parfaite ne tombe dans la
fenêtre (bord de bande σ sous −3 eV).

**États en excès** (w₂ > 0,0812 = 3 × 0,0271 ; parfaite : max 0,060, aucun état au-dessus du seuil) :

- 9×9 lacune non relaxée : paire paire (σ) dégénérée, bandes 324–325, **+0,101 eV** (+0,060 rel. E_F), w₂ = 0,713, w₁ = 0,281 ;
  impairs (π) : bandes 318–319 dégénérées à −1,759 (w₂ 0,114), bande 320 à **−0,737 eV** (w₂ 0,246), bande 326 à +0,269 (0,102).
- R1 nspin1 : un état pair, bande 324, **−0,002 eV** (+0,003 rel. E_F), w₂ 0,667 ; impairs : 317 à −1,820 (0,084), 320 à −0,806 (0,222), 325 à +0,252 (0,100).
- R1 nspin2 ↑ : un état pair, bande 321, **−0,850 eV**, w₂ 0,635 ; projwfc : s + p_x,y de l'atome 81 = 0,641, s + p_x,y des atomes 64 + 80 = 0,106,
  p_z des trois voisins 0,000 (|ψ|² projeté 0,986) ; impairs : 317 (−1,801, 0,097), 319 (−1,743, 0,136), 320 (−0,853, 0,237 ; p_z voisins 0,262), 325 (+0,225, 0,094).
- R1 nspin2 ↓ : aucun état pair dans la fenêtre ; impairs : 318 (−1,743, 0,113), 319 (−1,734, 0,145), 320 (−0,756, 0,218 ; p_z voisins 0,242), 324 (+0,255, 0,116).
- Pour les états π localisés de nspin2, le poids projwfc sur les p_z des seconds voisins est 0,056–0,122 (états à −1,8 / −1,7 eV) et ≤ 0,018 (états à −0,85 et +0,23).

Entrées `projwfc.in` préparées, non lancées : `d1/projwfc_inputs/{defect_9x9_d,defect_9x9_p,vac_9x9_relax_nspin1}/projwfc.in`
(mêmes paramètres que R1b).

#### D1 — alignement (Lu 2019) et fenêtres

| géométrie | atome le plus loin (1-based) / distance (Å) / déplacement vs parfaite (Å) | décalage 0,5 Å (meV) | décalage 1,0 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ 3D (meV) | points dans la sphère (d / p) |
|---|---|---|---|---|---|
| D | 161 / 18.508 / 0.00e+00 | -21.59 | -24.65 | +9.09 | 8611 / 8611 (1,0 Å) ; 1039 / 1039 (0,5 Å) |
| N1 | 161 / 18.509 / 1.29e-03 | +92.65 | +20.22 | +9.07 | 8599 / 8611 (1,0 Å) ; 1045 / 1039 (0,5 Å) |
| N2 | 161 / 18.510 / 2.57e-03 | +330.25 | +25.73 | -17.29 | 8607 / 8611 (1,0 Å) ; 1057 / 1039 (0,5 Å) |
| N2u | 161 / 18.510 / 2.57e-03 | +317.07 | +15.89 | -18.61 | 8607 / 8611 (1,0 Å) ; 1057 / 1039 (0,5 Å) |
| N2d | 161 / 18.510 / 2.57e-03 | +343.43 | +35.57 | -15.97 | 8607 / 8611 (1,0 Å) ; 1057 / 1039 (0,5 Å) |

| géométrie | E_F (eV) | décalage appliqué (1,0 Å, eV) | bandes de la fenêtre | états | pairs (σ) | impairs (π) | max\|1 − \|⟨σ_h⟩\|\| | max\|z − z₀\| (Å) |
|---|---|---|---|---|---|---|---|---|
| 9x9 parfaite | -4.21952 | +0.00000 | 299–326 | 28 | 0 | 28 | 8.5e-15 | 0.0e+00 |
| 9x9 lacune non relaxée | -4.22195 | -0.02465 | 297–326 | 30 | 2 | 28 | 9.6e-12 | 0.0e+00 |
| R1 nspin1 | -4.22394 | +0.02022 | 297–325 | 29 | 1 | 28 | 1.2e-11 | 5.2e-07 |
| R1 nspin2 ↑ | -4.22442 | +0.02573 | 297–325 | 29 | 1 | 28 | 6.1e-10 | 6.0e-08 |
| R1 nspin2 ↓ | -4.22442 | +0.02573 | 297–324 | 28 | 0 | 28 | 1.6e-13 | 6.0e-08 |

Seuil de localisation : w₂ > 0.0812 (3 × moyenne des 28 états de la fenêtre de la parfaite, 0.0271 ; max dans la parfaite 0.0602 ; fraction d'aire du disque 0.0296).

| géométrie | bande | ε − E_D (eV, aligné) | ε − E_F(défaut) (eV) | ⟨σ_h⟩ | w₂ (2 Å) | w₁ (1 Å) | projwfc : p_z 3 voisins | s+p_x,y atome 81 | s+p_x,y atomes 64+80 | p_z 6 seconds voisins | \|ψ\|² projeté |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 9x9 lacune non relaxée | 318 | -1.759 | -1.801 | -1 | 0.114 | 0.008 | — | — | — | — | — |
| 9x9 lacune non relaxée | 319 | -1.759 | -1.801 | -1 | 0.114 | 0.008 | — | — | — | — | — |
| 9x9 lacune non relaxée | 320 | -0.737 | -0.778 | -1 | 0.246 | 0.047 | — | — | — | — | — |
| 9x9 lacune non relaxée | 324 | +0.101 | +0.060 | +1 | 0.713 | 0.281 | — | — | — | — | — |
| 9x9 lacune non relaxée | 325 | +0.101 | +0.060 | +1 | 0.713 | 0.281 | — | — | — | — | — |
| 9x9 lacune non relaxée | 326 | +0.269 | +0.228 | -1 | 0.102 | 0.022 | — | — | — | — | — |
| R1 nspin1 | 317 | -1.820 | -1.814 | -1 | 0.084 | 0.007 | — | — | — | — | — |
| R1 nspin1 | 320 | -0.806 | -0.800 | -1 | 0.222 | 0.047 | — | — | — | — | — |
| R1 nspin1 | 324 | -0.002 | +0.003 | +1 | 0.667 | 0.219 | — | — | — | — | — |
| R1 nspin1 | 325 | +0.252 | +0.257 | -1 | 0.100 | 0.024 | — | — | — | — | — |
| R1 nspin2 ↑ | 317 | -1.801 | -1.789 | -1 | 0.097 | 0.008 | 0.097 | 0.000 | 0.000 | 0.056 | 0.996 |
| R1 nspin2 ↑ | 319 | -1.743 | -1.731 | -1 | 0.136 | 0.015 | 0.120 | 0.000 | 0.000 | 0.122 | 0.996 |
| R1 nspin2 ↑ | 320 | -0.853 | -0.842 | -1 | 0.237 | 0.059 | 0.262 | 0.000 | 0.000 | 0.002 | 0.994 |
| R1 nspin2 ↑ | 321 | -0.850 | -0.839 | +1 | 0.635 | 0.205 | 0.000 | 0.641 | 0.106 | 0.000 | 0.986 |
| R1 nspin2 ↑ | 325 | +0.225 | +0.237 | -1 | 0.094 | 0.026 | 0.105 | 0.000 | 0.000 | 0.014 | 0.991 |
| R1 nspin2 ↓ | 318 | -1.743 | -1.731 | -1 | 0.113 | 0.010 | 0.112 | 0.000 | 0.000 | 0.080 | 0.996 |
| R1 nspin2 ↓ | 319 | -1.734 | -1.722 | -1 | 0.145 | 0.016 | 0.132 | 0.000 | 0.000 | 0.108 | 0.996 |
| R1 nspin2 ↓ | 320 | -0.756 | -0.745 | -1 | 0.218 | 0.053 | 0.242 | 0.000 | 0.000 | 0.006 | 0.994 |
| R1 nspin2 ↓ | 324 | +0.255 | +0.267 | -1 | 0.116 | 0.031 | 0.130 | 0.000 | 0.000 | 0.018 | 0.991 |



Figure : `fig/d1_states` (ε − E_D contre w₂, parité en couleur, seuil en tirets, cinq géométries).

## D4 — Reconstruction de la super-cellule 9×9 non relaxée (J4, jobs 21796855 / 21797745 / 21797985 ; 3 min)

Bases : (a) états de Bloch de maille repliés |nk⟩, k sur la grille 9×9 (81 k, indices 0, 3, 6, … de la 27×27), H = diag(ε) + M/81,
diagonalisé par bloc de parité (parité de chaque |nk⟩ par ⟨σ_h⟩ de la maille) ; (b) |wk⟩ : H = V†εV ⊕ M_W(k,k′)/81 ; (c) |wr⟩ :
H_SC = H(R) replié + M_W(R,R′) replié modulo 9 (étiquettes recentrées, lacune dans la maille r = 0). Poids w₂ des états propres
reconstruits en ondes planes de la super-cellule (même disque de 2 Å et même seuil 0,0812 que D1) ; pour (b), (c) le poids sur les six
fonctions de Wannier « site + voisins » (sp² ×3 et p_z de A en r = 0, p_z de B dans les trois mailles voisines) est aussi donné.
Vérification du chemin des poids : un état de Bloch pur donne w₂ = 0,027–0,030 (fraction d'aire 0,0296).

**Portes.** Porte 1 : 4,9e-14 eV (seuil 1e-8) : OK. Porte 2, littérale, avec le M_W de production : 3,7e-7 eV (seuil 1e-9) :
**échec littéral** ; cause identifiée : les 729 points k du `.save` dense sont écrits avec un bruit de 1,3e-7 sur 27k (8 chiffres
dans le XML ; 7,8e-14 pour la grille 9×9), et la double TF de production (R jusqu'à 13) en hérite. Avec M_W recalculé à k = m/27
exacts, porte 2 = 2,4e-12 eV : OK, le repliement est exact. Effet du bruit sur la production : V_loc (145 × 145) diffère au plus de
1,6e-7 eV (p_z–p_z 6,616636914 contre 6,616636914), M_SC replié de 7,2e-7 eV. Les variantes (c) utilisent le M_W de production tel
quel. Résidu à un corps H(R) contre V†εV sur les 81 k : 7,5e-6 eV (valeurs propres 1,3e-5) ; écart littéral (b) [V†εV] contre
(c-all) [H(R)] : 1,3e-5 eV.

**E_D par base** : maille grossière 16 bandes −4,238470 (bandes 3–4 à K : −4,238470 / −4,238470) = parfaite 9×9 (−4,238471) ;
maille dense 20 bandes −4,238895 = H(R) −4,238896. Parité des états de maille : 20 bandes exacte à 1,6e-6 ; 16 bandes déviation
jusqu'à 6,5e-3 sur 24 états des bandes 15–16 (mélange σ/π de bandes quasi dégénérées) ; couplage résiduel pair–impair de H : (a1) 9,8e-4 eV,
(a2) 1e-9, (a3) 8e-6, (b), (c) ≤ 2,4e-9 eV.

**M grossier (16 bandes, juin) contre M dense restreint aux 81 k et tronqué à 16 bandes** (gauge-invariant seulement par les valeurs
singulières) : valeurs singulières, écart relatif max 1,4e-3 (production, bandes 1–8 : 5,5e-4) ; éléments |M| : jusqu'à 3,45 eV
d'écart (max |M| 7,15), diagonale bandes 1–8 jusqu'à 1,29 eV ; ε 16 bandes contre 20 bandes : 0,92 meV ; la 16ᵉ et la 17ᵉ bande de
la maille dense sont dégénérées (gap minimal 1,4e-7 eV) : le sous-espace à 16 bandes n'est pas unique.

**Résultats par variante** (ε − E_D en eV ; w₂ ; QE = 9×9 lacune de D1, aligné) — même fenêtre de 30 états (2 pairs + 28 impairs) pour
toutes les variantes :

| variante | pairs (σ) localisés | impairs (π) localisés |
|---|---|---|
| QE (D1) | +0,101 ×2 (w₂ 0,713) | −1,759 ×2 (0,114) ; **−0,737 (0,246)** ; +0,269 (0,102) |
| (a1) M_ed 16 bandes, total | **−2,812 ×2 (0,287)** | −1,803 ×2 (0,101) ; **−1,350 (0,226)** |
| (a1) M^L seul | aucun (max w₂ 0,056 à −2,779) | −1,803 ×2 (0,101) ; −1,780 (0,103) |
| (a1) M^NL seul | −2,847 ×2 (0,255) | −1,370 (0,222) |
| (a2) M dense ⊂ 81 k, 16 bandes | aucun (−2,812 ×2 : 0,058 / 0,074) | −2,654 (0,084) ; −1,803 ×2 (0,091 / 0,090) ; −1,350 (0,084) |
| (a3) idem 20 bandes | aucun (−2,848 ×2 : 0,067 / 0,055) | −2,656 (0,083) ; −1,803 ×2 (0,091 / 0,090) ; −1,359 (0,084) |
| (b) 5 WF, V†εV | −2,514 ×2 (0,106 / 0,088 ; poids WF site 0,448) | −2,647 (0,086) ; −1,803 ×2 (0,091 / 0,090) ; −1,315 (0,083 ; poids WF site 0,247) |
| (c-all) H(R) + M_W replié, toutes mailles | −2,514 ×2 (0,106 / 0,088 ; 0,448) | −2,647 (0,086) ; −1,803 ×2 ; −1,315 (0,083 ; 0,247) |
| (c-3) idem, R_cut = 3 | −2,516 / −2,515 (0,106 / 0,088) | −2,649 (0,086) ; −1,803 ×2 ; −1,312 (0,083) |

Marches (états localisés, même parité, correspondance au plus proche) :

- QE → (a1) : pair +0,101 (0,713) → −2,812 (0,287), Δ = −2,91 eV ; impair −0,737 (0,246) → −1,350 (0,226), Δ = −0,61 ; −1,759 (0,114) → −1,803 (0,101), Δ = −0,044 ;
  +0,269 (0,102) → −1,350, Δ = −1,62 (pas d'état localisé au-dessus de −1,3 dans (a1)). Médiane |Δε| des 28 impairs triés : 7,3 meV ; des 2 pairs : 2,91 eV.
- (a1) → (a2) (« remplissage / densité », en fait M grossier de juin contre M dense restreint, même 16 bandes) : énergies inchangées
  (médiane |Δε| 8e-5 pair, 1,4e-5 impair) mais la paire paire à −2,812 passe sous le seuil (w₂ 0,287 → 0,058 / 0,074) et l'impair −1,350
  passe de 0,226 à 0,084.
- (a2) → (a3) (16 → 20 bandes) : impairs identiques à 1e-8 (−1,350 → −1,359 : Δ −0,009) ; pairs : −2,812 → −2,848 (Δ −0,036), toujours sous le seuil.
- (a3) → (b) (20 bandes → 5 WF) : pairs −2,848 → −2,514 (Δ +0,334), w₂ 0,067 / 0,055 → 0,106 / 0,088 (poids WF site 0,448) ; impairs −1,359 → −1,315 (Δ +0,044), −2,656 → −2,647.
- (b) → (c-all) : identiques (Δε ≤ 3e-6 eV, Δw₂ ≤ 1,7e-5 pair, 1,1e-2 impair au sein des paires dégénérées).
- (c-all) → (c-3) (R_cut 3) : Δε = −1,2 meV (pairs), −1,6 / +0,1 / +2,6 meV (impairs) ; w₂ inchangés à 2e-3.

Décalage rigide résiduel entre l'alignement de D1 (Lu, 1,0 Å, + E_D de la parfaite) et l'alignement de chaque modèle sur son E_D,
ajusté sur les 269 états QE sous E_D − 4 eV (médiane) : +1,9 (a1), +2,0 (a2), +2,1 (a3), +1,5 (b), +1,5 (c-all), +1,2 meV (c-3) ;
résidu maximal des paires triées 1,4–2,0 eV (états profonds mal appariés par tri). Décalage Lu de D1 : −24,65 meV.

Le secteur σ de (b) et (c) est rapporté, mais la base à 5 fonctions de Wannier ne contient pas de σ* (bandes sp² anti-liantes) :
les 2 états pairs de la fenêtre y sont les seuls états σ disponibles au-dessus de −3 eV.

#### D4 — portes et références

Porte 1 (H(R) replié contre valeurs de Wannier aux 81 k) : 4.88e-14 eV (seuil 1e-8) → OK. Porte 2 (repliement de M_W, même partie à un corps) : 3.67e-07 eV avec M_W de production (k du XML dense, bruit 1.3e-07 sur 27k), 2.44e-12 eV avec M_W recalculé à k = m/27 exacts (seuil 1e-9) → OK ; V_loc production contre k exacts : max\|diff\| 1.57e-07 eV (p_z–p_z 6.616636914 contre 6.616636914). Résidu H(R) contre V†εV aux 81 k : 7.51e-06 eV (valeurs propres 1.33e-05) ; écart littéral (b) contre (c-all) avec H(R) : 1.26e-05 eV.

E_D par base : maille grossière 16 bandes à K -4.23847 (bandes 3–4 : -4.23847, -4.23847) ; maille dense 20 bandes -4.23890 ; H(R) -4.23890 ; super-cellule parfaite -4.23847 eV. Parité des états de maille : 16 bandes max\|1 − \|⟨σ⟩\|\| = 6.5e-03 (786 pairs, 510 impairs), 20 bandes 1.6e-06 (961 / 659).

M grossier (16 bandes, juin) contre M dense restreint aux 81 k et tronqué à 16 bandes : max\|\|M_c\| − \|M_d\|\| = 3.454 eV (max\|M\| 7.151) ; valeurs singulières : écart relatif max 1.39e-03 ; diagonale, bandes 1–8 : 1.29e+00 eV, toutes bandes 1.542 eV (bande 13, k 33) ; ε 16 bandes contre 20 bandes : 9.19e-04 eV ; gap minimal bandes 16–17 (dense) 1.4e-07 eV ; états de la maille 16 bandes à parité non entière (> 1e-4) : 24.

| variante | dim | couplage résiduel pair–impair (eV) | E_D (eV) | pairs (σ) dans la fenêtre | impairs (π) | états localisés (ε − E_D ; w₂) pairs | impairs | décalage rigide résiduel (meV, états < E_D − 4) / résidu max (eV) |
|---|---|---|---|---|---|---|---|---|
| QE_D1 | — | — | — | 2 | 28 | +0.101 ; 0.713; +0.101 ; 0.713 | -1.759 ; 0.114; -1.759 ; 0.114; -0.737 ; 0.246; +0.269 ; 0.102 | — |
| a1_tot | 1296 | 9.8e-04 | -4.23847 | 2 | 28 | -2.812 ; 0.287; -2.812 ; 0.287 | -1.803 ; 0.101; -1.803 ; 0.101; -1.350 ; 0.226 | +1.9 / 1.93 |
| a1_L | 1296 | 3.3e-05 | -4.23847 | 0 | 28 | aucun | -1.803 ; 0.101; -1.803 ; 0.101; -1.780 ; 0.103 | — |
| a1_NL | 1296 | 9.5e-04 | -4.23847 | 2 | 28 | -2.847 ; 0.255; -2.847 ; 0.255 | -1.370 ; 0.222 | — |
| a2_16 | 1296 | 1.0e-09 | -4.23890 | 2 | 28 | aucun | -2.654 ; 0.084; -1.803 ; 0.091; -1.803 ; 0.090; -1.350 ; 0.084 | +2.0 / 1.93 |
| a3_20 | 1620 | 8.0e-06 | -4.23890 | 2 | 28 | aucun | -2.656 ; 0.083; -1.803 ; 0.091; -1.803 ; 0.090; -1.359 ; 0.084 | +2.1 / 2.00 |
| b_5wf | 405 | 1.5e-09 | -4.23890 | 2 | 28 | -2.514 ; 0.106; -2.514 ; 0.088 | -2.647 ; 0.086; -1.803 ; 0.091; -1.803 ; 0.090; -1.315 ; 0.083 | +1.5 / 1.38 |
| c_all | 405 | 2.4e-09 | -4.23890 | 2 | 28 | -2.514 ; 0.106; -2.514 ; 0.088 | -2.647 ; 0.086; -1.803 ; 0.091; -1.803 ; 0.090; -1.315 ; 0.083 | +1.5 / 1.38 |
| c_3 | 405 | 2.4e-09 | -4.23890 | 2 | 28 | -2.516 ; 0.106; -2.515 ; 0.088 | -2.649 ; 0.086; -1.803 ; 0.089; -1.803 ; 0.089; -1.312 ; 0.083 | +1.2 / 1.38 |

Marches (états localisés de l'étape i → état localisé le plus proche de l'étape i + 1, même parité) :

| marche | parité | états dans la fenêtre (i → i+1) | médiane \|Δε\| des spectres triés (eV) | correspondances (ε_i ; w₂,i → ε_{i+1} ; w₂,i+1 ; Δε) |
|---|---|---|---|---|
| QE_D1->a1_tot | even | 2 → 2 | +2.9134 | +0.101 ; 0.713 → -2.812 ; 0.287 ; -2.913 (> 0,1); +0.101 ; 0.713 → -2.812 ; 0.287 ; -2.913 (> 0,1) |
| QE_D1->a1_tot | odd | 28 → 28 | +0.0073 | -1.759 ; 0.114 → -1.803 ; 0.101 ; -0.044; -1.759 ; 0.114 → -1.803 ; 0.101 ; -0.044; -0.737 ; 0.246 → -1.350 ; 0.226 ; -0.612 (> 0,1); +0.269 ; 0.102 → -1.350 ; 0.226 ; -1.619 (> 0,1) |
| a1_tot->a2_16 | even | 2 → 2 | +0.0001 | -2.812 ; 0.287 → absent; -2.812 ; 0.287 → absent |
| a1_tot->a2_16 | odd | 28 → 28 | +0.0000 | -1.803 ; 0.101 → -1.803 ; 0.091 ; +0.000; -1.803 ; 0.101 → -1.803 ; 0.091 ; +0.000; -1.350 ; 0.226 → -1.350 ; 0.084 ; +0.000 |
| a2_16->a3_20 | even | 2 → 2 | +0.0365 | aucun état localisé |
| a2_16->a3_20 | odd | 28 → 28 | +0.0000 | -2.654 ; 0.084 → -2.656 ; 0.083 ; -0.002; -1.803 ; 0.091 → -1.803 ; 0.090 ; -0.000; -1.803 ; 0.090 → -1.803 ; 0.090 ; -0.000; -1.350 ; 0.084 → -1.359 ; 0.084 ; -0.009 |
| a3_20->b_5wf | even | 2 → 2 | +0.3344 | aucun état localisé |
| a3_20->b_5wf | odd | 28 → 28 | +0.0000 | -2.656 ; 0.083 → -2.647 ; 0.086 ; +0.009; -1.803 ; 0.091 → -1.803 ; 0.091 ; +0.000; -1.803 ; 0.090 → -1.803 ; 0.091 ; +0.000; -1.359 ; 0.084 → -1.315 ; 0.083 ; +0.044 |
| c_all->c_3 | even | 2 → 2 | — | -2.514 ; 0.106 → -2.515 ; 0.088 ; -0.001; -2.514 ; 0.088 → -2.515 ; 0.088 ; -0.001 |
| c_all->c_3 | odd | 28 → 28 | — | -2.647 ; 0.086 → -2.649 ; 0.086 ; -0.002; -1.803 ; 0.091 → -1.803 ; 0.089 ; +0.000; -1.803 ; 0.090 → -1.803 ; 0.089 ; +0.000; -1.315 ; 0.083 → -1.312 ; 0.083 ; +0.003 |
| b_5wf->c_all | even | — → — | — | -2.514 ; 0.106 → -2.514 ; 0.106 ; +0.000; -2.514 ; 0.088 → -2.514 ; 0.106 ; +0.000 |
| b_5wf->c_all | odd | — → — | — | -2.647 ; 0.086 → -2.647 ; 0.086 ; +0.000; -1.803 ; 0.091 → -1.803 ; 0.090 ; -0.000; -1.803 ; 0.090 → -1.803 ; 0.090 ; -0.000; -1.315 ; 0.083 → -1.315 ; 0.083 ; +0.000 |

Tous les états de la fenêtre par variante et parité (ε − E_D ; w₂) :

- QE_D1, even : +0.101 (0.713), +0.101 (0.713)
- QE_D1, odd : -2.793 (0.009), -2.793 (0.009), -2.790 (0.001), -2.777 (0.043), -2.777 (0.043), -2.535 (0.058), -2.105 (0.001), -2.105 (0.001), -2.103 (0.003), -2.086 (0.043), -2.086 (0.043), -2.001 (0.010), -1.815 (0.017), -1.811 (0.000), -1.811 (0.000), -1.808 (0.000), -1.807 (0.003), -1.807 (0.003), -1.798 (0.005), -1.798 (0.005), -1.798 (0.000), -1.759 (0.114), -1.759 (0.114), -0.737 (0.246), +0.003 (0.011), +0.015 (0.033), +0.015 (0.033), +0.269 (0.102)
- a1_tot, even : -2.812 (0.287), -2.812 (0.287)
- a1_tot, odd : -2.790 (0.009), -2.790 (0.009), -2.790 (0.001), -2.790 (0.046), -2.790 (0.046), -2.654 (0.067), -2.100 (0.001), -2.100 (0.001), -2.100 (0.003), -2.099 (0.053), -2.099 (0.053), -2.021 (0.019), -1.804 (0.018), -1.804 (0.000), -1.804 (0.000), -1.804 (0.000), -1.804 (0.002), -1.804 (0.002), -1.804 (0.005), -1.804 (0.005), -1.804 (0.000), -1.803 (0.101), -1.803 (0.101), -1.350 (0.226), -0.000 (0.011), -0.000 (0.034), -0.000 (0.034), +0.108 (0.036)
- a1_L, even : 
- a1_L, odd : -2.790 (0.009), -2.790 (0.009), -2.790 (0.001), -2.790 (0.046), -2.790 (0.046), -2.779 (0.056), -2.100 (0.001), -2.100 (0.001), -2.100 (0.003), -2.099 (0.053), -2.099 (0.053), -2.089 (0.047), -1.804 (0.018), -1.804 (0.000), -1.804 (0.000), -1.804 (0.000), -1.804 (0.002), -1.804 (0.002), -1.804 (0.005), -1.804 (0.005), -1.804 (0.000), -1.803 (0.101), -1.803 (0.101), -1.780 (0.103), -0.000 (0.011), -0.000 (0.034), -0.000 (0.034), +0.007 (0.028)
- a1_NL, even : -2.847 (0.255), -2.847 (0.255)
- a1_NL, odd : -2.790 (0.011), -2.790 (0.038), -2.790 (0.011), -2.790 (0.038), -2.790 (0.011), -2.659 (0.067), -2.099 (0.004), -2.099 (0.050), -2.099 (0.004), -2.099 (0.050), -2.099 (0.004), -2.022 (0.019), -1.803 (0.023), -1.803 (0.008), -1.803 (0.023), -1.803 (0.027), -1.803 (0.023), -1.803 (0.030), -1.803 (0.023), -1.803 (0.026), -1.803 (0.023), -1.803 (0.008), -1.803 (0.023), -1.370 (0.222), -0.000 (0.021), +0.000 (0.034), +0.000 (0.024), +0.105 (0.035)
- a2_16, even : -2.812 (0.058), -2.812 (0.074)
- a2_16, odd : -2.790 (0.036), -2.790 (0.028), -2.790 (0.028), -2.790 (0.025), -2.790 (0.019), -2.654 (0.084), -2.100 (0.013), -2.100 (0.017), -2.100 (0.023), -2.099 (0.040), -2.099 (0.039), -2.021 (0.007), -1.804 (0.046), -1.804 (0.007), -1.804 (0.006), -1.804 (0.000), -1.804 (0.003), -1.804 (0.003), -1.804 (0.008), -1.804 (0.008), -1.804 (0.001), -1.803 (0.091), -1.803 (0.090), -1.350 (0.084), -0.000 (0.017), -0.000 (0.034), -0.000 (0.035), +0.108 (0.023)
- a3_20, even : -2.848 (0.067), -2.848 (0.055)
- a3_20, odd : -2.790 (0.036), -2.790 (0.028), -2.790 (0.028), -2.790 (0.025), -2.790 (0.019), -2.656 (0.083), -2.100 (0.013), -2.100 (0.017), -2.100 (0.023), -2.099 (0.040), -2.099 (0.039), -2.022 (0.007), -1.804 (0.046), -1.804 (0.007), -1.804 (0.006), -1.804 (0.000), -1.804 (0.003), -1.804 (0.003), -1.804 (0.008), -1.804 (0.008), -1.804 (0.001), -1.803 (0.091), -1.803 (0.090), -1.359 (0.084), -0.000 (0.017), -0.000 (0.034), -0.000 (0.035), +0.106 (0.023)
- b_5wf, even : -2.514 (0.106), -2.514 (0.088)
- b_5wf, odd : -2.790 (0.036), -2.790 (0.028), -2.790 (0.028), -2.790 (0.025), -2.790 (0.019), -2.647 (0.086), -2.100 (0.013), -2.100 (0.017), -2.100 (0.023), -2.099 (0.040), -2.099 (0.039), -2.019 (0.007), -1.804 (0.046), -1.804 (0.007), -1.804 (0.006), -1.804 (0.000), -1.804 (0.003), -1.804 (0.003), -1.804 (0.008), -1.804 (0.008), -1.804 (0.001), -1.803 (0.091), -1.803 (0.090), -1.315 (0.083), -0.000 (0.017), -0.000 (0.034), -0.000 (0.035), +0.117 (0.023)
- c_all, even : -2.514 (0.106), -2.514 (0.088) ; poids Wannier site+voisins : 0.448, 0.448
- c_all, odd : -2.790 (0.025), -2.790 (0.039), -2.790 (0.028), -2.790 (0.028), -2.790 (0.015), -2.647 (0.086), -2.100 (0.017), -2.100 (0.013), -2.100 (0.023), -2.099 (0.039), -2.099 (0.040), -2.019 (0.007), -1.804 (0.046), -1.804 (0.007), -1.804 (0.006), -1.804 (0.000), -1.804 (0.003), -1.804 (0.003), -1.804 (0.008), -1.804 (0.008), -1.804 (0.001), -1.803 (0.091), -1.803 (0.090), -1.315 (0.083), -0.000 (0.017), -0.000 (0.034), -0.000 (0.035), +0.117 (0.023) ; poids Wannier site+voisins : 0.003, 0.003, 0.000, 0.045, 0.045, 0.071, 0.000, 0.000, 0.000, 0.050, 0.050, 0.019, 0.000, 0.000, 0.000, 0.000, 0.001, 0.001, 0.000, 0.000, 0.000, 0.105, 0.105, 0.247, 0.000, 0.037, 0.037, 0.039
- c_3, even : -2.516 (0.106), -2.515 (0.088) ; poids Wannier site+voisins : 0.448, 0.448
- c_3, odd : -2.790 (0.038), -2.790 (0.026), -2.790 (0.026), -2.790 (0.030), -2.790 (0.016), -2.649 (0.086), -2.099 (0.022), -2.099 (0.013), -2.099 (0.018), -2.099 (0.039), -2.099 (0.040), -2.020 (0.007), -1.804 (0.046), -1.804 (0.007), -1.804 (0.008), -1.803 (0.010), -1.803 (0.006), -1.803 (0.000), -1.803 (0.003), -1.803 (0.003), -1.803 (0.002), -1.803 (0.089), -1.803 (0.089), -1.312 (0.083), -0.000 (0.017), +0.000 (0.034), +0.000 (0.035), +0.117 (0.023) ; poids Wannier site+voisins : 0.011, 0.003, 0.002, 0.043, 0.037, 0.071, 0.000, 0.001, 0.000, 0.050, 0.050, 0.019, 0.000, 0.002, 0.000, 0.000, 0.000, 0.000, 0.001, 0.000, 0.000, 0.104, 0.105, 0.247, 0.000, 0.037, 0.037, 0.039



Figure : `fig/d4_ladder` (niveaux de la fenêtre par variante et parité, taille ∝ w₂).

## D2 — Éléments de M (9×9 de production, M dense 27×27 → base de Wannier ; J1, job 21796852, 3 min 56)

Vérifications préalables du GO (J1) :

| vérification | valeur | règle | issue |
|---|---|---|---|
| max\|M_L_9x9 (juin) − M_L_dense_9x9_coarsecheck\| / max\|M_L\| sur (16, 81, 16, 81) | 1,3e-15 (max\|diff\| 1,1e-17 Ha ; max\|M_L\| 8,55e-3 Ha) ; moyenne diagonale 76,58 meV pour les deux | ≤ 2e-3 | `M_ed_9x9.npy` (juin) sert pour D4 (a1) ; L = `M_L_9x9.npy`, NL = M_ed − M_L |
| M^NL grossier de juin (M_ed − M_L) contre M^NL recalculé avec la chaîne actuelle (`compute_M_NL`, 172 s) | écart relatif 1,1e-16 (max\|M^NL\| 0,270 Ha) | — | chaîne NL inchangée depuis juin ; fichier `prep/M_NL_coarse_new_9x9.npy` |
| ‖M_W,tot − M_W,L − M_W,NL‖ / ‖M_W,tot‖ (linéarité de la rotation) | 1,4e-15 | ≤ 1e-12 | OK |
| recentrage | R_d = [4, 4, 0] pour tot, L et NL ; ‖M_W(R, 0)‖ : 9,365 (R = 0), 0,658 / 0,437 (1ʳᵉ couronne), 0,460 / 0,165 / 0,035 (2ᵉ) eV | production identique | OK |
| max\|hors bloc σ–π\| de M_loc (145 × 145, R_cut 3) | tot 2,07e-9 eV (max\|M_loc\| 6,617 ; ratio 3,1e-10) ; L 2,07e-9 (ratio 6,7e-9) ; NL 1,3e-10 (ratio 2,1e-11) ; résidu hermitien ≤ 2,4e-15 | > 1e-3 × max ⇒ D3 approximatif | D3 par bloc exact |
| rotation d'un fichier dense de 3,4 G | 7 s par fichier (V†MV + double TF), MaxRSS 17,9 G | — | — |

### Base de Wannier (eV ; R du voisin B relatif à la maille de la lacune ; H(R) de `wannier_tb.dat`)

| élément | total | local (M^L) | non local (M^NL) | H(R) |
|---|---|---|---|---|
| p_z–p_z du site lacunaire (R = 0) | **6,6166** | 0,3113 | 6,3053 | sur site −3,9125 |
| p_z(lacune)–p_z(B) pour R = (0,0,0), (−1,0,0), (0,−1,0) | −0,4910 (×3, égaux à 2e-9) | −0,0112 | −0,4797 | t = −2,9089 (×3) |
| diagonale p_z des trois voisins B | 0,04428 (×3) | 0,00713 | 0,03714 | — |
| diagonale des trois sp² du carbone retiré (R = 0) | 0,6310, 0,6310, 0,6310 | 0,1882 (×3) | 0,4428, 0,4428, 0,4427 | sur site −15,0847 |
| sp²–sp² du carbone retiré (R = 0) | −2,6533 | +0,0117 | −2,6650 | −2,1542 |
| ‖bloc sur site 5 × 5‖ | 9,3647 | 0,4520 | 9,1335 | — |
| parties imaginaires | ≤ 3e-11 | | | |

Comparaisons demandées : t (premiers voisins, H(R)) = −2,909 eV ; seuil U_c de D6.1 = 5,94 eV ; Kaasbjerg PRB 101, 045433
(2020), p. 5 : ~70 eV Å² → 70 / A_cell = 13,29 eV avec A_cell = 5,266 Å² (a = 2,4659 Å), 70 / (A_cell/2) = 26,59 eV
(la valeur « ≈ +27 eV » du prompt correspond à l'aire par atome ; avec 5,24 Å² : 13,36 / 26,72 eV).

### Base de Bloch (27×27 dense, 20 bandes, norme unit_cell × A_cell = convention intensive du mémoire, eV Å²)

Paire π/π* dégénérée : bandes 3–4 (0-based) à K (indice 495, k = (2/3, 1/3), poids p_z 1,000) et à K′ (indice 261,
(1/3, 2/3)) ; ε = −4,238895 eV aux deux points (écart 1,8e-7 eV). Demi-trace M̄ = ½ Tr du bloc 2 × 2, norme = Frobenius.

| bloc | partie | M̄ (eV) | M̄ (eV Å²) | ‖bloc‖ (eV) | ‖bloc‖ (eV Å²) |
|---|---|---|---|---|---|
| intravallée (K, K) | total | 3,2426 | 17,07 | 6,489 | 34,17 |
| | M^L | 0,1396 | 0,735 | 0,283 | 1,49 |
| | M^NL | 3,1030 | 16,34 | 6,206 | 32,68 |
| intervallée (K′, K) | total | 3,1405 + 0,8034 i | 16,54 + 4,23 i | 6,511 | 34,29 |
| | M^L | 0,1471 + 0,0376 i | 0,775 + 0,198 i | 0,305 | 1,61 |
| | M^NL | 2,9934 + 0,7658 i | 15,76 + 4,03 i | 6,206 | 32,68 |

Moyenne ⟨ΔV^L⟩ de production : la construction de M dense n'en soustrait aucune (`compute_M_dense_stages.py`, `subtract_mean=False`) ;
le contrôle C14 soustrait la moyenne diagonale de M^L, 67,00 meV (`resonance_9x9.npz`, clé `ML_diag_mean`) ; produit
⟨ΔV^L⟩ · A_sc = 0,06700 × 426,53 = 28,58 eV Å² (A_sc = 81 × 5,266). Moyenne diagonale de M^L grossier 16 bandes (J1) : 76,58 meV.

## D3 — Critère de pôle avec M_loc → α M_loc, par bloc de parité (J5a g⁽⁰⁾ 600² : job 21797930, 2 min 46 ; J5b : job 21797931, 7 min 32)

Première soumission (21796856) annulée après 20 min : la diagonalisation par lots de 2 401 matrices 145 × 145 avec 16 fils BLAS
tournait 13 fois trop lentement (61 ms contre 4,7 ms par matrice avec 1 fil BLAS et 16 fils Python) ; correctif dans le module
nouveau `pole_criterion.py` (`_batched`), résultats identiques à 1e-14, aucune routine de production touchée.

**Porte de régression** (α = 1, matrice complète, 300², fenêtre ±3 eV, mêmes lignes que `resonance_criteria.py`) : min |det|/max =
**2,0908e-4 à −2,530 eV**, λ_min = **+0,00029 + 0,01078 i** (production : 2,091e-4 à −2,530 ; +0,0003 + 0,0108 i) ; minima secondaires
−2,185 (2,61e-2), −2,132, −2,082, −2,027, −1,972 eV (4,05e-2) identiques au log de production : **OK**.

**Le minimum de production à −2,530 eV** : vecteur propre de λ_min à **100 % dans le bloc σ** (poids π 4e-20) ; 0,976 sur les trois sp²
de l'atome retiré (R = 0), 0,024 sur les autres sp², 0 sur les p_z. Il est **doublement dégénéré** : à −2,530 eV les deux plus petites valeurs
propres sont +0,0003 + 0,0108 i et −0,0001 + 0,0108 i, de vecteurs (poids sur les trois sp²) [0,488, 0,488, 0] avec phases (−1, +1, 0) et
[0,163, 0,163, 0,650] — les deux combinaisons orthogonales à la combinaison symétrique des trois sp² (doublet E) ; la troisième valeur
propre σ reste à 0,9995. Les « sauts de branche » signalés à ces racines (recouvrement 0,0 ou 0,2) sont ceux d'un doublet : le vecteur de
plus petit module bascule entre les deux membres dégénérés ; le seul recouvrement 1,0 (racine −2,531 à α = 1) correspond au même membre
suivi. Aucune valeur propre du bloc π n'approche zéro sur la fenêtre : min |λ| = 0,775 (300²) / 0,789 (600²) à −0,905 / −0,830 eV.

**Balayage α** (fenêtre [−3, +1], |det| normalisé sur la fenêtre ; 300² et 600² donnent les mêmes positions à ≤ 0,08 eV pour les
minima σ et ≤ 0,1 eV pour les extrema π) :

| α | matrice complète : min \|det\|/max (position) ; min \|λ\| (position), bloc | bloc σ : min \|λ\| (position) | bloc π : min \|det\|/max (position) ; min \|λ\| (position) | pic de −Im T̄(K) 300² / 600² (max, eV) | Re T̄(E_D) |
|---|---|---|---|---|---|
| 0,5 | 7,6e-2 (−3,0, bord) ; 0,215 (−3,0), σ | 0,215 (−3,0, bord) | 0,44 (−1,29) ; 0,81 (−1,29) | −2,035 / −1,957 (1,3) | 1,423 |
| 1 | 3,8e-4 (−2,530) ; 0,0108 (−2,530), σ | 0,0108 (−2,530) | 0,24 (−0,905) ; 0,78 (−0,905) | −1,292 / −1,242 (3,2) | 2,521 |
| 1,5 | 5,3e-6 (−1,380) ; 0,0051 (−1,380), σ | 0,0051 (−1,380) | 0,17 (−0,635) ; 0,78 (−0,635) | −0,905 / −0,882 (5,0) | 3,396 |
| 2 | 8,9e-7 (+0,013) ; 0,0034 (+0,013), σ | 0,0034 (+0,013) | 0,14 (−0,477) ; 0,82 (−0,477) | −0,685 / −0,720 (6,7) | 4,110 |
| 3 | 3,6e-3 (−0,317) ; 0,286 (+1,0, bord), σ | 0,286 (+1,0, bord) | 0,10 (−0,427) ; 0,90 (−0,427) | −0,525 / −0,535 (9,3) | 5,210 |
| 4 | 5,1e-3 (−0,317) ; 0,715 (+1,0), σ | 0,715 (+1,0) | 0,086 (−0,377) ; 0,96 (−0,125) | −0,427 / −0,435 (11,5) | 6,023 |
| 6 | 5,9e-3 (−0,265) ; 0,938 (−0,125), **π** (p_z lacune 0,99) | 0,996 (−3,0) | 0,070 (−0,267) ; 0,94 (−0,125) | −0,375 / −0,352 (14,1) | 7,162 |
| 10 | 5,7e-3 (−0,220) ; 0,898 (−0,125), π (p_z lacune 0,99) | 0,994 (−3,0) | 0,056 (−0,222) ; 0,90 (−0,125) | −0,267 / −0,280 (17,6) | 8,527 |

Racines de Re λ_min (matrice complète, 300²) : α = 1 : −2,531 ; 1,5 : −2,901 (saut), −1,379 ; 2 : −2,529 (saut), +0,013 ; 3 : −1,376 (saut,
λ ≈ −0,998 : passage par Re λ = −1, pas un zéro) ; 4 : +0,064 (saut) ; 6, 10 : aucune. Le doublet σ de plus petit module est à λ ≈ 0,0108 (α = 1),
0,0051 (1,5), 0,0034 (2) puis sort de la fenêtre par le haut pour α ≥ 3 (0,286 à +1,0 eV pour α = 3). Pour α ≥ 2, une valeur propre portée à 0,98–0,99
par la p_z de la lacune apparaît juste au-dessus de E_D (+0,01 à +0,02 eV) ; c'est elle qui devient λ_min de la matrice complète pour α ≥ 6, à −0,125 eV.

**Variantes** (300² / 600² identiques) :
- locale seule × α (V = αV_L + V_NL) : le doublet σ se déplace lentement, −2,575 (α 0,5), −2,530 (1), −2,482 (1,5), −2,432 (2), −2,330 (3),
  −2,222 (4), −1,990 (6), −1,482 (10) eV, |λ|min 0,0115 → 0,0053 ; pic de −Im T̄ −1,29 → −0,905 ; Re T̄(E_D) 2,47 → 3,49.
- non locale seule × α (V = V_L + αV_NL) : −2,530 (1), −1,445 (1,5), −0,137 (2), puis hors fenêtre par le haut (α = 3 : 0,242 à +1,0 ;
  4 : 0,649) ; pic de −Im T̄ −1,29 → −0,31 ; Re T̄(E_D) 2,52 → 8,04. Pour α_NL = 4 (300²) et 6, 10, λ_min de la matrice complète est de
  nouveau la valeur propre p_z de la lacune (0,98–0,99 à −0,12 eV).

T̄(K) : la production (`resonance_9x9.npz`) donne le pic de −Im T̄ à −1,2925 et Re T̄(E_D) = 2,5213 ; ici à α = 1 : −1,292 (300²) et 2,521.

#### D3 — porte de régression et balayage α

Porte (α = 1, matrice complète, 300², fenêtre ±3 eV) : min \|det\|/max = 2.0908e-04 à -2.530 eV (attendu 2,091e-4 à −2,530) ; min \|λ\| = 0.0108, λ = +0.0003 +0.0108 i à -2.530 eV (attendu +0,0003 + 0,0108 i) → OK. Minima locaux de \|det\|/max : -2.530 (2.091e-04, min|λ| 0.0108), -2.185 (2.609e-02, min|λ| 0.1446), -2.132 (2.901e-02, min|λ| 0.1615), -2.082 (3.262e-02, min|λ| 0.1768), -2.027 (3.658e-02, min|λ| 0.1926), -1.972 (4.051e-02, min|λ| 0.2077).

Vecteur propre de λ_min à −2,530 eV : poids bloc π 0.0000, bloc σ 1.0000 ; groupes : pz_vac 0.0000, pz_nn 0.0000, sp2_A0 0.9755, pz_other 0.0000, sp2_other 0.0245 ; λ = +0.0003 +0.0108 i.


##### Balayage α, 300² (fenêtre [−3, +1] eV ; \|det\| normalisé sur la fenêtre)

| variante | α | bloc | min \|det\|/max (position) | minima locaux de \|det\|/max | min \|λ\| (position) ; λ | poids de v_min : p_z lacune / p_z 3 voisins / sp² A(0) / autres p_z / autres sp² | bloc π / σ | racines de Re λ_min (recouvrement ; saut) | pic de −Im T̄(K) (max) | T̄(E_D) |
|---|---|---|---|---|---|---|---|---|---|---|
| tot | 0.5 | pi | 4.382e-01 (-1.292) | -1.407 (4.40e-01), -1.350 (4.39e-01), -1.292 (4.38e-01), -1.240 (4.39e-01), -1.185 (4.38e-01), -1.127 (4.39e-01) | 0.8115 (-1.292) ; +0.7154 +0.3829 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -2.035 (1.353 eV) | +1.4227 -0.0289 i |
| tot | 0.5 | sigma | 8.170e-02 (-3.000) |  | 0.2153 (-3.000) ; +0.2113 +0.0414 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | — | — |
| tot | 0.5 | full | 7.610e-02 (-3.000) | -2.085 (3.48e-01), -2.030 (3.46e-01), -1.975 (3.47e-01), -1.920 (3.47e-01), -1.862 (3.48e-01), -1.805 (3.49e-01) | 0.2153 (-3.000) ; +0.2113 +0.0414 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | -2.035 (1.353 eV) | +1.4227 -0.0289 i |
| tot | 1.0 | pi | 2.442e-01 (-0.905) | -1.015 (2.49e-01), -0.955 (2.47e-01), -0.905 (2.44e-01), -0.850 (2.46e-01), -0.797 (2.44e-01), -0.747 (2.47e-01) | 0.7753 (-0.905) ; +0.5811 +0.5132 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| tot | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| tot | 1.0 | full | 3.831e-04 (-2.530) | -2.530 (3.83e-04), -2.185 (4.78e-02), -2.132 (5.32e-02), -2.082 (5.98e-02), -2.027 (6.70e-02), -1.972 (7.42e-02) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| tot | 1.5 | pi | 1.707e-01 (-0.635) | -0.792 (1.76e-01), -0.742 (1.75e-01), -0.682 (1.71e-01), -0.635 (1.71e-01), -0.587 (1.73e-01), -0.527 (1.76e-01) | 0.7845 (-0.635) ; +0.5742 +0.5345 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.905 (5.148 eV) | +3.3961 -0.1649 i |
| tot | 1.5 | sigma | 1.254e-05 (-1.380) | -1.380 (1.25e-05) | 0.0051 (-1.380) ; -0.0001 +0.0051 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.901 (0.24; saut), -1.379 (0.00; saut) | — | — |
| tot | 1.5 | full | 5.271e-06 (-1.380) | -1.380 (5.27e-06), -1.022 (9.90e-04), -0.907 (1.46e-03), -0.800 (1.97e-03), -0.752 (2.23e-03), -0.642 (2.84e-03) | 0.0051 (-1.380) ; -0.0001 +0.0051 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.901 (0.24; saut), -1.379 (0.00; saut) | -0.905 (5.148 eV) | +3.3961 -0.1649 i |
| tot | 2.0 | pi | 1.357e-01 (-0.477) | -0.680 (1.42e-01), -0.632 (1.38e-01), -0.582 (1.37e-01), -0.525 (1.36e-01), -0.477 (1.36e-01), -0.430 (1.39e-01) | 0.8197 (-0.477) ; +0.6172 +0.5395 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.685 (6.784 eV) | +4.1103 -0.2418 i |
| tot | 2.0 | sigma | 2.111e-06 (+0.013) | +0.013 (2.11e-06) | 0.0034 (+0.013) ; -0.0000 +0.0034 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.529 (0.21; saut), +0.013 (0.00; saut) | — | — |
| tot | 2.0 | full | 8.864e-07 (+0.013) | -1.552 (1.54e-02), -1.495 (1.31e-02), -1.435 (1.11e-02), -1.215 (5.54e-03), -1.107 (3.87e-03), +0.013 (8.86e-07) | 0.0034 (+0.013) ; -0.0000 +0.0034 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.529 (0.21; saut), +0.013 (0.00; saut) | -0.685 (6.784 eV) | +4.1103 -0.2418 i |
| tot | 3.0 | pi | 1.003e-01 (-0.427) | -0.520 (1.06e-01), -0.475 (1.02e-01), -0.427 (1.00e-01), -0.382 (1.02e-01), -0.337 (1.06e-01), -0.287 (1.12e-01) | 0.8987 (-0.427) ; +0.5335 +0.7233 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.525 (9.470 eV) | +5.2104 -0.3892 i |
| tot | 3.0 | sigma | 4.482e-03 (+1.000) |  | 0.2862 (+1.000) ; -0.2862 +0.0037 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.376 (0.16; saut) | — | — |
| tot | 3.0 | full | 3.599e-03 (-0.317) | -0.425 (3.75e-03), -0.377 (3.64e-03), -0.332 (3.60e-03), -0.317 (3.60e-03), -0.270 (3.60e-03), -0.230 (3.67e-03) | 0.2862 (+1.000) ; -0.2862 +0.0037 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.376 (0.16; saut) | -0.525 (9.470 eV) | +5.2104 -0.3892 i |
| tot | 4.0 | pi | 8.547e-02 (-0.377) | -0.470 (9.14e-02), -0.425 (8.69e-02), -0.377 (8.55e-02), -0.335 (8.62e-02), -0.285 (8.96e-02), -0.232 (9.56e-02) | 0.9586 (-0.125) ; +0.9583 +0.0273 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.427 (11.649 eV) | +6.0234 -0.5210 i |
| tot | 4.0 | sigma | 1.332e-02 (+1.000) |  | 0.7149 (+1.000) ; -0.7149 +0.0049 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | +0.016 (0.13; saut) | — | — |
| tot | 4.0 | full | 5.134e-03 (-0.317) | -0.422 (5.56e-03), -0.375 (5.28e-03), -0.317 (5.13e-03), -0.270 (5.14e-03), -0.230 (5.29e-03), -0.182 (5.60e-03) | 0.7149 (+1.000) ; -0.7149 +0.0049 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | +0.064 (0.00; saut) | -0.427 (11.649 eV) | +6.0234 -0.5210 i |
| tot | 6.0 | pi | 6.957e-02 (-0.267) | -0.420 (8.03e-02), -0.370 (7.46e-02), -0.312 (7.04e-02), -0.267 (6.96e-02), -0.227 (7.24e-02), -0.182 (7.92e-02) | 0.9383 (-0.125) ; +0.9374 +0.0409 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.375 (14.155 eV) | +7.1619 -0.7386 i |
| tot | 6.0 | sigma | 2.406e-02 (+1.000) |  | 0.9964 (-3.000) ; +0.9964 +0.0000 i | 0.000 / 0.000 / 0.041 / 0.000 / 0.959 | 0.000 / 1.000 | aucune | — | — |
| tot | 6.0 | full | 5.846e-03 (-0.265) | -0.370 (6.68e-03), -0.312 (6.09e-03), -0.265 (5.85e-03), -0.227 (5.95e-03), -0.182 (6.33e-03), -0.137 (7.04e-03) | 0.9383 (-0.125) ; +0.9374 +0.0409 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.375 (14.155 eV) | +7.1619 -0.7386 i |
| tot | 10.0 | pi | 5.604e-02 (-0.222) | -0.365 (7.24e-02), -0.307 (6.28e-02), -0.262 (5.70e-02), -0.222 (5.60e-02), -0.180 (5.92e-02), -0.135 (6.81e-02) | 0.8982 (-0.125) ; +0.8956 +0.0682 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.267 (17.738 eV) | +8.5269 -1.0508 i |
| tot | 10.0 | sigma | 3.228e-02 (+1.000) |  | 0.9941 (-3.000) ; +0.9941 +0.0001 i | 0.000 / 0.000 / 0.041 / 0.000 / 0.959 | 0.000 / 1.000 | aucune | — | — |
| tot | 10.0 | full | 5.722e-03 (-0.220) | -0.465 (1.08e-02), -0.417 (9.28e-03), -0.262 (5.94e-03), -0.220 (5.72e-03), -0.180 (5.91e-03), -0.135 (6.65e-03) | 0.8982 (-0.125) ; +0.8956 +0.0682 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.267 (17.738 eV) | +8.5269 -1.0508 i |
| L | 0.5 | pi | 2.507e-01 (-0.905) | -1.015 (2.55e-01), -0.955 (2.54e-01), -0.905 (2.51e-01), -0.852 (2.52e-01), -0.797 (2.51e-01), -0.750 (2.54e-01) | 0.7792 (-0.905) ; +0.5964 +0.5016 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.292 (3.129 eV) | +2.4662 -0.0869 i |
| L | 0.5 | sigma | 4.334e-04 (-2.575) | -2.575 (4.33e-04) | 0.0115 (-2.575) ; +0.0004 +0.0115 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.00; saut) | — | — |
| L | 0.5 | full | 4.775e-04 (-2.575) | -2.575 (4.78e-04), -2.252 (5.93e-02), -2.185 (6.95e-02), -2.132 (7.49e-02), -2.082 (8.23e-02), -2.027 (9.04e-02) | 0.0115 (-2.575) ; +0.0004 +0.0115 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.00; saut) | -1.292 (3.129 eV) | +2.4662 -0.0869 i |
| L | 1.0 | pi | 2.442e-01 (-0.905) | -1.015 (2.49e-01), -0.955 (2.47e-01), -0.905 (2.44e-01), -0.850 (2.46e-01), -0.797 (2.44e-01), -0.747 (2.47e-01) | 0.7753 (-0.905) ; +0.5811 +0.5132 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| L | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| L | 1.0 | full | 3.831e-04 (-2.530) | -2.530 (3.83e-04), -2.185 (4.78e-02), -2.132 (5.32e-02), -2.082 (5.98e-02), -2.027 (6.70e-02), -1.972 (7.42e-02) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| L | 1.5 | pi | 2.376e-01 (-0.797) | -0.952 (2.41e-01), -0.905 (2.38e-01), -0.850 (2.39e-01), -0.797 (2.38e-01), -0.747 (2.40e-01), -0.687 (2.42e-01) | 0.7701 (-0.797) ; +0.6155 +0.4629 i | 0.982 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.185 (3.357 eV) | +2.5762 -0.0947 i |
| L | 1.5 | sigma | 2.595e-04 (-2.482) | -2.482 (2.59e-04) | 0.0102 (-2.482) ; +0.0002 +0.0102 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.483 (1.00) | — | — |
| L | 1.5 | full | 3.091e-04 (-2.482) | -2.482 (3.09e-04), -2.135 (3.55e-02), -2.082 (4.12e-02), -2.027 (4.74e-02), -1.972 (5.35e-02), -1.917 (5.94e-02) | 0.0102 (-2.482) ; +0.0002 +0.0102 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.483 (1.00) | -1.185 (3.357 eV) | +2.5762 -0.0947 i |
| L | 2.0 | pi | 2.312e-01 (-0.797) | -0.952 (2.36e-01), -0.905 (2.32e-01), -0.850 (2.33e-01), -0.797 (2.31e-01), -0.747 (2.34e-01), -0.685 (2.35e-01) | 0.7651 (-0.797) ; +0.6015 +0.4730 i | 0.983 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.185 (3.473 eV) | +2.6309 -0.0988 i |
| L | 2.0 | sigma | 2.051e-04 (-2.432) | -2.432 (2.05e-04) | 0.0096 (-2.432) ; +0.0003 +0.0096 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.433 (0.00; saut) | — | — |
| L | 2.0 | full | 2.615e-04 (-2.432) | -2.432 (2.62e-04), -2.135 (2.28e-02), -2.082 (2.75e-02), -2.030 (3.28e-02), -1.972 (3.80e-02), -1.917 (4.31e-02) | 0.0096 (-2.432) ; +0.0003 +0.0096 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.433 (0.00; saut) | -1.185 (3.473 eV) | +2.6309 -0.0988 i |
| L | 3.0 | pi | 2.196e-01 (-0.797) | -0.905 (2.22e-01), -0.850 (2.22e-01), -0.797 (2.20e-01), -0.747 (2.21e-01), -0.685 (2.22e-01), -0.637 (2.24e-01) | 0.7563 (-0.797) ; +0.5733 +0.4933 i | 0.983 / 0.014 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.125 (3.704 eV) | +2.7399 -0.1070 i |
| L | 3.0 | sigma | 1.331e-04 (-2.330) | -2.330 (1.33e-04) | 0.0087 (-2.330) ; +0.0002 +0.0087 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.330 (0.00; saut) | — | — |
| L | 3.0 | full | 1.656e-04 (-2.330) | -2.330 (1.66e-04), -2.140 (7.55e-03), -2.085 (1.07e-02), -2.030 (1.41e-02), -1.975 (1.78e-02), -1.917 (2.14e-02) | 0.0087 (-2.330) ; +0.0002 +0.0087 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.330 (0.00; saut) | -1.125 (3.704 eV) | +2.7399 -0.1070 i |
| L | 4.0 | pi | 2.092e-01 (-0.795) | -0.902 (2.13e-01), -0.850 (2.12e-01), -0.795 (2.09e-01), -0.747 (2.10e-01), -0.685 (2.11e-01), -0.637 (2.12e-01) | 0.7489 (-0.795) ; +0.5501 +0.5081 i | 0.984 / 0.013 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.125 (3.931 eV) | +2.8483 -0.1156 i |
| L | 4.0 | sigma | 9.010e-05 (-2.222) | -2.222 (9.01e-05) | 0.0079 (-2.222) ; +0.0001 +0.0079 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.223 (1.00) | — | — |
| L | 4.0 | full | 8.637e-05 (-2.222) | -2.222 (8.64e-05), -2.037 (4.57e-03), -1.977 (6.83e-03), -1.920 (9.26e-03), -1.862 (1.18e-02), -1.805 (1.42e-02) | 0.0079 (-2.222) ; +0.0001 +0.0079 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.223 (1.00) | -1.125 (3.931 eV) | +2.8483 -0.1156 i |
| L | 6.0 | pi | 1.907e-01 (-0.685) | -0.847 (1.96e-01), -0.795 (1.92e-01), -0.745 (1.92e-01), -0.685 (1.91e-01), -0.637 (1.91e-01), -0.587 (1.95e-01) | 0.7345 (-0.685) ; +0.5635 +0.4711 i | 0.986 / 0.012 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -1.015 (4.417 eV) | +3.0633 -0.1335 i |
| L | 6.0 | sigma | 4.557e-05 (-1.990) | -1.990 (4.56e-05) | 0.0068 (-1.990) ; +0.0004 +0.0068 i | 0.000 / 0.000 / 0.978 / 0.000 / 0.022 | 0.000 / 1.000 | -1.991 (0.00; saut) | — | — |
| L | 6.0 | full | 3.415e-05 (-1.990) | -1.990 (3.41e-05), -1.637 (5.26e-03), -1.580 (6.45e-03), -1.520 (7.62e-03), -1.462 (8.76e-03), -1.405 (9.83e-03) | 0.0068 (-1.990) ; +0.0004 +0.0068 i | 0.000 / 0.000 / 0.978 / 0.000 / 0.022 | 0.000 / 1.000 | -1.991 (0.00; saut) | -1.015 (4.417 eV) | +3.0633 -0.1335 i |
| L | 10.0 | pi | 1.601e-01 (-0.635) | -0.742 (1.65e-01), -0.682 (1.61e-01), -0.635 (1.60e-01), -0.587 (1.62e-01), -0.527 (1.64e-01), -0.480 (1.66e-01) | 0.7040 (-0.635) ; +0.5005 +0.4950 i | 0.988 / 0.010 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.905 (5.372 eV) | +3.4870 -0.1727 i |
| L | 10.0 | sigma | 1.530e-05 (-1.482) | -1.482 (1.53e-05) | 0.0053 (-1.482) ; +0.0002 +0.0053 i | 0.000 / 0.000 / 0.979 / 0.000 / 0.021 | 0.000 / 1.000 | -2.926 (0.22; saut), -1.484 (0.00; saut) | — | — |
| L | 10.0 | full | 7.304e-06 (-1.482) | -1.482 (7.30e-06), -1.017 (1.85e-03), -0.905 (2.43e-03), -0.855 (2.75e-03), -0.797 (3.03e-03), -0.750 (3.34e-03) | 0.0053 (-1.482) ; +0.0002 +0.0053 i | 0.000 / 0.000 / 0.979 / 0.000 / 0.021 | 0.000 / 1.000 | -2.926 (0.22; saut), -1.484 (0.00; saut) | -0.905 (5.372 eV) | +3.4870 -0.1727 i |
| NL | 0.5 | pi | 4.204e-01 (-1.185) | -1.350 (4.22e-01), -1.292 (4.21e-01), -1.237 (4.21e-01), -1.185 (4.20e-01), -1.125 (4.21e-01), -1.072 (4.23e-01) | 0.8026 (-1.185) ; +0.7144 +0.3658 i | 0.982 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -2.035 (1.442 eV) | +1.4866 -0.0315 i |
| NL | 0.5 | sigma | 5.743e-02 (-3.000) |  | 0.1778 (-3.000) ; +0.1724 +0.0433 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | — | — |
| NL | 0.5 | full | 5.366e-02 (-3.000) | -2.082 (3.33e-01), -2.030 (3.32e-01), -1.975 (3.32e-01), -1.917 (3.33e-01), -1.862 (3.33e-01), -1.805 (3.34e-01) | 0.1778 (-3.000) ; +0.1724 +0.0433 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | -2.035 (1.442 eV) | +1.4866 -0.0315 i |
| NL | 1.0 | pi | 2.442e-01 (-0.905) | -1.015 (2.49e-01), -0.955 (2.47e-01), -0.905 (2.44e-01), -0.850 (2.46e-01), -0.797 (2.44e-01), -0.747 (2.47e-01) | 0.7753 (-0.905) ; +0.5811 +0.5132 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| NL | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| NL | 1.0 | full | 3.831e-04 (-2.530) | -2.530 (3.83e-04), -2.185 (4.78e-02), -2.132 (5.32e-02), -2.082 (5.98e-02), -2.027 (6.70e-02), -1.972 (7.42e-02) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.292 (3.236 eV) | +2.5213 -0.0908 i |
| NL | 1.5 | pi | 1.741e-01 (-0.635) | -0.792 (1.79e-01), -0.742 (1.78e-01), -0.682 (1.75e-01), -0.635 (1.74e-01), -0.587 (1.77e-01), -0.527 (1.80e-01) | 0.7883 (-0.635) ; +0.5862 +0.5271 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.905 (5.037 eV) | +3.3474 -0.1603 i |
| NL | 1.5 | sigma | 1.404e-05 (-1.445) | -1.445 (1.40e-05) | 0.0052 (-1.445) ; +0.0001 +0.0052 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.914 (0.26; saut), -1.446 (0.00; saut) | — | — |
| NL | 1.5 | full | 6.081e-06 (-1.445) | -1.445 (6.08e-06), -1.020 (1.49e-03), -0.907 (2.05e-03), -0.857 (2.36e-03), -0.800 (2.65e-03), -0.752 (2.96e-03) | 0.0052 (-1.445) ; +0.0001 +0.0052 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.914 (0.26; saut), -1.446 (0.00; saut) | -0.905 (5.037 eV) | +3.3474 -0.1603 i |
| NL | 2.0 | pi | 1.400e-01 (-0.525) | -0.680 (1.45e-01), -0.632 (1.41e-01), -0.582 (1.41e-01), -0.525 (1.40e-01), -0.477 (1.40e-01), -0.432 (1.44e-01) | 0.8274 (-0.525) ; +0.5873 +0.5829 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.685 (6.540 eV) | +4.0222 -0.2317 i |
| NL | 2.0 | sigma | 2.461e-06 (-0.135) | -0.135 (2.46e-06) | 0.0036 (-0.135) ; +0.0002 +0.0036 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.23; saut), -0.136 (0.00; saut) | — | — |
| NL | 2.0 | full | 7.854e-07 (-0.137) | -1.730 (2.27e-02), -1.670 (1.92e-02), -1.610 (1.62e-02), -1.552 (1.37e-02), -1.495 (1.15e-02), -0.137 (7.85e-07) | 0.0036 (-0.135) ; +0.0002 +0.0036 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.23; saut), -0.136 (0.00; saut) | -0.685 (6.540 eV) | +4.0222 -0.2317 i |
| NL | 3.0 | pi | 1.045e-01 (-0.427) | -0.577 (1.15e-01), -0.522 (1.09e-01), -0.475 (1.05e-01), -0.427 (1.05e-01), -0.385 (1.07e-01), -0.340 (1.11e-01) | 0.9088 (-0.427) ; +0.5699 +0.7079 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.527 (9.057 eV) | +5.0588 -0.3672 i |
| NL | 3.0 | sigma | 3.525e-03 (+1.000) |  | 0.2419 (+1.000) ; -0.2419 +0.0035 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.511 (0.20; saut) | — | — |
| NL | 3.0 | full | 3.316e-03 (+0.990) | -0.377 (3.45e-03), -0.332 (3.41e-03), -0.317 (3.41e-03), -0.270 (3.40e-03), -0.230 (3.45e-03), +0.990 (3.32e-03) | 0.2419 (+1.000) ; -0.2419 +0.0035 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.511 (0.20; saut) | -0.527 (9.057 eV) | +5.0588 -0.3672 i |
| NL | 4.0 | pi | 8.936e-02 (-0.377) | -0.517 (1.01e-01), -0.472 (9.36e-02), -0.425 (8.99e-02), -0.377 (8.94e-02), -0.335 (9.09e-02), -0.285 (9.50e-02) | 0.9821 (-0.125) ; +0.9820 +0.0125 i | 0.993 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.475 (11.120 eV) | +5.8185 -0.4867 i |
| NL | 4.0 | sigma | 1.220e-02 (+1.000) |  | 0.6485 (+1.000) ; -0.6485 +0.0047 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -0.209 (0.19; saut) | — | — |
| NL | 4.0 | full | 5.148e-03 (-0.317) | -0.422 (5.48e-03), -0.375 (5.25e-03), -0.332 (5.16e-03), -0.317 (5.15e-03), -0.270 (5.17e-03), -0.230 (5.33e-03) | 0.6485 (+1.000) ; -0.6485 +0.0047 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -0.171 (0.00; saut) | -0.475 (11.120 eV) | +5.8185 -0.4867 i |
| NL | 6.0 | pi | 7.429e-02 (-0.315) | -0.422 (8.16e-02), -0.372 (7.70e-02), -0.315 (7.43e-02), -0.267 (7.46e-02), -0.230 (7.82e-02), -0.182 (8.57e-02) | 0.9771 (-0.125) ; +0.9770 +0.0162 i | 0.993 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.377 (13.564 eV) | +6.8610 -0.6786 i |
| NL | 6.0 | sigma | 2.358e-02 (+1.000) |  | 0.9993 (-3.000) ; +0.9993 +0.0000 i | 0.000 / 0.000 / 0.169 / 0.000 / 0.831 | 0.000 / 1.000 | aucune | — | — |
| NL | 6.0 | full | 6.167e-03 (-0.267) | -0.420 (7.43e-03), -0.370 (6.79e-03), -0.312 (6.32e-03), -0.267 (6.17e-03), -0.227 (6.32e-03), -0.182 (6.73e-03) | 0.9771 (-0.125) ; +0.9770 +0.0162 i | 0.993 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.377 (13.564 eV) | +6.8610 -0.6786 i |
| NL | 10.0 | pi | 6.107e-02 (-0.262) | -0.367 (7.30e-02), -0.310 (6.52e-02), -0.262 (6.11e-02), -0.225 (6.16e-02), -0.180 (6.61e-02), -0.137 (7.60e-02) | 0.9671 (-0.125) ; +0.9668 +0.0237 i | 0.993 / 0.006 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.312 (16.252 eV) | +8.0376 -0.9347 i |
| NL | 10.0 | sigma | 3.334e-02 (+1.000) |  | 0.9991 (-3.000) ; +0.9991 +0.0000 i | 0.000 / 0.000 / 0.298 / 0.000 / 0.702 | 0.000 / 1.000 | aucune | — | — |
| NL | 10.0 | full | 6.412e-03 (-0.222) | -0.365 (8.21e-03), -0.307 (7.10e-03), -0.262 (6.49e-03), -0.222 (6.41e-03), -0.180 (6.73e-03), -0.135 (7.57e-03) | 0.9671 (-0.125) ; +0.9668 +0.0237 i | 0.993 / 0.006 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.312 (16.252 eV) | +8.0376 -0.9347 i |

##### Balayage α, 600² (fenêtre [−3, +1] eV ; \|det\| normalisé sur la fenêtre)

| variante | α | bloc | min \|det\|/max (position) | minima locaux de \|det\|/max | min \|λ\| (position) ; λ | poids de v_min : p_z lacune / p_z 3 voisins / sp² A(0) / autres p_z / autres sp² | bloc π / σ | racines de Re λ_min (recouvrement ; saut) | pic de −Im T̄(K) (max) | T̄(E_D) |
|---|---|---|---|---|---|---|---|---|---|---|
| tot | 0.5 | pi | 4.471e-01 (-1.215) | -1.297 (4.48e-01), -1.270 (4.47e-01), -1.242 (4.47e-01), -1.215 (4.47e-01), -1.187 (4.47e-01), -1.160 (4.47e-01) | 0.8216 (-1.215) ; +0.7376 +0.3618 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.957 (1.271 eV) | +1.4227 -0.0274 i |
| tot | 0.5 | sigma | 8.170e-02 (-3.000) |  | 0.2153 (-3.000) ; +0.2113 +0.0414 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | — | — |
| tot | 0.5 | full | 7.715e-02 (-3.000) | -1.897 (3.61e-01), -1.867 (3.60e-01), -1.840 (3.60e-01), -1.812 (3.60e-01), -1.782 (3.60e-01), -1.755 (3.60e-01) | 0.2153 (-3.000) ; +0.2113 +0.0414 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | -1.957 (1.271 eV) | +1.4227 -0.0274 i |
| tot | 1.0 | pi | 2.479e-01 (-0.830) | -0.910 (2.49e-01), -0.882 (2.48e-01), -0.855 (2.48e-01), -0.830 (2.48e-01), -0.802 (2.48e-01), -0.777 (2.48e-01) | 0.7893 (-0.830) ; +0.6272 +0.4792 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| tot | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| tot | 1.0 | full | 3.810e-04 (-2.530) | -2.530 (3.81e-04), -1.472 (1.19e-01), -1.442 (1.21e-01), -1.415 (1.22e-01), -1.387 (1.24e-01), -1.357 (1.25e-01) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| tot | 1.5 | pi | 1.736e-01 (-0.642) | -0.692 (1.74e-01), -0.665 (1.74e-01), -0.642 (1.74e-01), -0.617 (1.74e-01), -0.590 (1.74e-01), -0.567 (1.75e-01) | 0.8030 (-0.642) ; +0.5803 +0.5551 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.882 (4.984 eV) | +3.3969 -0.1564 i |
| tot | 1.5 | sigma | 1.254e-05 (-1.380) | -1.380 (1.25e-05) | 0.0051 (-1.380) ; -0.0001 +0.0051 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.901 (0.24; saut), -1.379 (0.00; saut) | — | — |
| tot | 1.5 | full | 5.146e-06 (-1.380) | -1.380 (5.15e-06) | 0.0051 (-1.380) ; -0.0001 +0.0051 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.901 (0.24; saut), -1.379 (0.00; saut) | -0.882 (4.984 eV) | +3.3969 -0.1564 i |
| tot | 2.0 | pi | 1.373e-01 (-0.537) | -0.582 (1.38e-01), -0.560 (1.38e-01), -0.537 (1.37e-01), -0.512 (1.37e-01), -0.487 (1.38e-01), -0.465 (1.39e-01) | 0.8363 (-0.537) ; +0.5602 +0.6210 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.720 (6.607 eV) | +4.1117 -0.2293 i |
| tot | 2.0 | sigma | 2.111e-06 (+0.013) | +0.013 (2.11e-06) | 0.0034 (+0.013) ; -0.0000 +0.0034 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.529 (0.21; saut), +0.013 (0.00; saut) | — | — |
| tot | 2.0 | full | 8.766e-07 (+0.013) | -2.432 (2.25e-01), -2.412 (2.25e-01), +0.013 (8.77e-07) | 0.0034 (+0.013) ; -0.0000 +0.0034 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.529 (0.21; saut), +0.013 (0.00; saut) | -0.720 (6.607 eV) | +4.1117 -0.2293 i |
| tot | 3.0 | pi | 1.025e-01 (-0.410) | -2.030 (5.25e-01), -2.002 (5.10e-01), -1.972 (4.96e-01), -0.455 (1.03e-01), -0.432 (1.03e-01), -0.410 (1.02e-01) | 0.9267 (-0.410) ; +0.5874 +0.7167 i | 0.982 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.535 (9.225 eV) | +5.2133 -0.3693 i |
| tot | 3.0 | sigma | 4.482e-03 (+1.000) |  | 0.2862 (+1.000) ; -0.2862 +0.0037 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.376 (0.16; saut) | — | — |
| tot | 3.0 | full | 3.560e-03 (-0.307) | -0.330 (3.57e-03), -0.307 (3.56e-03), -0.287 (3.57e-03), +0.248 (5.65e-03), +0.275 (5.65e-03), +0.300 (5.64e-03) | 0.2862 (+1.000) ; -0.2862 +0.0037 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.376 (0.16; saut) | -0.535 (9.225 eV) | +5.2133 -0.3693 i |
| tot | 4.0 | pi | 8.581e-02 (-0.355) | -2.057 (5.51e-01), -2.030 (5.35e-01), -2.002 (5.20e-01), -1.972 (5.05e-01), -0.355 (8.58e-02), -0.332 (8.60e-02) | 0.9605 (-0.130) ; +0.9600 +0.0293 i | 0.992 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.435 (11.227 eV) | +6.0279 -0.4944 i |
| tot | 4.0 | sigma | 1.332e-02 (+1.000) |  | 0.7149 (+1.000) ; -0.7149 +0.0049 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | +0.016 (0.13; saut) | — | — |
| tot | 4.0 | full | 5.079e-03 (-0.307) | -0.327 (5.11e-03), -0.307 (5.08e-03), -0.287 (5.09e-03), +0.518 (1.19e-02), +0.543 (1.19e-02), +0.570 (1.19e-02) | 0.7149 (+1.000) ; -0.7149 +0.0049 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | +0.064 (0.00; saut) | -0.435 (11.227 eV) | +6.0279 -0.4944 i |
| tot | 6.0 | pi | 6.973e-02 (-0.285) | -2.030 (5.46e-01), -2.002 (5.31e-01), -1.972 (5.17e-01), -1.942 (5.03e-01), -0.302 (7.00e-02), -0.285 (6.97e-02) | 0.9410 (-0.127) ; +0.9401 +0.0425 i | 0.993 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.352 (14.029 eV) | +7.1695 -0.7011 i |
| tot | 6.0 | sigma | 2.406e-02 (+1.000) |  | 0.9964 (-3.000) ; +0.9964 +0.0000 i | 0.000 / 0.000 / 0.041 / 0.000 / 0.959 | 0.000 / 1.000 | aucune | — | — |
| tot | 6.0 | full | 5.868e-03 (-0.260) | -2.492 (4.00e-01), -2.470 (4.00e-01), -0.277 (5.89e-03), -0.260 (5.87e-03) | 0.9410 (-0.127) ; +0.9401 +0.0425 i | 0.993 / 0.006 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.352 (14.029 eV) | +7.1695 -0.7011 i |
| tot | 10.0 | pi | 5.678e-02 (-0.232) | -2.990 (5.70e-01), -2.030 (5.57e-01), -2.002 (5.42e-01), -1.972 (5.28e-01), -1.942 (5.15e-01), -0.232 (5.68e-02) | 0.9029 (-0.125) ; +0.9003 +0.0684 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.280 (17.386 eV) | +8.5398 -0.9980 i |
| tot | 10.0 | sigma | 3.228e-02 (+1.000) |  | 0.9941 (-3.000) ; +0.9941 +0.0001 i | 0.000 / 0.000 / 0.041 / 0.000 / 0.959 | 0.000 / 1.000 | aucune | — | — |
| tot | 10.0 | full | 5.799e-03 (-0.225) | -2.520 (4.25e-01), -2.495 (4.24e-01), -2.472 (4.25e-01), -0.225 (5.80e-03) | 0.9029 (-0.125) ; +0.9003 +0.0684 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.280 (17.386 eV) | +8.5398 -0.9980 i |
| L | 0.5 | pi | 2.548e-01 (-0.855) | -0.910 (2.55e-01), -0.882 (2.55e-01), -0.855 (2.55e-01), -0.830 (2.55e-01), -0.802 (2.55e-01), -0.777 (2.56e-01) | 0.7940 (-0.855) ; +0.6305 +0.4826 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.272 (3.042 eV) | +2.4665 -0.0824 i |
| L | 0.5 | sigma | 4.334e-04 (-2.575) | -2.575 (4.33e-04) | 0.0115 (-2.575) ; +0.0004 +0.0115 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.00; saut) | — | — |
| L | 0.5 | full | 4.751e-04 (-2.575) | -2.575 (4.75e-04), -1.615 (1.41e-01), -1.585 (1.43e-01), -1.557 (1.45e-01), -1.527 (1.47e-01), -1.500 (1.49e-01) | 0.0115 (-2.575) ; +0.0004 +0.0115 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.00; saut) | -1.272 (3.042 eV) | +2.4665 -0.0824 i |
| L | 1.0 | pi | 2.479e-01 (-0.830) | -0.910 (2.49e-01), -0.882 (2.48e-01), -0.855 (2.48e-01), -0.830 (2.48e-01), -0.802 (2.48e-01), -0.777 (2.48e-01) | 0.7893 (-0.830) ; +0.6272 +0.4792 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| L | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| L | 1.0 | full | 3.810e-04 (-2.530) | -2.530 (3.81e-04), -1.472 (1.19e-01), -1.442 (1.21e-01), -1.415 (1.22e-01), -1.387 (1.24e-01), -1.357 (1.25e-01) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| L | 1.5 | pi | 2.413e-01 (-0.802) | -0.880 (2.42e-01), -0.855 (2.41e-01), -0.827 (2.41e-01), -0.802 (2.41e-01), -0.777 (2.42e-01), -0.750 (2.42e-01) | 0.7847 (-0.802) ; +0.6260 +0.4732 i | 0.982 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.215 (3.263 eV) | +2.5765 -0.0899 i |
| L | 1.5 | sigma | 2.595e-04 (-2.482) | -2.482 (2.59e-04) | 0.0102 (-2.482) ; +0.0002 +0.0102 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.483 (1.00) | — | — |
| L | 1.5 | full | 3.098e-04 (-2.482) | -2.482 (3.10e-04), -1.302 (1.01e-01), -1.275 (1.02e-01), -1.162 (1.07e-01) | 0.0102 (-2.482) ; +0.0002 +0.0102 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.483 (1.00) | -1.215 (3.263 eV) | +2.5765 -0.0899 i |
| L | 2.0 | pi | 2.349e-01 (-0.802) | -0.880 (2.36e-01), -0.855 (2.35e-01), -0.827 (2.35e-01), -0.802 (2.35e-01), -0.777 (2.35e-01), -0.750 (2.35e-01) | 0.7801 (-0.802) ; +0.6122 +0.4835 i | 0.983 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.187 (3.375 eV) | +2.6313 -0.0937 i |
| L | 2.0 | sigma | 2.051e-04 (-2.432) | -2.432 (2.05e-04) | 0.0096 (-2.432) ; +0.0003 +0.0096 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.433 (0.00; saut) | — | — |
| L | 2.0 | full | 2.602e-04 (-2.432) | -2.432 (2.60e-04) | 0.0096 (-2.432) ; +0.0003 +0.0096 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.433 (0.00; saut) | -1.187 (3.375 eV) | +2.6313 -0.0937 i |
| L | 3.0 | pi | 2.231e-01 (-0.775) | -0.855 (2.24e-01), -0.827 (2.24e-01), -0.800 (2.23e-01), -0.775 (2.23e-01), -0.747 (2.23e-01), -0.722 (2.24e-01) | 0.7716 (-0.775) ; +0.5991 +0.4862 i | 0.983 / 0.014 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.132 (3.600 eV) | +2.7404 -0.1015 i |
| L | 3.0 | sigma | 1.331e-04 (-2.330) | -2.330 (1.33e-04) | 0.0087 (-2.330) ; +0.0002 +0.0087 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.330 (0.00; saut) | — | — |
| L | 3.0 | full | 1.609e-04 (-2.330) | -2.330 (1.61e-04) | 0.0087 (-2.330) ; +0.0002 +0.0087 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.330 (0.00; saut) | -1.132 (3.600 eV) | +2.7404 -0.1015 i |
| L | 4.0 | pi | 2.122e-01 (-0.747) | -0.825 (2.13e-01), -0.800 (2.13e-01), -0.775 (2.12e-01), -0.747 (2.12e-01), -0.722 (2.12e-01), -0.697 (2.13e-01) | 0.7634 (-0.747) ; +0.5874 +0.4876 i | 0.984 / 0.013 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.075 (3.827 eV) | +2.8488 -0.1096 i |
| L | 4.0 | sigma | 9.010e-05 (-2.222) | -2.222 (9.01e-05) | 0.0079 (-2.222) ; +0.0001 +0.0079 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.223 (1.00) | — | — |
| L | 4.0 | full | 8.363e-05 (-2.222) | -2.222 (8.36e-05) | 0.0079 (-2.222) ; +0.0001 +0.0079 i | 0.000 / 0.000 / 0.977 / 0.000 / 0.023 | 0.000 / 1.000 | -2.223 (1.00) | -1.075 (3.827 eV) | +2.8488 -0.1096 i |
| L | 6.0 | pi | 1.931e-01 (-0.695) | -0.772 (1.94e-01), -0.745 (1.94e-01), -0.720 (1.93e-01), -0.695 (1.93e-01), -0.667 (1.93e-01), -0.645 (1.94e-01) | 0.7482 (-0.695) ; +0.5675 +0.4876 i | 0.986 / 0.012 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.992 (4.288 eV) | +3.0639 -0.1266 i |
| L | 6.0 | sigma | 4.557e-05 (-1.990) | -1.990 (4.56e-05) | 0.0068 (-1.990) ; +0.0004 +0.0068 i | 0.000 / 0.000 / 0.978 / 0.000 / 0.022 | 0.000 / 1.000 | -1.991 (0.00; saut) | — | — |
| L | 6.0 | full | 3.204e-05 (-1.990) | -1.990 (3.20e-05) | 0.0068 (-1.990) ; +0.0004 +0.0068 i | 0.000 / 0.000 / 0.978 / 0.000 / 0.022 | 0.000 / 1.000 | -1.991 (0.00; saut) | -0.992 (4.288 eV) | +3.0639 -0.1266 i |
| L | 10.0 | pi | 1.628e-01 (-0.615) | -0.690 (1.64e-01), -0.665 (1.63e-01), -0.640 (1.63e-01), -0.615 (1.63e-01), -0.590 (1.63e-01), -0.565 (1.64e-01) | 0.7214 (-0.615) ; +0.5284 +0.4912 i | 0.988 / 0.010 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.855 (5.228 eV) | +3.4879 -0.1638 i |
| L | 10.0 | sigma | 1.530e-05 (-1.482) | -1.482 (1.53e-05) | 0.0053 (-1.482) ; +0.0002 +0.0053 i | 0.000 / 0.000 / 0.979 / 0.000 / 0.021 | 0.000 / 1.000 | -2.926 (0.22; saut), -1.484 (0.00; saut) | — | — |
| L | 10.0 | full | 6.980e-06 (-1.482) | -1.482 (6.98e-06) | 0.0053 (-1.482) ; +0.0002 +0.0053 i | 0.000 / 0.000 / 0.979 / 0.000 / 0.021 | 0.000 / 1.000 | -2.926 (0.22; saut), -1.484 (0.00; saut) | -0.855 (5.228 eV) | +3.4879 -0.1638 i |
| NL | 0.5 | pi | 4.289e-01 (-1.187) | -1.270 (4.30e-01), -1.242 (4.29e-01), -1.215 (4.29e-01), -1.187 (4.29e-01), -1.160 (4.29e-01), -1.132 (4.29e-01) | 0.8132 (-1.187) ; +0.7246 +0.3692 i | 0.982 / 0.015 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.930 (1.366 eV) | +1.4867 -0.0299 i |
| NL | 0.5 | sigma | 5.743e-02 (-3.000) |  | 0.1778 (-3.000) ; +0.1724 +0.0433 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | — | — |
| NL | 0.5 | full | 5.442e-02 (-3.000) | -1.867 (3.45e-01), -1.840 (3.45e-01), -1.810 (3.45e-01), -1.782 (3.45e-01), -1.755 (3.45e-01), -1.725 (3.45e-01) | 0.1778 (-3.000) ; +0.1724 +0.0433 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | aucune | -1.930 (1.366 eV) | +1.4867 -0.0299 i |
| NL | 1.0 | pi | 2.479e-01 (-0.830) | -0.910 (2.49e-01), -0.882 (2.48e-01), -0.855 (2.48e-01), -0.830 (2.48e-01), -0.802 (2.48e-01), -0.777 (2.48e-01) | 0.7893 (-0.830) ; +0.6272 +0.4792 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| NL | 1.0 | sigma | 3.328e-04 (-2.530) | -2.530 (3.33e-04) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | — | — |
| NL | 1.0 | full | 3.810e-04 (-2.530) | -2.530 (3.81e-04), -1.472 (1.19e-01), -1.442 (1.21e-01), -1.415 (1.22e-01), -1.387 (1.24e-01), -1.357 (1.25e-01) | 0.0108 (-2.530) ; +0.0003 +0.0108 i | 0.000 / 0.000 / 0.976 / 0.000 / 0.024 | 0.000 / 1.000 | -2.531 (1.00) | -1.242 (3.153 eV) | +2.5216 -0.0861 i |
| NL | 1.5 | pi | 1.770e-01 (-0.642) | -0.720 (1.78e-01), -0.692 (1.77e-01), -0.667 (1.77e-01), -0.642 (1.77e-01), -0.617 (1.77e-01), -0.592 (1.78e-01) | 0.8064 (-0.642) ; +0.5921 +0.5474 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.910 (4.877 eV) | +3.3482 -0.1520 i |
| NL | 1.5 | sigma | 1.404e-05 (-1.445) | -1.445 (1.40e-05) | 0.0052 (-1.445) ; +0.0001 +0.0052 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.914 (0.26; saut), -1.446 (0.00; saut) | — | — |
| NL | 1.5 | full | 6.088e-06 (-1.445) | -1.445 (6.09e-06) | 0.0052 (-1.445) ; +0.0001 +0.0052 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.914 (0.26; saut), -1.446 (0.00; saut) | -0.910 (4.877 eV) | +3.3482 -0.1520 i |
| NL | 2.0 | pi | 1.414e-01 (-0.537) | -1.972 (4.79e-01), -0.585 (1.42e-01), -0.562 (1.41e-01), -0.537 (1.41e-01), -0.515 (1.42e-01), -0.490 (1.42e-01) | 0.8418 (-0.537) ; +0.5811 +0.6090 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.722 (6.401 eV) | +4.0236 -0.2197 i |
| NL | 2.0 | sigma | 2.461e-06 (-0.135) | -0.135 (2.46e-06) | 0.0036 (-0.135) ; +0.0002 +0.0036 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.23; saut), -0.136 (0.00; saut) | — | — |
| NL | 2.0 | full | 7.799e-07 (-0.137) | -2.412 (2.17e-01), -0.137 (7.80e-07) | 0.0036 (-0.135) ; +0.0002 +0.0036 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -2.576 (0.23; saut), -0.136 (0.00; saut) | -0.722 (6.401 eV) | +4.0236 -0.2197 i |
| NL | 3.0 | pi | 1.067e-01 (-0.432) | -2.030 (5.22e-01), -2.002 (5.07e-01), -1.972 (4.93e-01), -0.455 (1.07e-01), -0.432 (1.07e-01), -0.412 (1.07e-01) | 0.9367 (-0.432) ; +0.5774 +0.7376 i | 0.981 / 0.016 / 0.000 / 0.003 / 0.000 | 1.000 / 0.000 | aucune | -0.537 (8.857 eV) | +5.0614 -0.3484 i |
| NL | 3.0 | sigma | 3.525e-03 (+1.000) |  | 0.2419 (+1.000) ; -0.2419 +0.0035 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.511 (0.20; saut) | — | — |
| NL | 3.0 | full | 3.266e-03 (+1.000) | -2.462 (3.16e-01), -2.442 (3.16e-01), -0.330 (3.38e-03), -0.307 (3.37e-03), -0.287 (3.37e-03) | 0.2419 (+1.000) ; -0.2419 +0.0035 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -1.511 (0.20; saut) | -0.537 (8.857 eV) | +5.0614 -0.3484 i |
| NL | 4.0 | pi | 9.008e-02 (-0.357) | -2.030 (5.31e-01), -2.002 (5.16e-01), -1.972 (5.02e-01), -0.377 (9.02e-02), -0.357 (9.01e-02), -0.335 (9.06e-02) | 0.9830 (-0.117) ; +0.9829 +0.0113 i | 0.993 / 0.005 / 0.000 / 0.002 / 0.000 | 1.000 / 0.000 | aucune | -0.457 (10.718 eV) | +5.8225 -0.4618 i |
| NL | 4.0 | sigma | 1.220e-02 (+1.000) |  | 0.6485 (+1.000) ; -0.6485 +0.0047 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -0.209 (0.19; saut) | — | — |
| NL | 4.0 | full | 5.096e-03 (-0.307) | -0.330 (5.11e-03), -0.307 (5.10e-03), -0.287 (5.12e-03), +0.945 (1.11e-02), +0.970 (1.11e-02), +0.998 (1.11e-02) | 0.6485 (+1.000) ; -0.6485 +0.0047 i | 0.000 / 0.000 / 0.975 / 0.000 / 0.025 | 0.000 / 1.000 | -0.176 (0.00; saut) | -0.457 (10.718 eV) | +5.8225 -0.4618 i |
| NL | 6.0 | pi | 7.399e-02 (-0.307) | -2.030 (5.40e-01), -2.002 (5.26e-01), -1.972 (5.12e-01), -0.325 (7.45e-02), -0.307 (7.40e-02), -0.287 (7.42e-02) | 0.9782 (-0.115) ; +0.9781 +0.0142 i | 0.993 / 0.005 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.357 (13.262 eV) | +6.8677 -0.6440 i |
| NL | 6.0 | sigma | 2.358e-02 (+1.000) |  | 0.9993 (-3.000) ; +0.9993 +0.0000 i | 0.000 / 0.000 / 0.169 / 0.000 / 0.831 | 0.000 / 1.000 | aucune | — | — |
| NL | 6.0 | full | 6.175e-03 (-0.282) | -2.492 (4.00e-01), -2.470 (4.01e-01), -0.282 (6.18e-03), -0.265 (6.19e-03) | 0.9782 (-0.115) ; +0.9781 +0.0142 i | 0.993 / 0.005 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.357 (13.262 eV) | +6.8677 -0.6440 i |
| NL | 10.0 | pi | 6.166e-02 (-0.255) | -2.057 (5.65e-01), -2.030 (5.49e-01), -2.002 (5.34e-01), -1.972 (5.21e-01), -1.942 (5.07e-01), -0.255 (6.17e-02) | 0.9686 (-0.115) ; +0.9683 +0.0206 i | 0.993 / 0.005 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.305 (16.180 eV) | +8.0484 -0.8875 i |
| NL | 10.0 | sigma | 3.334e-02 (+1.000) |  | 0.9991 (-3.000) ; +0.9991 +0.0000 i | 0.000 / 0.000 / 0.298 / 0.000 / 0.702 | 0.000 / 1.000 | aucune | — | — |
| NL | 10.0 | full | 6.457e-03 (-0.235) | -2.520 (4.31e-01), -2.495 (4.30e-01), -2.472 (4.31e-01), -0.235 (6.46e-03) | 0.9686 (-0.115) ; +0.9683 +0.0206 i | 0.993 / 0.005 / 0.000 / 0.001 / 0.000 | 1.000 / 0.000 | aucune | -0.305 (16.180 eV) | +8.0484 -0.8875 i |


Figures : `fig/d3_alpha_full`, `fig/d3_alpha_pi` (|det|/max, min|λ|, −Im T̄(K) contre ε − E_D pour les huit α, 600²).

## D5 — Alignement de ΔV^L, toutes tailles non relaxées (J3, `d5`, 35 s)

Composante G = 0 = moyenne 3D de ΔV^L = V_d − V_p (pp.x plot_num 1, sans soustraction) ; ligne le long de a₁ du site de la lacune à la
demi-boîte (interpolation linéaire sur la grille) ; plans frontière = plans i₁ = (s₁ + ½)·n₁ et i₂ = (s₂ + ½)·n₂ ; « aligné » =
ΔV − décalage Lu (1,0 Å). Soustraction en production : aucune dans la construction de M dense (`compute_M_dense_stages.py`,
`subtract_mean=False`) ; 67,0 meV (moyenne diagonale de M^L dense) dans le seul contrôle C14 ; moyenne diagonale de M^L grossier
16 bandes 76,58 meV (J1).

#### D5 — ΔV^L par taille (eV sauf mention ; ligne a₁ de la lacune à la demi-boîte ; plans frontière à la demi-boîte)

| N | A_sc (Å²) | atome loin / distance (Å) | ⟨ΔV⟩_3D = composante G = 0 (meV) | ΔV au site (eV) | bord ligne a₁ brut / − moyenne / aligné (meV) | plans frontière max\|ΔV − c\|, tous z : brut / − moyenne / aligné (meV) | idem, plan des atomes (meV) | décalage Lu 0,5 / 1,0 Å (meV) |
|---|---|---|---|---|---|---|---|---|
| 5 | 131.6 | 49 / 9.97 | +28.62 | +115.500 | -141.6 / -170.2 / -19.3 | 151.3 / 179.9 / 210.9 | 151.3 / 179.9 / 135.0 | -131.64 / -122.33 |
| 6 | 189.6 | 71 / 12.81 | +20.83 | +116.220 | -141.1 / -161.9 / -42.5 | 150.5 / 171.4 / 184.9 | 142.2 / 163.0 / 78.3 | -125.83 / -98.62 |
| 7 | 258.0 | 97 / 14.24 | +14.55 | +114.698 | -40.7 / -55.3 / -2.2 | 109.1 / 123.7 / 129.5 | 48.9 / 63.5 / 77.5 | -37.93 / -38.54 |
| 8 | 337.0 | 1 / 17.08 | +11.23 | +115.649 | -25.8 / -37.0 / -2.9 | 36.8 / 48.1 / 41.1 | 36.4 / 47.6 / 26.5 | -24.62 / -22.97 |
| 9 | 426.5 | 161 / 18.51 | +9.09 | +116.243 | -24.0 / -33.0 / +0.7 | 34.1 / 43.2 / 41.9 | 33.0 / 42.1 / 17.0 | -21.59 / -24.65 |
| 10 | 526.6 | 198 / 19.93 | +7.17 | +115.586 | -11.8 / -18.9 / -3.5 | 17.1 / 24.3 / 11.3 | 17.1 / 24.2 / 10.4 | -7.80 / -8.22 |
| 11 | 637.2 | 241 / 22.78 | +5.95 | +114.415 | -11.4 / -17.4 / +0.7 | 28.3 / 34.3 / 33.3 | 15.5 / 21.5 / 21.2 | -12.09 / -12.16 |
| 12 | 758.3 | 287 / 25.63 | +5.06 | +116.226 | -39.7 / -44.8 / -10.5 | 46.1 / 51.1 / 46.3 | 42.5 / 47.6 / 16.8 | -36.07 / -29.24 |



Dépendance en taille (meV) : ⟨ΔV⟩_3D 28,6 (5×5) → 20,8 → 14,6 → 11,2 → 9,1 (9×9) → 7,2 → 6,0 → 5,1 (12×12) ; bord de la ligne a₁ brut
−141,6, −141,1, −40,7, −25,8, −24,0, −11,8, −11,4, −39,7 ; aligné −19,3, −42,5, −2,2, −2,9, +0,7, −3,5, +0,7, −10,5 ; décalage Lu (1,0 Å)
−122,3, −98,6, −38,5, −23,0, −24,7, −8,2, −12,2, −29,2 (0,5 Å : −131,6, −125,8, −37,9, −24,6, −21,6, −7,8, −12,1, −36,1) ; max |ΔV| sur les
plans frontière (tous z) aligné 210,9, 184,9, 129,5, 41,1, 41,9, 11,3, 33,3, 46,3 ; dans le plan des atomes 135,0, 78,3, 77,5, 26,5, 17,0,
10,4, 21,2, 16,8. Les familles N = 3m (6, 9, 12) et N ≠ 3m ne suivent pas la même suite (CLAUDE.md, N mod 3).

Figure : `fig/d5_boundary` (bord contre N ; profils ΔV^L(a₁) des huit tailles).

## R — Préparation des géométries relaxées (J3, `r`, 6 s ; rien codé)

Lecture du code (fichiers et lignes de la chaîne actuelle) :

- **Partie non locale** (`defects/non_local.py`, `compute_M_NL`, l. 109–201) : les positions de **tous** les atomes sont lues
  séparément dans `sc_p` (`tau_s_p`) et `sc_d` (`tau_s_d`), les phases e^{∓iK·τ} sont calculées atome par atome et M^NL = M_d − M_p :
  une géométrie relaxée est traitée par construction (il suffit de passer le `.save` relaxé comme `sc_d`). Aucun indice de spin ;
  les projecteurs KB (PP_DIJ, β_ℓ) ne dépendent pas du spin : **ΔV^NL est identique pour ↑ et ↓** à géométrie donnée.
- **Partie locale** (`defects/local_R.py` l. 15–97, `local_G.py`) : ΔV = V_d − V_p lus par `qe_io.get_pot` (un fichier filplot,
  aucune notion de spin, l. 393–453) ; positions non utilisées. nspin = 2 exige deux potentiels pp.x (spin_component 1 et 2, fournis
  par J0) et **deux passages de la partie locale seulement** ; la somme M^L_σ + M^NL par spin reste à assembler par un pilote.
- **Lecteur de wfc** (`qe_io._read_all_wfc`, l. 342–391) : `npol = 1` exigé ; les fichiers `wfcup1`/`wfcdw1` d'un run lsda (même
  attribut `ik`) ne sont pas distingués ; gamma_only non géré. Sans effet sur la chaîne M (seules les wfc de maille sont lues).
- Manquant : pilote à deux passages locaux + un passage non local avec `sc_d` relaxé, somme par spin, sidecar `matrix_io` portant le
  spin, et deux potentiels pp.x de la parfaite si l'on voulait aussi une référence spin-polarisée.

Décroissance de ΔV au bord pour R1 (potentiel de la géométrie relaxée moins V_p de la 9×9 parfaite non relaxée, même grille) :

#### R — ΔV = V(géométrie) − V_p(9×9 parfaite non relaxée) (eV)

| géométrie | décalage Lu utilisé (eV) | ⟨ΔV⟩_3D (meV) | ΔV au site (eV) | ligne a₁ : t = 0 / 0,125 / 0,25 / 0,375 / 0,5 (eV) | bord ligne brut / aligné (meV) | plans frontière max\|ΔV\| brut / aligné (meV) | idem plan des atomes (meV) |
|---|---|---|---|---|---|---|---|
| D_unrelaxed | -0.0247 | +9.09 | +116.243 | +116.243 / -0.208 / -0.054 / -0.043 / -0.024 | -24.0 / +0.7 | 34.1 / 41.9 | 33.0 / 17.0 |
| N1 | +0.0202 | +9.07 | +116.079 | +116.079 / -4.337 / -1.858 / +0.221 / +0.089 | +89.0 / +68.8 | 2587.8 / 2567.6 | 2587.8 / 2567.6 |
| N2u | +0.0159 | -18.61 | +115.619 | +115.619 / -8.130 / -3.069 / +0.262 / +0.092 | +91.7 / +75.8 | 3673.8 / 3657.9 | 3673.8 / 3657.9 |
| N2d | +0.0356 | -15.97 | +114.915 | +114.915 / -7.953 / -3.039 / +0.262 / +0.089 | +89.2 / +53.6 | 3668.8 / 3633.3 | 3668.8 / 3633.3 |

V_↑ − V_↓ (R1 nspin2) : max\|·\| = 10.130 eV, moyenne 3D -2.64 meV, au site +0.704 eV, plans frontière max 56.99 meV.




Sur les plans frontière, max |ΔV| passe de 34 meV (non relaxée) à 2,59 (nspin1) et 3,67 eV (nspin2) : les plans coupent des cœurs
atomiques déplacés de ~1e-3 Å entre les deux géométries. Le bord de la ligne a₁ vaut +89 / +92 / +89 meV (nspin1, ↑, ↓) contre
−24 meV (non relaxée) ; la ligne passe de −4,3 eV (t = 0,125) à +0,22 eV (t = 0,375) pour nspin1 contre −0,21 / −0,04 eV pour la non
relaxée. V_↑ − V_↓ : max 10,13 eV (près des cœurs), moyenne 3D −2,64 meV, au site de la lacune +0,704 eV, plans frontière 57 meV.

## Fichiers de la campagne

Répertoire de travail `graphene/qe/defects/R4_quasi_lie/` : `README.md`, `R4_rapport.md` (ce fichier), `r4_driver.py`, `submit_r4.sh`,
`JOBID`, `r4_log.txt`, `sections/` (rédaction par section, `tables.md` généré), `prep/` (json, `M_coarse_*.npy`, `M_NL_coarse_new_9x9.npy`
+ sidecar), `d6/`, `d1/` (json, npz, `projwfc_inputs/`), `d4/`, `d3/`, `d5/`, `r/` (json + npz), `cache/` (M_W tot/L/NL 638 Mo, V_loc,
g⁽⁰⁾ jouet et réel à 300², 600², 1200² : 3.4G au total), `fig/` (pdf + png), `slurm-r4-*`, `failed_runs/`. Copie versionnée
`graphene-raman/article/R4_quasi_lie/` : tout sauf `cache/`, `slurm-*`, `JOBID`, `failed_runs/`, `__pycache__`, fichiers `.npy`
(`M_coarse_*.npy` 27 Mo, `M_NL_coarse_new_9x9.npy` ; `d3_curves.npz` 4,8 Mio est copié) ; le marqueur `D4_FAIL` de la première exécution de D4 est rangé dans `failed_runs/`.
Modules nouveaux (non commités) : `src/electron_defect_interaction/{io/qe_gamma_io.py, io/projwfc_io.py, defects/alignment.py,
defects/many_body/pole_criterion.py, defects/many_body/tb_models.py, wannier/supercell_fold.py}`. `tb_models.py` n'est utilisé que par
D6 (modèle synthétique, V_loc « U » et « site retiré ») : matériel de test, à déplacer sous `tests/` si le dépôt le préfère ;
les cinq autres sont des lecteurs et des outils réutilisables. Fichier de production modifié : aucun (`article/R1_vacancy_relaxed/README.md`
et `super_cell_relaxed/9x9/R1_rapport.md` : une ligne pour J0). Git : rien ajouté, rien commité (Greg fait add, commit, push).
Aucune suppression ; les caches (3.4G) et les potentiels pp.x de R1 (3 × 241 Mo + 3 × 184 Mo de `pp*.out`) sont les seuls fichiers
nouveaux volumineux ; leur nettoyage éventuel passe par un manifeste et un GO séparé.

Figures : `d6_toy_1200`, `d6_real_1200`, `d6_toy_grid`, `d1_states`, `d4_ladder`, `d3_alpha_full`, `d3_alpha_pi`, `d5_boundary`.

**STOP — rapport R4 terminé le 2026-09-25.**

## Erratum (R5, 2026-09-25)

Ajouté le 2026-09-25 par la campagne R5 (`graphene/qe/defects/R5_base_vs_M/R5_rapport.md`, §0.3 et §A.4) ; rien n'est modifié au-dessus.

**Cause.** Les poids w₂ des variantes (a2), (a3), (b), (c-all), (c-3) du §D4 sont faux. Le pilote R4 construisait l'index d'ondes planes
de la super-cellule avec `sc_planewave_index(k81, G20, …)`, k81 étant les 81 k de la maille grossière (écrits dans [−4/9, 4/9]),
alors que les coefficients C20 lus dans le `.save` dense sont référés aux k du dense, écrits dans [0, 26/27] : 56 des 81 k diffèrent
d'un vecteur entier ΔG ∈ {(0,1), (1,0), (1,1)}. Chaque ψ_nk du dense était donc multiplié par une phase périodique e^{2πi·9ΔG·x} de
module 1 : |ψ_nk|² intact (états de Bloch purs à w₂ = 0,027–0,030, ce qui a masqué l'erreur), interférences entre k fausses pour toute
superposition. (a1) utilise C16 et k81 du même `.save` : correct. Les énergies, les portes 1 et 2, les « poids WF site + voisins »
et les marches en énergie du §D4 ne dépendent pas de cet index et restent valables ; seuls les w₂ et les étiquettes « localisé » des
cinq variantes sont à remplacer. Vérification R5 (A.4) : w₂(a2) = w₂(a1) à 2,7e-4 sur les cinq états localisés de (a1).

**Tableau D4 corrigé** (ε − E_D en eV ; w₂ corrigé ; entre parenthèses la valeur fausse du §D4 ; même seuil 0,0812) :

| variante | pairs (σ) localisés | impairs (π) localisés |
|---|---|---|
| QE (D1) | +0,101 ×2 (0,713) | −1,759 ×2 (0,114) ; −0,737 (0,246) ; +0,269 (0,102) |
| (a1) M_ed 16 bandes, total | −2,812 ×2 (0,287) | −1,803 ×2 (0,101) ; −1,350 (0,226) |
| (a1) M^L seul | aucun | −1,803 ×2 (0,101) ; −1,780 (0,103) |
| (a1) M^NL seul | −2,847 ×2 (0,255) | −1,370 (0,222) |
| (a2) M dense ⊂ 81 k, 16 bandes | −2,812 ×2 (**0,287** ; était 0,058 / 0,074, « sous le seuil ») | −1,803 ×2 (**0,101** ; 0,091 / 0,090) ; −1,350 (**0,226** ; 0,084) |
| (a3) idem 20 bandes | −2,848 ×2 (**0,254** ; 0,067 / 0,055, « sous le seuil ») | −1,803 ×2 (**0,101** ; 0,091 / 0,090) ; −1,359 (**0,225** ; 0,084) |
| (b) 5 WF, V†εV | −2,514 ×2 (**0,489** ; 0,106 / 0,088) | −1,803 ×2 (**0,101** ; 0,091 / 0,090) ; −1,315 (**0,232** ; 0,083) |
| (c-all) H(R) + M_W replié, toutes mailles | −2,514 ×2 (**0,489** ; 0,106 / 0,088) | −1,803 ×2 (**0,101**) ; −1,315 (**0,232** ; 0,083) |
| (c-3) idem, R_cut = 3 | −2,516 / −2,515 (**0,490** ; 0,106 / 0,088) | −1,803 ×2 (**0,100 / 0,101**) ; −1,312 (**0,232** ; 0,083) |

Conséquences sur le texte du §D4 : la phrase « (a1) → (a2) : … la paire paire à −2,812 passe sous le seuil (w₂ 0,287 → 0,058 / 0,074) et
l'impair −1,350 passe de 0,226 à 0,084 » est caduque (w₂ identiques à 3e-4) ; « (a2) → (a3) : … toujours sous le seuil » et
« (a3) → (b) : … w₂ 0,067 / 0,055 → 0,106 / 0,088 » sont à lire avec les valeurs corrigées (0,254 → 0,489 pour la paire σ, 0,225 → 0,232
pour l'impair) ; les lignes « aucun état localisé » des marches (a1)→(a2), (a2)→(a3), (a3)→(b) pour les pairs sont fausses : la paire σ
est localisée dans toutes les variantes. Les six états localisés par variante sont donc les mêmes de (a1) à (c-3), à ≤ 0,3 eV près en
énergie. Un second constat de R5 (A.2) concerne la normalisation relative de M^L et M^NL dans les fichiers M de production ; il est
rapporté dans `R5_rapport.md` et n'est pas repris ici.
