# R6 — Production corrigée : normalisation de la partie locale de M — rapport de campagne

Statut PRODUCTION. Prompt R6 (Greg, 2026-09-25). Ordre : phase 0 → STOP → étape G (Greg corrige le noyau, commit) → (GO 1)
assemblage et portes → STOP → (GO 2) validation super-cellule → STOP → (GO 3) production du chapitre 4 → rapport → STOP.
Dépôt `graphene-raman`, HEAD f4b7ec3 au moment de la phase 0 ; racine `$GRAPHENE_RAMAN`. Chiffres bruts, sans interprétation.

## Phase 0 — Inventaire (2026-09-25 ; rien n'est calculé, rien n'est modifié dans le dépôt ni dans `results/`)

Seules opérations faites sur les données : lectures de sidecars, de logs, de XML, et trois différences numériques sur les M
grossiers (nœud de connexion, < 200 Mo, § 0.1 c). Aucun job, aucun fichier écrit hors de ce répertoire.

### 0.1 Inventaire des M

#### a) Grossiers (maille N×N, 16 bandes, N² k ; `results/M/`, Hartree, sidecar `bloch_norm = unit_cell`)

| N | M_ed (mtime) | M_L | M_NL | sidecar (mtime, note) | chaîne qui l'a produit | L + NL vérifiable ? |
|---|---|---|---|---|---|---|
| 5×5 | `M_ed_5x5.npy` (16, 25)² 2,56 Mo, 2026-08-31 16:01 | absent | absent | 09-07 16:11, « raw, pre-migration » | juin (job 14411452, 18 juin 15:50, pilote d'époque `backend=qe, method=real`, M^L + M^NL dans un seul processus, seule la somme sauvée) | **non** (aucune partie sur disque) |
| 6×6 | `M_ed_6x6.npy` (16, 36)² 5,3 Mo, 08-31 16:01 | absent | absent | 09-07 16:45, « raw, pre-migration » | juin (même pilote, log absent de `logs/` ; `M_ed_1_14411649`) | **non** |
| 7×7 | `M_ed_7x7.npy` (16, 49)² 9,8 Mo, 09-07 16:04 | `M_L_7x7.npy` 09-07 15:55 | `M_NL_7x7.npy` 09-07 15:56 | 09-07 16:57, « M^L + M^NL, unit-cell Bloch norm ; M^L from the resampled-grid kernel (217 = 7×31) », `N_cells = 49` | `compute_M.py` (job 20435732 : ml 60,6 s sur 32 rangs, nl 64,3 s, combine 0,3 s). Le combine du 07 15:56 a écrit (L+NL)/49 en `supercell` ; le fichier a été réécrit à 16:04 comme L + NL (unit_cell) | oui : max\|M_ed − (M_L + M_NL)\| / max\|M_ed\| = **0,0** ; (L+NL)/49 : écart 0,98 |
| 8×8 | `M_ed_8x8.npy` (16, 64)² 16,8 Mo, 08-31 16:01 | absent | absent | 09-07 16:11, « raw, pre-migration » | juin (`M_ed_3_144…`) | **non** |
| 9×9 | `M_ed_9x9.npy` (16, 81)² 26,9 Mo, 08-31 16:01 | `M_L_9x9.npy` 06-18 17:29 | `_test_mnl/M_NL_serial.npy` (06-19 16:52, (16, 81)², test MPI vs série PASS 2,8e-14) ; `M_NL_mpi.npy` | 09-07 16:11, « raw, pre-migration » | juin en deux étapes (job 14415762, `M_ed2_0` : M^L MPI 32 rangs 3 min → `M_L_9x9.npy` ; M^NL série + combine 3 min) | oui : max\|M_ed − (M_L + M_NL_serial)\| / max = **1,35e-13** (M_NL_mpi : 2,4e-13) ; M_NL = M_ed − M_L exact |
| 10×10 | `M_ed_10x10.npy` (16, 100)² 41 Mo, 08-31 16:01 | `M_L_10x10.npy` 06-18 17:31 | absent | 09-07 16:11 | juin deux étapes (`M_ed2_1_14415762`) | par différence : M_ed − M_L hermitien à 7,9e-15 |
| 11×11 | `M_ed_11x11.npy` (16, 121)² 60 Mo, 08-31 16:01 | `M_L_11x11.npy` 06-18 17:35 | absent | 09-07 16:11 | juin deux étapes (`M_ed2_2_14415762`) ; **grille FFT 360 non multiple de 11 (11 × 30 = 330), noyau de juin sans rééchantillonnage** (même défaut que l'ancien 7×7 sur 216, corrigé le 07-09 par `prep_realspace_inputs`) | par différence (hermitien 1,4e-14) ; M_L suspect |
| 12×12 | `M_ed_12x12.npy` (16, 144)² 85 Mo, 08-31 16:01 | `M_L_12x12.npy` 06-18 17:39 | absent | 09-07 18:12 | juin deux étapes (`M_ed2_3_14415762`, M^L 13 min, total 21 min) | par différence (hermitien 1,4e-14) |

Grandeurs brutes (Ha) : max\|M_ed\| 5×5 0,2820, 6×6 0,2539, 7×7 0,2589, 8×8 0,2613, 9×9 0,2628, 10×10 0,2644, 11×11 0,2667,
12×12 0,2659 ; max\|M_L\| 7×7 1,306e-2, 9×9 8,547e-3, 10×10 6,460e-3, 11×11 9,546e-3, 12×12 4,942e-3 ; max\|M_NL\| 7×7 0,2698,
9×9 0,2698, et max\|M_ed − M_L\| 0,2698 pour 10, 11, 12. Produits max\|M_L\|·N_cells : 7×7 0,640, 9×9 0,692, 10×10 0,646,
11×11 1,155, 12×12 0,712. Hermiticité max\|M − M†\| : M_ed 2,8e-17 (5, 6, 8), 3,7e-15 (7), 5,4e-15 (9), 7,9e-15 (10), 1,4e-14 (11, 12).
Le `M_L_dense_9x9_coarsecheck.npy` (noyau `compute_ML_R_mpi_shared` sur la grille grossière, job 20203627) reproduit `M_L_9x9.npy` (juin,
`compute_ML_R_mpi`) à 1,3e-15 relatif : les deux noyaux MPI portent la même normalisation.

Fichiers `M_ed_*_norm.npy` et `M_dense_*_norm.npy` (copies /N_cells en `supercell`, `migrate_M_norm.py`) : **absents** (supprimés
au ménage de septembre) ; `compute_spectral.py`, `compute_tmatrix.py`, `_eta_scan.py` qui les lisent sont hors chaîne.

#### b) Denses (maille dense D×D = p·N, 20 bandes ; `results/M/`, Hartree, sidecars `unit_cell`, `part`, `p`, `D`, `kernel = R`)

| N | p, D, N_kd | `M_L_dense` (mtime, job) | `M_NL_dense` | `M_dense` = L + NL (mtime) | taille par fichier |
|---|---|---|---|---|---|
| 5×5 | 5, 25, 625 | 09-04 12:36 (20212556, ml 32 rangs × 6 fils) | 09-04 12:44 | 09-04 12:45 | 2,50 Go |
| 6×6 | 4, 24, 576 | 09-07 16:31 (20435993 : u_nk 31,9 Go/nœud, 14,2 min de blocs, 27 min total) | 09-07 16:42 | 09-07 16:42 | 2,12 Go |
| 7×7 | 4, 28, 784 | 09-07 16:35 (20435537, 1 h 02 ; grille 216 rééchantillonnée 217) | 09-07 16:56 | 09-07 16:56 | 3,93 Go |
| 8×8 | 4, 32, 1024 | 09-04 13:56 (20212560) | 09-04 14:32 | 09-04 14:32 | 6,71 Go |
| 9×9 | 3, 27, 729 | 09-04 13:15 (20212562) | 09-04 13:47 | 09-04 13:47 | 3,40 Go |
| 12×12 | 2, 24, 576 | 09-07 17:27 (20435997, 1 h 41) | 09-07 18:09 | 09-07 18:09 | 2,12 Go |

Total dense : 3 × 20,8 Go = 62,4 Go. Pas de dense pour 10×10 ni 11×11 (pas de `.save` dense, pas de wannierisation). Chaîne :
`scripts/compute_M_dense_stages.py` (`submit_M_dense.sh`, 1 nœud exclusif 32 × 6) : stage `ml` = `local_R.compute_ML_R_mpi_shared`
(l. 47–51, `uc_dense`, potentiel de la super-cellule N×N, `Omega_sc` de la N×N), stage `nl` = `non_local.compute_M_NL` sur le `.save`
dense (l. 73–79), stage `combine` (l. 82–92) : `M = M_L + M_NL`, sidecar `unit_cell` + `N_kd` ; `--out-norm` (M/N_kd, `supercell`)
écrit puis supprimé au ménage.

**Facteur de correction du dense** : le noyau reçoit `Omega_sc` de la super-cellule N×N (`P["scp"]`), donc M^L_dense est en norme
« super-cellule N×N » quel que soit p ; le facteur est **N_cells = N²** (81 pour la 9×9 à D = 27, vérifié par R5-A.2 sur
`M_dense_9x9`, base 20 : rapport 81,0000), et non N_kd = D².

#### c) Autres M

- R5, 128 bandes × 81 k (9×9 grossière, `.save` neuf `qe_tmp/R5_uc9x9_nb128/`) : `R5_base_vs_M/b/M_{L,NL,ed}_9x9_nb128.npy`
  3 × 1,72 Go, sidecars `unit_cell`, chaîne de production grossière (`compute_M.py` ml/nl, combine unit_cell) ; même facteur 81 sur
  M_L (R5 B.2). Entrée de l'étape 2 (c-all à 128 bandes).
- `_ML_7x7_ref_b4.npy`, `M_L_dense_7x7_coarsecheck.npy` (4 bandes, validation grille), `_test_mnl/` : contrôles, pas des entrées.

#### d) Où l'assemblage L + NL a lieu

| endroit | ligne | ce qu'il écrit | convention |
|---|---|---|---|
| `scripts/compute_M.py` `stage_combine` | l. 144–159 | `M = (M_L + M_NL) / N_cells`, sidecar **`supercell`** | ≠ production (les `M_ed_*` de production sont L + NL en `unit_cell`) ; le 7×7 a été réécrit à la main le 07-09 |
| `scripts/compute_M_dense_stages.py` `stage_combine` | l. 82–92 | `M = M_L + M_NL`, sidecar `unit_cell`, `N_kd` ; `--out-norm` = M/N_kd `supercell` | production |
| `local_G.compute_M_dense` | l. 479–505 | `M_L + M_NL` (monolithique, non utilisé : corruption de tas nk ≥ 81) | — |
| pilote de juin (`backend=qe, method=real`) | absent du dépôt (`compute_M_cluster.py` n'existe plus) | somme dans un processus (5, 6, 8) ou M_L sauvé puis somme (9–12) | `unit_cell` |
| R5 `r5_driver.py` (b2sum) | — | `M_ed_9x9_nb128 = M_L + M_NL` | `unit_cell` |

Facteur 1/N_cells en aval (jamais dans les fichiers) : `wannier/supercell_fold.py:34–40` (`bloch_folded_hamiltonian`, H = diag ε + M/N_cells),
l. 59–62 (`kbasis_matrix`), `scripts/test_local_tmatrix_real.py:42–48` (`compute_T(…, Mbk5 / Nc)` puis × Nc) ; la matrice t locale
(`local_tmatrix.py`) travaille avec V_loc intensif tiré de M `unit_cell` (contrat commenté dans `compute_spectral_wannier.py:59–64`).

### 0.2 Le noyau local

#### Lignes où la normalisation entre (`src/electron_defect_interaction/defects/local_R.py`, HEAD f4b7ec3)

| noyau | lignes | contenu |
|---|---|---|
| `compute_ML_R` (série) | 49, 65, 72 | `A_sc, Omega_sc = io.get_A_volume(sc_wfk_path)` ; `psi = compute_psi_nk_fold_sc(…, Omega_sc, Ndiag, ngfft, …)` ; `dV = Omega_sc / prod(ngfft)`. La division par √Ω est dans `wavefunctions/fold_wfk_to_sc.py:93` : `psi = (u * phase_k) / np.sqrt(float(Omega_sc))` |
| `compute_ML_R_mpi` | 187, 198–199, 220, 224 | `Omega_sc = float(prep["Omega_sc"])` ; `dV = Omega_sc / Ntot` ; **`inv_sqrtO = 1.0 / np.sqrt(Omega_sc)`** ; `Psi = (u[:, :, pu] * phase) * inv_sqrtO` ; `M_local += dV * (Psi.conj() * Vb) @ Psi.T` |
| `compute_ML_R_mpi_shared` | 294, 321–322, 339, 342 | idem : `Omega_sc = meta["Omega_sc"]` ; `dV = Omega_sc / Ntot` ; **`inv_sqrtO = 1.0 / np.sqrt(Omega_sc)`** ; `Psi = … * inv_sqrtO` |
| `_build_u_uc` | 155–170 | u_nk = ifftn(C) · N_uc sur la grille de maille : Σ_{r∈maille} \|u\|²/N_uc = Σ_G \|C\|² = 1 (« no k-phase, no 1/sqrt(Omega) ») |
| `local_G.prep_reciprocal_inputs` (noyau réciproque, `--kernel G`, **non utilisé en production**, sidecars `kernel = R`) | 40, 58 | `A_sc, Omega_sc = …` ; **`Ved_G = np.fft.fftn(Ved) * (1 / np.prod(ngfft))`** = (1/Ω_sc) ∫ V_ed e^{−iG·r} : même norme super-cellule ; `zero_pad_potential` (l. 294–335) et `_reconstruct_Ved` (l. 338–340) propagent ce coefficient |

Bilan par élément (grille N_tot = N_cells·N_uc points, dV = Ω_sc/N_tot = Ω_uc/N_uc) : M^L_actuel = dV Σ_r ψ* ΔV ψ avec ψ = u e^{ik·r}/√Ω_sc
= (1/N_tot) Σ_{r∈SC} u* ΔV u, et ∫_SC \|ψ\|² = 1 (norme super-cellule). En norme maille, ψ = u e^{ik·r}/√Ω_uc, ∫_maille \|ψ\|² = 1,
et ⟨ψ\|ΔV\|ψ′⟩_SC = (1/Ω_uc) ∫_SC u* ΔV u′ = N_cells × M^L_actuel. C'est le rapport 81,0000 mesuré par R5-A.2 (toutes les paires,
bases 16 et 20).

Noyau non local (`src/electron_defect_interaction/defects/non_local.py`) : l. 148 `A_uc, Omega_uc = io.get_A_volume(uc_wfk_path)`,
l. 186–187 `pref = (4*np.pi / np.sqrt(Omega_uc))` ; MPI l. 257, 273 `pref = 4 * np.pi / np.sqrt(Omega_uc)` ; B = pref Σ_G C* F Y e^{−iK·τ},
M^NL = Σ D_li B B* : (4π)²/Ω_uc, norme maille. R5-A.2 : direct = M^NL/81 à 3,4e-15 eV.

#### Changement minimal proposé (diff `phase0/local_R.py.diff`, 17 lignes, NON appliqué)

- `compute_ML_R` : `N_cells = prod(Ndiag)` ; `Omega_uc = Omega_sc / N_cells` ; appel `compute_psi_nk_fold_sc(…, Omega_uc, …)` (la
  routine ne se sert de son argument Ω que pour le 1/√Ω, l. 93, et pour `check_normalize`, désactivé par défaut ; `dV` reste Ω_sc/N_tot).
- `compute_ML_R_mpi` et `compute_ML_R_mpi_shared` : `inv_sqrtO = 1.0 / np.sqrt(Omega_sc / N_cells)` avec `N_cells = prod(Ndiag)`
  (resp. `prod(meta["Ndiag"])`) ; `dV` inchangé. Docstrings : « / sqrt(Omega_uc), unit-cell norm ».
- Effet exact : M^L_nouveau = N_cells × M^L_ancien pour tout potentiel (√N_cells est exact en flottant pour N entier :
  √81 = 9 ; écart attendu ≤ 1e-15 relatif, cible 1e-12 de l'étape 1.1).
- Noyau G (`phase0/local_G.py.diff`, 6 lignes) : `Ved_G = fftn(Ved) * (N_cells / prod(ngfft))` (Ndiag calculé avant), pour que
  `test_zero_pad_dense.py` / `test_pad_vs_full_supercell.py` (G contre G, invariants) et un éventuel `--kernel G` restent cohérents
  avec R. Optionnel : ce noyau n'a produit aucun M de production.
- Inchangés : `fold_wfk_to_sc.py` (paramètre renommé sémantiquement, pas de code), `wfk.compute_psi_nk` (maille), `non_local.py`.

#### Pourquoi la reconstruction KS ne pouvait pas voir le facteur

1. `scripts/test_ks_reconstruction.py:79–137` (`reconstruct_ks_hamiltonian`) et `scripts/ks_reconstruction_all.py:49–73` (`reconstruct`)
   n'appellent pas `compute_ML_R` : le terme local est calculé **sur la grille de la maille**, avec `compute_psi_nk(…, Om)` où
   `Om = Ω_uc` et `dV = Ω_uc/N_r` (l. 123–130 et l. 67–71). Le potentiel de la super-cellule parfaite y est **restreint** à la maille
   (`sc_pot_on_uc_grid`, restriction directe ou de Fourier). Ω_sc n'y entre jamais : la norme maille est la seule présente.
2. Même un test sur la super-cellule avec un potentiel périodique V_p ne peut pas le voir : pour V_p périodique et k = k′,
   Σ_{r∈SC} u* V_p u = N_cells Σ_{r∈maille} u* V_p u, donc M^L_actuel(V_p) = (N_cells/N_tot) Σ_maille u* V_p u = (1/N_uc) Σ_maille u* V_p u =
   l'élément de matrice de la maille : la normalisation 1/Ω_sc est compensée exactement par les N_cells copies identiques. Le facteur
   n'apparaît que pour un ΔV **localisé** dans une seule maille (une seule copie : (1/N_tot) Σ = (1/N_cells)(1/N_uc) Σ_maille).
   Le test nul (`null_test`, l. 138–149 : `compute_ML_R(uc, sc_p, pot_p, pot_p)`) donne 0 quelle que soit la norme.
3. Conséquence pour le test 1.1 « M^L de V_p inchangé par la correction » : il est vrai pour la reconstruction KS sur la maille
   (Ndiag = (1, 1, 1), N_cells = 1 : la correction est l'identité) — c'est le test proposé (`ks_reconstruction_all.py` rejoué, écarts
   attendus identiques à `results/M/ks_reconstruction.npz` à 1e-12). Il est **faux** pour un `compute_ML_R(uc, sc_p, pot_p, ·, pristine=True)`
   sur la super-cellule : le noyau corrigé y donne N_cells × l'élément de maille (M^L_nouveau(V_p) = N_cells × M^L_ancien(V_p), par
   construction, comme M^NL_p de la super-cellule = N_cells × celui de la maille). À préciser par Greg avant GO 1 (§ 0.5, décision 3).

### 0.3 Tout ce qui lit M en aval

Convention : les M2 gardent `bloch_norm = unit_cell` (désormais vrai pour les deux parties) et reçoivent la clé sidecar
`M_normalization = "v2 …"`. Pour chaque lecteur, le seul changement de physique est le fichier lu ; les changements de code sont
(i) le chemin `results/M` → `cfg["results_dir"]` (= `results/M2`), (ii) la porte de version `require_normalization = "v2"` de
`load_M_checked` (diff `phase0/matrix_io.py.diff`), (iii) rien d'autre. `config.dense_paths` (`config.py:22–29`) porte le chemin dense
pour dix scripts (diff `phase0/config.py.diff` : `results_dir(cfg)` + clés obligatoires).

| script | entrée M / dérivée (chemin actuel) | sortie | changement R6 |
|---|---|---|---|
| `compute_spectral_wannier.py` (niveau 1 : cartes R_cut × grille × η ; `submit_spectral_wannier_dense.sh`) | `dense_paths → M_dense_{S}` ; coarse : `results/M/M_ed_{S}` (l. 55–66) | `specwd_{S}_prod.npz` | chemin + porte v2 ; sorties dans `results/M2/` |
| `resonance_metrics.py` (Γ_T, Γ_Born, δρ, ρ, T̄) | `dense_paths` + `M_L_dense_{S}` (l. 24–27, moyenne diagonale de M^L = décalage G≈0, 67,00 meV en v1 → × N_cells en v2) | `resonance_{S}.npz` | chemin + porte v2 ; la variante « M − ⟨M^L⟩·1 » (l. 42, 61) suit automatiquement ; 3.5 (C = décalage Lu ajouté à ΔV^L) est une nouvelle option, pas un remplacement |
| `resonance_criteria.py` (det, λ, Friedel, Lloyd, Γ à c = 0,1 %) | `dense_paths` (l. 25) ; lit `resonance_{S}.npz` (l. 81) | `resonance_criteria_{S}.npz` | chemin + porte v2 |
| `rcut_resigma.py` (Re Σ par R_cut ; `--nk-int`) | `dense_paths` (l. 24) | `resigma_{S}_rc*.npz`, `..._nk*.npz` | chemin + porte v2 |
| `m_rcut_convergence.py` (tab:rcut_M) | `dense_paths` (l. 20) | `m_rcut_convergence.csv` (l. 56) | chemin + porte v2 |
| `mwr_locality_coarse_vs_dense.py` (localité de M_W, coarse vs dense) | `dense_paths` + `results/M/M_ed_{S}` (l. 30) pour `wannier/{S}` ∈ {5x5, 7x7, 8x8} | `mwr_locality.npz` | chemin + porte v2 ; **exige les M2 grossiers 5×5, 7×7, 8×8** |
| `analyze_M.py` (§4.1.5 : cartes M̃, L/NL, échelle, tests) | `dense_paths` ; littéraux `M_ed_{S}`, `M_dense_{S}_nb16` (optionnel, absents), `M_L_dense_{S}_coarsecheck`, `M_L_{S}` (l. 61, 70, 109, 115, 134, 141–144) | `M_analysis.npz`, `M_tests_summary.csv` | chemin + porte v2 ; **exige les M2 grossiers de toutes les tailles denses** (test « convention intensive », l. 115–119) et un `coarsecheck` v2 pour 9×9 (recalcul 1.1) |
| `lnl_frobenius_all.py` (tab:L_NL) | `dense_paths` + `M_L_dense`/`M_NL_dense` (l. 31) | `lnl_frobenius.csv` | chemin + porte v2 ; ratio ⟨\|M^NL\|_F⟩/⟨\|M^L\|_F⟩ divisé par N_cells (v1 : 6,4 → 36,3 pour N = 5 → 12) |
| `_bz_ratio_LNL.py`, `check_onsite_and_NL.py`, `check_M_dense_vs_coarse.py`, `check_ML_coarse_kernel.py`, `check_M_dense_nb20_vs_nb16.py`, `validate_ML_grid_7x7.py`, `_old_vs_new_7x7.py`, `_mcheck.py`, `_normtest.py` | littéraux `results/M/…` | stdout / csv | contrôles ponctuels : chemin (`_nb16`, `obsolete_grid_7x7` : fichiers supprimés) |
| `test_local_tmatrix_real.py` (**test d'or**, bloquant ; `submit_golden_dense.sh` 5×5 dense, 3 h 36) | `results/M/M_dense_{N}` ou `M_ed_{N}` (l. 25, 27) | stdout PASS/FAIL | chemin + porte v2 ; à rejouer avec M2 avant tout GO 3 (compare local et dense sur le même M2 : pas de valeur de référence stockée) |
| `level2_families.py` | `mwr_locality.npz`, `M_analysis.npz`, `specwd_{S}_prod.npz` | `level2_families.csv` | chemin des dérivées |
| `make_figures.py`, `make_figures_memoire.py` (figures du mémoire ch. 4 : `fig_spectral`, `fig_plateau`, `fig_rcut`, `fig_level2`, `fig_locality`, `fig_M_map`, `fig_M_scaling`, `fig_Ved_*`, `fig_ks_reconstruction`, versions `_final`) | `specwd_*`, `resonance_{REF}`, `resonance_criteria_{REF}`, `mwr_locality`, `M_analysis`, `ks_reconstruction`, `ved_analysis` | `figures/*.pdf, .png`, `level1_summary.csv`, `level2_summary.csv` | chemin des dérivées (`--analysis`, `--resonance`, `--locality` ont déjà des options ; défauts à pointer sur `results/M2`) ; `ks_reconstruction.npz` et `ved_analysis.npz` inchangés (ne lisent pas M) |
| `make_figures_epw.py` (**chapitre 5**, `fig_epw_vs_ed` : Γ^ed contre Γ^ep) et `epw_ed_vs_ep.py` | `results/M/resonance_9x9.npz` (l. 71–72 ; `--resonance`) | `fig_epw_vs_ed*.pdf` | chemin : la figure du chapitre 5 dépend de M2 |
| `nkint_check_post.py` | `resigma_9x9_rc3_nk*.npz` (`--pattern`) | `nkint_check_9x9.csv` | chemin |
| `summarize_level1_maps.py` | `results/M/logs/specwd_{J}.out` (logs !) | stdout | chemin des logs (`results/M2/logs/`) |
| `sampling_table.py`, `ks_reconstruction_all.py`, `analyze_Ved.py`, `compute_convergence.py` | ne lisent pas M (XML, `.save`, `Vks`) ; `compute_convergence` lit `gamma_*.npz` (obsolètes) | — | aucun (sorties recopiées ou liées dans `results/M2/` pour les figures) |
| `compute_spectral.py`, `compute_tmatrix.py`, `_eta_scan.py`, `migrate_M_norm.py` | `M_ed_{S}_norm.npy` (absents) | — | hors chaîne (à laisser) |
| `submit_M.sh`, `submit_M_dense.sh`, `submit_spectral_wannier*.sh`, `submit_rcut_resigma.sh`, `submit_nkint_check.sh`, `submit_golden_dense.sh`, `submit_lnl_frobenius.sh`, `submit_tmatrix.sh`, `submit_spectral.sh` | `--out results/M/…`, `#SBATCH --output=results/M/logs/…` | | chemins → `results/M2/` (variable `RES=$(python -c 'print(cfg["results_dir"])')` ou littéral) |
| `src/…/config.py:27` `dense_paths` | `results/M/M_dense_{size}.npy` codé en dur | | `results_dir(cfg)` (diff) |
| `src/…/io/matrix_io.py` `load_M_checked` | sidecar | | clé `M_normalization`, refus des v1 quand `require_normalization="v2"` (diff) |
| `scripts/compute_M.py` combine | — | `M_ed` en `supercell` /N_cells | passer en `unit_cell` L + NL + tag v2, comme la chaîne dense (diff `phase0/compute_M.py.diff`) |
| R4 `r4_driver.py`, R5 `r5_driver.py` (étape 2) | `results/M/M_ed_9x9`, `M_dense_9x9`, `M_L_*`, `M_NL_*`, `b/M_*_nb128` | tableaux D4 | copies R6 pointant sur `results/M2/` et sur les nb128 réassemblés |
| `NOTES_TGAMMA.md` §2, `results/M/*.csv` | valeurs v1 citées | | table de correspondance ancien → nouveau (3.7) ; **`défauts.tex` n'est pas sur Rorqual** (`memoire/` du dépôt ne contient que `EM1_tb/` ; aucun `.tex` du mémoire dans `$PROJECTS/gregb26`) : la colonne « figure ou tableau du mémoire » viendra des labels cités (`tab:rcut_M`, `tab:L_NL`, `tab:tests_M`) et de NOTES_TGAMMA, sauf copie du tex fournie par Greg |

Fichiers dérivés suivis par git (`.gitignore` l. 33–43 : `results/M/specwd_*_prod.npz`, `resonance_*.npz`, `mwr_locality.npz`, `M_analysis.npz`,
`ks_reconstruction.npz`, `ved_analysis.npz`, `*.csv` ; 49 fichiers) : `results/M2/*` est ignoré par la règle `results/*` ; il faudra les mêmes
exceptions pour `results/M2/` (diff `.gitignore` à ajouter à l'étape 3, Greg).

### 0.4 Fichiers de super-cellule pour la porte A.2, par taille, et coût

Tous les `.save` de super-cellule sont sur le scratch (`qe_tmp/defect_{N}_{p,d}/defect_{N}_{p,d}.save/`, `gamma_only = true`, un seul
`wfc1.dat`, ecutwfc 50 Ha) et miroités dans `graphene/qe/qe_tmp_backup/` (`MD5SUMS_2026-09-17.txt`). Les `.save` denses de maille
(`defect_uc_dense_{24,25,27,28,32}`, 20 bandes, 1,7–2,9 Go) et grossiers (16 bandes, grille 30×30×192) y sont aussi.

| N | grille FFT SC | N_grid | npwx (demi-sphère) / sphère complète 2·npwx − 1 | `.save` (p, d) | grille = N × 30 ? | états purs : C de |
|---|---|---|---|---|---|---|
| 5×5 | 150×150×192 | 4 320 000 | 119 010 / 238 019 | 253 Mo × 2 | oui | grossier 5×5 (25 k) ; dense 25×25 (625 k, 25 k coïncidents) |
| 6×6 | 180×180×192 | 6 220 800 | 171 510 / 343 019 | 480 / 474 Mo | oui | grossier ; dense 24×24 (36 coïncidents) |
| 7×7 | 216×216×192 | 8 957 952 | 233 409 / 466 817 | 838 / 831 Mo | **non** (noyau sur 217×217×192, ΔV rééchantillonné) → la porte doit poser Ψ et ΔV sur la même grille 217 (`fourier_resample`) | grossier ; dense 28×28 (49 coïncidents) |
| 8×8 | 240×240×192 | 11 059 200 | 304 854 / 609 707 | 1,4 Go × 2 | oui | grossier ; dense 32×32 (64 coïncidents) |
| 9×9 | 270×270×192 | 14 002 200 | 385 710 / 771 419 | 2,1 Go × 2 | oui | grossier ; dense 27×27 (81) ; nb128 (81) ; 30 états QE (R5) |
| 10×10 | 300×300×192 | 17 280 000 | 476 289 / 952 577 | 3,2 Go × 2 | oui | grossier seulement |
| 11×11 | 360×360×192 | 24 883 200 | 576 360 / 1 152 719 | 4,5 Go × 2 | **non** (11 × 30 = 330 ; noyau actuel rééchantillonnerait sur 363) | grossier seulement (M_L de juin suspect, § 0.1) |
| 12×12 | 360×360×192 | 24 883 200 | 685 905 / 1 371 809 | 6,3 Go × 2 | oui | grossier ; dense 24×24 (144 coïncidents) |

Porte sur états purs (c = δ) : elle n'a besoin d'aucune fonction d'onde de super-cellule (seulement Ω_sc, la grille, ΔV = V_d − V_p des
`Vks_*`, la position τ de l'atome retiré via `alignment.vacancy_site`, et la sphère complète de la super-cellule pour le projecteur KB) ;
les `wfc1.dat` servent (i) aux 30 états QE de la 9×9 (1.2) et (ii) à la parité. Pour le M2 dense, les états purs testables sur la grille
N×N sont ceux des k **coïncidents** avec la grille grossière (k = i/N, N² parmi D²), comme dans R5 (bug 0.3 : prendre les k du `.save`
qui fournit C) ; une paire à k non coïncident demanderait la cellule zéro-padded p·N (grille p² fois plus grande : 9×9 126 M points,
8×8 177 M points, 2,8 Go par tableau complexe) — option, pas nécessaire pour le facteur (global).

Coût mesuré (R5, J1 21811069, 16 cœurs, 9×9) : grille 14,0 M points, ΔV lu et projecteur KB (771 419 g) en 1 s, 6 paires base 16 en 2 s,
7 paires base 20 en 1 s, 30 états en 10 s. Mémoire : 3 tableaux complexes N_grid (12×12 : 3 × 398 Mo), W du projecteur
(nfull × 12 complexes : 12×12 263 Mo), C dense ≤ 2,9 Go. Estimation par taille : < 5 min (lecture des `.save` denses comprise), < 16 Go.
Une seule tâche pour les 8 tailles × (grossier, dense) : 16 cœurs, 64 Go, ≤ 1 h.

### 0.5 Diffs proposés, plan, jobs, coûts (rien n'est appliqué)

- `phase0/production.json.diff` : clés `M_normalization = "v2 : L et NL en norme unit_cell, 2026-09-25"`, `M_normalization_note`,
  `results_dir = "results/M2"`, `results_dir_frozen = "results/M"` (insérées après `_comment`, rien d'autre ne change).
- `phase0/config.py.diff` : `load_production` exige `M_normalization` et `results_dir` (refus sinon), les affiche ; `results_dir(cfg)` ;
  `dense_paths` lit `results_dir(cfg)`.
- `phase0/matrix_io.py.diff` : constantes `M_NORM_V2`/`M_NORM_V1`, paramètre `require_normalization` de `load_M_checked` (un sidecar
  sans la clé est v1 → refus quand v2 est exigé ; règle « refuser, jamais avertir »).
- `phase0/compute_M.py.diff` : stage combine en `unit_cell` L + NL + tag v2 (convention de production) ; tags v2 sur les parties.
- `phase0/local_R.py.diff`, `phase0/local_G.py.diff` : noyau (étape G, Greg).
- `phase0/README_results_M_gele.md` : texte du README à déposer dans `results/M/` au GO 1 (« résultats obtenus avec M^L non normalisé
  (facteur N_cells manquant), remplacés par R6 »).
- Modules R5 promus au GO 1 (rôles nommés) : `r5_sc_projection.py` → `src/electron_defect_interaction/wavefunctions/sc_projection.py`
  (états Γ de super-cellule sur la grille FFT, sphère complétée, projection sur la base de Bloch repliée) ; `r5_deltav_pw.py` →
  `src/electron_defect_interaction/defects/deltav_pw.py` (application directe de ΔV^L sur la grille et de ΔV^NL par le projecteur KB de
  l'atome retiré ; `check_pure_bloch` = **porte A.2**). Nouveaux scripts : `scripts/assemble_M2.py` (M2 = N_cells·M^L + M^NL, sidecar
  `unit_cell`, `hartree`, `assembled_from`, `N_cells`, `M_normalization`, `date`, hermiticité, contrôle M2 − N_cells·M^L − M^NL = 0),
  `scripts/gate_M_normalization.py` (porte A.2 par taille, tableau N × (L, NL), refus si > 1e-6 eV), `scripts/submit_r6_*.sh`.

Plan et coûts (temps `sacct` des jobs de septembre, rrg-cotemich-ac) :

| étape | job | contenu | ressources, durée estimée | précédent mesuré |
|---|---|---|---|---|
| G | — | Greg : `local_R.py` (+ `local_G.py`), `production.json`, `config.py`, `matrix_io.py`, `compute_M.py` (ou Code pour les quatre derniers au GO 1, à dire) ; commit | — | — |
| 1.1 | J1 `r6_kernel` | M^L grossier 9×9 avec le noyau corrigé (`compute_M.py --stage ml`, sortie dans ce répertoire, pas dans `results/`) ; comparaison à 81 × `M_L_9x9.npy` (cible 1e-12) ; porte A.2 ; `ks_reconstruction_all.py` rejoué (V_p, maille) ; variante V_p super-cellule 4 bandes (attendu × 81) | 1 nœud exclusif 32 rangs, 20 min | M^L 9×9 grossier 3 min (juin, 32 rangs) ; 7×7 61 s |
| 1.1 bis | J1b `r6_coarse_L` (si décision 1 = oui) | M^L grossier 5×5, 6×6, 8×8 (noyau corrigé) ; M^NL = M_ed − M^L_nouveau/N_cells ; porte A.2 | même nœud, 10 min | 7×7 : 61 s |
| 1.2 | J2 `r6_assemble` | M2 grossiers (7, 9, 10, 12 ; + 5, 6, 8 selon J1b ; 11 selon décision 2) et denses (6 tailles) → `results/M2/` : `M_L_dense_*` = N_cells·M^L (nouveau fichier), `M_NL_dense_*` = lien symbolique vers `results/M/` (inchangé) ou copie, `M_dense_*` = somme ; sidecars ; hermiticité ; M2 − N_cells·M^L − M^NL = 0 ; nb128 de R5 ; 30 états QE (identité c†Hc, s = −Lu) | 16 cœurs, 96 Go, 1 h (lecture 41,6 Go, écriture 41,6 Go par mmap et blocs) | — |
| 1.3 | J3 `r6_gate` | porte A.2 sur chaque M2 (grossier et dense, k coïncidents, paires intra-k et inter-k) → tableau N × (L, NL) | 16 cœurs, 64 Go, ≤ 1 h | R5 A.2 : secondes par base |
| 2 | J4 `r6_d4` | pilote R5 `a` adapté à M2 : (a1) 16 b, (a3) 20 b, (b) 5 WF, (c-all), (c-3), (c-all) 128 b projeté ; portes de repliement R4 ; états localisés par parité contre QE | 16 cœurs, 96 Go, 30 min | R5 J1 2 min 10, J5 21 min |
| 2 | J5 `r6_crit9` | `resonance_criteria.py` 9×9 (det, λ par bloc, −Im T̄, pic) et `resonance_metrics.py` avec M2 | 16 cœurs, 64 Go, 2 h | crit9 1 h 34 ; reson9 7 min 39 |
| 2 | J6 `r6_golden` | test d'or 5×5 dense avec M2 (bloquant avant GO 3) | 16 cœurs, 120 Go, 4 h | 3 h 36 |
| 3.1 | J7 | `mwr_locality_coarse_vs_dense.py`, `m_rcut_convergence.py` (9×9, 12×12, R_cut 0…6) | 8 cœurs, 10 min | 45 s ; 5 min |
| 3.2 | J8 | `analyze_M.py`, `lnl_frobenius_all.py`, `_bz_ratio_LNL.py` | 16 cœurs, 15 min | 1 min 28 ; 2 min |
| 3.3 | J9 × 6 | niveau 1 : `compute_spectral_wannier.py --dense` par taille (36 combinaisons 5, 7, 8, 9 ; 16 pour 6, 12) ; `rcut_resigma.py` R_cut 0…4 et N_k^int 150…600 (9×9) ; niveau 2 (`level2_families.py`, `make_figures.py` → csv) | 16 cœurs, 64 Go ; 1 h 15 (9×9), 10–40 min (6, 12) ; resigma 13 + 9 min ; nkint 38 min | specwd 20294202 1 h 15 ; nkint 21337627 37 min 40 |
| 3.4 | J10 | `resonance_metrics.py` + `resonance_criteria.py` 6×6, 12×12 (3 h dépassées le 2026-09-07 : `reson_6x6` et `reson_12x12` TIMEOUT ; les npz présents datent du même jour) | 16 cœurs, 6 h | 3 h 00 (TIMEOUT) |
| 3.5 | J11 | C14 : `resonance_metrics.py --shift-L +25 meV` (option nouvelle : constante ajoutée à ΔV^L, R_cut 3) | 16 cœurs, 30 min | 7 min 39 |
| 3.6 | — | tests de la chaîne (tab:tests_M) : hermiticité, k coïncidents, nb16→20 (fichiers absents : ligne « non rejouable »), noyau p = 1, fermetures Fourier, convention intensive, test d'or, g₀ par lots ; **+ ligne porte A.2** | dans J8 / J6 | — |
| 3.7 | — | table ancien → nouveau (NOTES_TGAMMA §2, `level1_summary.csv`, `level2_*.csv`, `m_rcut_*.csv`, `nkint_check_9x9.csv`, `lnl_frobenius.csv`, `M_tests_summary.csv`, `M_analysis.npz`, `mwr_locality.npz`, `resonance_*`) | nœud de connexion | — |
| fin | — | figures (`make_figures*.py`, style du mémoire) ; manifestes : aucune suppression ; md5 + miroir des M2 denses (41,6 Go nouveaux sur `/project` : 39 To / 180 To utilisés) ; `.save` 128 b de R5 (1,5 Go, `qe_tmp/R5_uc9x9_nb128/`, pas encore miroité) → `qe_tmp_backup/R5_uc9x9_nb128/` avec md5 | | |

Total : ≈ 15–20 h de tâches à 16 cœurs, 1 nœud exclusif 30 min ; une journée civile avec soumissions parallèles.

Décisions attendues de Greg avant GO 1 :

1. **M2 grossiers 5×5, 6×6, 8×8** : aucune partie L/NL sur disque (juin, somme seule). Options : (a) recalculer M^L grossier avec le noyau
   corrigé (≈ 1 min par taille, J1b ; M^NL = M_ed − M^L_nouveau/N_cells, exact à 1e-13) ; (b) pas de M2 grossier pour ces tailles (alors
   `mwr_locality` coarse 5×5/8×8 et le test « convention intensive » d'`analyze_M` perdent ces tailles). Recommandation : (a).
2. **10×10 et 11×11** : grossiers seulement, aucun consommateur de production (pas de dense, pas de wannierisation, `sampling_table`
   ne lit pas M). 10×10 : réassemblage possible (parties par différence). 11×11 : `M_L_11x11.npy` de juin calculé sur une grille 360 non
   commensurable (11 × 30 = 330) par le noyau d'avant la correction du 07-09 (même défaut que l'ancien 7×7) ; max\|M_L\|·N_cells 1,155
   contre 0,64–0,71 pour les autres tailles. Options : recalcul de M^L 11×11 avec le noyau corrigé (rééchantillonnage 363, ≈ 5 min) ou
   exclusion. Recommandation : 10×10 réassemblé, 11×11 exclu et noté.
3. **Test « M^L de V_p inchangé »** (1.1) : sur la maille (reconstruction KS rejouée, identité exacte) ou sur la super-cellule (attendu
   × N_cells, § 0.2). Recommandation : les deux, avec l'attendu de chacun.
4. **Noyau G** (`local_G.py`) : appliquer le diff (cohérence R/G) ou laisser avec une note.
5. **Qui applique** `production.json`, `config.py`, `matrix_io.py`, `compute_M.py` : Greg (avec le noyau) ou Code au GO 1.
6. **Disque** : `M_NL_dense_*` en liens symboliques vers `results/M/` (+41,6 Go) ou copies (+62,4 Go).
7. **`.gitignore`** : exceptions `results/M2/` identiques à `results/M/` (l. 35–43) pour suivre les npz/csv de production.

**STOP — phase 0 terminée le 2026-09-25.** Rien n'est calculé, rien n'est appliqué ; attente de l'étape G puis du GO 1.

## Étape G — lecture du commit de Greg (ee051ec, 2026-09-25 16:00, « correction de la normalisation du noyau local »)

Le clone local était en retard d'un commit (Greg a poussé depuis une autre machine) : `git pull --ff-only` fait par Code (avance
rapide, aucune fusion, arbre propre ; seul acte git de la campagne). Diff du commit : `local_R.py` (+5/−2) et `fold_wfk_to_sc.py` (+1/−1).

| noyau | changement commité | lecture |
|---|---|---|
| `compute_ML_R_mpi` (l. 194, 201) | `N_cells = np.prod(Ndiag)` ; `inv_sqrtO = 1.0 / np.sqrt(Omega_sc/N_cells)` | conforme au diff proposé (1/√Ω_uc) ; c'est le noyau de `compute_M.py --stage ml` (tous les M^L grossiers) |
| `compute_ML_R_mpi_shared` (l. 322, 325) | mêmes deux lignes, placées après le `node.bcast(meta)` | **écart d'exécution attendu** : dans cette fonction `Ndiag` n'est défini que sur le rang 0 du nœud (l. 284, sous `if nrank == 0:`) ; les autres rangs ont `meta["Ndiag"]` seulement → `NameError: name 'Ndiag' is not defined` hors rang 0. Noyau de la chaîne dense (`compute_M_dense_stages.py`), non utilisé par l'étape 1 (aucun M dense n'est recalculé) ; testé dans J1 (`--coarse 5x5`, 2 bandes). Correction d'une ligne : `N_cells = np.prod(meta["Ndiag"])` |
| `compute_psi_nk_fold_sc` (`fold_wfk_to_sc.py:93`) | `psi = (u * phase_k) / np.sqrt(float(Omega_sc)/N)` | **écart** : dans cette routine `N = np.prod(ngfft)` (l. 51) est le nombre de points de grille de la super-cellule, pas N_cells ; ψ y est multipliée par √(N_grid/Ω_sc) au lieu de √(N_cells/Ω_sc), soit un facteur √N_uc = √(N_grid/N_cells) de trop (5×5 : N_uc = 30·30·192 = 172 800 → M^L série × 172 800 par rapport au noyau MPI). Seul le noyau **série** `compute_ML_R` en dépend (`scripts/run.py`, `test_ks_reconstruction.null_test` (0 quelle que soit la norme), `validate_ML_grid_7x7.py`) ; aucun M de production. Testé dans J1 (`serial`, 5×5, 2 bandes). Correction : `np.sqrt(float(Omega_sc) / float(np.prod(Ndiag)))` |

Aucun M de l'étape 1 ne passe par ces deux noyaux ; ils sont rapportés (écart → rapporter), non corrigés par Code.

## Étape 1 (GO 1 reçu le 2026-09-25) — Assemblage et portes

Décisions prises par Code faute de réponse explicite aux sept points de la phase 0 (toutes réversibles, tout est nouveau fichier) :
(1) M^L grossier recalculé avec le noyau corrigé pour 5×5, 6×6, 8×8 (parties absentes) ; (2) 10×10 réassemblé (parties par différence)
et 11×11 réassemblé avec un M^L recalculé (grille 363) et M^NL = M_ed − M_L(juin) ; (3) test V_p fait des deux façons (maille : reconstruction
KS rejouée ; super-cellule : attendu × N_cells) ; (4) noyau G non touché ; (5) Code a appliqué les diffs `matrix_io.py` (porte v2, compatible)
et `compute_M.py` (tag v2 sur les parties, combine en `unit_cell`) — `config.py` et `production.json` restent à appliquer par Greg avant l'étape 3
(couplés) ; (6) `M_NL_dense_*` en liens symboliques vers `results/M/` ; (7) `.gitignore` à Greg. Modules promus : `src/electron_defect_interaction/
wavefunctions/sc_projection.py` et `defects/deltav_pw.py` (en-tête de rôle, code inchangé). Scripts versionnés nouveaux : `scripts/assemble_M2.py`,
`scripts/gate_M_normalization.py`. `results/M/README.md` déposé (gel), `results/M2/` créé. Le paquet n'a pas d'installation éditable dans
`.venv` (import impossible sans `sys.path`) : les jobs préfixent `PYTHONPATH=$PROJ/src:$PYTHONPATH` (le préfixe garde h5py de scipy-stack).

Jobs (répertoire de travail, `JOBID`) : J1 `r6kernel` 21819922 (nœud exclusif 32 rangs : M^L grossier 9, 7, 5, 6, 8, 11 ; V_p ; série ; KS ;
partagé) → J2 `r6assemble` 21819923 (afterok J1) → J3 `r6gate` 21819938 et J4 `r6states` 21819939 (afterok J2).
Vérifications préalables sur le nœud de connexion : `assemble_M2.py` reproduit bit à bit 49·M_L + M_NL du 7×7 et M_NL(9×9) = M_ed − M_L égale
`_test_mnl/M_NL_serial.npy` à 1,3e-13 ; la porte sur les fichiers v1 gelés donne direct/(M^L/N_cells) = 49,000000 (7×7 grossier, grille 217)
et 25,000000 (5×5 dense, 25 k coïncidents), NL exact (≤ 2e-15 eV), verdict REFUSÉ comme attendu.

### Résultats de l'étape 1 (jobs du 2026-09-25 : J1 21819922 échoué au 11×11 puis J1b 21820489 COMPLETED 24 min ; J2b 21820490 COMPLETED 10 min 53 ; J3b 21820491 COMPLETED 3 min 03 ; J4b 21820492 COMPLETED 43 s)

Le premier J1 a calculé les M^L grossiers 9×9 (143 s, 32 rangs), 7×7 (62 s, grille rééchantillonnée 217), 5×5 (26 s), 6×6 (33 s), 8×8 (89 s) puis a
échoué sur le 11×11 : `ValueError: operands could not be broadcast together with shapes (16,121,50000) (1,181,50000)` — le `.save` de maille
11×11 (`data/graphene/unit_cell/qe/defect_11x11.save`) a été écrasé par un run `bands` (XML : `<calculation>bands`, 181 k du chemin
Γ-K-M-Γ ; 121 fichiers wfc de la grille encore présents). Aucun M^L 11×11 n'est recalculable sans nscf QE → **11×11 exclu de R6**
(aucun consommateur de production ; son M_L de juin est de toute façon sur une grille non commensurable). Chaîne resoumise sans 11×11.

#### 1.1 Noyau corrigé (J1b, `kernel/*.json`, journal `r6_log.txt`)

| test | résultat | attendu |
|---|---|---|
| M^L 9×9 grossier (16 b, 81 k), noyau `compute_ML_R_mpi` ee051ec, contre 81 × `M_L_9x9.npy` (juin) | max\|v2 − 81 v1\| / max\|81 v1\| = **8,97e-16** ; rapport v2/(81 v1) médian 1,000000 (min 1,0000, max 1,0000) ; hermiticité 2,3e-16 ; max\|v2\| 0,6923 Ha | ≤ 1e-12 : **OK** |
| idem contre 81 × `M_L_dense_9x9_coarsecheck.npy` (noyau partagé de septembre) | 1,36e-15 | — |
| M^L 7×7 grossier contre 49 × `M_L_7x7.npy` (sept., grille 217) | **8,72e-16** ; rapport médian 1,000000 ; hermiticité 2,2e-16 ; max 0,6400 Ha | OK |
| M^L 5×5, 6×6, 8×8 (parties absentes en v1) | max\|v2\| 0,6710 / 0,8025 / 0,6519 Ha (max/N_cells 2,68e-2 / 2,23e-2 / 1,02e-2 Ha) ; hermiticité 2,2e-16 | — |
| Porte A.2 sur le M^L 9×9 recalculé (NL = `_test_mnl/M_NL_serial.npy` de juin), 8 paires dont (3,K,3,K) L direct +0,155642 = M^L/81 | max\|Δ_L\| = **6,79e-15 eV**, max\|Δ_NL\| = 3,44e-15 eV ; rapports médians 1,000000 / 1,000000 | ≤ 1e-6 eV : **OK** |
| Potentiel périodique, maille : `ks_reconstruction_all.py --sizes 9x9` rejoué après ee051ec (sortie hors `results/`) contre `results/M/ks_reconstruction.npz` | 24 clés (9×9 grossier et dense), écart **0,00e+00** ; grossier max 0,0015 meV, dense max 0,4606 meV (identiques aux anciens) | inchangé : **OK** (la reconstruction ne passe pas par le noyau) |
| Potentiel périodique, super-cellule : `compute_ML_R_mpi` avec `pristine=True` (V_p, 9×9, 4 bandes, `subtract_mean=False`) contre ⟨u\|V_p\|u⟩ sur la grille de la maille (V_p restreint, périodicité 0,018 meV) | blocs k = k′ : max\|M^L(V_p) − 81 ⟨u\|V_p\|u⟩\| / max = **1,07e-7** (1,6e-5 Ha absolu ; sans le facteur 81 : 80,0) ; blocs k ≠ k′ : max 3,6e-6 Ha = 2,5e-8 × max\|M\| ; moyenne diagonale −3241,92 eV = 81 × (−40,02 eV) | × N_cells par construction (§0.2) : **conforme** |
| Noyau série `compute_ML_R` (`fold_wfk_to_sc.py:93` ee051ec), 5×5, bandes 0–1, contre le noyau MPI corrigé | rapport série/MPI = **172 800** (min = max = médiane) = N_uc = 30·30·192 ; max\|M\| 9,23e4 Ha | **écart confirmé** : diviseur √(Ω_sc/N_grid) au lieu de √(Ω_sc/N_cells) ; correction `np.sqrt(float(Omega_sc)/np.prod(Ndiag))` (Greg) |
| Noyau partagé `compute_ML_R_mpi_shared` (`--coarse 5x5`, 2 bandes, 32 rangs) | `UnboundLocalError: cannot access local variable 'Ndiag'` sur 31 rangs (tous sauf le rang 0 du nœud) ; le rang 0 attend au `Reduce`, étape tuée par le `timeout 20m` ; aucun fichier produit | **écart confirmé** : `N_cells = np.prod(meta["Ndiag"])` (Greg) ; ce noyau sera nécessaire pour tout futur M dense |

#### 1.2 Réassemblage M2 (J2b, `scripts/assemble_M2.py`, `assemble_summary.jsonl`, `results/M2/`)

M2 = N_cells·M^L + M^NL, N_cells = N² (dense compris), par blocs (memmap), vérifications relues sur disque. Pour les 15 assemblages :
max\|out_L − N_cells·M^L\| = 0, max\|out_NL − source\| = 0, max\|M2 − (out_L + out_NL)\| = 0 (exact au bit), hermiticité max\|M2 − M2†\|/max\|M2\| ≤ 2,0e-14.

| fichier (`results/M2/`) | N_cells | source de M^L | source de M^NL | max\|N·M^L\| (Ha) | max\|M^NL\| | max\|M2\| | herm. rel. | durée |
|---|---|---|---|---|---|---|---|---|
| `M_ed_5x5` (16, 25) | 25 | recalculé ee051ec (`ml/M_L_5x5_v2`) | M_ed(juin) − M^L_v2/25 | 0,6710 | 0,2698 | 0,9260 | 2,4e-16 | 2 s |
| `M_ed_6x6` (16, 36) | 36 | recalculé | M_ed(juin) − M^L_v2/36 | 0,8025 | 0,2698 | 1,0034 | 2,2e-16 | 1 s |
| `M_ed_7x7` (16, 49) | 49 | `M_L_7x7` (sept.) × 49 | `M_NL_7x7` (sept.) | 0,6400 | 0,2698 | 0,8544 | 4,3e-15 | 4 s |
| `M_ed_8x8` (16, 64) | 64 | recalculé | M_ed(juin) − M^L_v2/64 | 0,6519 | 0,2698 | 0,8934 | 3,7e-16 | 2 s |
| `M_ed_9x9` (16, 81) | 81 | `M_L_9x9` (juin) × 81 | M_ed(juin) − M_L(juin) (= `M_NL_serial` à 1,3e-13) | 0,6923 | 0,2698 | 0,9114 | 5,9e-15 | 1 s |
| `M_ed_10x10` (16, 100) | 100 | `M_L_10x10` (juin) × 100 | M_ed − M_L (juin) | 0,6460 | 0,2698 | 0,8557 | 9,3e-15 | 79 s |
| `M_ed_12x12` (16, 144) | 144 | `M_L_12x12` (juin) × 144 | M_ed − M_L (juin) | 0,7116 | 0,2698 | 0,9355 | 1,5e-14 | 3 s |
| `M_dense_5x5` (20, 625) | 25 | `M_L_dense_5x5` × 25 | lien → `results/M/M_NL_dense_5x5` | 0,6879 | 0,2893 | 0,9632 | 3,7e-15 | 56 s |
| `M_dense_6x6` (20, 576) | 36 | × 36 | lien | 0,7409 | 0,2913 | 1,0235 | 5,2e-15 | 50 s |
| `M_dense_7x7` (20, 784) | 49 | × 49 | lien | 0,6952 | 0,2913 | 0,9733 | 7,5e-15 | 83 s |
| `M_dense_8x8` (20, 1024) | 64 | × 64 | lien | 0,7052 | 0,2913 | 0,9826 | 1,1e-14 | 200 s |
| `M_dense_9x9` (20, 729) | 81 | × 81 | lien | 0,7591 | 0,2896 | 1,0425 | 1,0e-14 | 72 s |
| `M_dense_12x12` (20, 576) | 144 | × 144 | lien | 0,7781 | 0,2913 | 1,0612 | 2,0e-14 | 46 s |
| `M_L_dense_9x9_coarsecheck` (16, 81) | 81 | × 81 (L seul) | — | 0,6923 | — | — | — | 1 s |
| `M_ed_9x9_nb128` (128, 81) | 81 | R5 `b/M_L_9x9_nb128` × 81 | R5 `b/M_NL_9x9_nb128` (copie) | 0,8315 | 0,3784 | 1,2099 | 7,9e-15 | 43 s |

Sidecars : `bloch_norm = unit_cell`, `units = hartree`, `M_normalization = "v2 : L et NL en norme unit_cell, 2026-09-25"`, `N_cells`, `part`,
`assembled_from` (chemins, facteur, mtime et taille de la source, mode), `date`, `md5`, `checks`, `campaign = R6`, clés v1 reprises (`p`, `D`,
`N_kd`, `vacancy_*`). `results/M2/MD5SUMS_2026-09-25.txt` : 43 fichiers, 47,3 Go réels + 20,8 Go de liens vers `results/M/` (inchangés).
Grandeur brute : max\|N_cells·M^L\| vaut 0,64–0,80 Ha pour toutes les tailles (contre 0,26–0,29 pour M^NL), alors que max\|M^L\| v1 valait
4,9e-3–1,3e-2 Ha.

9×9, porte A.2 sur états purs et identité des 30 états QE (J4b, `states/`) : base 16 (M2 grossier) états purs max\|Δ\| = 6,80e-15 eV (L 6,8e-15,
NL 3,4e-15), rapport L direct/(M2_L/81) = 1,000000 ; base 20 (M2 dense restreint aux 81 k coïncidents) 3,54e-9 eV (L 3,5e-9, NL 1,5e-9) ;
30 états : (i) c†(M2/81)c = (ii) application directe à 2,7e-15 (L) / 3,6e-14 eV (NL) sur la base 16, 6,7e-9 / 8,6e-9 sur la base 20 ; constante
d'alignement s = **+24,66 meV** (min +24,66, max +24,67 ; base 20 : +25,06 à +25,17) pour Lu = −24,65 meV. Les valeurs par état
(σ 324/325 : c†(M2/81)c = +5,849 eV dont L +5,052 ; π 320 : +0,367 dont L +0,338 ; etc., tableau `states/states_tables.md`) sont identiques
à la variante diagnostique « M^L × 81 » de R5 (§A.2), qui est maintenant la convention.

#### 1.3 Porte A.2 pour chaque M2 (J3b, `gate/gate_*.json`, `gate/gate_table.md`)

Huit paires par fichier (intra-k à K ou au k le plus proche, inter-k K→K′, Γ, bandes hautes, paires mixtes) ; pour le dense, les N² k coïncidents
avec la grille grossière, coefficients et k du `.save` dense ; union des ondes planes repliées = sphère de la super-cellule pour les 14 cas
(238 019 à 1 371 809 g) ; 7×7 sur la grille 217 rééchantillonnée.

| N | grossier : max\|Δ_L\| / max\|Δ_NL\| (eV) | dense : max\|Δ_L\| / max\|Δ_NL\| (eV) | rapports direct/(M/N_cells) médians (L ; NL) | verdict |
|---|---|---|---|---|
| 5×5 | 4,16e-15 / 2,91e-15 | 3,55e-15 / 9,44e-16 | 1,000000 ; 1,000000 | OK |
| 6×6 | 3,30e-14 / 6,40e-15 | 3,45e-10 / 5,20e-11 | 1,000000 ; 1,000000 | OK |
| 7×7 | 1,51e-14 / 8,08e-15 | 1,21e-8 / 4,04e-9 | 1,000000 ; 1,000000 | OK |
| 8×8 | 1,32e-14 / 6,60e-15 | 4,66e-15 / 1,56e-15 | 1,000000 ; 1,000000 | OK |
| 9×9 | 6,80e-15 / 3,41e-15 | 3,54e-9 / 1,54e-9 | 1,000000 ; 1,000000 | OK |
| 9×9, 128 bandes | 3,52e-15 / 5,76e-16 | — | 1,000000 ; 1,000000 | OK |
| 10×10 | 4,34e-15 / 5,55e-16 | — | 1,000000 ; 1,000000 | OK |
| 11×11 | — (exclu) | — | — | — |
| 12×12 | 6,63e-15 / 3,10e-15 | 1,10e-10 / 3,19e-11 | 1,000000 ; 1,000000 | OK |
| contrôle : v1 gelé `results/M/M_L_dense_9x9` | — | 1,47e-1 / 1,54e-9 | **81,000000** ; 1,000000 | REFUSÉ (attendu) |

Seuil 1e-6 eV ; les écarts denses de 1e-10–1e-8 eV sont ceux des k du XML dense (bruit 1e-7 sur k, R5 A.2) ; les grossiers sont à la précision machine.
Tous les M2 entrent en production ; aucun fichier de `results/M/` n'a été modifié (README de gel déposé).

#### Fichiers, état du dépôt, points pour Greg

- Répertoire de travail : `ml/` (5 M^L recalculés, 61 Mo), `kernel/` (json + `Mp_9x9_4b.npy` 1,7 Mo), `gate/` (15 json + table), `states/`,
  `assemble_summary.jsonl`, `assemble_ls.txt`, `ksrec/` (liens + npz rejoué), `r6_log.txt`, `JOBID`, `slurm-r6-*`. Copie versionnée
  `article/R6_production_corrigee/` : rapport, scripts R6, json, tables (pas les npy ni les slurm).
- Non commité dans le dépôt (Greg) : `scripts/assemble_M2.py`, `scripts/gate_M_normalization.py`, `src/electron_defect_interaction/wavefunctions/
  sc_projection.py`, `src/electron_defect_interaction/defects/deltav_pw.py` (nouveaux) ; `scripts/compute_M.py`, `src/electron_defect_interaction/
  io/matrix_io.py` (diffs de la phase 0 appliqués) ; `article/R6_production_corrigee/`.
- Avant l'étape 3 (Greg) : `config/production.json` + `config.py` (diffs `phase0/`, couplés : `results_dir`, `M_normalization`) ; corrections des deux
  noyaux non utilisés ici (`fold_wfk_to_sc.py:93`, `local_R.py:322`) ; exceptions `.gitignore` pour `results/M2/` ; le `.save` 128 b de R5 reste non miroité.
- Disque : `results/M2/` 47,3 Go réels sur `/project` (manifeste md5 écrit ; miroir à faire en fin de campagne).

**STOP — étape 1 terminée le 2026-09-25 (17 h 10). Attente du GO 2 (validation contre la super-cellule 9×9 : escalier D4 avec M2).**
