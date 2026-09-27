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

## Étape 2 (GO 2 reçu le 2026-09-25, 22 h 50) — Validation contre la super-cellule 9×9 : escalier D4 avec M2

Git : `pull --ff-only` (Code) vers 09280bb « fix normalisation v2 » (Greg : `fold_wfk_to_sc.py:93` → √(Ω_sc/prod(Ndiag)), `local_R.py:322` →
`meta["Ndiag"]` ; les deux écarts de l'étape 1 sont clos) ; b9f36b0 « R6-R7 » contient les fichiers de l'étape 1. Pilote `r6_d4.py` (prep, d4,
d3, b3) réutilisant les pilotes R4/R5 ; jobs J5 `r6d4` 21833649 (prep + d4 en 42 s, puis d3 : TIMEOUT à 4 h, voir 2.3), J6 `r6b3` 21833650
(COMPLETED 10 min 49), J7 `r6d3` (relance de d3 à 1 fil BLAS).

### 2.0 M_W et V_loc de M2 (prep, `d4/prep_results.json`, `cache/`)

Rotation V†M2V des trois parties denses (8 s chacune), double TF, recentrage R_d = [4, 4, 0], R_loc = 29 mailles (R_cut 3, dim 145) : mêmes
R, R_d et R_loc que R4. Linéarité ‖M_W(tot) − M_W(L) − M_W(NL)‖/‖M_W(tot)‖ = 1,3e-15. Contrôle contre les rotations v1 de R4 :
max|V_loc(M2) − (81 V_L(v1) + V_NL(v1))| = 8,8e-15 eV (L 1,1e-14, NL 0,0), M_W idem 8,8e-15 : la rotation est linéaire, M2 = 81 L + NL
exactement en base de Wannier. Éléments sur le site de la lacune (eV, base de Wannier, maille R = 0) :

| élément | M2 tot | M2 L | M2 NL | v1 tot (R4) | v1 L |
|---|---|---|---|---|---|
| p_z–p_z lacune (A, R = 0) | **31,5208** | 25,2154 | 6,3053 | 6,6166 | 0,3113 |
| p_z(B) diagonal, R = 0 et trois voisins | 0,6151 | 0,5779 | 0,0371 | — | — |
| \|p_z(A)–p_z(B)\| voisins | 1,3905 | 0,9108 | 0,4797 | — | — |
| sp² diagonaux (×3) | 15,689 | 15,246 | 0,4428 | — | — |
| ‖M_W(0, 0)‖ | 41,879 | 36,614 | 9,134 | 9,365 | 0,452 |
| localité : ‖M_W(R, 0)‖ à \|R\| = 1 (a) | 1,963 / 1,963 / 1,388 / 1,388 | | | 0,658 / 0,437 | |
| hors bloc σ–π de V_loc | 0,0 | 0,0 | 0,0 | | |

### 2.1 Escalier D4 avec M2 (d4, `d4/d4_results.json`, `d4/d4_tables.md`, `fig/d4_ladder_M2`)

Portes de R4 : porte 1 (H(R) replié contre Wannier aux 81 k) 4,88e-14 eV (seuil 1e-8) → OK ; porte 2 (repliement de M_W(M2), même partie à
un corps) 1,93e-6 eV avec les k du XML dense, **2,26e-12 eV avec k = m/27** (seuil 1e-9) → OK ; V_loc production contre k exacts 7,5e-7 eV
(p_z–p_z 31,520784 / 31,520784) ; résidu H(R) contre V†εV aux 81 k 7,5e-6 eV. E_D : maille 16 b −4,238470, dense 20 b −4,238895,
Wannier −4,238896, parfaite −4,238471 eV. M2 : max|tot − L − NL| 4,5e-15 (16 b), 4,9e-15 eV (20 b).

États localisés (ε − E_D en eV ; w₂ disque 2 Å, seuil 0,0812 ; QE = 9×9 lacune de D1, aligné Lu ; même fenêtre de 30 états, 2 pairs + 28 impairs) :

| variante | dim | couplage pair–impair (eV) | pairs (σ) localisés | impairs (π) localisés | décalage rigide résiduel (meV, 269 états < E_D − 4) |
|---|---|---|---|---|---|
| QE (D1) | — | — | +0,101 ×2 (0,713) | −1,759 ×2 (0,114) ; **−0,737 (0,246)** ; +0,269 (0,102) | — |
| (a1) M2 16 bandes, total | 1296 | 3,7e-3 | aucun (2 états pairs de la fenêtre à w₂ < seuil) | −1,775 ×2 (0,116) ; **−0,727 (0,236)** ; +0,263 (0,107) | +24,3 |
| (a1) M2^L seul | 1296 | — | −2,241 (0,618) | −1,775 ×2 (0,116) ; −0,742 (0,239) ; +0,255 (0,105) | — |
| (a1) M^NL seul | 1296 | — | −2,847 ×2 (0,255) | −1,370 (0,222) | — |
| (a2) M2 dense ⊂ 81 k, 16 bandes | 1296 | 3,2e-8 | aucun | −1,775 ×2 (0,116) ; −0,727 (0,236) ; +0,263 (0,107) | +24,4 |
| (a3) idem 20 bandes | 1620 | 2,9e-5 | aucun | −1,776 ×2 (0,116) ; **−0,734 (0,238)** ; +0,259 (0,106) | +24,4 |
| (b) 5 WF, V†εV ⊕ M2_W/81 | 405 | 1,7e-8 | **−0,812 (0,821)** | −1,774 ×2 (0,116) ; **−0,689 (0,225)** ; +0,286 (0,115) | +24,0 |
| (c-all) H(R) + M2_W replié, toutes mailles | 405 | 1,7e-7 | −0,812 (0,821 ; poids WF site + voisins 0,791) | −1,774 ×2 (0,116) ; −0,689 (0,225) ; +0,286 (0,115) | +24,0 |
| (c-3) idem, R_cut 3 | 405 | 1,7e-7 | −0,813 (0,821 ; 0,791) | −1,765 / −1,763 (0,115) ; −0,677 (0,227) ; +0,302 (0,113) | +6,8 |

Attendu R5 (variante M^L × 81) : (a1) −0,727 (0,236), (a3) −0,734 (0,238) — **reproduits à l'identique** (R5 : −0,72741 / −0,73436 ; ici −0,72741 /
−0,73436, w₂ 0,2358 / 0,2380). Pour mémoire, v1 (R4, w₂ corrigés) : (a1) −1,350, (a3) −1,359, (b)/(c-all) −1,315, paire σ −2,812 / −2,848 / −2,514.

Marches (états localisés, même parité, plus proche voisin ; Δε en eV) :
- QE → (a1) : π −0,737 → −0,727 (+0,010) ; −1,759 ×2 → −1,775 (−0,016) ; +0,269 → +0,263 (−0,007) ; σ +0,101 ×2 → absent (aucun état pair localisé
  dans (a1), (a2), (a3) : les deux états pairs de la fenêtre y sont à −2,818/−2,814 avec w₂ ≤ 0,04). Médiane |Δε| des 28 impairs triés : 24 meV (= le décalage
  rigide résiduel de 24,3 meV ; QE aligné Lu contre modèle aligné sur son E_D).
- (a1) → (a2) : identiques (Δε ≤ 1e-6, Δw₂ ≤ 1e-6) ; (a2) → (a3) (16 → 20 bandes) : π −0,727 → −0,734 (−0,007), +0,263 → +0,259 ; −1,775 → −1,776.
- (a3) → (b) (20 bandes → 5 WF) : π −0,734 → −0,689 (+0,046), +0,259 → +0,286 (+0,028), −1,776 → −1,774 ; **paire σ : apparaît à −0,812 (w₂ 0,821)**,
  un seul état pair localisé (l'autre membre de la paire QE +0,101 ×2 n'a pas d'équivalent : la base à 5 WF n'a que 3 sp² sur A par maille).
- (b) → (c-all) : identiques (Δε ≤ 1e-6). (c-all) → (c-3) : Δε = −1 meV (σ), +8 / +11 / +16 meV (π) ; décalage rigide résiduel +6,8 meV contre +24,0.
- QE → (c-3) (chaîne complète de production, R_cut 3) : π −0,737 → −0,677 (+0,060) ; +0,269 → +0,302 (+0,033) ; −1,759 → −1,763 (−0,003) ;
  σ +0,101 → −0,813 (−0,915).

Tous les états de la fenêtre par variante : `d4/d4_tables.md`. Figure `fig/d4_ladder_M2` (niveaux par variante et parité, taille ∝ w₂).

### 2.2 (a1) à n bandes avec le M2 à 128 bandes (b3, J6, `b3/b3_tables_M2.md`, `fig/b3_ladder_M2`)

Le « (c-all) avec le M2 à 128 bandes projeté sur 5 WF » n'apporte rien de plus que (c-all) à 20 bandes : les 5 fonctions de Wannier sont
construites dans les 20 bandes du `.save` dense (V = U_dis·U de dimension 20 × 5, jauge du `.save` dense), donc P M128 P se réduit au bloc 20 × 20 de
M128, dans une jauge différente (`.save` nb128) pour laquelle il n'y a pas de U. Remplacé par l'échelle gauge-invariante (a1) à n bandes sur le M2
nb128 réassemblé (`results/M2/M_ed_9x9_nb128.npy`, 1,7 Go, max|M2| 32,9 eV), contre R5 B.3 (variante M^L × 81) :

| n | dim | π localisé (ε − E_D ; w₂) | Δε_π vs n précédent (meV) | σ localisés (×2) | δ(π_320) / δ(σ_324) | médiane \|Δε\| vs QE σ / π (eV) | R5 (M^L × 81) π |
|---|---|---|---|---|---|---|---|
| 16 | 1 296 | −0,727 (0,236) | — | aucun | 0,0009 / 0,0363 | — / 0,024 | −0,727 (0,236) |
| 24 | 1 944 | −0,741 (0,240) | −13,2 | +0,843 (0,674) | 0,0005 / 0,0176 | 0,742 / 0,025 | −0,741 (0,240) |
| 32 | 2 592 | −0,747 (0,242) | −6,3 | +0,428 (0,697) | 0,0003 / 0,0057 | 0,326 / 0,025 | −0,747 (0,242) |
| 48 | 3 888 | −0,754 (0,244) | −7,5 | +0,345 (0,700) | 0,0001 / 0,0039 | 0,244 / 0,025 | −0,754 (0,244) |
| 64 | 5 184 | −0,757 (0,245) | −2,9 | +0,211 (0,708) | 0,0001 / 0,0017 | 0,109 / 0,025 | −0,757 (0,245) |
| 96 | 7 776 | −0,759 (0,245) | −2,2 | +0,137 (0,713) | 0,0000 / 0,0006 | 0,035 / 0,025 | −0,759 (0,245) |
| 128 | 10 368 | **−0,760 (0,245)** | −0,9 | **+0,116 (0,714)** | 0,0000 / 0,0003 | 0,014 / 0,025 | −0,760 (0,245) |

QE : π −0,737 (0,246), σ +0,101 ×2 (0,713). Premier n avec une paire σ à moins de 0,3 eV de +0,101 : 48. Identique à R5 B.3 (variante M^L × 81)
à toutes les décimales rapportées : le fichier M2 nb128 réassemblé est bien le M de cette variante.

### 2.3 Critère de pôle et −Im T̄(K) avec V_loc(M2) (d3, J7 21850495, 1 min 57 ; `d3/d3_results_M2.json`, `d3/d3_tables_M2.md`, `d3/d3_curves_M2.npz`, `fig/d3_pole_M2`)

Même construction que R4 D3 et que `resonance_criteria.py` : V_loc = P M_W P sur les 29 mailles (dim 145), α = 1, η 0,02 eV, g₀ du vrai H(R)
(caches R4 réutilisés : mêmes R_loc), fenêtre [−3, +3] à 300² et [−3, +1] à 600², |det[1 − g₀V]| normalisé sur la fenêtre, λ_min de 1 − g₀V,
T̄(K) = trace/2 de la paire π ; blocs par parité (complet, π, σ). Exécution : J5 avait lancé d3 avec 16 fils BLAS (61 ms par matrice au lieu de
4,7 ms, piège consigné dans R4) : 50 min par variante, mur de 4 h atteint avant l'écriture des fichiers ; J7 (1 fil BLAS, 16 fils Python) : 2 min.

| V_loc | grille | bloc | min \|det\|/max (ε − E_D, eV) | min \|λ\| ; λ (ε − E_D) | vecteur propre (poids) | racines de Re λ_min | pic de −Im T̄(K) (max, eV) | Re T̄(E_D) | zéros de Re T̄ |
|---|---|---|---|---|---|---|---|---|---|
| **M2 tot** | 300² | complet | **1,32e-4 (−0,812)** | **0,0019 ; +0,0000 + 0,0019 i (−0,812)** | σ : sp² de l'atome retiré 0,989, autres sp² 0,011, π 0 | −0,813 | **−0,177 (24,45)** | **+11,12** | −2,446 ; −0,222 ; +1,639 |
| M2 tot | 300² | σ | 3,48e-4 (−0,812) | 0,0019 (−0,812) | idem | −0,813 | — | — | — |
| M2 tot | 300² | π | 2,33e-2 (−0,172) | 0,3967 ; +0,219 + 0,331 i (−0,170) | p_z lacune 0,997 | −0,371 ; −0,254 (λ ≈ 0,4–0,6 : pas de zéro) | −0,177 (24,45) | +11,12 | idem |
| M2 tot | 600² | complet | 1,27e-4 (−0,812) | 0,0019 (−0,812) | idem | −0,813 | −0,190 (23,64) | +11,15 | −2,443 ; −0,225 |
| M2 tot | 600² | π | 3,52e-2 (−0,162) | 0,4066 ; +0,272 + 0,302 i (−0,152) | p_z lacune 0,997 | −0,366 ; −0,247 | −0,190 (23,64) | +11,15 | idem |
| M2 L seul | 300² | complet | 1,50e-4 (+3,0, bord) | 0,0577 ; −0,058 + 0,001 i (+3,0, bord) | σ | aucune | −0,217 (22,51) | +10,32 | −2,446 ; −0,255 ; +1,639 |
| M2 L seul | 300² | π | 2,54e-2 (−0,175) | 0,3721 (−0,172) | p_z lacune | −0,431 ; −0,265 | −0,217 | +10,32 | |
| NL seul | 300² | complet | 2,48e-4 (−2,620) | 0,0123 ; +0,0003 + 0,0123 i (−2,620) | σ (sp² 0,975) | −2,620 | −1,292 (3,02) | +2,41 | aucun |
| NL seul | 300² | π | 1,83e-1 (−0,905) | 0,7837 (−0,905) | p_z lacune | aucune | −1,292 | +2,41 | |
| v1 tot (R4, pour mémoire) | 300² | complet | 2,09e-4 (−2,530) | 0,0108 ; +0,0003 + 0,0108 i (−2,530) | σ (sp² 0,976) | −2,531 | −1,292 (3,24) | +2,52 | aucun |
| v1 tot | 300² | π | 1,74e-1 (−0,905) | 0,7753 (−0,905) | p_z lacune | aucune | −1,292 | +2,52 | |
| v1 tot | 600² | complet | 3,81e-4 (−2,530) | 0,0108 (−2,530) | | −2,531 | −1,242 (3,15) | +2,52 | |

Minima locaux secondaires de |det|/max (M2 tot, complet, 300²) : −0,630 (7,5e-4), −0,255 (5,7e-4), −0,210 (5,0e-4), −0,172 (4,7e-4), −0,132 (5,2e-4) ;
bloc π : −0,462 (6,3e-2), −0,412, −0,255, −0,210, −0,172 (2,3e-2), −0,132. Bloc π (M2, 300²) : min |λ| par minimum local 0,574 (−0,252), 0,456 (−0,207),
0,397 (−0,170), 0,399 (−0,130) ; racines de Re λ_min à −0,371 et −0,254 avec |λ| ≈ 0,4–0,6 (passage de Re λ par zéro, pas de zéro de λ).
v1 et M2 se retrouvent en regard de −0,73 (π QE −0,737 ; (c-3) −0,677) : v1 : σ à −2,530, pic −Im T̄ à −1,292, bloc π sans zéro (min |λ| 0,78 à −0,905) ;
M2 : σ à **−0,812** (= la paire σ de (b)/(c) du §2.1, −0,812/−0,813), pic −Im T̄ à **−0,177** (300²) / −0,190 (600²) avec un maximum 7,5 fois plus haut
(24,5 contre 3,2 eV), bloc π sans zéro (min |λ| 0,40 à −0,17). Le pic de −Im T̄ suit la partie locale (L seul : −0,217 ; NL seul : −1,292 comme v1).
Figure `fig/d3_pole_M2` : |det|/max, min |λ|, −Im T̄(K) à 300², M2 (plein) contre v1 (tirets), trait à −0,737.

### Fichiers et état à la fin de l'étape 2

- Répertoire de travail : `cache/Mwr_M2_9x9.npz` (638 Mo), `cache/Vloc_M2_9x9.npz`, `d4/` (json, npz, tables), `d3/` (json, npz 1,3 Mo, tables),
  `b3/`, `fig/{d4_ladder_M2, b3_ladder_M2, d3_pole_M2}`, `r6_d4.py`, `submit_r6.sh` (tâches d4, b3, d3), `slurm-r6-*`, `JOBID`.
  Copie versionnée `article/R6_production_corrigee/etape2/` (pilote, json, tables, figures ; pas les npz > 5 Mo ni les slurm).
- Aucun fichier de production modifié ; `results/M/` intact ; `results/M2/` inchangé depuis l'étape 1. Git : rien fait par Code depuis le pull ff-only.
- Reste à Greg avant GO 3 : `config/production.json` + `config.py` (diffs `phase0/`, couplés) ; exceptions `.gitignore` pour `results/M2/` ;
  commit de `article/R6_production_corrigee/` (étape 2) et de `scripts/`… si retouchés (aucun script de `scripts/` ni de `src/` modifié à l'étape 2).

**STOP — étape 2 terminée le 2026-09-26 (9 h 10). Attente du GO 3 (production du chapitre 4 avec M2, config v2).**

Mise à jour du 2026-09-26 (Code, à la demande de Greg) : diffs `phase0/production.json.diff` et `phase0/config.py.diff` appliqués (`results_dir = results/M2`, `M_normalization = v2`, `results_dir(cfg)`, clés obligatoires) ; exceptions `.gitignore` pour `results/M2/` (mêmes règles que `results/M/` + `MD5SUMS_*.txt`). Vérifié : `load_production()` exige et affiche les nouvelles clés, `dense_paths` pointe sur `results/M2/M_dense_<S>.npy` (six fichiers présents, sidecar v2). Non commité.

## Étape 3 (GO 3 reçu le 2026-09-26 avec les ajouts 3.0, 3.0b, 3.3 porte, 3.4 blocs)

### 3.0 Lecteurs de M (aucun calcul ; `etape3/diff_3.0.patch`, `etape3/diff_3.0.stat.txt`, `etape3/grep_3.0.txt`)

Refonte systématique (script `refactor_3_0.py`, conservé dans `etape3/`) puis corrections à la main ; rien n'est exécuté sur les données.

Règles appliquées :
- **Chemins** : tout littéral `results/M/` des scripts Python devient `f"{RES}/…"` avec `RES = results_dir(load_production(verbose=False))`
  (import `config.results_dir` ajouté ; `RES` défini après le bloc d'imports, donc disponible pour les valeurs par défaut d'argparse).
  Dans les docstrings, commentaires et messages, `results/M/` devient `<results_dir>/`. Dans les `submit_*.sh` : `RES=$(python -c … results_dir …)`
  après le chargement des modules, `mkdir -p "$RES/logs"`, toutes les sorties `--out "$RES/…"` ; les directives `#SBATCH --output/--error` ne
  peuvent pas être calculées : littéral `results/M2/logs/` (= `results_dir` de la config ; c'est le seul littéral `M2` des scripts).
- **Porte v2** : `load_M_checked(…, require_normalization=matrix_io.M_NORM_V2)` dans les dix scripts qui chargent M ainsi ; pour les lectures par
  `np.load(…, mmap_mode="r")` (gros fichiers denses), nouvelle fonction `matrix_io.check_manifest(path, require_normalization=…)` (même refus que
  `load_M_checked`, sans charger) appelée avant chaque mmap (`analyze_M.mmap_M`, `lnl_frobenius_all.mmap_M`, `resonance_metrics` (M^L dense),
  `_bz_ratio_LNL`, `check_onsite_and_NL`, `check_M_dense_vs_coarse`, `check_ML_coarse_kernel`, `check_M_dense_nb20_vs_nb16`, `_mcheck`, `_normtest`).
  `analyze_M` : le repli « fichier grossier sans sidecar » (`np.load(pb)`) est supprimé, les M2 grossiers ont un sidecar.
- **Fichiers absents → « non rejouable », pas de plantage** : `analyze_M` (bloc nbnd 16 → 20 : ligne `non rejouable / fichier absent` dans
  `M_tests_summary.csv`), `check_M_dense_nb20_vs_nb16.py`, `_old_vs_new_7x7.py` (`obsolete_grid_7x7` supprimé au ménage), `validate_ML_grid_7x7.py`
  (référence série à recalculer avec le noyau corrigé), `compute_convergence.py` (aucun `gamma_*.npz`).
- **Chaîne obsolète** (`compute_spectral.py`, `compute_tmatrix.py`, `_eta_scan.py`, `migrate_M_norm.py` : lisent les `M_ed_*_norm.npy` supprimés,
  écrivaient `gamma_*`, `dos_*` sous `results/M`) : chemins convertis comme les autres **et** refus explicite à l'exécution (`SystemExit "[obsolète, R6 …]"`).
- **Producteurs** : `compute_M_dense_stages.py` tague désormais `M_normalization = v2` sur ses trois sorties (ml, nl, combine, + `N_cells`) et son
  combine exige v2 sur les parties (comme `compute_M.py` depuis l'étape 1) ; `submit_M_dense.sh` n'écrit plus la copie `_norm` (supprimée au ménage).
- **PYTHONPATH** : le paquet n'a pas d'installation éditable dans `.venv` (constat de l'étape 1) ; chaque `submit_*.sh` exporte
  `PYTHONPATH="$PROJ/src:$PYTHONPATH"` (préfixe : h5py de scipy-stack conservé). Sans cela aucun script de production ne s'importe.

Fichiers touchés : 36 scripts Python, 11 `submit_*.sh` (dont `_test_mnl_mpi.sh`), `src/…/io/matrix_io.py` (`check_manifest`) ; avec la config et le
`.gitignore` de la veille : 50 fichiers, +360/−177 lignes (`etape3/diff_3.0.stat.txt`).

Contrôles :
- `grep -n "results/M/" scripts/*.py scripts/*.sh src/ tests/` (hors `results/M2/`) : **aucune ligne** (`etape3/grep_3.0.txt`) ; occurrences restantes
  de `results/M2` : les `#SBATCH --output/--error` de 11 scripts ; 36 scripts définissent `RES = results_dir(…)`, 11 submit définissent `RES=$(…)` ;
  aucun `load_M_checked` sans porte v2 dans la chaîne.
- `py_compile` des 80 scripts Python : OK ; `bash -n` des 15 scripts shell : OK.
- `--help` (argparse, donc `RES` et la config chargés) : `make_figures`, `make_figures_epw`, `nkint_check_post` (après correction du gabarit
  `--pattern`, devenu f-string par la conversion : `{{S}}`, `{{N}}`), `epw_ed_vs_ep`, `compute_spectral_wannier`, `resonance_metrics`, `rcut_resigma`,
  `resonance_criteria`, `m_rcut_convergence`, `ks_reconstruction_all` : OK.
- Gardes : `check_M_dense_nb20_vs_nb16.py 9x9`, `_old_vs_new_7x7.py`, `validate_ML_grid_7x7.py cmp` → `[non rejouable]` ; `compute_spectral.py` → `[obsolète]`.
- Exécution réelle sur M2 : `check_ML_coarse_kernel.py 9x9` (M_L_9x9 contre coarsecheck, tous deux v2) : rel 1,28e-15, hermiticité 2,3e-16, PASS.
- Porte v2 : `check_manifest("results/M/M_dense_9x9.npy", v2)` refuse (`M_normalization='v1'`), le même appel sur `results/M2/` accepte.
- **Garde-fou** : `chmod -R a-w results/M` appliqué (répertoire, `logs/`, fichiers : `r--r-----`) ; `touch results/M/_test_write` → « Permission denied » ;
  lecture des `M_NL_dense_*` de `results/M2/` (liens) vérifiée. À rétablir (`chmod -R u+w results/M`) seulement sur GO.

Non touchés (à signaler) : `CLAUDE.md` (« Données : lues uniquement dans results/M/*.npz », section figures) et `NOTES_TGAMMA.md` citent encore
`results/M/` : documentation, pas des scripts ; à mettre à jour par Greg avec le rapport final. Un `.gitignore.swp` (vim) traîne à la racine du dépôt.

**STOP 3.0 — relecture et commit par Greg (`git diff` = `etape3/diff_3.0.patch`), puis GO 3b (test d'or 5×5 dense avec M2, bloquant).**

### 3.0b Test d'or (GO 3b ; `scripts/submit_golden_dense.sh 5x5`, job `r6golden` 21852238, COMPLETED 1 h 52 min 27 s, 16 cœurs, 120 Go)

Commit de Greg après 3.0 : 2e49ffa « r6 stage 3 » (clone à jour, arbre propre). `test_local_tmatrix_real.py 5x5 --dense` sur `results/M2/M_dense_5x5.npy`
(sidecar v2 exigé) : R_d = [2, 2, 0], ‖V_loc(0, 0)‖ = **40,8800 eV** (v1 : 9,9183), résidu d'hermiticité 1,3e-14 ; Γ dense × N_cells et Γ local sur
[0,2646, 14,93] eV ; **max|Γ_loc − Γ_dense N_cells| / max = 1,80e-13 (seuil 1e-8) : PASS** ; positivité min Γ_loc = +0,2646 : OK. Code de sortie 0 →
3.1 autorisé.

Second contrôle du même job (`check_M_dense_vs_coarse.py 5x5`, M2 complet dense contre grossier aux 25 k coïncidents, 16 bandes) : écart relatif
maximal des valeurs singulières 3,48e-1, des normes de Frobenius 2,20e-1, max|M| 0,9632 (dense) / 0,9260 (grossier) → « CHECK ». Valeur v1 du même
contrôle : voir ci-dessous (§3.6) ; le test de la chaîne `analyze_M` restreint ce contrôle aux bandes 1–8 et 1–15 (la 16ᵉ bande est dégénérée avec la
17ᵉ dans le `.save` dense, R4 §D4), c'est cette version qui entre dans `tab:tests_M`.

### 3.1–3.3 Jobs A/B/C (après le test d'or ; `etape3/runbook_3.sh`, `JOBID`)

| job | script | durée | produit (`results/M2/`) |
|---|---|---|---|
| specwd_5x5 21857272 / 6x6 21857273 / 7x7 21857274 / 8x8 21857275 / 9x9 21857276 / 12x12 21857277 | `submit_spectral_wannier_dense.sh S prod` (grilles 60/120/240, η 0,05/0,02/0,01, R_cut 0…4 : 45 combinaisons par taille) | 1 h 38 / 1 h 47 / 52 min / 1 h 23 / 53 min / 2 h 10 | `specwd_<S>_prod.npz` |
| nkint 21857278 | `submit_nkint_check.sh 9x9` (N_k^int 150, 300, 450, 600) | 7 min | `resigma_9x9_rc3_nk*.npz` |
| resigma 21857279 / 21857280 | `submit_rcut_resigma.sh 9x9 0,1,2,3` / `9x9 4` | 3 min / 2 min | `resigma_9x9_rc0123.npz`, `resigma_9x9_rc4.npz` |
| post_locality 21857281 | `submit_post.sh locality` (`mwr_locality_coarse_vs_dense.py` ; `m_rcut_convergence.py` 9×9 et 12×12, R_cut 0…6) | 9 min | `mwr_locality.npz`, `m_rcut_convergence.csv` |
| post_ksrec 21857284 | `submit_post.sh ksrec` (`ks_reconstruction_all.py` six tailles ; copie de `ved_analysis.npz` v1, indépendant de M ; `sampling_table.py`) | 16 min | `ks_reconstruction.npz`, `ved_analysis.npz`, `sampling_table.csv` |
| post_analyze 21857283, 21857345 | `submit_post.sh analyze` | FAILED | `analyze_M.py` : `NameError: ks` — l'échantillon de k (graine 0) avait été déplacé dans le bloc « nb16 » par la refonte 3.0 et sert aussi à la fermeture Wannier ; corrigé (défini hors du bloc, même graine), 3ᵉ soumission (voir `JOBID`). Ligne « test d'or » d'`analyze_M` paramétrée par `GOLDEN_RESULT="1.80e-13,r6golden_21852238"` (elle était codée en dur avec la valeur v1) |

### 3.3 Porte de niveau 1 sur les paramètres gelés (`r6_level1_gate.py`, `etape3/level1_gate.md`)

Même critère que celui qui a figé `config/production.json` le 2026-09-05 (NOTES_TGAMMA §3 : C11 ≤ 5 % quand la grille double et quand η est
divisé par 2 ; C10 ≤ 5 % de R_cut 3 à 4 ; N_k^int rapporté avec le même seuil sur la médiane de la fenêtre, Γ(E_D) à part). Médianes de |Γ|·N_cells
(meV, états de la fenêtre ±3 eV) :

| taille | médiane (R_cut 3, 240², η 0,02) | E_res (eV) | C11 grille 120 → 240 | C11 η 0,02 → 0,01 | C10 R_cut 3 → 4 | verdict |
|---|---|---|---|---|---|---|
| 5×5 | 3 288,90 | −0,284 | 0,34 % | 1,06 % | 0,42 % | OK |
| 6×6 | 3 155,98 | −0,227 | 0,01 % | 0,52 % | 1,58 % | OK |
| 7×7 | 3 052,03 | −0,269 | 0,09 % | 0,53 % | 2,99 % | OK |
| 8×8 | 3 020,42 | −0,269 | 0,10 % | 0,45 % | 2,08 % | OK |
| **9×9** | **3 132,60** | **−0,175** | 0,49 % | 0,10 % | 0,47 % | **OK** |
| 12×12 | 3 162,78 | −0,175 | 0,41 % | 0,44 % | 0,04 % | OK |

N_k^int (9×9, R_cut 3) : 150 → 3 155,57 (+0,73 % vs 300), 300 → 3 132,60, 450 → 3 135,87, 600 → 3 135,65 (300 vs 600 : 0,10 %) ; E_res −0,238 / −0,175 /
−0,202 / −0,202 ; Γ_T(E_D) (états à K) 4 801,9 / 3 630,7 / 3 489,5 / 3 463,9 meV (300 → 600 : +4,81 % ; v1 : +5,36 %). R_cut par `rcut_resigma` (9×9) :
0 → 3 009,73 (+4,38 % vs R_cut 4), 1 → 3 223,74, 2 → 3 170,21, 3 → 3 132,60 (+0,47 %), 4 → 3 147,47 ; E_res −0,175 pour R_cut 0, 2, 3, 4 (−0,227 à R_cut 1) ;
Re Σ médian à R_cut 3 : −26,65 meV (v1 : +722,3), |Re Σ|/Γ médian 0,125 (v1 : 0,227). Pour mémoire, v1 : médiane 9×9 2 473,55 meV, E_res −1,238 eV.

**Verdict : les paramètres gelés satisfont le critère du 2026-09-05 avec M2 ; la config n'est pas modifiée ; 3.4 lancé.**

Niveau 2 (post_fig 21872955 ; `level2_summary.csv`, `level2_families.csv` ; R_cut 3, 240², η 0,02) :

| taille | famille | sous-réseau de la lacune | médiane Γ·N_cells v2 (meV) | v1 | p_z–p_z lacune (eV) | Re M^L(K) / Re M^NL(K), paire π (eV) | écart Δk à K (Å⁻¹) |
|---|---|---|---|---|---|---|---|
| 6×6 | 3m | B | 3 155,98 | 2 509,20 | 31,502 | +11,356 / +3,103 | 0,000 |
| 9×9 | 3m | A | 3 132,60 | 2 473,55 | 31,521 | +11,305 / +3,103 | 0,000 |
| 12×12 | 3m | B | 3 162,78 | 2 458,74 | 31,528 | +11,222 / +3,103 | 0,000 |
| 5×5 | non-3m | A | 3 288,90 | 2 524,30 | 31,002 | +10,534 / +3,100 | 0,055 |
| 7×7 | non-3m | A | 3 052,03 | 2 487,46 | 30,936 | +10,486 / +3,101 | 0,050 |
| 8×8 | non-3m | A | 3 020,42 | 2 478,47 | 31,092 | +10,734 / +3,101 | 0,043 |

Familles : 3m moyenne 3 150,5 meV, écart (max − min)/moyenne 0,96 % (v1 : 2 480,5 meV ; 2,03 %) ; non-3m 3 120,5 meV, 8,60 % (v1 : 2 496,7 meV ;
1,84 %). Re M^L(K) v1 : +0,316 (6×6), +0,140 (9×9), +0,078 (12×12), +0,421 (5×5), +0,214 (7×7), +0,168 (8×8) ; Re M^NL(K) inchangé.

### 3.1 Rotation de Wannier, recentrage, localité de M_W, tableau R_cut (post_locality 21857281 ; `mwr_locality.npz`, `m_rcut_convergence.csv`)

Recentrage inchangé (R_d = [4, 4, 0] pour la 9×9 dense, [2, 2, 0] pour la 5×5 et la 6×6). Localité de M_W (eV ; ‖M_W(0, 0)‖, p_z–p_z du site
de la lacune, p_z de l'autre sous-réseau, premières couronnes ‖M_W(R, 0)‖ à |R| = 1) :

| taille (dense) | ‖M_W(0,0)‖ v2 | v1 | p_z–p_z lacune v2 | v1 | p_z autre sous-réseau v2 | couronne \|R\| = 1, v2 (v1) |
|---|---|---|---|---|---|---|
| 5×5 | 40,880 | 9,918 | 31,002 (A) | 7,296 | 0,316 | 1,965 / 1,385 / 1,495 (0,681 / 0,441 / 0,468) |
| 6×6 | 35,290 | — | 31,502 (B) | 7,007 | 0,569 | 2,327 / 0,757 / 0,832 |
| 7×7 | 40,777 | — | 30,937 (A) | 6,807 | 0,319 | 1,962 / 1,385 / 1,498 |
| 8×8 | 41,072 | — | 31,092 (A) | 6,691 | 0,393 | 1,962 / 1,386 / 1,498 |
| **9×9** | **41,879** | 9,365 | **31,521 (A)** | 6,617 | 0,615 | 1,963 / 1,388 / 1,501 (0,658 / 0,437 / 0,460) |
| 12×12 | 35,329 | — | 31,528 (B) | 6,483 | 0,635 | 2,328 / 0,756 / 0,832 |

Grossiers (wannierisations N×N, 5×5, 7×7, 8×8 ; aliasés) : ‖M_W(0,0)‖ 4,010 / 2,235 / 1,958 eV (v1 : 0,924 pour la 5×5).

Tableau R_cut (tab:rcut_M ; `m_rcut_convergence.py`, grille fine 60², max|ΔM|/max|M| de la paire π reconstruite avec V_loc tronqué) :

| R_cut | n_sites | 9×9 v2 | 9×9 v1 | 12×12 v2 | 12×12 v1 |
|---|---|---|---|---|---|
| 0 | 1 | 3,16e-1 | 4,09e-1 | 8,67e-1 | 7,31e-1 |
| 1 | 5 | 2,13e-1 | 1,74e-1 | 2,17e-1 | 1,86e-1 |
| 2 | 13 | 1,14e-1 | 4,97e-2 | 1,25e-1 | 1,16e-1 |
| **3** | 29 | **6,68e-2** | 2,27e-2 | **9,02e-2** | 2,58e-2 |
| 4 | 49 | 4,42e-2 | 1,57e-2 | 7,57e-2 | 1,57e-2 |
| 5 | 81 | 1,75e-2 | 1,25e-2 | 5,36e-2 | 1,24e-2 |
| 6 | 113 | 1,21e-2 | 1,08e-2 | 3,29e-2 | 1,15e-2 |

max|M_π| (eV) : 25,44 (9×9), 25,94 (12×12) (v1 : 6,46 / 6,37).

### 3.2 Cartes de couplage, comparaison Kaasbjerg, tableau L / NL (post_analyze 21862371 ; `M_analysis.npz`, `lnl_frobenius.csv`)

9×9 dense 27×27, K à l'indice 495, paire π/π* = bandes (3, 4) (poids p_z 1,0), A_cell = 5,266 Å² : M̃ = |M|·A_cell (eV Å²) max 116,46 (π ; à K 107,22) et
146,02 (π* ; à K 109,79) — v1 : 23,95 (23,87) et 34,16. Kaasbjerg : V₀ = 70 eV Å², soit 13,3 eV par maille (27 eV par atome pour A_cell 5,24 Å²).
Échelle intensive (max|M|, bandes 1–16, eV) : dense 23,85 / 25,15 / 23,82 / 24,31 / 25,35 / 25,72 et grossier 25,20 / 27,30 / 23,25 / 24,31 / 24,80 / 25,46
pour N = 5, 6, 7, 8, 9, 12 (v1 dense : 0,26–0,27 eV, maxima dominés par M^NL).

Tableau L / NL (tab:L_NL ; `lnl_frobenius_all.py`, moyennes des normes de Frobenius des blocs 2 × 2 π/π* sur tous les couples (k′, k)) :

| taille | ⟨‖M^NL‖_F⟩ / ⟨‖M^L‖_F⟩ v2 | v1 | ⟨‖M^L‖_F⟩ v2 (eV) | v1 | ⟨‖M^NL‖_F⟩ (eV, inchangé) |
|---|---|---|---|---|---|
| 5×5 | 0,258 | 6,45 | 24,44 | 0,977 | 6,300 |
| 6×6 | 0,252 | 9,09 | 24,83 | 0,690 | 6,268 |
| 7×7 | 0,258 | 12,66 | 24,36 | 0,497 | 6,291 |
| 8×8 | 0,257 | 16,43 | 24,52 | 0,383 | 6,292 |
| 9×9 | 0,252 | 20,44 | 25,01 | 0,309 | 6,311 |
| 12×12 | 0,252 | 36,31 | 24,86 | 0,173 | 6,268 |

(v1 : le rapport croissait comme N_cells parce que M^L portait 1/N_cells de trop ; v2 : rapport 0,25 indépendant de N.)

### 3.4 Résonance et critères (post_res9x9 21862372 [6 h 35 : métriques 7 min 49, critères 6 h 27], post_res6x6 21862373 [2 h 39], post_res12x12 21862374 [1 h 38] ; `resonance_<S>.npz`, `resonance_criteria_<S>.npz`)

`resonance_metrics.py` (R_cut 3, 240², η 0,02, N_k^int 300, ρ₀ 900², c = 1 %) ; énergies en eV relatives à E_D ; médianes des courbes sur ±3 eV :

| grandeur | 9×9 v2 | 9×9 v1 | 6×6 v2 (v1) | 12×12 v2 (v1) |
|---|---|---|---|---|
| médiane courbe Γ_T (meV) | **3 740,2** | 2 343,6 | 3 780,1 (2 453,6) | 3 749,5 (2 304,4) |
| médiane états \|Γ_T\|·N_cells (meV) ; E_res états ±1,5 eV | 3 132,9 ; −0,175 | 2 473,5 ; −1,238 | 3 156,5 ; −0,227 | 3 166,6 ; −0,175 |
| médiane Γ_Born (meV) ; Born/T médian (min, max) | 177 204,7 ; 46,5 (1,11, 211,7) | 7 405,8 ; 3,30 (0,68, 16,5) | 174 332,9 ; 46,5 | 179 060,4 ; 45,9 |
| pic Γ_T ; pic Γ_Born | −0,180 ; +1,695 | −1,240 ; +1,695 | −0,240 ; +1,695 | −0,180 ; +1,695 |
| pic Γ_T/ρ₀ ; pic δρ ; pic ρ_dis | −0,175 ; −0,815 ; +1,710 | −0,015 ; −2,530 ; −2,530 | −0,175 ; −0,847 ; +1,705 | −0,170 ; −0,820 ; +1,710 |
| pic \|T̄\| ; pic −Im T̄ ; min \|Re T̄\| | −0,172 ; −0,177 ; −2,447 | −0,905 ; −1,292 ; −2,145 | −0,177 ; −0,220 ; +1,435 | −0,170 ; −0,175 ; −2,207 |
| Re T̄(E_D) ; Im T̄(E_D) | +11,122 ; −1,796 | +2,521 ; −0,091 | +9,218 ; −1,394 | +12,429 ; −2,103 |
| zéros de Re T̄ | −2,447 ; −0,222 ; +1,638 | aucun | −0,262 ; +1,433 | −2,207 ; −2,170 ; −2,160 ; −0,210 ; puis 9 entre +2,69 et +2,94 |
| moyenne diagonale de M^L (meV, C14 ancien) ; effet sur Γ_T (max rel / médian) | 5 427,3 ; 5,2e-2 / 1,1e-2 | 67,0 ; 6,4e-4 / 1,8e-4 | 5 498,5 ; 6,5e-2 / 1,8e-2 | 5 417,4 ; 6,3e-2 / 1,3e-2 |

Critères par bloc (`resonance_criteria.py --blocks full,pi,sigma`, même g₀, fenêtre ±3 eV pour det/λ, toute la bande pour Friedel et Lloyd) — 9×9, 6×6 et 12×12 :

| taille | bloc | min \|det\|/max (ε − E_D) | λ_min (ε − E_D) | ∫δρ bande (Tr[t g₀′]) / Lloyd / ±3 eV (états) |
|---|---|---|---|---|
| **9×9** | complet | 1,32e-4 (−0,812) | +0,0000 + 0,0019 i (−0,812) | −1,0007 / −1,0007 / +0,669 |
| 9×9 | π | 2,33e-2 (−0,172) | +0,219 + 0,331 i, \|λ\| 0,397 (−0,170) | −0,9980 / −0,9980 / −0,285 |
| 9×9 | σ | 3,48e-4 (−0,812) | +0,0000 + 0,0019 i (−0,812) | −0,0028 / −0,0028 / +0,955 |
| 6×6 | complet | 1,43e-4 (−0,847) | +0,0001 + 0,0019 i (−0,847) | −1,0009 / −1,0009 / +0,635 |
| 6×6 | π | 3,30e-2 (−0,175) | +0,340 + 0,361 i, \|λ\| 0,496 (−0,175) | −0,9978 / −0,9978 / −0,317 |
| 6×6 | σ | 3,03e-4 (−0,847) | +0,0001 + 0,0019 i (−0,847) | −0,0030 / −0,0030 / +0,952 |
| 12×12 | complet | 1,44e-4 (−0,817) | −0,0000 + 0,0019 i (−0,817) | −1,0005 / −1,0005 / +0,697 |
| 12×12 | π | 1,83e-2 (−0,170) | +0,251 + 0,236 i, \|λ\| 0,344 (−0,127) | −0,9981 / −0,9981 / −0,261 |
| 12×12 | σ | 4,44e-4 (−0,817) | −0,0000 + 0,0019 i (−0,817) | −0,0024 / −0,0024 / +0,958 |
| 6×6 v1 | complet | 2,01e-4 (−2,415) | +0,0003 + 0,0094 i (−2,410) | −0,0341 / −0,0341 / +1,785 |
| 9×9 v1 | complet | 2,09e-4 (−2,530) | +0,0003 + 0,0108 i (−2,530) | −0,0569 / −0,0569 / +1,782 |

Pôle σ signalé (§2.3) : minimum du bloc σ le plus proche de −0,81 eV : 9×9 à −0,812 (|λ| 0,0019, |det|/max 3,5e-4 ; même énergie que le d3
de l'étape 2), 6×6 à −0,847 (|λ| 0,0019, |det|/max 3,0e-4), 12×12 à −0,817 (0,0019 ; 4,4e-4) ; c'est le minimum global de la matrice complète dans
les trois cas, et le bloc π n'a pas de zéro (|λ| ≥ 0,34). Autres minima locaux de |det|/max de la matrice complète 9×9 : −0,630 (7,5e-4 ; |λ| 0,017),
−0,255 (5,7e-4), −0,210 (5,0e-4), −0,172 (4,7e-4), −0,132 (5,2e-4). Écart ponctuel maximal entre les deux formules de δρ (Tr[t g₀′] et Lloyd) sur la
bande, 9×9 : 0,630 (complet), 0,263 (π), 0,594 (σ) états/eV ; leurs intégrales sont égales aux quatre décimales imprimées.

Γ_T à c = 0,1 % sur ±1 eV (lacune ici, azote substitutionnel dans Kaasbjerg Fig. 17 : ordre de grandeur seulement) : 9×9 min 2,25 meV (+1,00), max 37,84
(−0,18), 7,79 à E_D, 13,86 à −0,5 eV, 2,62 à +0,5 eV, ħ/Γ à ∓0,3 eV 24 / 230 fs ; 6×6 : 2,03 (+1,00), 33,47 (−0,24), 6,85 à E_D, 23 / 259 fs ; 12×12 : 2,45
(+1,00), 40,52 (−0,18), 8,43 à E_D, 26 / 213 fs (v1 9×9 : 0,63 (+0,24) / 5,55 (−1,00) / 1,27 à E_D ; 419 / 1 025 fs ; v1 6×6 : 0,68 / 6,19 / 1,37 ; 371 / 961 fs).

### 3.5 C14 redéfini (post_c14 21862375, 10 min ; `resonance_9x9_shiftL.npz`)

C ajouté uniformément à ΔV^L sur la super-cellule = C·1 sur V_loc en base de Wannier (29 mailles × 5 fonctions), R_cut 3, 240², η 0,02, 9×9 :

| variante | médiane états \|Γ_T\|·N_cells (meV) | E_res états ±1,5 eV (eV) | pic de la courbe Γ_T (eV) | max rel \|ΔΓ_T\| / médian |
|---|---|---|---|---|
| sans décalage | 3 132,87 | −0,175 | −0,180 | — |
| C = +25 meV (= −Lu, R5) | 3 187,96 | −0,175 | −0,180 | 1,92e-1 / — |
| C = −25 meV | 3 158,54 | −0,227 | −0,235 | 1,63e-1 / — |
| C14 ancien : M − ⟨M^L⟩_diag·1 (5 427 meV) | Γ_T_noshift | — | — | 5,2e-2 / 1,1e-2 |

### 3.6 Tests de la chaîne (tab:tests_M ; `M_tests_summary.csv` v2, + ligne porte A.2 ajoutée par `r6_tests_gate_row.py` dans le job final)

| test | valeur v2 | seuil | verdict | v1 |
|---|---|---|---|---|
| hermiticité M2 dense 5 / 6 / 7 / 8 / 9 / 12 | 4,1e-15 / 5,8e-15 / 8,4e-15 / 1,2e-14 / 1,2e-14 / 2,3e-14 | 1e-12 | OK | 1,4e-14 … 8,1e-14 |
| padding, k coïncidents, bandes 1–8 / 1–15, 9×9 | 5,0e-4 / 1,3e-3 | 2e-3 | OK | 5,5e-4 / 1,2e-3 |
| non-régression nbnd 16 → 20 | non rejouable (fichiers `_nb16` supprimés au ménage) | 1e-5 | fichier absent | 1,2e-7 |
| noyau dense à p = 1 vs noyau N×N, 9×9 | 1,3e-15 | 1e-10 | OK | 1,3e-15 |
| fermeture Fourier (k → R → k) ; Bloch → Wannier → Bloch | 2,4e-14 ; 1,3e-15 | 1e-12 | OK | 2,4e-14 ; 1,3e-15 |
| convention intensive : (max−min)/moyenne de max\|M\| sur N = 5, 7, 8, 9 (bandes 1–16) | **7,7e-2** | 5e-2 | **À VOIR** | 2,7e-2 (OK) |
| test d'or 5×5 dense (local vs compute_T) | 1,80e-13 | 1e-8 | OK | 2,3e-14 |
| g₀ par lots vs référence | 1,2e-14 (inchangé) | 1e-12 | OK | 1,2e-14 |
| porte A.2 (ΔV appliqué aux états purs contre M2/N_cells), 8 grossiers (5, 6, 7, 8, 9, 9 à 128 b, 10, 12) + 6 denses (5, 6, 7, 8, 9, 12) ; ajoutée le 2026-09-27 (problème 1 ci-dessous) | L ≤ 3,3e-14 / 1,2e-8 eV, NL ≤ 8,1e-15 / 4,0e-9 eV | 1e-6 eV | OK | (v1 : L = N_cells × direct, refusé) |
| reconstruction KS (fig_ks_reconstruction) | identique à v1 au bit (0,00) | — | OK | — |

La ligne « convention intensive » sort du seuil avec M2 (7,7 % contre 5 %) : max|M| des bandes 1–16 vaut 23,85 (5×5), 23,82 (7×7), 24,31 (8×8), 25,35 eV
(9×9) ; en v1 le maximum était porté par M^NL (inchangé, 6,3 eV de norme moyenne), en v2 par M^L. Rapporté, sans jugement.

### 3.7 Table de correspondance ancien → nouveau (`etape3/table_v1_v2.md`, 126 lignes ; `r6_compare_v1_v2.py`)

Colonnes : grandeur ; v1 (`results/M/`, et `results/epw/` pour le chapitre 5) ; v2 (`results/M2/`) ; fichier source ; figure, tableau ou section de
NOTES_TGAMMA / NOTES_EPW. Blocs : niveau 1 (médianes et E_res des six tailles, variantes 9×9 de grille, η et R_cut), niveau 2 (p_z–p_z, Re M^L et Re M^NL
à K), tab:L_NL, localité de M_W, tab:rcut_M, N_k^int, résonance 9×9, critères 9×9 par bloc, Friedel et Lloyd, Γ_T à c = 0,1 %, tab:tests_M (dont les deux
lignes « porte A.2 », nouvelles), C14 redéfini, Γ^ed/Γ^ep du chapitre 5. `défauts.tex` n'est pas sur Rorqual : la dernière colonne donne les étiquettes
(tab:…, fig_…) et les sections des notes ; le texte du mémoire est à reprendre par Greg avec cette table.

Principaux changements (9×9 sauf mention) :

| grandeur | v1 | v2 |
|---|---|---|
| médiane Γ·N_cells (R_cut 3, 240², η 0,02) | 2 473,55 meV | 3 132,60 meV |
| E_res (argmax des états ±1,5 eV) | −1,238 eV | −0,175 eV |
| pic −Im T̄(K) ; Re T̄(E_D) | −1,292 eV ; +2,521 eV | −0,177 eV ; +11,122 eV |
| minimum de det et λ (matrice complète) | −2,530 eV (λ = 0,0108 i) | −0,812 eV (λ = 0,0019 i), porté par le bloc σ |
| ∫δρ sur toute la bande (Friedel) | −0,057 | −1,0007 (π −0,998 ; σ −0,003) |
| Born/T médian | 3,30 | 46,5 |
| p_z–p_z sur le site de la lacune | 6,617 eV | 31,521 eV |
| ⟨‖M^NL‖_F⟩/⟨‖M^L‖_F⟩ | 20,44 (croît comme N_cells) | 0,252 (0,25 pour les six tailles) |
| tab:rcut_M, R_cut 3 | 2,27e-2 | 6,68e-2 |
| Re Σ médian (R_cut 3) | +722,3 meV | −26,65 meV |
| Γ^ed/Γ^ep médian (ch. 5 ; c = 1 %, 300 K, ±3 eV) | 0,361 | 0,537 |

### 3.7b Chapitre 5 : Γ^ed/Γ^ep et fig_epw_vs_ed (2026-09-27, nœud de connexion, accord de Greg)

`epw_ed_vs_ep.py --selfen results/epw/selfen_240_dg0.02_mv0.02_T300.npz --tag 24k24q_mv0.02` (chaîne mv0.02 de production, NOTES_EPW) avec la résonance
9×9 de M2 : médiane Γ^ed 37,462 meV (±1,2 eV : 59,298), Γ^ep 56,956 meV (±1,2 eV : 21,653 ; inchangé), rapport médian 0,537, min 0,123 à +1,875 eV,
max 32,442 à −0,180 eV, croisements (rapport = 1) à −1,504 et +0,690 eV. v1 : 0,361 ; 0,088 à +1,875 ; 2,011 à −0,755 ; −1,619 et −0,215 ; Γ^ed 23,497
(±1,2 eV : 10,737) meV. Sortie `results/M2/ed_vs_ep_24k24q_mv0.02.npz` ; les `results/epw/ed_vs_ep_*.npz` (v1, commités) sont intacts. Figure :
`etape3/figures_fix/fig_epw_vs_ed.{pdf,png}`. Non refaits : `ed_vs_ep_24k24q.npz` et `fig_epw_vs_ed_mv0.002` (ancienne chaîne EPW, Γ^ed v1).

### Figures (post_fig 21872955, puis relecture et corrections du 2026-09-27)

Le job a régénéré dans `figures/` (versions v1 sauvegardées dans `etape3/figures_v1/`, identiques au bit aux blobs de HEAD 2e49ffa) :
- qui dépendent de M2 : fig_convergence, fig_locality_final, fig_spectral_final, fig_M_map_final, fig_M_scaling_final (`make_figures_memoire.py`) ;
  fig_rcut, fig_plateau, fig_level2, fig_locality, fig_spectral, fig_M_map, fig_M_scaling (`make_figures.py`) ;
- qui ne dépendent pas de M : fig_Ved*, fig_ks_reconstruction, fig_epw_kohn_degauss (PNG identiques au bit à HEAD ; PDF différents seulement par
  /CreationDate, 6 à 7 octets).

Relecture visuelle : fig_locality_final, fig_M_map_final, fig_convergence (b), fig_level2, fig_plateau et fig_M_scaling sont lisibles. Défauts trouvés :

| figure | défaut avec M2 | suite |
|---|---|---|
| fig_M_scaling_final | vide : `set_ylim(5, 9)` calé sur v1, alors que max\|M\| vaut 23,3 à 27,3 eV avec M2 | `set_ylim(0, 1,25 max)`, comme fig_M_scaling de `make_figures.py` |
| fig_spectral_final (a) | Born hors cadre (`set_ylim(4e2, 8e4)` calé sur v1 ; Born M2 de 4,2e4 à 4,7e5 meV) ; légende sur le pic de T | bornes tirées des données (0,5 × min, 10 × max), légende au-dessus des courbes |
| fig_epw_vs_ed (ch. 5) | non régénérée par le job (problème 2) ; avec M2, la légende couvre le pic de Γ^ed | M2 et marge au-dessus des courbes (4 × max) |
| fig_convergence (a) | la légende semi-transparente couvre les points R_cut 3–4 de 5×5, 6×6 et 12×12 | non corrigé (choix de mise en page) |
| fig_spectral (figure de travail) | (a) en échelle linéaire écrasé par Born ; titres (a) et (b) qui se chevauchent ; titre (f) coupé | non corrigé |

Versions corrigées : `etape3/figures_fix/` (les 6 figures du mémoire et les 8 figures EPW, régénérées par les scripts corrigés). Seules fig_spectral_final,
fig_M_scaling_final et fig_epw_vs_ed diffèrent de `figures/` ; les sept figures EPW sans M sont identiques au bit à HEAD. Installation par Greg :
`bash etape3/install_figures_R6.sh [--pdf-dates]` (copie des trois figures ; avec l'option, PDF à date seule rétablis depuis HEAD ; puis copie versionnée
dans `article/R6_production_corrigee/etape3/figures_v2/`). Le classificateur du mode automatique a refusé à Code l'écrasement de fichiers de `figures/`.

### Problèmes d'exécution (job du 2026-09-26 au soir) et corrections (2026-09-27)

1. `r6_tests_gate_row.py` (job r6final) cherchait les json de la porte dans `etape1/gate/` (chemin de la copie `article/`) au lieu de `gate/` : aucune ligne
   ajoutée, message « porte A.2 ajoutée » quand même. Corrigé : chemin ; refus (sortie ≠ 0) si aucun json ; message avec le nombre de lignes. Relancé :
   deux lignes (8 grossiers, 6 denses) dans `M_tests_summary.csv` ; table 3.7 et copies `article/` régénérées.
2. `submit_post.sh figures` lançait `make_figures_epw.py` sans les options de production de NOTES_EPW. Conséquences : fig_epw_validation, fig_epw_gamma,
   fig_epw_phonselfen et fig_epw_decay réécrites avec la chaîne mv0.002 (rétablies par Greg par `git checkout`, 2026-09-27 14:39) ; fig_epw_vs_ed non produite
   (`selfen_prod_T300.npz` absent). Corrigé avec l'accord de Greg :
   - `submit_post.sh` passe les options (`--prod-tag 240_dg0.02_mv0.02 --phself-tag path_1200_dg0.02_mv0.02 --val-tag 24k24q_mv0.02 --sel-suffix _mv0.02 --control
     --kohn-val-tags 24k24q,24k24q_mv0.02 --dfpt-tag 24k24q`) et lance `epw_ed_vs_ep.py` ;
   - `epw_ed_vs_ep.py` écrit dans `results_dir` (sinon il aurait écrasé `results/epw/ed_vs_ep_*.npz`, qui sont v1) ;
   - `.gitignore` : `!results/M2/ed_vs_ep_*.npz`.

   Contrôle : avec ces options, les sept figures EPW sans M sortent identiques au bit à HEAD.
3. Limites d'axe calées sur v1 dans `make_figures_memoire.py` (deux) et place de la légende de fig_epw_vs_ed dans `make_figures_epw.py` : voir Figures.

### Manifestes

- Aucune suppression pendant R6. `results/M/` intact : README déposé (« résultats obtenus avec M^L non normalisé (facteur N_cells manquant), remplacés par
  R6 »), en lecture seule (`chmod -R a-w`) depuis 3.0. Le droit d'écriture (`chmod -R u+w results/M`) sera rétabli sur GO.
- `results/M2/` (sur `/project`) : 43 fichiers M (47,3 Go réels, plus 20,8 Go de liens `M_NL_dense_*` vers `results/M/`). `MD5SUMS_2026-09-25.txt` relu le
  2026-09-26 par le job r6final : 43 OK sur 43. Produits de l'étape 3 :
  - suivis par les exceptions `.gitignore` : `specwd_*_prod.npz` (6), `resonance_*.npz` (6x6, 9x9, 12x12, 9x9_shiftL, criteria ×3), `mwr_locality.npz`,
    `M_analysis.npz`, `ks_reconstruction.npz`, `ved_analysis.npz` (copie v1, indépendant de M), `ed_vs_ep_24k24q_mv0.02.npz`, `*.csv` (9) ;
  - non suivis : `resigma_9x9_*.npz` (55 Mo), `logs/`.

  Le miroir des M2 denses prévu en phase 0 n'est pas fait : décision à Greg.
- `.save` à 128 bandes de R5 (`qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.save`, scratch, 1,5 Go) miroité dans `qe_tmp_backup/R5_uc9x9_nb128/`
  (`MD5SUMS_2026-09-26.txt`, `md5sum -c` : 0 écart). La copie du scratch n'est pas supprimée (GO séparé).

### État du dépôt (rien n'est commité par Code)

- Modifiés depuis 2e49ffa : `scripts/analyze_M.py` (échantillon `ks` hors du bloc nb16, `GOLDEN_RESULT`), `scripts/resonance_criteria.py` (`--blocks`,
  `--flag-eV`), `scripts/resonance_metrics.py` (`--shift-L-meV`), `scripts/epw_ed_vs_ep.py` (sortie dans `results_dir`), `scripts/make_figures_memoire.py` et
  `scripts/make_figures_epw.py` (limites d'axe), `.gitignore` (`ed_vs_ep`), `article/R6_production_corrigee/R6_rapport.md`.
- Nouveaux : `scripts/submit_post.sh` ; `results/M2/` (produits ci-dessus) ; `article/R6_production_corrigee/etape3/` (csv_v2, figures_v2, table_v1_v2.md,
  level1_gate.md, scripts R6).
- `figures/` : 33 fichiers modifiés par le job ; après `install_figures_R6.sh --pdf-dates`, restent les figures qui dépendent de M2.
- Documentation qui cite encore v1 : `NOTES_TGAMMA.md` §2 et §6, `CLAUDE.md` (« Données : results/M/ »), NOTES_EPW (Γ^ed/Γ^ep 0,361) : à mettre à jour
  avec la table 3.7 (par Greg, ou par Code sur demande).
- Les modifications de `article/R7_tailles_3m/` visibles dans `git status` viennent de R7, pas de R6.

**STOP — étape 3 terminée le 2026-09-27.** Reste à Greg : installer les figures corrigées (`etape3/install_figures_R6.sh`), commit, GO pour rétablir
l'écriture sur `results/M/`, décision sur le miroir des M2 denses, mise à jour de NOTES_TGAMMA, CLAUDE.md, NOTES_EPW et `défauts.tex` avec la table 3.7.

