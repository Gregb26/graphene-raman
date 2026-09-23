# Inventaire du projet sur Rorqual — 2026-09-16 (P12, lecture seule)

Session strictement en lecture seule : aucune modification, aucun déplacement.
Seule écriture : ce fichier. Non commité.

Périmètre : `~/links/projects/rrg-cotemich-ac/gregb26/` (= `/project/rrg-cotemich-ac/gregb26/`)
et `$SCRATCH` = `/scratch/gregb26` (= `/lustre10/scratch/gregb26`). Aucun répertoire
d'un autre utilisateur n'a été consulté. Les liens `~/links/{projects,scratch,nearlines}`
pointent vers ces volumes ; le nearline (269 Go, niveau groupe) n'a pas été inspecté.

Tailles = tailles apparentes (`du --apparent-size`, Gio). Dates = mtime / atime des
fichiers. Attention : Lustre met l'atime à jour paresseusement (relatime), une
date d'accès ancienne n'exclut pas une lecture récente via `mmap`, mais c'est
l'atime que la purge regarde.

---

## 0. En-tête : occupation, quotas, urgence scratch

### Quotas (`diskusage_report`, 2026-09-16)

| Espace | Occupé | Quota | Fichiers |
|---|---|---|---|
| /home (gregb26) | 4,4 Go | 50 Go | 62 k / 500 k |
| /scratch (gregb26) | 107 Go (113,3 Gio apparent) | 20 To | 15 775 / 1 M |
| /project rrg-cotemich-ac (groupe) | 40 To | 180 To | 4,7 M / 15 M |
| /project def-cotemich (groupe) | 12 Go | 10 To | 353 k / 500 k |
| /nearline rrg-cotemich-ac (groupe) | 269 Go | 100 To | 25 / 5000 |

Part de gregb26 dans /project rrg-cotemich-ac (somme des `du`) : **≈ 1 931 Gio ≈ 1,93 Tio,
≈ 115 600 fichiers** (hors contenu des `.git`). Détail : `graphene/` 1 763,5 G,
`ab-initio-defects/` 160,8 G, `codes/` 2,8 G, `abinit_assignment/` 2,4 G, `dft/` 1,2 G,
`abinit_processing/` 0,02 G, `ai/` 0.

### Politique de purge du scratch

- Pas d'information dans le motd de Rorqual. La page `docs.alliancecan.ca/wiki/Scratch_purging_policy`
  n'est pas accessible depuis le nœud de connexion (bloquée par un anti-robot Anubis) ;
  non vérifiée dans cette session.
- Règle Alliance telle que je la connais (à confirmer avec les courriels de notification) :
  fichiers dont **atime ET ctime datent de plus de 60 jours** ; liste des candidats
  envoyée en début de mois, suppression le 15 du mois.
- `/scratch/.purged/lustre10/scratch/` existe (mécanisme en place) mais est **vide** au
  2026-09-16 : aucune liste de candidats visible. Les fichiers ci-dessous ont déjà plus de
  60 jours et sont toujours là au lendemain d'un 15 : soit la purge n'est pas encore
  active sur Rorqual, soit la liste n'est pas exposée. **Ne pas compter dessus.**

### URGENT — sur le scratch et NULLE PART AILLEURS, déjà éligibles à la purge (atime = ctime ≤ 2026-06-18)

Les fonctions d'onde de toute la chaîne M vivent uniquement sur le scratch : le dépôt
n'en a que des **liens symboliques** (`ab-initio-defects/data/graphene/**` → 1 189 liens vers
`/home/gregb26/links/scratch/qe_tmp/…`) et les `outdir` des entrées QE dans
`graphene/qe/defects/**` pointent tous vers `qe_tmp/`. Le projet ne contient que
`scf.in/out`, `pp.out` et les potentiels `Vks_*` (pp.x).

| Répertoire scratch | Taille | atime des wfc | Contenu | Coût de recalcul |
|---|---|---|---|---|
| `qe_tmp/defect_{5..12}x{5..12}_{p,d}/*.save` (16 dirs) | **43,2 G** | 2026-06-17/18 | SCF supercellules pristine/défaut : `wfc1.hdf5`, `charge-density.hdf5`, `data-file-schema.xml`, `C.upf` | 16 SCF supercellule (jusqu'à 12×12 : heures-cœur ×centaines) |
| `qe_tmp/defect_unit_cell_{10x10,11x11}` (.save) | 1,36 G | 2026-08-31 | SCF/NSCF maille unitaire grossière 10×10, 11×11 | modéré |
| `qe_tmp/graphene_scf` (.save 1,95 G) | 4,0 G | 2026-05-20 (accédé 06-18) | SCF graphène pristine utilisé par `graphene/qe/{bands,nscf,dos}.in` | faible (SCF maille unitaire) |

Fichiers > 1 Go concernés (tous atime = 2026-06-17/18) : `defect_12x12_{p,d}.save/wfc1.hdf5`
(6,2 G ×2), `11x11` (4,4 G ×2), `10x10` (3,05 G ×2), `9x9` (2,04 G ×2), `8x8` (1,3 G ×2).

Encore accédés récemment (donc pas encore éligibles, mais **uniques eux aussi**) :

| Répertoire scratch | Taille (.save) | atime | Contenu | Coût |
|---|---|---|---|---|
| `qe_tmp/defect_uc_dense_{24,25,27,28,32}/*.save` | 11,0 G (+12,2 G de `.wfcN` en vrac, doublons) | 2026-09-07/14 | NSCF denses nbnd=20 (chaîne de production M, `config.dense_paths`) | 5 NSCF denses (heures) |
| `qe_tmp/defect_unit_cell_{5..9}x…, 12x12` (.save) | ≈ 1,0 G | 2026-09-07/09 | mailles unitaires grossières (chaîne grossière, golden test) | faible |
| `qe_tmp/defect_unit_cell_{5x5,7x7,8x8}_wann` | 0,68 G | 2026-09-03 | pw2wannier ; les sorties Wannier sont dans le projet | faible |

Total « irremplaçable, uniquement sur scratch » : **≈ 57 G** (43,2 + 11,0 + ≈ 3,1).
Le reste du scratch (≈ 56 G : `qe_conv/` 44,9 G et fichiers `.wfcN` en vrac) est
reproductible ou redondant (voir §2 et §3c).

---

## 1. Arborescence (≈ 3 niveaux)

### 1.1 Premier niveau de `~/links/projects/rrg-cotemich-ac/gregb26/`

| Répertoire | Taille | Fichiers | Contenu le plus récent (mtime) | Dernier accès |
|---|---|---|---|---|
| `graphene/` | 1 763,5 G | 79 629 | 2026-09-14 | 2026-09-16 |
| `ab-initio-defects/` (dépôt git, `main`, dernier commit ab98f6b 2026-09-15) | 160,8 G | 16 885 | 2026-09-16 | 2026-09-16 |
| `codes/` | 2,8 G | 16 700 | 2026-01-19 | 2026-09-16 |
| `abinit_assignment/` | 2,4 G | 208 | 2026-05-28 | 2026-05-29 |
| `dft/` | 1,2 G | 701 | 2025-11-24 | 2025-11-26 |
| `abinit_processing/` | 0,02 G | 1 467 | 2026-01-19 | 2026-09-14 |
| `ai/` | 0,0 G | 10 | 2025-10-17 | 2026-02-04 |
| `submit.sbatch` | 732 o | 1 | 2025-11-24 | — |

### 1.2 `graphene/` (1 763,5 G)

| Sous-répertoire | Taille | Fichiers | mtime max | Note |
|---|---|---|---|---|
| `graphene/qe/epw/` | 1 289,6 G | 57 332 | 2026-09-14 | EPW ; 24k-24q = production |
| ├ `24k-24q/` | 179,3 G | ≈ 12 000 | 2026-09-14 | production (§3a) |
| ├ `36k-30q/` | 581,4 G | 11 775 wfc | 2026-04-30 | test avril ; `phonons/_ph0` |
| ├ `30k-24q/` | 263,6 G | 7 899 wfc | 2026-04-30 | test avril ; `phonons/_ph0` |
| ├ `16k-8q/` | 151,7 G | 9 360 wfc | 2026-04-29 | dfpt_g_G/K 72 G ×2, epw_g_*, validation_16k8q |
| ├ `24k-12q/` | 45,3 G | | 2026-04-28 | `phonons/_ph0` |
| ├ `16k-16q/` | 35,4 G | | 2026-04-28 | `phonons/_ph0` |
| ├ `16k-12q/` | 20,4 G | | 2026-04-29 | `phonons/_ph0` |
| ├ `12k-12q/` | 11,7 G | | 2026-04-28 | `phonons/_ph0` |
| ├ `30k-30q/` | 0,45 G | | 2026-04-29 | `_ph0` supprimé 2026-09-10 |
| ├ `36k-36q/`, `8k-8q/` | 0,31 G, 0,03 G | | 2026-04-29 | |
| └ `NOTES_EPW.md`, `CLAUDE.md` | | | 2026-09-14 | |
| `graphene/abinit/` | 420,7 G | ≈ 9 250 | 2026-02-02 | ère ABINIT (nov. 2025 – févr. 2026), remplacée par QE |
| ├ `pristine/unit_cell/` | 396,6 G | 8 780 | | `wannier/bands` 375 G (WFK.nc), `wannier/pi_bands` 11,7 G |
| ├ `pristine/supercell/` | 6,5 G | | | |
| └ `defective/` | 17,6 G | 436 | 2025-12-21 | 5x5…10x10, vacancies |
| `graphene/qe/test/` | 31,5 G | 11 910 | 2026-04-21 | `_ph0` 30,8 G (test ph.x 16×16 q, dvscf 0,46 G), `tutorial01` 0,6 G |
| `graphene/qe/defects/` | 21,7 G | 821 | 2026-09-07 | entrées/sorties QE + `Vks_*` (3,6 G) + Wannier ; `unit_cell/` 15,4 G, `super_cell/` 6,3 G |
| `graphene/qe/convergence/` | 3,6 M | | 2026-05-27 | ecut, kpoint, smearing, phonons (sorties seulement ; `outdir` sur scratch `qe_conv/`) |
| `graphene/qe/` (racine) | ≈ 40 M | | 2026-05-20 | scf/nscf/bands/dos, dyn*, ifc, modes du graphène pristine (janv. 2026) |

### 1.3 `ab-initio-defects/` (160,8 G)

| Sous-répertoire | Taille | Fichiers | mtime max | Note |
|---|---|---|---|---|
| `results/M/` | 132,4 G | 540 | 2026-09-15 | matrices M (§3) ; `obsolete_grid_7x7/` 15 G, `_test_mnl/` 52 M, `logs/` 2,3 M |
| `jobs/` | 27,9 G | ≈ 150 | 2026-02-03 | ère ABINIT (`wfk_*.nc`, `pot_*.nc`, `M_test.npy`) ; ignoré par git |
| `.venv/` | 0,29 G | 15 838 | 2026-09-07 | |
| `.git/` | 0,11 G | | 2026-09-16 | |
| `wannier/` | 60 M | | 2026-09-07 | 12 grilles (5x5…32x32, `_nb16`) ; `24x24`, `*_nb16` non suivis par git |
| `figures/` | 15 M | 53 | 2026-09-15 | 27 figures (pdf+png) |
| `results/epw/`, `results/test_A/` | 2,3 M, 0,3 M | | 2026-09-14 | npz commités |
| `data/` | (liens) | 1 189 liens | 2026-06-18 | → scratch `qe_tmp/` |
| `scripts/` (71), `src/`, `tests/`, `config/`, `notebooks/` (1,6 M) | < 5 M | | | |

### 1.4 Autres

| Répertoire | Taille | Note |
|---|---|---|
| `codes/qe/` (1,06 G dont `.git` 0,80 G), `codes/abinit/` (1,14 G dont `.git` 0,62 G), `codes/wannier90/` (0,57 G, `libwannier90.a` compilé) | 2,8 G | sources ; 1 lien cassé : `codes/wannier90/test-suite/library-mode-test/ref/gaas.win` |
| `dft/{compmatphys,copper,diamond,iron}` | 1,2 G | tutoriels nov. 2025 |
| `abinit_assignment/` | 2,4 G | devoir GaAs mai 2026 (WFK, DEN, GSR) |
| `abinit_processing/` | 0,02 G | `pseudo/` (**référencé par ≥ 42 entrées QE comme `pseudo_dir`**), `scripts/`, `venv/` (août 2025), `cif/` vide |
| `ai/IonicConductivityFlow/` | 149 K | dépôt git, oct. 2025 |

### 1.5 Scratch `/scratch/gregb26` (113,3 G, 15 775 fichiers)

| Répertoire | Taille | Fichiers | mtime max | atime max |
|---|---|---|---|---|
| `qe_tmp/` | 68,4 G | 6 279 | 2026-09-07 | 2026-09-14 |
| ├ `defect_{N}x{N}_{p,d}/` (16) | 43,2 G | 5 chacun | 2026-06-18 | 2026-09-09 (xml seulement ; wfc : 06-18) |
| ├ `defect_uc_dense_{24,25,27,28,32}/` | 23,3 G (dont .save 11,0 G) | 643–1 091 | 2026-09-07 | 2026-09-14 |
| ├ `defect_unit_cell_{5..12}x…/` (8) | 2,3 G | 29–248 | 2026-06-19 | 2026-09-09 |
| ├ `defect_unit_cell_{5x5,7x7,8x8}_wann/` | 0,68 G | | 2026-09-02 | 2026-09-03 |
| └ `graphene_scf/` | 4,0 G | 949 | 2026-05-20 | 2026-06-18 |
| `qe_conv/` | 44,9 G | 9 492 | 2026-05-29 | 2026-05-29 |
| ├ `phonons/` | 34,8 G | 1 043 | 2026-05-29 | 2026-05-29 |
| │  └ `_ph0/` | 32,5 G | 128 `.wfcN` + 233 xml phsave | | pas de dvscf, pas de q_* |
| ├ `smearing/` (24 runs) | 8,0 G | 6 600 | 2026-05-19 | 2026-05-19 |
| ├ `kpoint/` (5 runs) | 1,2 G | 969 | 2026-05-19 | 2026-05-19 |
| └ `ecut_{40..130}_k{12,24}/` (20 runs) | 0,9 G | | 2026-05-19 | 2026-05-19 |
| `compmatphys/basic/` | 0 | 4 | 2025-10-14 | 2025-10-17 |
| `jobs/`, `graphene/` | 0 | 0 | 2026-02-03, 2026-07-16 | répertoires vides |

---

## 2. Scratch — détail et diagnostic

Éligible à la purge (atime **et** ctime < 2026-07-18, 60 jours) : **≈ 88 G, 10 500 fichiers**
(`find -atime +60 -ctime +60`). Ventilation :

| Bloc | Taille | Statut |
|---|---|---|
| `qe_tmp/defect_{N}_{p,d}` (SCF supercellules) | 43,2 G | **IRREMPLAÇABLE, unique** — rapatrier |
| `qe_tmp/defect_unit_cell_{10x10,11x11}` | 1,4 G | unique ; utilisé par la chaîne grossière (liens `data/`) |
| `qe_tmp/graphene_scf` | 4,0 G | reproductible (SCF maille unitaire ; entrées `graphene/qe/scf.in`) ; 2,07 G de `.wfcN` en vrac redondants |
| `qe_conv/phonons` | 34,8 G | reproductible ; résultats (dyn, q2r, matdyn) déjà dans `graphene/qe/convergence/phonons/` ; `_ph0` = 32,5 G de wfc de reprise ph.x |
| `qe_conv/smearing`, `kpoint`, `ecut_*` | 10,1 G | reproductible ; sorties de convergence dans `graphene/qe/convergence/` |

Ce qui est sur le scratch et existe aussi dans le projet : rien de volumineux. Les seuls
« doublons » sont internes au scratch : les `prefix.wfcN` en vrac (copies par processus MPI)
à côté de `prefix.save/wfcN.hdf5` — 12,2 G dans `defect_uc_dense_*`, 2,1 G dans
`graphene_scf`, 1,0 G dans `defect_unit_cell_*`, 33,7 G dans `qe_conv/phonons(+_ph0)`.

Qui lit le scratch (voir §4) : `link_data.sh` (construit `data/`), `config.dense_paths()`,
`check_M_dense_vs_coarse*.py`, `check_onsite_and_NL.py`, `test_local_tmatrix_real.py`,
`compute_M_dense_stages.py`, et tous les `outdir` des entrées QE de `graphene/qe/defects/`
et `graphene/qe/{bands,nscf,dos}.in`.

---

## 3. Classification

### 3a. IRREMPLAÇABLE (heures-cœur à refaire)

| Emplacement | Taille | Contenu |
|---|---|---|
| scratch `qe_tmp/defect_{N}_{p,d}/*.save` | 43,2 G | SCF supercellules 5×5…12×12 (voir §0) |
| scratch `qe_tmp/defect_uc_dense_D/*.save` | 11,0 G | NSCF denses nbnd=20, D = 24/25/27/28/32 |
| scratch `qe_tmp/defect_unit_cell_N/*.save` | ≈ 2,3 G | mailles unitaires grossières (nb faible mais chaîne grossière + tests) |
| `results/M/M_{,L_,NL_}dense_{5,6,7,8,9,12}x…npy` (nb20) | 58 G | matrices M denses de production (Hartree, sidecars .json) ; 8×8 = 4 × 6,25 G |
| `results/M/M_ed_*.npy`, `M_L_*`, `M_NL_*` (grossières) | 0,5 G | M grossières (compute_M.py ; recalcul minutes à ~1 h selon N) |
| `graphene/qe/epw/24k-24q/phonons/save/` + `graphene.dyn*` | 0,97 G | DFPT 24×24 q : 61 dvscf + dyn (`epw_pp_save.py`) |
| `graphene/qe/epw/24k-24q/epw1/` : `epwdata.fmt`, `graphene.epmatwp` (0,81 G), `wigner/vmedata/dmedata/crystal.fmt`, `.chk/.mmn/.ukk`, `epw1.out` | ≈ 1,1 G | produits de P1 ; les runs epw2/epw6 y pointent par liens relatifs `../epw1/*.fmt` |
| `graphene/qe/epw/24k-24q/epw1/graphene.epb{1..16}` | 30,4 G | matrice e-ph Bloch (intermédiaire ; regénérable en relançant epw1 depuis `save/`, plusieurs heures) |
| `graphene/qe/epw/24k-24q/dfpt_g_{G,K}/` : `ph.out` (prt), dvscf (0,48 G), dyn xml, `nscf.out` | ≈ 2,3 G ×2 | vérification D² DFPT (§5.4) — **hors** les 66,3 G de wfc de chaque dir |
| `graphene/qe/epw/24k-24q/epw2_selfen_*`, `epw6_phself_*`, `epw_g_{G,K,ring}`, `bands/` (sorties) | ≈ 0,1 G | sorties EPW (elecselfen 240², phonselfen 1200², prtgkk) — petites mais des heures |
| `graphene/qe/defects/**` : `scf.in/out`, `nscf.in/out`, `Vks_*` (24 fichiers, 3,6 G), `wannier.*` (chk, mmn, amn, tb.dat, hr.dat) | ≈ 15 G | provenance + potentiels pp.x + Wannier ; les `pp.out` (184 M ×16 ≈ 3 G) sont le stdout de pp.x (reproductibles depuis les .save) |
| `graphene/qe/` racine : `graphene.dyn*`, `.ifc`, `.modes`, `.freq`, `.dos`, `bands.dat` | ≈ 40 M | phonons/bandes du graphène pristine (janv./mai 2026) |
| `ab-initio-defects/wannier/*` | 60 M | tb.dat, hr.dat, chk par grille (5 non suivis : `24x24`, `*_nb16`) |
| `results/epw/*.npz`, `results/M/specwd_*_prod.npz`, `resonance_*.npz`, csv | < 5 M | commités (produits de production) |
| `graphene/qe/epw/16k-8q/` : `epw_g_{G,K}` (1,3 G ×2), `dfpt_g_*/ph.out+dvscf`, `validation_16k8q` | ≈ 5 G | grille de validation (results/epw/validation_16k8q.npz) — **hors** les 145 G de wfc |

### 3b. REPRODUCTIBLE en < 1 h par un script du dépôt

| Emplacement | Taille | Script régénérateur |
|---|---|---|
| `results/M/M_dense_*_norm.npy` (6 fichiers nb20 + 4 nb16) | 19,4 G + 9,9 G | `scripts/migrate_M_norm.py` (copie renormalisée `supercell` de `M_dense_*`) |
| `results/M/M_ed_*_norm.npy` | 0,25 G | `scripts/migrate_M_norm.py` |
| `results/M/specwd_*_prod.npz` | | `scripts/compute_spectral_wannier.py` via `submit_spectral_wannier_dense.sh` (~1 h pour R_cut 3, cf. mémoire) |
| `results/M/resonance_*.npz`, `resonance_criteria_*.npz` | | `resonance_metrics.py`, `resonance_criteria.py` |
| `results/M/M_analysis.npz`, `ved_analysis.npz`, `mwr_locality.npz`, `ks_reconstruction.npz` | | `analyze_M.py`, `analyze_Ved.py`, `mwr_locality_coarse_vs_dense.py`, `ks_reconstruction_all.py` |
| `results/M/level1_summary.csv`, `level2_*.csv`, `lnl_frobenius.csv`, `sampling_table.csv`, `m_rcut_*.csv`, `convergence.npz/png` | | `summarize_level1_maps.py`, `level2_families.py`, `lnl_frobenius_all.py`, `sampling_table.py`, `m_rcut_convergence.py`, `compute_convergence.py` |
| `results/M/dos_*.npz/png`, `gamma_*.npz/png` | | `compute_tmatrix.py`, `compute_spectral.py` |
| `results/M/resigma_9x9_*.npz` | | `rcut_resigma.py` (`submit_rcut_resigma.sh`) |
| `results/M/*_coarsecheck.npy`, `_ML_7x7_ref_b4.npy`, `_test_mnl/` | 0,1 G | `check_ML_coarse_kernel.py`, `validate_ML_grid_7x7.py`, `_test_mnl_mpi.sh` |
| `results/epw/*.npz` | 2,3 M | `epw_selfen_post.py`, `epw_phself_post.py`, `epw_validate.py`, `epw_d2_extract.py`, `epw_ring_check.py` |
| `results/test_A/` | 0,3 M | `run_test_A_batch.py` (`submit_test_A.sh`) |
| `figures/` | 15 M | `make_figures.py`, `make_figures_memoire.py`, `make_figures_epw.py` |
| `ab-initio-defects/data/` (liens) | — | `scripts/link_data.sh NxN` |
| `graphene/qe/defects/**/pp.out` | ≈ 3 G | pp.x (`submit.pp`) — nécessite les `.save` du scratch |
| `graphene/qe/epw/24k-24q/phonons/save/` | 0,97 G | `epw_pp_save.py` à partir de `phonons/_ph0` — **mais `_ph0` a été supprimé le 2026-09-10** : `save/` est donc irremplaçable (classé 3a) |
| scratch `qe_tmp/graphene_scf`, `qe_conv/*` | 49 G | `graphene/qe/submit.scf`, `graphene/qe/convergence/*/submit.*`, `phonons.sh` (résultats déjà dans le projet) |

### 3c. JETABLE (candidats — rien n'a été supprimé)

| Emplacement | Taille | Pourquoi |
|---|---|---|
| `graphene/qe/epw/36k-30q/phonons/_ph0` | 579,6 G (wfc) | grille test avril 2026, jamais accédée depuis 2026-05-05 ; dvscf 1,4 G et dyn à part |
| `graphene/qe/epw/30k-24q/phonons/_ph0` | 262,4 G (wfc) | idem |
| `graphene/qe/epw/16k-8q/dfpt_g_{G,K}/` wfc | 144,7 G | wfc des runs prt ; garder `ph.out`, dvscf, dyn, `epw_g_*` (≈ 5 G) |
| `graphene/qe/epw/24k-12q`, `16k-16q`, `16k-12q`, `12k-12q` (`phonons/_ph0`) | 45 + 35 + 20 + 12 = 112 G | grilles test avril (non arbitrées, cf. NOTES_EPW) |
| `graphene/qe/epw/24k-24q/dfpt_g_{G,K}/` wfc + `_ph0` wfc | 66,3 G ×2 | fonctions d'onde des runs prt terminés (garder les sorties, ≈ 2,3 G ×2) |
| `graphene/qe/epw/24k-24q/epw1/*.wfc*`, `graphene.save/wfc*.hdf5` | 5,3 + 2,7 G | wfc NSCF de epw1 (regénérables ; les .epb et .fmt suffisent à epw2/epw6) |
| `graphene/abinit/pristine/` | 403,1 G | ère ABINIT (nov. 2025 – févr. 2026), méthode abandonnée ; 375 G de `graphene_w90o_DS4_WFK.nc` dont 212 G de doublons exacts (§5) |
| `graphene/abinit/defective/` | 17,6 G | idem ; `10x10/defect_not_relaxed/test{,/test2}` dupliqués |
| `ab-initio-defects/jobs/` | 27,9 G | wfk/pot ABINIT déc. 2025 – févr. 2026 (`M_test.npy`, `wfk_*.nc`) ; ignoré par git |
| `results/M/obsolete_grid_7x7/` | 15,0 G | déjà étiqueté obsolète ; comparé par `_old_vs_new_7x7.py` |
| `results/M/*_nb16.npy` (16 fichiers) | 39,5 G | artefacts nbnd=16 conservés pour trace ; `check_M_dense_nb20_vs_nb16.py` les lit |
| `results/M/M_dense_*_norm.npy` | 29,3 G | copies renormalisées, regénérables (3b) |
| `graphene/qe/test/` | 31,5 G | `_ph0` d'un test ph.x 16×16 q (avril) + tutorial01 ; dvscf 0,46 G |
| `abinit_assignment/` | 2,4 G | devoir de cours mai 2026 (WFK, DEN) |
| `dft/` | 1,2 G | tutoriels nov. 2025 |
| `codes/{qe,abinit,wannier90}` | 2,8 G | sources + `.git` (1,4 G) ; recompilables — **vérifier** que `codes/wannier90/libwannier90.a` n'est pas lié par un module utilisé |
| `abinit_processing/venv/` | ~0,02 G | venv août 2025 ; **ne pas toucher `abinit_processing/pseudo/`** (pseudo_dir de ≥ 42 entrées QE) |
| `ab-initio-defects/hamiltonian_reconstruction.png` (racine), `scripts/_test_mnl_mpi.sh` (non suivi), `.pytest_cache` | < 1 M | vrac |
| scratch `qe_conv/` | 44,9 G | convergences mai 2026, résultats dans le projet |
| scratch `.wfcN` en vrac (hors `.save`) dans `qe_tmp/defect_uc_dense_*`, `graphene_scf`, `defect_unit_cell_*` | 15,3 G | copies par processus, redondantes avec `.save/wfcN.hdf5` |
| scratch `compmatphys/`, `jobs/`, `graphene/` | 0 | vides (oct. 2025 / févr. 2026) |
| `/home/gregb26/.vscode-server` (info) | 5,2 G | cache VS Code (hors périmètre) |

Total « jetable » projet ≈ 1 550 G ; scratch ≈ 60 G.

---

## 4. Chemins codés en dur (carte de ce qu'un déplacement casserait)

### 4.1 Racines absolues

| Racine | Où | Effet d'un déplacement |
|---|---|---|
| `/home/gregb26/links/scratch/qe_tmp/` | `src/electron_defect_interaction/config.py:22` (`dense_paths`, défaut), `scripts/link_data.sh`, `check_M_dense_vs_coarse{,_bands}.py`, `check_onsite_and_NL.py`, `test_local_tmatrix_real.py`, `compute_M_dense_stages.py`, 1 189 liens `data/**`, tous les `outdir` de `graphene/qe/defects/**/*.in` (16 supercellules, 5 denses, 12 mailles unitaires), `graphene/qe/{bands,nscf,dos}.in` (`qe_tmp/graphene_scf`) | casse la chaîne M dense et grossière, les tests réels, `data/` |
| `/home/gregb26/links/scratch/qe_conv/` | `graphene/qe/convergence/**/*.in` (phonons, smearing ×72, kpoint, ecut) | seulement une relance des convergences |
| `/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/ab-initio-defects` | tous les `scripts/submit_*.sh` (12), `_test_mnl_mpi.sh`, `link_data.sh` | `cd` en tête des jobs SLURM |
| `…/gregb26/graphene/qe/epw/24k-24q` | `epw_d2_extract.py`, `epw_ring_check.py` (`epw_g_ring`), `submit_epw_p1_post.sh` | post-traitement EPW |
| `…/gregb26/graphene/qe/defects/unit_cell`, `…/ab-initio-defects/wannier` | `finalize_wannier.py` | import des Wannier |
| `…/gregb26/graphene/qe/defects` | `link_data.sh` (`Vks_*`) | |
| `…/gregb26/abinit_processing/pseudo` | `pseudo_dir` de 42 entrées `graphene/qe/epw/**/*.in` (scf, nscf, bands) + `graphene/qe/convergence/phonons/scf.in` (+ probablement `graphene/qe/defects/**`, non compté) | déplacer `abinit_processing/` casse toute relance QE |
| `EDI_DATA` (env, défaut `data/graphene`) | `scripts/_paths.py` | |

### 4.2 Chemins relatifs (depuis la racine du dépôt) par script

| Script | Chemins référencés |
|---|---|
| `analyze_M.py` | `data/graphene`, `results/M/M_analysis.npz`, `M_dense_{N}[_nb16].npy`, `M_ed_{N}.npy`, `M_L_{N}.npy`, `M_L_dense_{N}_coarsecheck.npy`, `M_tests_summary.csv` |
| `analyze_Ved.py` | `data/graphene`, `results/M/ved_analysis.npz` |
| `check_M_dense_nb20_vs_nb16.py` | `results/M/M_dense_{N}[_nb16].npy` |
| `check_M_dense_vs_coarse.py`, `_bands.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, scratch `defect_uc_dense_{D}.save`, `results/M/M_L[_dense]_{N}.npy` |
| `check_ML_coarse_kernel.py` | `results/M/M_L_{N}.npy`, `M_L_dense_{N}_coarsecheck.npy` |
| `check_onsite_and_NL.py` | scratch dense `.save`, `results/M/M_{,L_,NL_}dense_{N}.npy`, `wannier/{D}x{D}` |
| `compare_bands_qe.py`, `compare_bands_w90_qe.py` | `results/compare_bands_*.png` |
| `compute_convergence.py` | `results/M/{prefix}_*x*.npz`, `convergence.npz/png` |
| `compute_M_dense_stages.py` | `data/graphene/supercell/qe/defect_{N}_{p,d}.save[/Vks_{N}_{p,d}]`, `data/graphene/unit_cell/qe/defect_{N}.save[/C.upf]`, scratch `qe_tmp` |
| `compute_M.py` | `results/M/M_ed_NxN.npy`, `M_L_NxN.npy`, `M_NL_NxN.npy` |
| `compute_spectral.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, `results/M/gamma_{N}.npz/png`, `M_ed_{N}_norm.npy` |
| `compute_spectral_wannier.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, `results/M/M_ed_{N}.npy` |
| `compute_tmatrix.py` | idem + `results/M/dos_{N}.npz/png`, `M_ed_{N}_norm.npy` |
| `epw_d2_extract.py` | abs. `24k-24q`, `results/epw/{d2_extract_24k24q,phself_path_1200_dg0.02,validation_24k24q}.npz` |
| `epw_phself_post.py`, `epw_selfen_post.py`, `epw_validate.py` | `results/epw/{phself,selfen,validation}_*.npz` |
| `epw_ring_check.py` | abs. `24k-24q/epw_g_ring`, `results/epw/ring_check_24k24q.npz` |
| `_eta_scan.py`, `_mcheck.py`, `_normtest.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, `results/M/M_ed_{N}[_norm].npy` |
| `finalize_wannier.py` | abs. `ab-initio-defects/wannier`, abs. `graphene/qe/defects/unit_cell` |
| `ks_reconstruction_all.py` | `data/graphene`, `results/M/ks_reconstruction.npz` |
| `level2_families.py` | `results/M/{level2_families.csv,M_analysis.npz,mwr_locality.npz,specwd_{N}_prod.npz}` |
| `link_data.sh` | abs. ROOT/PROJ/SCRATCH/GRAPHENE ; `data/graphene/{unit_cell,supercell}/qe/defect_*` |
| `lnl_frobenius_all.py` | `results/M/lnl_frobenius.csv` |
| `make_figures.py` | `figures/memoire.mplstyle`, `results/M/{ks_reconstruction.npz,level1_summary.csv,level2_summary.csv,M_analysis.npz,mwr_locality.npz,resonance_NxN.npz,specwd_*.npz,ved_analysis.npz}` |
| `make_figures_memoire.py` | `figures/memoire.mplstyle`, `results/M/{M_analysis,mwr_locality,resonance_criteria_*,resonance_*,specwd_*_prod,ved_analysis}.npz` |
| `make_figures_epw.py` | `figures/memoire.mplstyle`, `results/epw/{phself_*,selfen_*_T*,validation_24k24q}.npz`, `results/M/resonance_NxN.npz` |
| `migrate_M_norm.py` | `results/M/M_ed_*`, `M_ed_{N}_norm.npy` |
| `m_rcut_convergence.py` | `results/M/m_rcut_convergence.csv` |
| `mwr_locality_coarse_vs_dense.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, `results/M/M_ed_{N}.npy`, `mwr_locality.npz`, `wannier/{…}` |
| `_old_vs_new_7x7.py` | `results/M/{,obsolete_grid_7x7/}M_{,L_}dense_7x7.npy` |
| `_paths.py` | `data/graphene` (ou `$EDI_DATA`) |
| `rcut_resigma.py` | `results/M/resigma_9x9_rc0123.npz` |
| `resonance_criteria.py`, `resonance_metrics.py` | `results/M/resonance[_criteria]_{N}.npz` |
| `run_test_A_batch.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, `results/test_A` |
| `sampling_table.py` | `data/graphene/supercell/qe/defect_{N}_{p,d}.save/data-file-schema.xml`, `results/M/{ks_reconstruction.npz,sampling_table.csv,ved_analysis.npz}` |
| `submit_epw_p1_post.sh` | abs. dépôt, abs. `graphene/qe/epw/`, `results/epw/logs` |
| `submit_golden_dense.sh`, `submit_M_dense.sh`, `submit_M.sh`, `submit_lnl_frobenius.sh`, `submit_rcut_resigma.sh`, `submit_spectral*.sh`, `submit_tmatrix.sh`, `submit_test_A.sh`, `_test_mnl_mpi.sh` | abs. dépôt ; `results/M/logs/`, `results/M/M_{dense,ed,L,NL}_*`, `specwd_*`, `resigma_*`, `wannier/`, `data/graphene/**`, `results/test_A/`, `results/M/_test_mnl/` |
| `summarize_level1_maps.py` | `results/M/logs/specwd_{…}.out` |
| `tag_vacancy_sublattice.py` | `data/graphene/supercell/qe/defect_{N}_{p,d}.save`, `results/M/M_*` (sidecars) |
| `test_local_green_batch.py` | `wannier/NxN/wannier_tb.dat` |
| `test_local_tmatrix_real.py` | `data/graphene/unit_cell/qe/defect_{N}.save`, scratch dense `.save`, `results/M/M_{dense,ed}_{N}.npy`, `wannier/{…}` |
| `test_pad_vs_full_supercell.py`, `test_wannier.py` | `data/graphene/{supercell,unit_cell}/qe/…` |
| `validate_ML_grid_7x7.py` | `data/graphene`, `results/M/M_L_{dense_NxN_coarsecheck,NxN}.npy`, `_ML_NxN_ref_b4.npy` |

Aucun `Makefile` dans le dépôt ni dans `graphene/`. `jobs/` ne contient pas de scripts
SLURM (seulement `run.py`, `run_mpi.py`, `NL.py`, `Lr.py`, `LG_*.py` de l'ère ABINIT).
Dans `graphene/qe/epw/24k-24q/`, les runs `epw2_selfen_*` et `epw6_phself_*` pointent
vers `../epw1/{epwdata,vmedata,wigner}.fmt` par **liens relatifs** (35 liens chacun) :
déplacer `24k-24q/` en bloc est sûr, séparer `epw1/` des autres ne l'est pas.

---

## 5. Doublons et anomalies

### 5.1 Fichiers > 5 Go

| Taille | Fichier | mtime / atime |
|---|---|---|
| 6,25 G ×4 | `results/M/M_{,L_,NL_}dense_8x8.npy`, `M_dense_8x8_norm.npy` | 2026-09-04 / 09-05..09 |
| 6,20 G | scratch `qe_tmp/defect_12x12_p/defect_12x12_p.save/wfc1.hdf5` | 2026-06-18 / 06-18 |
| 6,18 G | scratch `qe_tmp/defect_12x12_d/defect_12x12_d.save/wfc1.hdf5` | 2026-06-18 / 06-18 |

Fichiers de 1 à 5 Go : 92 dans `graphene/abinit/pristine` (197,9 G), 40 dans `results/M`
(106,6 G), 16 `.epb` dans `epw/24k-24q/epw1` (30,4 G), 8 sur le scratch (21,6 G),
3 dans `graphene/abinit/defective` (6,5 G), 6 dans `ab-initio-defects/jobs` (11,2 G),
1 dans `abinit_assignment` (1,2 G).

### 5.2 Doublons volumineux (même nom + même taille, > 100 Mo) — 52 groupes, ≈ 212 G de copies excédentaires

| Groupe | Copies | Total | Où |
|---|---|---|---|
| `graphene_w90o_DS4_WFK.nc` (7 tailles : 3,18 / 2,79 / 2,39 / 2,00 / 1,60 / 1,21 / 0,82 G) | 14 / 9 / 20 / 10 / 13 / 8 / 15 | **≈ 180 G** | `graphene/abinit/pristine/unit_cell/wannier/bands/{5..20}x…x1/{5wann,8wann}/{16..64}bands/` — même WFK recopié pour chaque nombre de bandes |
| `graphene_w90o_DS3_WFK.nc`, `DS2_WFK.nc`, `graphene_w90_bandso_DS{1,2,3}_WFK` | 2 à 10 | ≈ 22 G | idem |
| `M_{,L_,NL_}dense_7x7.npy`, `M_dense_7x7_norm.npy` (3,66 G) | 2 | 14,6 G (7,3 G excédentaires) | `results/M/` vs `results/M/obsolete_grid_7x7/` — même nom/taille mais grilles différentes (comparés par `_old_vs_new_7x7.py`), pas forcément identiques |
| `defect_nro_{WFK,DEN,POT}.nc` | 2 | 5,3 G | `graphene/abinit/defective/10x10/defect_not_relaxed/test/` et `test/test2/` |
| `graphene_pristine_sco_WFK.nc` (0,98 G) | 2 | 2,0 G | `graphene/abinit/defective/8x8/pristine/` et `abinit/pristine/supercell/8x8/64G/` |
| `graphene_p_uc_5x5x1o_DS2_WFK.nc` | 2 | 1,5 G | `abinit/pristine/unit_cell/ground/5x5x1/calc{1,2}` |
| `graphene_vac1_a1o_WFK.nc`, `…relaxo_WFK` | 3, 2 | 0,9 G | `abinit/defective/5x5/` et `defective/vacancies/1/atom1/…` |
| `wfk_{p,d}[_sc].nc`, `wfk_uc.nc` (170,7 M, etc.) | 2–4 | ≈ 3 G | `ab-initio-defects/jobs/{6x6,8x8,16x16}` et `jobs/10x10x1_5x5_Sc/*/…` |
| `charge-density.hdf5` (0,10–0,14 G) | 2 | 0,7 G | scratch `qe_conv/ecut_*` vs `qe_conv/kpoint/*` (même maille) |

Sans même nom mais redondants : scratch `prefix.wfcN` en vrac vs `prefix.save/wfcN.hdf5`
(≈ 49 G, §2) ; `results/M/*_norm.npy` dérivables (29 G, §3b).

### 5.3 Répertoires vides (42)

- `graphene/qe/epw/{24k-24q,8k-8q,36k-30q,24k-12q,30k-24q,16k-16q,36k-36q,16k-12q,12k-12q,30k-30q}/epw3` (10)
- `graphene/qe/epw/16k-8q/phonon_line_coupling_strength`
- `dft/compmatphys/basic/{NaCl, AlFe/relaxation, AlFe/atompos, Al/band, Al/kpt}`
- `abinit_processing/cif`, `abinit_processing/venv/abienv/include/python3.11`, `ab-initio-defects/.venv/include/python3.11`
- `codes/qe/external/{devxlib,qe-gipaw,fox,pw2qmcpack,wannier90,lapack,d3q,mbd}`, `codes/wannier90/build`
- `codes/abinit/tests/modules_with_data/{MgO_eph_zpr,LiF_eph_varpeq,MgB2_eph4isotc,diamond_eph_gwpt}`
- `.git/refs/tags`, `.git/objects/info` (×4 dépôts)
- scratch : `compmatphys/basic/**` (9 sous-dirs vides), `jobs/`, `graphene/`

### 5.4 Antérieur à juin 2026 et non accédé depuis (mtime et atime < 2026-06-01)

Projet : **1 566 G, 83 701 fichiers**.

| Répertoire | Taille | Fichiers | mtime max | atime max |
|---|---|---|---|---|
| `graphene/qe/epw/` (grilles test, hors 24k-24q) | 1 110,1 G | 45 222 | 2026-04-30 | 2026-05-05 |
| `graphene/abinit/pristine/` | 403,1 G | 8 807 | 2026-02-02 | 2026-05-27 |
| `graphene/qe/test/` | 31,5 G | 11 882 | 2026-04-21 | 2026-04-22 |
| `graphene/abinit/defective/` | 17,6 G | 436 | 2025-12-21 | 2026-02-03 |
| `abinit_assignment/` | 2,4 G | 208 | 2026-05-28 | 2026-05-29 |
| `dft/{copper,diamond,compmatphys}` | 1,2 G | 695 | 2025-11-24 | 2025-11-26 |
| `ab-initio-defects/.venv/lib` (partiel) | 0,25 G | 14 709 | 2025-11-26 | 2026-01-23 |

(`ab-initio-defects/jobs/` date de déc. 2025 – févr. 2026 mais a été accédé après juin ;
il figure en 3c.)

Scratch : **48,9 G, 10 442 fichiers** (`qe_conv/` 44,9 G + `graphene_scf` 4,0 G),
plus les 43,2 G de supercellules dont les wfc n'ont pas été lus depuis le 2026-06-18 (§0).

### 5.5 Autres anomalies

- 1 lien symbolique cassé : `codes/wannier90/test-suite/library-mode-test/ref/gaas.win → ../gaas.win`.
- Dépôt git : 6 entrées non suivies (`scripts/_test_mnl_mpi.sh`, `wannier/24x24/`, `wannier/{25,27,28,32}x…_nb16/`).
- `graphene/qe/epw/24k-24q/epw1/.epw1.in.swp` (vim, avril) et `epw1_oom_20637122.out` (run OOM).
- `abinit_assignment/__ABI_MPIABORTFILE__.lock`, `fort.98`.
- `ab-initio-defects/hamiltonian_reconstruction.png` à la racine du dépôt (août), hors `figures/`.
- Le `_ph0` de `graphene/qe/epw/24k-24q/phonons` a été supprimé le 2026-09-10 : `phonons/save/`
  (dvscf + dyn, 0,97 G) est désormais la seule source DFPT de la grille de production.

---

## Fichiers de travail de la session (scratchpad, hors projet, effacés avec la session)

`du_*.txt`, `project_files.tsv` (115 601 lignes : taille, mtime, atime, chemin),
`scratch_files.tsv` (15 775 lignes), `project_symlinks.txt`, `project_emptydirs.txt`,
`dups.txt`, `grep_{scripts,epw,src}.txt`, `paths_by_script.tsv`
dans `/tmp/claude-3139538/-home-gregb26/04514fb5-014b-45e7-93dd-ade179002f69/scratchpad/`.
