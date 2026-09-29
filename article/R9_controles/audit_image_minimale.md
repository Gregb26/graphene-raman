# Audit — image minimale dans une cellule à 60° (écart 2 de R9) et contour de zone d'analyze_M.py (écart 3)

- Prompt : « Audit — image minimale dans une cellule à 60° (écart 2 de R9) et contour de zone de analyze_M.py (écart 3) » (Greg, 2026-09-28),
  à la suite du STOP de la clôture de R9.
- Statut : **TEST**, lecture seule. Aucune campagne, aucun job. Rien de modifié dans `src/`, `scripts/`, `config/`, `results/`, ni dans les
  pilotes et rapports existants ; seule la sous-commande `audit` a été ajoutée à `r9_driver.py`. Git en lecture (HEAD 1883861).
- Exécution : `python r9_driver.py audit` sur le nœud de connexion, géométrie seule. Sortie : `audit/audit_results.json`.
  L'option `--parts` permet de rejouer une partie seule ; le résultat est alors fusionné dans le json.
- Chiffres bruts. Les corrections sont proposées au §7 ; aucune n'a été appliquée.

## 0. Méthode

- **Réduction axe par axe.** C'est celle du code audité, recopiée à l'identique dans le pilote (la routine de production n'est pas appelée) :
  chaque composante réduite est ramenée dans [−0,5 ; 0,5) par `d → (d + 0,5) mod 1 − 0,5`, puis on prend la norme cartésienne.
- **Vraie image minimale.** On part de la réduction axe par axe, puis on garde la plus courte des 9 images `d + s`, avec
  s ∈ {−1, 0, 1}² dans le plan. L'axe a₃ est perpendiculaire au plan ; c'est vérifié sur chaque XML.
  - Pourquoi 9 images suffisent (a₁, a₂ à 60°) : la cellule de Wigner-Seitz a ses coordonnées réduites dans [−2/3 ; 2/3],
    donc `s = v* − d` tombe dans (−7/6 ; 7/6].
  - Les égalités entre images sont comptées à 1e-6 près.
- **Deux propriétés utiles pour lire les tableaux.** On note L la longueur d'un vecteur de la cellule, en Å, ou en mailles pour les réseaux R.
  1. Tout vecteur qui sort de la réduction axe par axe avec une longueur < L/2 est la vraie image. Raison : toute translation non nulle
     mesure au moins L.
  2. Toute distance vraie < (√3/4) L = 0,433 L est rendue exactement. C'est le rayon inscrit du parallélogramme [−½ ; ½)².
  - Au-delà, la réduction axe par axe peut rendre une image plus longue que la vraie : jusqu'à (√3/2) L pour un point à L/2 en
    distance vraie. La plus grande distance vraie dans la cellule est L/√3.
- **Géométries lues :**
  - super-cellules non relaxées 5…12, 15, 18, 21, 24, 27 : XML des `.save` de `qe_tmp`, les mêmes que R4, R5, R7 et R9 ;
  - 9×9 relaxée de R1 : nspin1 et nspin2 ;
  - grilles FFT : XML de chaque `.save` ;
  - étiquettes R de Wannier : boîte MP D×D de `config/production.json`, et N×N pour les grilles grossières.
- **Portes de régression** (toutes passées) :
  - Nombre de points des sphères de Lu (0,5 et 1,0 Å) autour de l'atome axe par axe : égal à celui de R7 D1
    (6, 9, 12, 15, 18, 21, 24, 27 ; par ex. 1 039 et 8 611 points ; 1 047 et 8 496 pour la 21×21).
  - Indice de l'atome de Lu : égal à celui enregistré par R4 D5 (5…12), R5 C (5…12) et R7 D1 (6…27).
  - Distances vraies d'A.1 (`a1_profiles_*.npz`) : écart 0,0 Å.
  - Valeur de Lu publiée (1,0 Å) : redonnée à 0,005 meV près.
  - Profils radiaux d'`analyze_Ved.py` : redonnés à 0,0 eV près, à partir des plans de `results/M2/ved_analysis.npz`.
  - Abscisses et poids de `results/M2/mwr_locality.npz` : redonnés, 0,0 sur les distances, 5e-18 relatif sur les poids.
  - Moyennes des disques de la carte D.1 de R9 : redonnées à 8,5e-14 eV Å² avec la convention d'étiquettes du pilote.

## 1. Inventaire

La colonne « Méthode » vaut :

- **axe** : réduction axe par axe ;
- **vraie** : recherche sur les images ;
- **arith.** : réduction d'indices entiers ou test d'égalité modulo 1, sans distance.

### 1.1 `src/electron_defect_interaction/`

| Fichier:ligne | Fonction | Méthode | Appelants (sites mesurés au §3) |
|---|---|---|---|
| `defects/alignment.py:22` | `min_image_dist` | axe, cartésienne | `vacancy_site`, `far_atom`, `far_atom_alignment` l. 76, 79, `atom_sphere_shifts` (`dist_axis`) ; `r4_driver.py:417` ; `r7_driver.py:99` |
| `defects/alignment.py:28` | `vacancy_site` | axe (min sur les atomes de d) | `far_atom_alignment`, `atom_sphere_shifts` ; `gate_M_normalization.py:157` ; `r4_driver.py:425, 608` ; `r5_driver.py:237, 633` ; `r7_driver.py:98` |
| `defects/alignment.py:35` | `far_atom` | axe (argmax) | `far_atom_alignment` l. 74 ; `atom_sphere_shifts` l. 177 (`i_far_axis`) |
| `defects/alignment.py:42` | `sphere_average` | axe, grille 3D | `far_atom_alignment` l. 81–82 |
| `defects/alignment.py:65` | `far_atom_alignment` (Lu) | axe | `r4_driver.py:454` (D1 : D, N1, N2), `:571` (D5) ; `r5_driver.py:672` (C) ; `r7_driver.py:120` (D1) ; `r9_driver.py:751` (A.1) |
| `defects/alignment.py:104` | `true_min_image_dist` (P1) | vraie | `atom_sphere_shifts` ; `r9_driver.py:775–776` |
| `defects/alignment.py:116` | `_sphere_points` (P1) | axe (critère de `sphere_average`) | `atom_sphere_shifts` |
| `defects/alignment.py:144` | `atom_sphere_shifts` (P1) | vraie pour `dist`, axe pour `dist_axis` et les sphères | `r9_driver.py:750` (A.1) |
| `io/qe_gamma_io.py:163` | `inplane_disc_mask` | axe, grille du plan | `r4_driver.py:462, 663` ; `r5_driver.py:254` et C (w₂ 2 Å, w₁ 1 Å) ; `r7_driver.py` (D1) |
| `defects/many_body/local_tmatrix.py:95` | `recenter_mwr` : `((R − R_d + D//2) % D) − D//2` | axe, entiers | toute la chaîne V_loc (`compute_spectral_wannier`, `rcut_resigma`, `resonance_metrics`, `resonance_criteria`, `m_rcut_convergence`, `mwr_locality_coarse_vs_dense`, R4…R9) |
| `defects/many_body/local_tmatrix.py:48` | `mwr_locality` : `|R − R0|` | norme réduite des étiquettes axe par axe | journal de `compute_spectral_wannier.py:84` |
| `wannier/wannier_interpolation.py:76–79` | `Mwk_to_Mwr`, étiquettes `arange(N) − N//2` | axe (boîte centrée sur l'origine) | toutes les étiquettes R ; interpolation hors grille par `Mwr_to_Mwk` (l. 91) et `Mwr_to_Mwk_pairs` (P2, l. 115) |
| `wannier/wannier_interpolation.py:166` | `_match_kpoint_order` | arith. (égalité ≤ 1e-5) | lecture des U |
| `wannier/wannier_interpolation.py:188` | `_infer_mp_grid` | arith. | — |
| `wannier/supercell_fold.py:30, 85, 99, 136` | `cell_of`, repliements | arith. (indices mod N) | R9 C |
| `utils/lattice.py:63` ; `io/qe_io.py:309` ; `wavefunctions/sc_projection.py:34` ; `utils/fft_utils.py:137` ; `utils/interpolation.py:26, 72` ; `defects/local_R.py:61, 133` ; `defects/local_G.py:56–445` ; `electron_photon/tb_model.py:79` ; `electron_phonon/selfen.py:32–43` ; `electron_phonon/phself.py:75` | réductions de k, d'indices de grille ou de Miller, multiplicités | arith. | — |

### 1.2 `scripts/`

| Fichier:ligne | Appel | Méthode | Ce qui est sélectionné ou calculé |
|---|---|---|---|
| `analyze_Ved.py:22, 63` ; `analyze_M.py:89–90` ; `tag_vacancy_sublattice.py:10` | lacune = argmax des min_d des normes réduites | axe, norme réduite | site de la lacune |
| `analyze_Ved.py:28–30, 40–42` | profil radial (anneaux de 0,05 Å jusqu'à \|a₁\|/2), non masqué | axe | `{S}_rad` → `fig_Ved_radial` (`make_figures.py:214–223`) |
| `analyze_Ved.py:79–86` | profil radial masqué (r ≥ 0,5 Å de chaque atome de d) | axe | `{S}_rad_masked` → `fig_Ved` (c) (`make_figures_memoire.py:127–133`), `fig_Ved_radial_masked` (`make_figures.py:226–247`) ; journal (valeur à 1, 2, 3 Å) |
| `analyze_Ved.py:69, 72–73` | voisins à 0,5–1,6 Å, déplacements | axe | journal d'anomalie |
| `analyze_Ved.py:32–35, 44` ; `analyze_M.py:96` ; `r4_driver.py:550–551` | plans frontière `(s_vac + ½)·n`, recentrage de la carte | aucune distance (plan défini par une coordonnée) | `ved_*_b1/b2`, `fig_Ved_boundary`, cartes |
| `analyze_Ved.py:90` | positions des atomes relatives à la lacune | axe | `{S}_atoms_xy` (tracé) |
| `analyze_M.py:37` (`kdist`), `:45, 108` ; `resonance_metrics.py:97` ; `check_onsite_and_NL.py:30` ; `gate_M_normalization.py:49, 168` ; `epw_d2_extract.py:66` | point de grille le plus proche de K, K′, Γ | axe, norme **réduite** | indice iK |
| `analyze_M.py:53–56` | `kfold` (9 images de G) | vraie | carte `fig_M_map*` |
| `analyze_M.py:56–57, 60` | `map_corners` | — | §5 (écart 3) |
| `mwr_locality_coarse_vs_dense.py:17–20` | abscisse \|R\| cartésienne des étiquettes recentrées | axe (via `recenter_mwr`) | `mwr_locality.npz` → `fig_locality_final` (`make_figures_memoire.py:61–72`), `fig_locality` (`make_figures.py:91–103`) |
| `compute_spectral_wannier.py:102` ; `rcut_resigma.py:46` ; `resonance_metrics.py:41` ; `resonance_criteria.py:34` ; `test_local_green_batch.py:10` | amas R_cut : `|Rn|_réduite ≤ R_cut` | axe (étiquettes) + norme réduite | V_loc de production |
| `m_rcut_convergence.py:30–44` | idem, et référence « toutes les R » sur la grille fine 60² | axe (étiquettes) | tab:rcut_M |
| `make_figures.py:131–138` ; `make_figures_memoire.py:92–96` | `bz_vertices` (Voronoi de `map_B`) | vraie | contour de `fig_M_map`, `fig_M_map_final` |
| `_bands.py:14–17, 39, 50` | étiquettes des sommets de chemins (`SYM_POINTS`, tolérance 2e-2) | arith. (égalité modulo 1) | étiquettes des figures de bandes |
| `finalize_wannier.py:64–65` ; `test_wannier.py:83` ; `check_M_dense_vs_coarse*.py` ; `m_rcut_convergence.py:22` ; `analyze_M.py:124, 158` ; `ks_reconstruction_all.py:34` ; `test_ks_reconstruction.py:57` | appariement de k, arrondis MP, multiplicités | arith. | — |

### 1.3 Pilotes R4…R9

| Fichier:ligne | Appel | Méthode | Ce qui est sélectionné ou calculé |
|---|---|---|---|
| `R4_quasi_lie/r4_driver.py:415–418` | premiers (< 1,8 Å) et seconds (2,2–2,7 Å) voisins | axe | indices D1 |
| `r4_driver.py:454` (D1 : D, N1, N2 relaxées), `:571` (D5 : 5…12) | Lu | axe | décalage appliqué aux états D1 ; table D5 |
| `r4_driver.py:462, 663` | disques w₂ (2 Å), w₁ (1 Å) | axe | seuil de localisation, poids |
| `r4_driver.py:266, 672` | K sur les grilles 27² et 9² | axe, norme réduite | iK |
| `R5_base_vs_M/r5_driver.py:108, 647` | K sur la grille N×N de la maille | axe, norme réduite | utilisé seulement si \|dk\| < 1e-4 |
| `r5_driver.py:672` ; `:254` et C | Lu (5…12) ; disques w₂ / w₁ | axe | niveaux π, σ alignés (R5 C) |
| `r5_driver.py:243, 674` ; `r7_driver.py:124` | ΔV au point de grille du site | arith. | valeur au site |
| `R7_tailles_3m/r7_driver.py:98–99, 120` | lacune, premiers voisins, Lu | axe | niveaux π, σ (R7 D1, 6…27), ajustements |
| `R7_tailles_3m/make_inputs_r7_relax.py:55–58` (`minimg_cart`), `:62–76`, `:105`, `:108–138` | champ u(r) de la 12×12 R2 transplanté pour \|r\| < R_CUT = 14,5 Å ; `rcut_max = L/2` (source relaxée) | axe (`d − rint(d)`) ; contrôle d'unicité l. 72 : 9 images | géométries initiales R7c 15…27 |
| `R7_tailles_3m/phase0_inventory.py:57, 67–68` | atome manquant, sous-réseau | axe / arith. (distance 0) | inventaire |
| `R9_controles/r9_driver.py:154` (`box_geometry`) | boîte N×N dans la boîte MP (c = R mod D, bord en mailles) | arith. (parallélogramme par définition) | A.2, A.3 |
| `r9_driver.py:771–776` (`orbital_sites`), `:757` (plateau) | distances vraies | vraie (P1) | plateau (i), (ii) d'A.1 |
| `r9_driver.py:750–756` | Lu (A.1) et C_N « Lu » | axe | colonne Lu ; variantes « aligné Lu » (A.3, B, C, D) |
| `r9_driver.py:534` (`crown_table`) | couronnes x autour de K, K′ (240²) | axe | E_res par couronne |
| `r9_driver.py:373, 1436` | iK (240² ; 27²) | axe, norme réduite | — |
| `r9_driver.py:1487–1493` (D.1) | disques autour de K, K′ : vraie (9 images) ; interpolation `Mwr_to_Mwk_pairs` avec les étiquettes brutes de `Mwk_to_Mwr` | vraie / axe | carte de Kaasbjerg |
| `r9_driver.py:1128–1136` (C), `NN_CELLS` l. 304 | voisins de la lacune en mailles fixes, repliement mod N | arith. | C |

Hors périmètre, cités pour mémoire :

- `super_cell_relaxed/9x9/analyze_relax.py:54` et `super_cell_relaxed/series/make_inputs_series.py:26` réduisent axe par axe,
  mais pour des distances entre voisins (< 3 Å) seulement ;
- le carnet `notebooks/defect_matrix.ipynb` n'a pas été lu.

## 2. Distances en jeu, par appel

| Appel | Distances en jeu | Seuil 0,433 L |
|---|---|---|
| `vacancy_site` (et variantes à norme réduite) | min sur les atomes de d : 0 (non relaxé), ≤ 0,25 Å (R1), 1,26–1,42 Å à la lacune | L ≥ 12,3 Å → seuil ≥ 5,3 Å : petites |
| voisins (< 1,8 ; 2,2–2,7 ; 0,5–1,6 Å) ; sphères de Lu (0,5 ; 1,0 Å) ; disques w₁, w₂ (1 ; 2 Å) ; masque de cœur (0,5 Å) | ≤ 2,7 Å | petites |
| partenaire de la parfaite (`far_atom_alignment` l. 76, 79) | 0 (non relaxé) ; 0,001–0,018 Å (R1) | petites |
| `far_atom` (atome de Lu) | argmax sur toute la cellule : jusqu'à L/√3 = 0,577 L en vrai, jusqu'à 0,87 L axe par axe | **grandes** |
| profil radial d'`analyze_Ved` | jusqu'à L/2 | **grandes** au-delà de 0,433 L |
| étiquettes R : amas R_cut ≤ 4 | \|Rn\| ≤ 4 mailles ; D ≥ 24 (dense), N ≥ 5 (grossier) | petites pour le dense ; grossier 5×5 : l'amas R_cut 3 couvre les 25 mailles |
| étiquettes R : abscisse de `fig_locality*`, interpolation hors grille (D.1, référence de tab:rcut_M) | toute la boîte D×D | **grandes** |
| k le plus proche de K | ≤ 1/(√3 n) \|b\| ; grilles avec K hors grille : 5, 7, 8, 10, 11 (grossières), 25, 28, 32 (denses) | petites (égalités de C₃) |
| couronnes 240² (x ≤ 400) | ≤ 20 pas de grille | petites |
| transplant R7c | \|r\| < 14,5 Å dans la 12×12 (L₁₂ = 29,59 Å, 0,433 L₁₂ = 12,81 Å) | **grandes** entre 12,81 et 14,5 Å |

## 3. Mesures sur les géométries réelles

### 3.1 Lacune

Pour les 13 tailles, la lacune (indice, min_d) est la même avec les trois méthodes : axe par axe cartésienne, vraie image, norme réduite
(celle d'`analyze_M`). min_d vaut 1,4237 Å et le suivant 0,0 Å.

Pour R1 nspin1 et nspin2, c'est aussi le même indice, 81. min_d vaut 1,3352 et 1,2627 Å ; le suivant vaut 0,1461 et 0,2481 Å.

### 3.2 Atome de Lu (`far_atom`)

Colonnes :

- **atome axe** : atome retenu par `far_atom`, avec sa distance axe par axe, sa distance vraie et son rang en distance vraie ;
- **vrais plus loin** : atomes à r_max (égalités à 1e-6 près) ;
- **r_max vrai** : vaut a_sc/√3 pour toutes les tailles ;
- **n (axe ≠ vraie)** : nombre d'atomes dont la distance axe par axe diffère de la vraie ;
- **d vraie min où ≠** : plus petite distance vraie parmi ces atomes.

Le partenaire de la parfaite est le même avec les deux méthodes. Son déplacement vaut 0,0 Å pour les géométries non relaxées, 0,001 Å
(atome axe par axe) et 0,012 à 0,018 Å (atome 140) pour R1.

| N | a_sc (Å) | atome axe | d axe (Å) | d vraie (Å) | rang | vrais plus loin | r_max vrai (Å) | n (axe ≠ vraie) | d vraie min où ≠ (Å) |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 12,33 | 49 | 9,97 | 6,21 | 8 | 2 | 7,12 | 5 / 49 | 5,70 |
| 6 | 14,80 | 71 | 12,81 | 7,40 | 9 | 2, 57 | 8,54 | 9 / 71 | 6,52 |
| 7 | 17,26 | 97 | 14,24 | 8,66 | 11 | 81 | 9,97 | 12 / 97 | 7,93 |
| 8 | 19,73 | 1 | 17,08 | 9,86 | 17 | 20 | 11,39 | 17 / 127 | 8,54 |
| 9 | 22,19 | 161 | 18,51 | 11,12 | 18 | 21, 140 | 12,81 | 21 / 161 | 9,97 |
| 10 | 24,66 | 198 | 19,93 | 12,41 | 20 | 23 | 14,24 | 25 / 199 | 10,75 |
| 11 | 27,12 | 241 | 22,78 | 13,58 | 23 | 26 | 15,66 | 33 / 241 | 12,16 |
| 12 | 29,59 | 287 | 25,63 | 14,80 | 36 | 28, 235 | 17,08 | 43 / 287 | 12,81 |
| 15 | 36,99 | 449 | 31,32 | 18,51 | 48 | 65, 384 | 21,35 | 65 / 449 | 16,42 |
| 18 | 44,39 | 1 | 38,44 | 22,19 | 63 | 115, 570 | 25,63 | 97 / 647 | 19,26 |
| 21 | 51,78 | 881 | 44,13 | 25,90 | 87 | 133, 748 | 29,90 | 133 / 881 | 22,78 |
| 24 | 59,18 | 2 | 49,83 | 29,62 | 114 | 201, 1000 | 34,17 | 178 / 1151 | 25,63 |
| 27 | 66,58 | 1457 | 56,95 | 33,30 | 141 | 225, 1232 | 38,44 | 225 / 1457 | 29,21 |
| 9 R1 nspin1 | 22,19 | 161 | 18,51 | 11,12 | 18 | 140 | 12,81 | 33 / 161 | 9,96 |
| 9 R1 nspin2 | 22,19 | 161 | 18,51 | 11,12 | 18 | 140 | 12,80 | 33 / 161 | 9,95 |

- L'atome axe par axe est à une distance vraie comprise entre 0,866 et 0,872 r_max.
- Il n'appartient à l'ensemble des vrais plus loin pour aucune taille.

### 3.3 Voisins, sphères, disques

- **Voisins.** Pour les 13 tailles et R1, les ensembles < 1,8 Å (3 atomes) et 2,2–2,7 Å (6 atomes) de la cellule d, et 0,5–1,6 Å de la
  parfaite (3 atomes), sont identiques axe par axe et en vraie image.
- **Sphères de Lu (0,5 et 1,0 Å).** Autour de l'atome axe par axe et de chaque atome vraiment le plus loin, pour les 13 tailles et R1, les
  ensembles de points sont identiques. Par exemple, autour de l'atome axe par axe : 1 039 et 8 611 points pour 5, 6, 8, 9, 10, 12, 15, 18, 24, 27 ; 1 129 et 9 154 pour 7 ;
  1 287 et 10 311 pour 11 ; 1 047 et 8 496 pour 21 ; 1 045 et 8 599 pour R1 nspin1.
- **Toutes les sphères de P1 (A.1).** Pour les 6 tailles de R9, il y a 1 584 sphères (tous les atomes, 0,5 et 1,0 Å). Aucune ne diffère
  entre axe par axe et vraie image, et aucun point n'est à moins de 1e-9 Å du rayon.
- **Disques w₁ (1 Å) et w₂ (2 Å).** Autour de la lacune, pour les 13 tailles, ils sont identiques : 547 et 2 161 points (567 et 2 271 pour
  7×7 ; 638 et 2 560 pour 11×11 ; 532 et 2 116 pour 21×21).

### 3.4 Profil radial d'`analyze_Ved.py`

Le plan z = z_C est reconstruit depuis `ved_analysis.npz`. La porte est passée : le profil axe par axe redonne `{S}_rad` et
`{S}_rad_masked` à 0,0 eV.

Les anneaux dont le nombre de points change commencent à celui qui contient (√3/4) a_sc. Dans ces anneaux, la réduction axe par axe exclut une partie
des points à distance vraie < a_sc/2 : par exemple, 490 points au lieu de 534 dans l'anneau à 9,625 Å de la 9×9.

| N | anneaux changés / total | premier anneau changé (Å) | 0,433 a_sc (Å) | max\|Δ\| non masqué (meV) à r (Å) | max\|Δ\| masqué (meV) à r (Å) |
|---|---|---|---|---|---|
| 5 | 18 / 124 (masqué 17) | 5,325 | 5,34 | 16,5 à 5,775 | 19,9 à 5,775 |
| 6 | 20 / 148 | 6,425 | 6,41 | 14,2 à 7,175 | 10,9 à 6,925 |
| 7 | 24 / 173 | 7,475 | 7,47 | 6,7 à 7,575 | 10,2 à 7,775 |
| 8 | 28 / 198 (masqué 27) | 8,525 | 8,54 | 1,8 à 9,525 | 1,9 à 8,775 |
| 9 | 30 / 222 | 9,625 | 9,61 | 1,4 à 9,675 | 1,5 à 10,075 |
| 10 | 34 / 247 | 10,675 | 10,68 | 0,9 à 10,825 | 1,0 à 10,675 |
| 11 | 37 / 272 | 11,775 | 11,75 | 1,2 à 11,925 | 1,7 à 11,925 |
| 12 | 40 / 296 | 12,825 | 12,81 | 2,2 à 14,325 | 1,8 à 14,325 |

Le masque de cœur (0,5 Å autour de chaque atome) est identique axe par axe et en vraie image pour les 8 tailles. Le détail anneau par anneau
(nombre de points, valeurs avant et après) est dans `audit_results.json` → `sizes.<S>.radial_analyze_Ved`.

### 3.5 Étiquettes R de Wannier

La boîte MP D×D est réduite axe par axe par `recenter_mwr`, ce qui donne les étiquettes centrées `arange(D) − D//2`. Distances en unités
de a.

| Taille | D (dense) | étiquettes à distance changée | x axe min où ≠ | x axe max → vrai max | décalage max | N (grossier) | changées |
|---|---|---|---|---|---|---|---|
| 5×5 | 25 | 96 / 625 | 13,11 | 20,78 → 13,89 | 8,25 | 5 | 2 / 25 |
| 6×6 | 24 | 89 / 576 | 12,53 | 20,78 → 13,86 | 8,78 | 6 | 5 / 36 |
| 7×7 | 28 | 123 / 784 | 14,53 | 24,25 → 15,62 | 10,25 | 7 | 6 / 49 |
| 8×8 | 32 | 161 / 1024 | 16,52 | 27,71 → 18,19 | 11,71 | 8 | 9 / 64 |
| 9×9 | 27 | 112 / 729 | 14,11 | 22,52 → 15,59 | 8,99 | 9 | 10 / 81 |
| 12×12 | 24 | 89 / 576 | 12,53 | 20,78 → 13,86 | 8,78 | 12 | 21 / 144 |

**Amas R_cut (norme réduite ≤ R_cut, R_cut 0…4, toutes les grilles denses et grossières).** Les ensembles sont identiques entre étiquettes
axe par axe et vraie image (en norme réduite). Il y a 1, 5, 13, 29 et 49 mailles, sauf sur les grilles grossières où la boîte est plus
petite que le disque :

- 5×5 : 25 et 25 mailles pour R_cut 3 et 4 ;
- 6×6 : 27 et 35 ;
- 7×7 : 45 pour R_cut 4 ;
- 8×8 : 47 pour R_cut 4.

**Abscisses de `fig_locality_final` et `fig_locality`.** Porte : distances à 0,0 près, poids à 5e-18 près (relatif).

- Grille élargie : tous les points dont l'abscisse change sont dans la fenêtre tracée (2e-5 ≤ w ≤ 60 eV).

  | Taille | points déplacés | abscisse axe par axe | w max parmi eux (eV) |
  |---|---|---|---|
  | 5×5 | 96 | ≥ 13,11 a | 3,5e-3 |
  | 6×6 | 89 | ≥ 12,53 a | 4,9e-3 |
  | 7×7 | 123 | ≥ 14,53 a | 2,3e-3 |
  | 8×8 | 161 | ≥ 16,52 a | 2,1e-3 |
  | 9×9 | 112 | ≥ 14,11 a | 2,2e-3 |
  | 12×12 | 89 | ≥ 12,53 a | 4,8e-3 |

- Grille grossière (croix) : tous les points déplacés sont tracés.

  | Taille | déplacements |
  |---|---|
  | 5×5 | 2 points, 3,46 → 2,65 a (w 2,18 et 1,97 eV) |
  | 7×7 | 4 points 4,36 → 3,46 a et 2 points 5,20 → 3,61 a (w 1,11–1,57 eV) |
  | 8×8 | 9 points : 4,58 → 3,61 ; 5,20 → 4,36 ; 5,29 → 3,46 ; 6,08 → 3,61 ; 6,93 → 4,00 a (w 0,47–1,29 eV) |

### 3.6 Interpolation hors grille : carte D.1 de R9 (9×9, D = 27, R_d = (4, 4))

- **Pilote D.1.** Il interpole avec les étiquettes brutes de `Mwk_to_Mwr` (boîte centrée sur l'origine, défaut en R_d).
- **Méthode de l'audit.** La carte est refaite dans le pilote avec les mêmes sommes de Fourier, recopiées, pour trois conventions
  d'étiquettes (bra et ket) :
  - **brut** : convention du pilote ;
  - **axe_Rd** : réduction axe par axe autour de R_d, c'est-à-dire les étiquettes de `recenter_mwr` ;
  - **vraie** : vraie image autour de R_d, avec un poids 1/n sur les n images à égalité, à la manière des dégénérescences de Wigner-Seitz.
- **Étiquettes changées.** 226 étiquettes changent de brut à vraie, dont 26 à égalité ; 130 de axe_Rd à vraie. Le plus grand bloc
  ‖M_W(R, R′)‖ concerné vaut 2,8e-3 eV, contre 41,88 eV pour le bloc sur site.

| Passage | variante | max\|ΔṼ\| sur la carte (eV Å²) | max\|ΔṼ\| dans les disques K, K′ (eV Å²) | max\|Δ moyenne de disque\| (eV Å²) |
|---|---|---|---|---|
| brut → vraie | brut | 0,443 | 0,094 | 0,009 |
| brut → vraie | exact | 0,443 | 0,094 | 0,009 |
| axe_Rd → vraie | brut | 0,366 | 0,070 | 0,026 |
| brut → axe_Rd | brut | 0,360 | 0,087 | 0,035 |

Moyennes des disques (valence | K, conduction | K, valence | K′, conduction | K′), en eV Å² :

| Convention | brut | exact |
|---|---|---|
| brut (R9) | 76,955 / 78,159 / 81,176 / 82,297 | 81,512 / 82,706 / 81,168 / 82,289 |
| vraie | 76,956 / 78,168 / 81,177 / 82,305 | 81,513 / 82,715 / 81,168 / 82,297 |

### 3.7 Grilles k et couronnes

- **K sur la grille** (grossières 6, 9, 12 ; denses 24, 27 ; sortie 240 ; internes 300, 450, 600, 900 ; q d'EPW 24) : le point choisi est
  exact.
- **K hors de la grille** (grossières 5, 7, 8, 10, 11 ; denses 25, 28, 32) :
  - l'argmin de la norme réduite axe par axe est unique ;
  - il appartient à l'ensemble des trois points équidistants en distance cartésienne vraie (images C₃ autour de K). Par exemple, en 25² :
    (17, 8), dans l'ensemble {(16, 8), (17, 8), (17, 9)}, à 0,0680 Å⁻¹ ;
  - même chose pour K′ ; Γ est sur toutes les grilles.
- **Couronnes 240² de `crown_table`** (x ≤ 400, 2 914 points k) : aucune différence entre axe par axe et vraie image.

### 3.8 Transplant R7c (`make_inputs_r7_relax.py`)

- **Côté source (12×12 de R2, lacune B).**
  - 225 sites ont |r| < 14,5 Å axe par axe, soit le « 225 atomes déplacés » de `transplant.log`.
  - 249 sites ont une distance vraie < 14,5 Å.
  - 24 sites manquent (0 en trop). Leur distance vraie va de 12,81 à 14,45 Å, leur distance axe par axe de 16,17 à 23,52 Å, et leur
    déplacement relaxé |u| de 0,0042 à 0,0128 Å.
  - Pour comparaison, le plus grand |u| parmi les sites retenus entre 12 et 14,5 Å vaut 0,0141 Å.
- **Côté cible (15…27).** L'ensemble à < 14,5 Å de la lacune est identique axe par axe et en vraie image (249 sites).
- **Les cinq `transplant.log`.** Ils donnent tous « source 12×12, R_cut 14.5 A » ; `rcut_max = L/2 − 0,3` (source relaxée N) n'a pas été
  employé.

## 4. Classement

### (a) Sans effet : résultat identique

- `vacancy_site` et les trois variantes à norme réduite (13 tailles, R1).
- Voisins de R4 D1 et R7 D1, voisins d'`analyze_Ved`.
- Partenaire de la parfaite et `far_displacement`.
- Sphères de Lu (0,5 et 1,0 Å), sphères de P1 (1 584), `sphere_average` et `_sphere_points`.
- Disques w₂ et w₁ d'`inplane_disc_mask` (R4, R5, R7).
- Masque de cœur d'`analyze_Ved`.
- Amas R_cut et `recenter_mwr` pour V_loc : tout le V_loc de production, T, Γ, niveau 1 et 2, R9 B et C.
- Point le plus proche de K. `analyze_M` §3 (`lnl_*`, repris par `level2_families.py`) et `gate_M_normalization` prennent un des trois
  points équidistants ; R5 C ne l'emploie que si |dk| < 1e-4.
- Couronnes de `crown_table`.
- `kfold` et `bz_vertices` (vraie image).
- Journal de `compute_spectral_wannier.py:84` (`mwr_locality`, norme réduite des étiquettes axe par axe) : seules les douze premières
  distances sont écrites, toutes dans l'amas, et aucune ne change.
- `box_geometry` et `orbital_sites` de R9 ; plateau (i) d'A.1 (distances vraies de P1).
- Toutes les réductions arithmétiques du §1.
- Plans frontière : pas de distance.

### (b) Effet sur une distance seulement affichée

1. **Distance de l'atome de Lu** (`dist_far`, colonne « atome / distance »). Elle apparaît dans :
   - la table D5 de R4 (`r4_driver.py:1013`, `d5/d5_results.json`) ;
   - le journal de R4 D1 ;
   - `alignment.dist_far_A` de R5 C (`c_results.json`) et de R7 D1 (`d1_results.json`, D1_tables) ;
   - R9 A.1, qui affiche déjà les deux distances.

   Valeurs affichées → vraies : 9,97 → 6,21 Å (5×5) … 18,51 → 11,12 Å (9×9) … 56,95 → 33,30 Å (27×27), comme au §3.2.
2. **Abscisse de `fig_locality_final` et `fig_locality`** (§3.5).
   - Grille élargie : points à x ≥ 12,5–16,5 a, w ≤ 4,9e-3 eV, déplacés jusqu'à 11,7 a.
   - Grille grossière : 2 (5×5), 6 (7×7) et 9 (8×8) croix déplacées, w entre 0,47 et 2,18 eV.
   - Le mémoire cite « Distance |R| (a) » sur l'axe (`make_figures_memoire.py:70`). `défauts.tex` n'est pas sur Rorqual. Les fichiers
     `figures/fig_locality_final.{pdf,png}` et `figures/fig_locality.{pdf,png}` existent dans le dépôt.

### (c) Effet sur une sélection qui entre dans un nombre

**c1. Atome de Lu (`far_atom`) → valeur de Lu.**

Valeurs de Lu (1,0 Å), en meV. Pour 0,5 Å : Δ = +386,96, +6,26, +3,15, +86,08, −36,88 et −1,59 meV (mêmes tailles, même ordre).

| N | Lu publiée (atome axe) | atomes vraiment les plus loin | moyenne | ΔLu = moyenne − publiée | plateau (i), pour mémoire |
|---|---|---|---|---|---|
| 5 | −122,33 | +176,14 | +176,14 | **+298,47** | −57,56 |
| 6 | −98,62 | −100,55 ; −86,57 | −93,56 | +5,05 | −50,51 |
| 7 | −38,54 | −37,04 | −37,04 | +1,50 | −26,91 |
| 8 | −22,97 | +44,78 | +44,78 | **+67,75** | −14,41 |
| 9 | −24,65 | −42,52 ; −49,24 | −45,88 | **−21,23** | −25,14 |
| 12 | −29,24 | −32,57 ; −28,66 | −30,61 | −1,38 | −18,69 |

Nombres touchés (Lu entre tel quel) :

- **R4.** D1 (9×9 non relaxée D, relaxées N1, N2 : décalage retiré des états de la fenêtre) et D5 (Lu 0,5 et 1,0 Å, 5…12).
- **R5 C (5…12).** Niveaux alignés x = ε − Lu − E_D des états π et σ.
- **R6 étape 2, D4.** Colonne QE π −0,737, σ +0,101 : c'est R5 C 9×9.
- **R7 D1 (6…27).** Niveaux π et σ, et ajustements (π, 1/N, ε∞ −0,052 ; σ plateau +0,102).
- **R9.**
  - Colonne Lu d'A.1 ; C_N = Lu, retenu pour les six tailles (« sans plateau »).
  - Variantes « aligné Lu » et « exact Lu » d'A.3, B, C et D ; colonne Lu de la clôture.
  - Colonne QE de C.1 (niveaux de R7).
  - Carte D.1 « exact » (C₉ = Lu 9×9).
- Si la table 3.7 du mémoire reprend les niveaux QE alignés, elle est touchée aussi (`défauts.tex` absent de Rorqual).

Écart sans recalcul :

- **Niveaux alignés.** x se déplace de −ΔLu (même état, w₂ inchangé, fenêtre [−3 ; +1] eV loin des états π et σ). Pour R5 C et R7 D1 :

  | N | π : publié → vraie | σ : publié → vraie |
  |---|---|---|
  | 5 | −0,328 → −0,627 | −0,131 → −0,429 |
  | 6 | −1,065 → −1,070 | +0,156 → +0,151 |
  | 7 | −0,495 → −0,496 | −0,282 → −0,284 |
  | 8 | −0,370 → −0,438 | −0,195 → −0,263 |
  | 9 | −0,737 → −0,716 | +0,101 → +0,123 |
  | 12 | −0,551 → −0,549 | +0,114 → +0,115 |

- **Non évalué** (il faudrait les moyennes de sphère autour des atomes vraiment les plus loin, donc lire V_ks) : 10×10 et 11×11
  (R4 D5, R5 C), 15…27 (R7 D1 et ses ajustements), R1 N1 et N2 (R4 D1).
- **Variantes « Lu » de R9 A.3, B, C et D.** L'entrée change de ΔLu, par exemple −21,23 meV en 9×9. Les sorties ne sont pas recalculées.
  Pour échelle, la clôture de R9 a déjà les variantes « plateau », avec C₉ = −25,14 contre −24,65 meV (Lu) et C₁₂ = −18,69 contre
  −29,24 meV.

**c2. Profil radial d'`analyze_Ved.py` (anneaux au-delà de 0,433 a_sc).**

- Nombres touchés :
  - `ved_analysis.npz` (`{S}_rad`, `{S}_rad_masked`) ;
  - `fig_Ved` (c), du mémoire, via `make_figures_memoire.py:127–133` ;
  - `fig_Ved_radial` et `fig_Ved_radial_masked` (`make_figures.py`) ;
  - la valeur « r = half box » du journal d'`analyze_Ved` (moyenne des trois derniers anneaux).
- Écart (§3.4) : jusqu'à 16,5 et 19,9 meV (5×5), 14,2 et 10,9 meV (6×6), 6,7 et 10,2 meV (7×7), ≤ 2,2 meV pour 8…12. Ce sont les
  valeurs non masquées puis masquées, sur 17 à 40 anneaux par taille.

**c3. Transplant R7c.**

- Nombres touchés : géométries initiales des relaxations R7c 15…27 (`relax.in`, `transplant.log`). 24 sites ne reçoivent pas leur
  déplacement ; |u| va de 0,004 à 0,013 Å.
- Les géométries relaxées et ΔE (−365,6 et −363,0 meV pour 15 et 18) ne dépendent de la géométrie initiale que par le chemin BFGS
  (forc_conv_thr 1e-4). L'écart n'est pas évaluable sans refaire une relaxation.

**c4. Interpolation hors grille avec les étiquettes axe par axe.**

- **D.1 de R9 (§3.6).** Max|ΔṼ| = 0,094 eV Å² dans les disques et 0,443 eV Å² sur toute la carte. Les moyennes de disque bougent de
  ≤ 0,009 eV Å² (brut → vraie), sur 76,96–82,71 eV Å².
- **Référence « toutes les R » de tab:rcut_M** (`m_rcut_convergence.py:30–44`, et `rcut_table` de R9 A.3).
  - La grille fine 60² n'est pas un multiple de D (27, 24…) : la référence dépend des étiquettes lointaines.
  - Les amas R_cut ≤ 6 ne changent pas.
  - Le max sur (k′, k) de 60² × 60² paires n'a pas été recalculé : c'est le coût d'un passage de `rcut_table`, hors du nœud de connexion.
  - Ordre de grandeur, par le même mécanisme (D.1, 9×9, un k) : 0,084 eV sur |M_π|, soit 3,3e-3 de max|M_π| = 25,44 eV. Pour
    comparaison, les valeurs de tab:rcut_M en 9×9 sont 6,68e-2 (R_cut 3), 1,75e-2 (5) et 1,21e-2 (6).

### Hors question (métrique, pas image) : constats consignés

- **Norme de R_cut.** R_cut est une norme sur les indices réduits. Il y a 29 mailles pour R_cut 3, contre 37 pour un disque cartésien de 3a.
  C'est déjà consigné dans `NOTES_TGAMMA.md` §4.2.
- **Point le plus proche de K.** Il est choisi en norme réduite ; les égalités C₃ cartésiennes sont départagées par cette norme (§3.7).
- **Commentaires de `make_inputs_r7_relax.py`** (l. 27, 105 ; docstring). Ils donnent L₁₂/2 = 14,795 Å, puis L/2, comme borne
  d'« unicité de l'image minimale ». La borne de la réduction `d − rint(d)` est (√3/4) L₁₂ = 12,81 Å (§0).

## 5. Contour de zone (écart 3)

- **Origine.** `scripts/analyze_M.py:56–57`, repris ci-dessous :

  ```python
  corners = [B·c for c in (1/3, 1/3), (1/3, −2/3), (−2/3, 1/3), (−1/3, −1/3), (−1/3, 2/3), (2/3, −1/3)]
  ```

  C'est rangé par angle et écrit sous la clé `map_corners` de `results/M2/M_analysis.npz` (l. 60). Le tableau stocké est égal à la formule
  (vérifié à 1e-12).
- **Utilisateurs.** Aucune figure du dépôt.
  - `fig_M_map` (`make_figures.py:131–138`) et `fig_M_map_final` (`make_figures_memoire.py:92–96`) tracent le contour par
    `bz_vertices(map_B)` : Voronoi du réseau réciproque, vraie zone.
  - Les coins d'`analyze_M` n'ont servi qu'au premier passage de la figure D de R9 (`kaasbjerg_fig3_map`). Ils ont été remplacés par les
    points K au passage `dfig` (R9, partie D).
  - Côté mémoire : `défauts.tex` n'est pas sur Rorqual. `figures/fig_M_map.{pdf,png}` et `figures/fig_M_map_final.{pdf,png}` existent.
- **Chiffres** (cellule QE : a₁, a₂ à 60°, b₁·b₂ = −½|b|², |b| = 2,9423 Å⁻¹ ; K = (2/3, 1/3) de `production.json` = (1,4711 ; −0,8494) Å⁻¹,
  |K| = 1,6987 Å⁻¹ ; `map_K` replié = (0 ; 1,6987) Å⁻¹) :

  | Coins | réduits | \|p\| (Å⁻¹) | sommet de Wigner-Seitz |
  |---|---|---|---|
  | actuels | (1/3, 1/3), (−1/3, −1/3) | 0,9808 | non (0 G équidistant) |
  | actuels | (1/3, −2/3), (−2/3, 1/3), (−1/3, 2/3), (2/3, −1/3) | 2,5948 | non (1 G équidistant) |
  | proposés | (2/3, 1/3), (1/3, 2/3), (−1/3, 1/3), (−2/3, −1/3), (−1/3, −2/3), (1/3, −1/3) | 1,6987 = \|K\| | oui (2 G équidistants) |

  Coins proposés en cartésien (Å⁻¹) : (±1,4711 ; ±0,8494) et (0 ; ±1,6987).
- **Remplacement proposé** pour `analyze_M.py:56` (non appliqué) : les six points K et K′ de la cellule QE, soit
  `[(2/3, 1/3), (1/3, 2/3), (−1/3, 1/3), (−2/3, −1/3), (−1/3, −2/3), (1/3, −1/3)]`. C'est la liste employée par R9 D (`r9_driver.py:1552`).
- **Constat voisin**, même nature, hors question : `scripts/_bands.py:16` (`SYM_POINTS["K"]`). Les images de K listées sont (2/3, 1/3),
  (1/3, 2/3), (−1/3, 1/3), (1/3, −1/3), (−1/3, −2/3), (2/3, −1/3). La norme de la vraie image de chacune est donnée dans
  `audit_results.json` → `zone_analyze_M.bands_SYM_POINTS_K` : 1,6987 Å⁻¹ (= |K|) pour les cinq premières, 0,9808 Å⁻¹ pour (2/3, −1/3).
  (2/3, −1/3) équivaut à (−1/3, −1/3) modulo G ; ce n'est pas un point K de cette cellule.

## 6. Écart 1 : libellé « convention intensive »

- **Ligne exacte** : `scripts/analyze_M.py:171`.

  ```python
  tests.append(("convention intensive (cellule unitaire)", "(max−min)/moyenne de max|M| sur N = 5,7,8,9 (bandes 1–16)", f"{spread:.1e}", "5e-2", "OK" if spread < 5e-2 else "À VOIR", "analyze_M.py"))
  ```

- **Ce que calcule la ligne** : `spread` est calculé l. 170 sur `SIZES`, l. 21, c'est-à-dire les tailles dont `M_dense` et `M_ed` sont
  présents : 5, 6, 7, 8, 9, 12.
- **Valeurs** (`M_analysis.npz`, max|M| des bandes 1–16, eV) : 23,850 (5×5) ; 25,148 (6×6) ; 23,818 (7×7) ; 24,311 (8×8) ; 25,348 (9×9) ;
  25,721 (12×12).
- **Étendues** :
  - sur les six tailles : 7,706e-2. C'est la valeur écrite, 7,7e-02 dans `results/M2/M_tests_summary.csv`, ligne « convention intensive
    (cellule unitaire) », verdict « À VOIR » ;
  - sur les tailles du libellé (5, 7, 8, 9) : 6,291e-2.
- **Table du mémoire concernée** : tab:tests_M. Elle est citée par `NOTES_TGAMMA.md` §7, constat 1 : « 7,7e-2 pour un seuil de 5e-2 ».
- **Pour correction dans le projet mémoire** : soit le libellé « N = 5, 6, 7, 8, 9, 12 », soit la valeur sur 5, 7, 8, 9 (6,3e-2).
  Rien n'est modifié.

## 7. Corrections proposées (non appliquées)

**P-c1. `alignment.far_atom` (`defects/alignment.py:35–39`).**

- Remplacer `d = min_image_dist(...)` par `d = true_min_image_dist(x_red_d, s_vac, A_cols)`. La fonction existe déjà (P1, l. 104).
- Atome de Lu : argmax_i min_{s∈{−1,0,1}³} |A (red(x_i − s_vac) + s)|.
- Règle d'égalité à fixer : plus petit indice, ou moyenne des décalages sur l'ensemble à r_max (2 atomes pour N = 3m, 1 sinon).
  C'est un choix de définition de Lu.
- `min_image_dist` (l. 22) peut rester pour les petites distances (voisins, `vacancy_site`, partenaire), qui sont exactes ;
  alternative : la renommer et la documenter.
- **Test qui montre le défaut.** Géométrie 9×9 de production : `far_atom` doit rendre une distance vraie de 12,813 Å (= a_sc/√3), atome 21
  ou 140. Il rend aujourd'hui l'atome 161 à 18,51 Å axe par axe, soit 11,12 Å en vrai.
- **Test synthétique.** Cellule à 60° de côté L, point en coordonnées réduites (0,45 ; 0,45) :
  - `min_image_dist` donne 0,779 L ;
  - `true_min_image_dist` donne 0,507 L.

**P-c2. `analyze_Ved.py:28–30` et `:79–80`.** Calculer R par la vraie image, avec r_true = min_s |(S1 + s₁) a₁ + (S2 + s₂) a₂|.

- **Test.** Le nombre de points par anneau à r ∈ [0,433 a_sc ; 0,5 a_sc) doit valoir, à l'aire près, celui d'un anneau complet. En 9×9,
  l'anneau à 9,625 Å a 490 points au lieu de 534, et celui à 11,075 Å en a 178 au lieu de 518.

**P-c3. `make_inputs_r7_relax.py:55–58`.** Faire de `minimg_cart` la vraie image (9 images, comme le contrôle de la l. 72).

- **Test.** Sur la 12×12 R2, les sites sous R_CUT = 14,5 Å doivent être 249, pas 225.
- Géométries R7c existantes : ne rien refaire sans décision (§4 c3).

**P-c4. Étiquettes R pour l'interpolation hors grille** (`wannier_interpolation.Mwr_to_Mwk` l. 91, `Mwr_to_Mwk_pairs` l. 115, et leurs
appelants hors grille : D.1 de R9, `m_rcut_convergence.py:30–44`).

- Recentrer sur R_d, puis prendre les images de Wigner-Seitz de R − R_d dans le super-réseau D, avec un poids 1/n sur les n images à
  égalité : `M(k′, k) = Σ_{R,R′} w_R w_R′ Σ_{images} e^{−2πi k′·R̃} M(R, R′) e^{2πi k·R̃′}`.
- **Test.** Sur la grille MP D×D, le résultat est inchangé (à 1e-12). Sur la carte D.1, la convention vraie redonne celle du §3.6
  (0,094 eV Å² dans les disques).

**P-b2. `mwr_locality_coarse_vs_dense.py:19–20`.** Remplacer `dcart = |Rn·A|` par min_s |(Rn + D s)·A|, avec D = MP.

- **Test.** L'abscisse maximale doit valoir D/√3 a au plus : 15,59 a pour la 9×9 dense, alors qu'elle vaut 22,52 a aujourd'hui.

**Écart 3 : `analyze_M.py:56`** (voir §5) et, à part, `_bands.py:16`.

**Écart 1 : `analyze_M.py:171`** (voir §6), pour le projet mémoire.

## Fichiers

- `r9_driver.py` : sous-commande `audit` (docstring d'en-tête) ; fonctions `red_axis`, `cart_axis`, `cart_true`, `aud_*`, `cmd_audit`.
  Rien d'autre n'est changé dans le pilote.
- `audit/audit_results.json` : tous les chiffres de ce rapport, par partie (`sizes`, `R1`, `R_grids`, `interp_D1`, `k_grids`,
  `crowns_240`, `transplant_R7c`, `zone_analyze_M`, `label_analyze_M`, `runs`).
- Copie : `article/R9_controles/audit_image_minimale.md`, `article/R9_controles/audit/audit_results.json` et
  `article/R9_controles/r9_driver.py`.

STOP — audit terminé le 2026-09-28. Aucune correction appliquée.
