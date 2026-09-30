# R10 — Base cohérente pour les chiffres du chapitre 4 — rapport de campagne

Statut **TEST** (post-traitement seul, aucun calcul QE, aucun `.save` produit) ; les sorties de C deviendront la production du ch. 4 dans un
nouveau répertoire de résultats (à confirmer au GO 2). Prompt R10 (Greg, 2026-09-29). Ordre : phase 0 → STOP → étape G (Greg) → GO 1 : A, B →
rapport intermédiaire → STOP → commit du code (Greg) → GO 2 : C → rapport → STOP. Dépôt `graphene-raman`, HEAD 7a64154 (arbre propre) ; racine
`$GRAPHENE_RAMAN`. Répertoire de travail `graphene/qe/defects/R10_plateau/`, copie versionnée `article/R10_plateau/`. Chiffres bruts, sans
interprétation. Les lignes marquées « arithmétique » sont des produits de chiffres déjà publiés ou du cône de Dirac (ħv_F = 5,459 eV Å, R5 C ;
|b| = 2,9422 Å⁻¹ ; ħv_F|b| = 16,06 eV) : elles dimensionnent le plan, ce ne sont pas des mesures.

## Phase 0 (2026-09-29 ; rien n'est calculé)

Opérations : lectures (rapports R4–R9, audit de l'image minimale ; `src/`, `scripts/`, `config/` ; pilotes `r5_driver.py`, `r7_driver.py`,
`r9_driver.py` ; json et npz de R5 C, R7 D1, R9 A/B/C/D, audit ; `results/M2/` : liste, en-têtes csv, clés des npz), `stat` des fichiers,
`sacct` des jobs de R6 étape 3 et de R9. Aucun job, aucune écriture hors de ce répertoire et de sa copie. Rien dans `src/`, `scripts/`,
`config/`, `results/`, `figures/`.

### 0.1 Potentiels et géométries

**Chemins.** Pour toutes les tailles, les pilotes R4/R5/R7/R9 lisent les mêmes fichiers :
- potentiels pp.x (`plot_num = 1`, Ry) : `graphene/qe/defects/super_cell/NxN/{defective,pristine}/Vks_NxN_{d,p}`, sur **/project**
  (`/lustre09/project/…`) : pas de purge ;
- XML (positions, cellule, grille FFT, valeurs propres, occupations) : scratch `/lustre10/scratch/gregb26/qe_tmp/defect_NxN_{d,p}/defect_NxN_{d,p}.save/data-file-schema.xml`,
  miroir `graphene/qe/qe_tmp_backup/defect_NxN_{d,p}/…/data-file-schema.xml` présent pour les 26 fichiers (même taille) ;
- pour 5…12, `data/graphene/supercell/qe/defect_NxN_{d,p}.save/` (chemins de R9 A.1) ne contient que des liens vers ces mêmes fichiers.

| N | atomes (d) | grille FFT | Vks_d (Mo) | Vks : dernier accès | XML d / p (ko) | XML : dernier accès | source des niveaux |
|---|---|---|---|---|---|---|---|
| 5 | 49 | 150 × 150 × 192 | 74 | 2026-09-28 | 43 / 42 | 2026-09-28 | R5 C |
| 6 | 71 | 180 × 180 × 192 | 107 | 2026-09-28 | 54 / 60 | 2026-09-28 | R5 C, R7 D1 |
| 7 | 97 | 216 × 216 × 192 | 154 | 2026-09-28 | 59 / 71 | 2026-09-28 | R5 C |
| 8 | 127 | 240 × 240 × 192 | 190 | 2026-09-28 | 68 / 84 | 2026-09-28 | R5 C |
| 9 | 161 | 270 × 270 × 192 | 241 | 2026-09-28 | 80 / 99 | 2026-09-28 | R5 C, R7 D1 |
| **10** | 199 | 300 × 300 × 192 | 297 | 2026-09-25 | 92 / 115 | 2026-09-28 | R4 D5, R5 C |
| **11** | 241 | 360 × 360 × 192 | 428 | 2026-09-25 | 106 / 134 | 2026-09-28 | R4 D5, R5 C |
| 12 | 287 | 360 × 360 × 192 | 428 | 2026-09-28 | 121 / 154 | 2026-09-28 | R5 C, R7 D1 |
| **15** | 449 | 450 × 450 × 192 | 669 | 2026-09-26 | 175 / 224 | 2026-09-28 | R7 D1 |
| **18** | 647 | 540 × 540 × 192 | 963 | 2026-09-26 | 240 / 310 | 2026-09-28 | R7 D1 |
| **21** | 881 | 625 × 625 × 192 | 1 290 | 2026-09-26 | 320 / 415 | 2026-09-28 | R7 D1 |
| **24** | 1 151 | 720 × 720 × 192 | 1 712 | 2026-09-26 | 412 / 538 | 2026-09-28 | R7 D1 |
| **27** | 1 457 | 810 × 810 × 192 | 2 167 | 2026-09-27 | 516 / 680 | 2026-09-28 | R7 D1 |

(en gras : tailles nouvelles pour A.1.) Grilles 11 (360) et 21 (625) non multiples de N : sans effet sur P1 (sphères en espace réel). Purge du scratch :
dernier accès des XML le 2026-09-28 (audit) → pas avant le 2026-11-27 (arithmétique, 60 jours) ; miroir complet de toute façon.

**Mémoire et durée de P1 par taille** (`alignment.atom_sphere_shifts`, lecture par `qe_io.get_pot`) :
- lecture : `get_pot` fait `readlines` + `" ".join` + `np.fromstring` ; pic ≈ 3 × la taille du fichier pendant la lecture, plus les tableaux déjà
  chargés (un V en float64 = 8 o × n₁n₂n₃ : 138 Mo (10), 199 Mo (11), 311 Mo (15), 448 Mo (18), 600 Mo (21), 796 Mo (24), 1 008 Mo (27)) ; estimation
  du pic (arithmétique) : 1,2 / 1,7 / 2,6 / 3,8 / 5,1 / 6,7 / **8,5 Go** pour 10 / 11 / 15 / 18 / 21 / 24 / 27 ;
- durée de lecture mesurée (R7 D1, `d1_results.json` : durée par taille moins lectures des fonctions d'onde = deux Vks + `far_atom_alignment` + valeurs
  propres) : 12 / 17 / 22 / 32 / **46 s** pour 15 / 18 / 21 / 24 / 27 ; R9 A.1 complet (P1 deux rayons + M_W) : 4–39 s par taille pour 5…12 ;
- boucle de P1 : une boîte d'indices locale par atome (≈ 31 × 31 × 23 points à 27×27 pour 1,0 Å) et une recherche du partenaire sur 27 images ;
  ordre de grandeur (arithmétique) ≈ 3 ms par atome et par rayon → ≈ 9 s pour 27×27 (1 457 atomes, deux rayons) ;
- total A.0 + A.1 (13 tailles) : ≈ 5–10 min, pic ≈ 10 Go. Proposé : 1 tâche, 8 cœurs, 64 Go, 1 h.

### 0.2 Niveaux QE alignés sur Lu : sources, E_D, Lu, relation avec le plateau

**Sources** (lecture seule) :
- 5…12 : `article/R5_base_vs_M/c/c_results.json` → `sizes[N]` : `pi_state` (x, w₂, w₁, bande, relEF), `sigma_doublet`, `alignment.shift_10` (Lu 1,0 Å),
  `E_D`, `E_D_source`, `E_F_p`, `E_F_d`, `gamma_gap_p`, `window_P`, `window_D`, `localized` ;
- 6…27 : `article/R7_tailles_3m/d1_8pts/d1_results.json` (mêmes clés) et `window_NxN.npz` (`D_x`, `D_e`, `D_e_aligned` = ε − Lu, `shift_Lu_1A`, `E_D`,
  `e_d_all`, `e_p_all`) ; `d1/` (cinq points) est la version antérieure, mêmes valeurs pour 6…18 ; porte de R7 : 6, 9, 12 = R5 C à 0 près.

**E_D par taille** :
- N = 3m (6, 9, 12, 15, 18, 21, 24, 27) : milieu du quadruplet de Dirac de la parfaite à Γ (`dirac_quadruplet`, étalement ≤ 3,2e-7 eV) ;
- N = 5, 7, 8, 10 : route « maille » E_K(maille 9×9) + Δ_Γ + Δ_fold (Δ_Γ = médiane sur les 8 états les plus bas à Γ de ε(maille N×N) − ε(maille 9×9) ;
  Δ_fold = médiane sur les états sous E_D − 4 eV de ε(parfaite N×N) − ε(maille N×N repliée)) ;
- N = 11 : K du run bands de la maille 11×11 (indice 61) + Δ_fold par l'état le plus bas du même run (la grille 121 k a été écrasée ; R5 écart 4).

| N | E_D (eV) | définition | Lu 1,0 Å utilisé (meV) | atome de Lu (axe par axe) | π : bande ; x ; w₂ | paire σ : bandes ; x ; w₂ |
|---|---|---|---|---|---|---|
| 5 | −4,23138 | maille | −122,33 | 49 | 98 ; −0,3284 ; 0,267 | 99, 100 ; −0,1308 ; 0,702 |
| 6 | −4,23745 | quadruplet | −98,62 | 71 | 140 ; −1,0654 ; 0,284 | 144, 145 ; +0,1560 ; 0,711 |
| 7 | −4,23904 | maille | −38,54 | 97 | 194 ; −0,4948 ; 0,301 | 195, 196 ; −0,2823 ; 0,703 |
| 8 | −4,23768 | maille | −22,97 | 1 | 254 ; −0,3698 ; 0,223 | 255, 256 ; −0,1954 ; 0,708 |
| 9 | −4,23847 | quadruplet | −24,65 | 161 | 320 ; −0,7373 ; 0,246 | 324, 325 ; +0,1014 ; 0,713 |
| 10 | −4,23889 | maille | −8,22 | 198 | 398 ; −0,4224 ; 0,241 | 399, 400 ; −0,2383 ; 0,708 |
| 11 | −4,23851 | K du même run | −12,16 | 241 | 482 ; −0,3319 ; 0,191 | 483, 484 ; −0,1641 ; 0,707 |
| 12 | −4,23869 | quadruplet | −29,24 | 287 | 572 ; −0,5508 ; 0,223 | 576, 577 ; +0,1140 ; 0,713 |
| 15 | −4,23878 | quadruplet | −15,64 | 449 | 896 ; −0,4593 ; 0,207 | 900, 901 ; +0,1033 ; 0,713 |
| 18 | −4,23884 | quadruplet | −18,18 | 1 | 1292 ; −0,3881 ; 0,194 | 1297, 1298 ; +0,1070 ; 0,712 |
| 21 | −4,23886 | quadruplet | −13,77 | 881 | 1760 ; −0,3437 ; 0,183 | 1765, 1766 ; +0,1030 ; 0,711 |
| 24 | −4,23888 | quadruplet | −13,12 | 2 | 2300 ; −0,3074 ; 0,174 | 2305, 2306 ; +0,1026 ; 0,712 |
| 27 | −4,23890 | quadruplet | −12,77 | 1457 | 2912 ; −0,2785 ; 0,166 | 2917, 2918 ; +0,1024 ; 0,712 |

(x = ε − Lu − E_D, eV ; bandes 1-based de la cellule avec lacune.)

**Vérification sur le code** : `r5_driver.window_states_gamma` (l. 603–605) calcule `x = e - shift - E_D` ; R5 C (l. 679) et R7 D1 (`analyze_size`,
appel de `r5.window_states_gamma`) passent `shift = Lu (1,0 Å)` pour la cellule avec lacune et 0 pour la parfaite. Donc
**x_plateau = ε − C_N − E_D = x_Lu + Lu − C_N** (vérifié). Ne dépendent pas du décalage : parité, w₂, w₁, seuil de localisation (3⟨w₂⟩ de la parfaite,
décalage 0), choix de l'état π (impair de plus grand w₂) et de la paire σ (deux pairs de plus grand w₂). Seule la borne de la fenêtre [−3 ; +1] eV bouge
de Lu − C_N ; les états π et σ sont à x ∈ [−1,07 ; +0,16].

### 0.3 E_res (B.1) : fonctions, coût, fenêtre

**Fonctions existantes pour Γ sur une grille de sortie quelconque** :

| fonction | où | grille d'énergie de g₀ | E_res | remarque |
|---|---|---|---|---|
| `lt.scattering_rate_fast` | `src/…/local_tmatrix.py` (niveau 1 de production, `compute_spectral_wannier.py`) | « sig » : [min ε_fenêtre − η ; max ε_fenêtre + η], pas η/8 ; **dépend de la grille de sortie** | argmax \|Γ\| des états à ±1,5 eV (dans le script) | refait g₀ à chaque appel : 8,5 min par appel à 900² (R9) |
| `resonance_setup` / `resonance_observables` | `r9_driver.py` (R9 B) | « res » : E_D ± (3 eV + η), 2 417 énergies, pas 2,5 meV ; **indépendante de la grille de sortie** | argmax \|Γ\| des états à ±1,5 eV (`rint` sur la grille « res ») ; courbes Γ_T à 5 et 2,5 meV ; T̄(K) | grille de sortie lue dans la config (240) : paramètre à ajouter dans le pilote R10 |
| `level1_stats` | `r9_driver.py` | « sig » (argmin, comme `rcut_resigma.py`) | argmax Γ des états à ±1,5 eV, avec l'état | même dépendance que `scattering_rate_fast` |
| `crown_table`, `crown_of_state` | `r9_driver.py` | — | couronne x = \|i b₁ + j b₂\|²/\|b₁\|² autour de K, K′ | tout N ; x entier exigé (K sur la grille) : 120, 240, 480, 960 sont multiples de 3 |

Proposé pour B.1 : chaîne « res » de R9 B (un seul g₀ 900² pour les quatre grilles et les deux variantes ; porte B.0 à 240² contre R9 B).

**Coût à 480² et 960²** : job R9 B aligné 21955646 = 15 min 07 pour 16 passages (2 tailles × 4 N_k^int × 2 variantes, g₀ en cache) → ≈ 55 s par passage à 240²
(41 266 états dans ±3 eV). Le coût suit le nombre d'états de la fenêtre (arithmétique) : ×4 à 480², ×16 à 960² (≈ 660 000 états) → ≈ 15 min par variante à 960².
Mémoire à 960² : amplitudes φ sur l'amas (921 600 × 5 × 145, complexe) 10,7 Go ; matrice des lorentziennes de la courbe 2 417 × 660 000 × 8 o = 12,8 Go par
tableau temporaire (à découper en tranches d'énergie dans le pilote : les sommes par ligne ne changent pas) → pic ≈ 40–50 Go. g₀ 900² « res » : 8,5 min (R9,
16 fils). Total B.1 en fenêtre complète : ≈ 1 h, 128 Go demandés. Restreint à [−1,0 ; +0,5] eV (marge de 0,5 eV = 25 η pour la courbe) : ≈ 32 000 états à 960²
(arithmétique, cône de Dirac), coût ≈ celui de 240².

**Fenêtre [−0,5 ; 0] eV sur 240²** (lecture des sorties de R9) :
- les états (Γ par état) ne sont pas stockés dans les sorties de R9 (`b_curves_*.npz` : courbes seulement) ; toutes les valeurs d'E_res de R9 (A.3, B, R.1, R.3 :
  six tailles, trois variantes, R_cut 3 et 4, N_k^int 300–900) sont dans **[−0,284 ; −0,134] eV** ⊂ [−0,5 ; 0] ;
- courbes Γ_T 9×9 (N_k^int 300 et 900, tel quel et plateau) : maximum hors de [−0,5 ; 0] (dans ±1,5 eV) = 36,5 / 37,2 / 32,0 / 32,6 % du pic, toujours à −0,5025 eV
  (bord de la fenêtre).

**Couronnes du cône de Dirac près de −0,2 eV** (arithmétique, E = ħv_F|b|√x/N, gauchissement trigonal ignoré) :

| grille | couronnes dans [0,12 ; 0,26] eV | exemples (x : E, eV) | écart entre couronnes voisines |
|---|---|---|---|
| 120² | 2 | 1 : 0,134 ; 3 : 0,232 | 98 meV |
| 240² | 5 | 7 : 0,177 ; 9 : 0,201 ; 12 : 0,232 | 9,5–43 meV |
| 480² | 16 | 28 : 0,177 ; 31 : 0,186 ; 36 : 0,201 ; 37 : 0,204 | 2,4–14,5 meV |
| 960² | 57 | 144 : 0,201 ; 147 : 0,203 ; 148 : 0,203 | 0,6–6,3 meV |

À 480² et 960², des couronnes voisines sont plus proches que le pas de la grille d'énergie de g₀ (2,5 meV = η/8, `ne_per_eta` gelé) : les états d'un même pas partagent
le même t(ε). Fait rapporté pour la lecture de B.1 (décision D8).

### 0.4 Familles (B.2) : où lire

| grandeur | source | remarque |
|---|---|---|
| E_F(d), E_F(p) | `c_results.json` (`E_F_d`, `E_F_p`) = XML `<band_structure><fermi_energy>` (Ha) | ΔE_F = E_F(d) − E_F(p) aussi dans `results/M2/sampling_table.csv` (`dE_F_meV` : −235,8 / −20,5 / −866,1 / −583,3 / −2,4 / −508,7 / −209,3 / +2,5 pour 5…12) |
| valeurs propres à Γ de la parfaite près de E_D | XML de la parfaite (`<ks_energies><eigenvalues>`) ; `c_results.json` : `gamma_gap_p` (homo, lumo, gap, par comptage d'électrons), `window_P` (x, parité, w₂ des états de la fenêtre) | gap 0 pour 3m (quadruplet à demi rempli) |
| états π et σ de la lacune | `c_results.json` : bandes du tableau de 0.2, `relEF` = ε − E_F(d) (sans alignement), x | — |
| occupations | XML de la cellule avec lacune : `<ks_energies><occupations size="nbnd">`, un seul k (Γ, poids 2), valeurs dans [0 ; 1] par état ; `<smearing degauss="5.0e-3">mv</smearing>` (Ha) = MV 0,01 Ry | **aucun lecteur dans `src/`** (grep « occupations » vide) : lecture ElementTree dans le pilote, comme `nelec` en R5/R7. Contrôle de format (9×9) : bande 320 (π) 1,000, bandes 324–325 (σ) 0,158 |
| C_N plateau | A.1 (10, 11) ; R9 `a/a1_results.json` → `C_i_eV` (5…9, 12), redonné par A.0 | — |
| max\|M\| dense (bandes 1–16) | `results/M2/M_analysis.npz` → `scale_<S>_dense` (5…9, 12) | 23,850 / 25,148 / 23,818 / 24,311 / 25,348 / 25,721 eV |
| pic de la courbe Γ_T aligné plateau | `R9_controles/a/a3_results_plateau.json` → `level1[S]` (5…9, 12 ; pas 2,5 meV) | 10×10 et 11×11 : pas de M dense |
| sampling_table.csv | `results/M2/sampling_table.csv` : N, N mod 3, dE_F, KS dense max/moyen, V_ed au bord, V_ed radial à 1,42 Å, masqué à 2 et 3 Å | — |

### 0.5 Code pour C (rien n'est écrit)

#### (a) Réservé à Greg — propositions

**(a1) Fonction de Wigner-Seitz partagée** — module `wannier/wannier_interpolation.py` (à côté de `Mwk_to_Mwr`) :

```
ws_images(R, R_center, MP, A_cols, tol=1e-6) -> dict
    R        (nR, 3) int   étiquettes de Mwk_to_Mwr (boîte MP D₁×D₂×1, arange(D) − D//2)
    R_center (3,) int      centre (R_d de recenter_mwr, étiquettes brutes)
    MP       (3,) int      grille MP (D₁, D₂, 1)
    A_cols   (3, 3) ou (2, 2) float   vecteurs de la maille en colonnes (unité de longueur libre : a ou Å)
  -> R_img (n_img, 3) int  images R_r + (D₁s₁, D₂s₂, 0), s ∈ {−1, 0, 1}², les plus proches de R_center
     idx   (n_img,)  int   indice r de l'étiquette d'origine
     w     (n_img,)  float 1/n_r sur les n_r images à égalité (|d − d_min| ≤ tol)
     dist  (nR,)     float distance vraie |A (R_img − R_center)| (unité de A_cols)
     n_tie (nR,)     int
ws_phase(k, ws, nR, sign) -> (nk, nR) complex      Φ(k, r) = Σ_{a : idx_a = r} w_a exp(sign · 2πi k·R_img,a)
Mwr_to_Mwk(Mwr, R, k, ws=None) ; Mwr_to_Mwk_pairs(Mwr, R, k_bra, k_ket, ws=None)
    avec ws : Mwk[w, k′, W, k] = Σ_{r, r′} Φ⁻(k′, r) Mwr[w, r, W, r′] Φ⁺(k, r′)   (Φ± = ws_phase(·, ws, nR, ±1)) ; ws=None : inchangé au bit
```

Recherche : R_r − R_center réduit d'abord axe par axe dans [−D/2 ; D/2) (9 images suffisent, audit §0), puis vraie distance dans la métrique de la cellule à 60° ;
même convention que la variante « vraie » de l'audit §3.6 (`r9_driver.aud_interp_D1`). `ws_phase` est exposée parce que `m_rcut_convergence.py` construit ses
matrices de phases lui-même (l. 33–34) au lieu d'appeler `Mwr_to_Mwk`.

Tests proposés (`tests/test_r10_functions.py`, synthétiques, rapides) : (T1) grille MP D×D (k = m/D), Mwr aléatoire, cellule à 60° : `Mwr_to_Mwk(ws)` et
`Mwr_to_Mwk_pairs(ws)` = sans ws à 1e-12 ; (T2) Σ_{a : idx_a = r} w_a = 1 pour tout r, chaque image ≡ R_r modulo D ; (T3) cas de l'audit : D = 20, R_center = 0,
R = (9, 9, 0) → images (−11, 9, 0) et (9, −11, 0), w = ½, dist = √103 a = 0,507 × (20 a) ; norme axe par axe 9√3 a = 0,779 × (20 a) ; (T4) étiquette déjà dans
la cellule de Wigner-Seitz : une image, w = 1, R_img = R. Contrôle d'intégration (porte C.0 c) : carte D.1 9×9 avec ws = « vraie | brut » de l'audit
(76,955553 / 78,168368 / 81,176724 / 82,305361 eV Å², `audit_results.json`), sans ws = R9 D.1 (76,955098 …).

**(a2) Fonction unique rotation → TF → alignement → recentrage** — module `defects/many_body/local_tmatrix.py` (à côté de `recenter_mwr`, `extract_V_loc`) :

```
defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N=0.0) -> dict(Mwr, R, Rn, R_d, in_box)
    Mwk = Mbk_to_Mwk(Mbk, U, U_dis) ; Mwr, R = Mwk_to_Mwr(Mwk, k, MP)
    in_box[r] = (R_r mod D)_a ∈ [0, n_box) pour a = 1, 2      (étiquettes brutes ; = box_geometry de r9_driver.py)
    Mwr[w, r, w, r] −= C_N  pour tout w et tout r de la boîte    (approximation (i))
    Rn, R_d = recenter_mwr(Mwr, R, MP)
    Mbk (nb, nk, nb, nk) eV ; U (nk, nw, nw) ; U_dis (nk, nb, nw) ; k (nk, 3) (tel que passé par chaque script : XML ou arrondi MP) ;
    MP (3,) ; n_box = N ; C_N en eV (unité de Mbk) ; sorties : Mwr (nw, nR, nw, nR), R, Rn (nR, 3) int, R_d (3,), in_box (nR,) bool
```

Grille grossière (D = N) : la boîte couvre toutes les mailles, l'identité M^L[C] = N_cells·C·𝕀 est exacte (R9 0.2, A.2 1,56e-13). Grille dense (D = pN) : les N² mailles
de la super-cellule dans la cellule zéro-paddée. Le recentrage vient après l'alignement : R_d = argmax ‖M_W(R, R)‖ (41,9 eV en 9×9 contre |C_N| ≤ 58 meV) ; test d'égalité
de R_d avec C_N = 0 inclus.

Tests proposés : (T5) M = 0, C_N = 1, grossière N×N (D = N = 4, nw = 2, U unitaire aléatoire) : M_W = −𝕀 sur tous les (R, R) et `Mwr_to_Mwk` = −N_cells·δ_kk′·𝕀 à 1e-12
(signe : ΔV − C_N) ; (T6) M = 0, C_N = 1, dense D = 6, n_box = 3 : `Mwr_to_Mwk` = −D_N(k − k′)·𝕀 (D_N(q) = Σ_{i,j<N} e^{2πi(q₁i+q₂j)}, R9 0.2) à 1e-12, in_box.sum() = N² ;
(T7) C_N = 0 : sorties identiques au bit à la chaîne actuelle (Mbk_to_Mwk → Mwk_to_Mwr → recenter_mwr). Contrôle d'intégration (porte C.0 b) : 9×9 M2 dense,
`extract_V_loc(defect_mwr(…, C_9))` = V_loc(C = 0) − C_9·diag(in_box répété 5 fois) **au bit** (V_variant « aligne » de R9 ; herm V_loc 4,2e-15 < 1e-10 : pas de
symétrisation, donc même opération flottante).

#### (b) Pour Code, après accord de Greg sur (a) — changements proposés par fichier (diffs écrits à l'étape G)

| fichier : lignes | changement | test ou contrôle |
|---|---|---|
| `scripts/analyze_Ved.py:28–30` et `:79–80` (P-c2) | R (distance au site, plan z_C) par la vraie image : `al.true_min_image_dist` sur les points de la grille du plan (i/n₁, j/n₂, z_vac), centre s_vac, cellule en Å ; X, Y (carte) inchangés ; anneaux jusqu'à \|a₁\|/2 inchangés | audit : anneau à 9,625 Å de la 9×9 : 534 points (490 aujourd'hui) ; 11,075 Å : 518 (178) ; anneaux < 0,433 a_sc redonnés à 0,0 |
| `scripts/mwr_locality_coarse_vs_dense.py:19–20` (P-b2) | `dcart = ws_images(Rn, 0, MP, A)["dist"]` (A de la l. 19, unités de a) ; tri et poids inchangés | abscisse max ≤ D/√3 a : 15,59 a pour la 9×9 dense (22,52 aujourd'hui) ; points déplacés = audit §3.5 |
| `scripts/analyze_M.py:169–171` | deux lignes « convention intensive (cellule unitaire), famille 3m (N = 6, 9, 12), M denses, bandes 1–16 » et « … famille non-3m (N = 5, 7, 8) … », chacune (max − min)/moyenne sur sa famille, seuil 5e-2 ; libellés exacts | arithmétique sur R9 : 3m 2,26e-2, non-3m 2,05e-2 (M2 brut) |
| `scripts/analyze_M.py:56–57, 60` | suppression de `corners` et de la clé `map_corners` (aucun lecteur, audit §5) | `M_analysis.npz` sans la clé ; figures inchangées (contour par `bz_vertices`) |
| `scripts/analyze_M.py:160` (fermeture Fourier) | `defect_mwr(…, C_N=0.0)` explicite (test de la chaîne k → R → k, pas de l'alignement) | valeur = M2 (2,4e-14) |
| `scripts/compute_spectral_wannier.py:79–82` | `defect_mwr(M, U, U_dis, k_coarse, MP, n_box=N, C_N=alignment_C(cfg, S))` (dense et grossier) | porte C.0 a/b |
| `scripts/rcut_resigma.py:30–31` | idem | idem |
| `scripts/resonance_metrics.py:39` (`vloc_from`) | idem pour V ; V_sub (C14 ancien) : décision D10 | porte C.0 |
| `scripts/resonance_criteria.py:33` | idem | — |
| `scripts/m_rcut_convergence.py:25` et `:33–34` | `defect_mwr` ; phases Pm, Pp par `ws_phase(kf, ws, nR, ∓1)` avec `ws = ws_images(Rn, 0, MP, A)` (référence « toutes les R », P-c4) ; masques R_cut sur Rn (amas ≤ 6 inchangés, audit) | C.2 : colonne « plateau, étiquettes actuelles » = R9 (2,470e-2 / 3,097e-2) |
| `scripts/mwr_locality_coarse_vs_dense.py:16` | `defect_mwr` (dense et grossier) | — |
| `scripts/make_figures_memoire.py` (fig_Ved (c), l. 127–133) | trait horizontal pointillé à C_N (config) pour chaque taille tracée, même couleur (C.3) | figure |
| `scripts/check_onsite_and_NL.py:48–50` (diagnostic, hors chaîne) ; `local_tmatrix.scattering_rate_from_wannier` l. 254–257 (aucun appelant) | branchement facultatif (décision D11) | — |
| lecture de C_N (tous les scripts ci-dessus) | `alignment_C(cfg, S)` (config.py, (c)) ; refus si la taille manque | — |
| chemins des matrices (option A de D4) | `MAT = matrices_dir(cfg)` : `compute_spectral_wannier.py:60` (M_ed grossier), `mwr_locality_coarse_vs_dense.py:31` (M_ed), `analyze_M.py:21, 116, 135, 148, 151` ; `dense_paths` couvre M_dense, M_L_dense, M_NL_dense | `grep RES}/M_` vide après changement |

Restent inchangés : `wannier_interpolation.wannier_interpolate` (générique, seul appelant `test_wannier.py`), scripts de test (`test_local_tmatrix_real.py`,
`test_wannier.py`, `test_local_rcut.py` : chaîne brute), `submit_*.sh` (voir 0.6, contraintes du lanceur). Non demandés ici : P-c1 (`far_atom`, Lu abandonné), P-c3
(transplant R7c), `_bands.py:16`.

#### (c) Config (Greg applique)

`config/production.json` (proposé) :

```
-  "results_dir": "results/M2",
-  "results_dir_frozen": "results/M",
+  "results_dir": "results/M2_plateau",
+  "matrices_dir": "results/M2",
+  "results_dir_frozen": ["results/M", "results/M2"],
+  "alignment": {
+    "version": "plateau (i), R10, 2026-09-xx",
+    "method": "C_N = moyenne des décalages de sphère de 1,0 Å (alignment.atom_sphere_shifts) sur les atomes à distance vraie >= 0,75 r_max
+               de la lacune ; M_W(R,R) - C_N·I_5 sur les mailles R mod D dans [0,N)^2 (étiquettes brutes de Mwk_to_Mwr), approximation (i),
+               rien d'autre ; matrices M2 brutes",
+    "C_N_eV": {"5x5": …, "6x6": …, …, "27x27": …},          (13 tailles, valeurs de A.1)
+    "source": "article/R10_plateau/a/a1_results.json", "source_md5": "…",
+    "labels_offgrid": "wigner_seitz"
+  }
```

`src/electron_defect_interaction/config.py` (proposé) : clés obligatoires + `matrices_dir`, `alignment` ; `matrices_dir(cfg, root=ROOT)` ; `alignment_C(cfg, size) -> float`
(eV ; `KeyError` si la taille manque) ; `dense_paths(…)["mfile"]` → `matrices_dir` ; docstring de `results_dir` (« produits » seulement) ; ligne `verbose` avec
`matrices=` et `C_N=`. `.gitignore` : mêmes exceptions pour `results/M2_plateau/` que pour `results/M2/` (npz de production, csv, README, MD5SUMS).

### 0.6 Inventaire de C

Sources : `article/R6_production_corrigee/etape3/table_v1_v2.md`, `scripts/make_figures_memoire.py`, `scripts/make_figures.py`, `scripts/submit_post.sh`,
`etape3/runbook_3.sh`. Coûts = `sacct` de R6 étape 3 (16 cœurs). « Dépend » = dépend de M_W ou de V_loc (donc de l'alignement).

| sortie (`results_dir`) | script (tâche) | figure / table | dépend | action en C | job R6, durée, MaxRSS |
|---|---|---|---|---|---|
| `specwd_<S>_prod.npz` ×6 | `compute_spectral_wannier.py` (`submit_spectral_wannier_dense.sh`, grilles 60/120/240, η 0,05/0,02/0,01, R_cut 0…4) | fig_convergence (a, b), fig_rcut, fig_plateau, fig_level2 ; `level1_summary.csv`, `level2_summary.csv` (make_figures.py), `level2_families.csv` ; §2 NOTES_TGAMMA | oui | refait | 21857272–77 : 1 h 38 / 1 h 47 / 52 min / 1 h 23 / 52 min / 2 h 10 (5, 6, 7, 8, 9, 12) ; 16–21 Go |
| `resigma_9x9_rc0123.npz`, `resigma_9x9_rc4.npz` → `m_rcut_resigma.csv` | `rcut_resigma.py` (`submit_rcut_resigma.sh`) ; csv par l'aide R6 `r6_m_rcut_resigma.py` (hors `scripts/`) → reprise dans le pilote (porte : sur les npz de M2, redonne le csv de M2 à 1e-9) | C10, C11 | oui | refait | 21857279, 21857280 : 3 min 05, 2 min 09 |
| `resigma_9x9_rc3_nk{150…600}.npz` → `nkint_check_9x9.csv` | `submit_nkint_check.sh` ; `nkint_check_post.py` | §6 NOTES_TGAMMA | oui | refait | 21857278 : 7 min 20 |
| `mwr_locality.npz` | `mwr_locality_coarse_vs_dense.py` (P-b2) | fig_locality_final, fig_locality ; p_z–p_z sur site (tab. familles) | oui (sur-site) | refait | 21857281 (avec la ligne suivante) : 9 min 32 ; 28 Go |
| `m_rcut_convergence.csv` (9×9, 12×12, R_cut 0…6) | `m_rcut_convergence.py` (P-c4) | tab:rcut_M | oui | refait | idem |
| `resonance_9x9.npz`, `resonance_criteria_9x9.npz` | `resonance_metrics.py`, `resonance_criteria.py` (`submit_post.sh resonance 9x9`) | fig_spectral_final (a–d), fig_spectral ; §2, C15, C16 ; entrée de fig_epw_vs_ed | oui | refait | 21862372 : 6 h 35 (critère à 16 fils BLAS) ; 31 Go |
| `resonance_6x6.npz`, `resonance_12x12.npz` (+ critères) | idem | R6 3.4 (Friedel 12×12) ; absents de table_v1_v2 | oui | décision D9 | 21862373 : 2 h 39 ; 21862374 : 1 h 38 |
| `resonance_9x9_shiftL.npz` | `resonance_metrics.py --shift-L-meV 25,-25` (`c14`) | R6 3.5 (C14 redéfini) | oui | décision D10 | 21862375 : 10 min 22 |
| `ed_vs_ep_24k24q_mv0.02.npz` | `epw_ed_vs_ep.py` (resonance_9x9 + `results/epw/selfen_240_dg0.02_mv0.02_T300.npz`) | fig_epw_vs_ed, NOTES_EPW (ch. 5, R6 3.7b) | oui | refait | secondes (post_fig 55 s) |
| `M_analysis.npz`, `M_tests_summary.csv` | `analyze_M.py` (M2 brut en base de Bloch ; M_W seulement dans la fermeture, C_N = 0) | fig_M_map(_final), fig_M_scaling(_final), fig_Ved_boundary, tab:tests_M ; Re M^L / Re M^NL à K (tab. familles) | non (M2 brut) ; décision D6 | relancé (changement de tab:tests_M, C.3) | 21862371 : 3 min 28 (avec la ligne suivante) |
| `lnl_frobenius.csv` | `lnl_frobenius_all.py` (M^L, M^NL bruts) | tab:L_NL | non ; D6 | relancé (même tâche, contrôle = M2 à 0) | idem |
| `level2_families.csv` | `level2_families.py` (mwr_locality, M_analysis, specwd) | tab. familles | oui | refait | secondes |
| `ved_analysis.npz` | `analyze_Ved.py` (P-c2), 8 tailles 5…12 | fig_Ved (c ; a, b inchangées), fig_Ved_radial, fig_Ved_radial_masked, fig_Ved_zoom, fig_Ved_map, fig_Ved_profile_mean | non (ΔV brut), mais P-c2 | refait | R6 : copie de v1 ; estimation ≈ 10 min |
| `sampling_table.csv` | `sampling_table.py` (ks_reconstruction, ved_analysis) | table d'échantillonnage | non | relancé après ved | secondes |
| `ks_reconstruction.npz` | `ks_reconstruction_all.py` (potentiel de la parfaite) | fig_ks_reconstruction | non | **copié** de `results/M2` (md5), pas refait | 21857284 : 16 min |
| test d'or 5×5 (ligne de tab:tests_M, `GOLDEN_RESULT`) | `r6golden` | tab:tests_M | cohérence de méthode | non rejoué (D12) | 21852238 : 3 h 36 |
| matrices `M_*.npy` | — | — | — | inchangées (M2 brutes) | — |
| figures du mémoire | `make_figures_memoire.py` (6), `make_figures.py` (14), `make_figures_epw.py` (fig_epw_vs_ed et les autres EPW, indépendantes) | ch. 4, ch. 5 | selon les npz | régénérées dans `R10_plateau/fig/` avec `--outdir` (C.4) | post_fig : 55 s |

Sorties du pilote en plus (C.2–C.5) : tab:rcut_M final à trois colonnes, carte de Kaasbjerg plateau (i) + Wigner-Seitz, planches avant/après, table de correspondance.

**Contraintes relevées pour le lanceur de C (aucun changement de `scripts/` nécessaire si elles sont respectées)** :
- `config.results_dir` sert aujourd'hui à la fois de répertoire des **matrices** et des **produits** (`dense_paths` et 8 chemins d'entrée `f"{RES}/M_…"` dans la chaîne : `analyze_M.py` ×6, `compute_spectral_wannier.py` ×1, `mwr_locality_coarse_vs_dense.py` ×1 ; plus les scripts de contrôle `check_*`) : changer
  `results_dir` seul ferait chercher les M dans le nouveau répertoire (décision D4) ;
- les `#SBATCH --output/--error` des `submit_*.sh` pointent sur `results/M2/logs/` (littéral) : le lanceur les remplace par `sbatch --output=… --error=…` vers
  `<nouveau>/logs/` ;
- `submit_post.sh ksrec` recopie `ved_analysis.npz` de `results/M` s'il est absent : non utilisé tel quel (le lanceur appelle `analyze_Ved.py` puis
  `sampling_table.py`) ;
- `submit_post.sh figures` écrit dans `figures/` (`make_figures.py`, `make_figures_memoire.py`, `make_figures_epw.py` sans `--outdir`) : non utilisé tel quel ;
  le lanceur appelle les mêmes scripts avec `--outdir graphene/qe/defects/R10_plateau/fig` et les options NOTES_EPW de `submit_post.sh` ;
- `resonance_metrics.py` : V_sub (« C14 ancien » = M − ⟨diag M^L⟩·𝕀 puis `vloc_from`) serait décalé deux fois si `vloc_from` applique C_N (D10).

**Nom proposé du nouveau répertoire** : `results/M2_plateau` (matrices M2, alignement plateau (i)). Taille attendue ≈ 115 Mo (produits de `results/M2` hors
matrices).

### 0.7 Plan, jobs, coûts ; décisions attendues

#### GO 1 (A, B) — pilote `r10_driver.py`, lanceur `submit_r10.sh` (diff archivé dans `submitted/<jobid>/`, md5 du pilote dans le journal)

| job | sous-commandes | contenu | ressources | durée estimée | dépend de |
|---|---|---|---|---|---|
| J1 `r10a` | `a0`, `a1`, `a2`, `a3` | A.0 : P1 sur 5, 6, 7, 8, 9, 12 → `C_i_eV` de R9 au bit ; valeurs « un atome » vraies de l'audit (`far_true_meV`, moyenne) à ≤ 1e-6 meV (échec : STOP sur A) ; A.1 : 10, 11, 15…27 (1,0 et 0,5 Å ; moyenne, rms, max\|écart\|, n ; un atome vrai : indices et valeurs ; Lu publié ; P1 à l'atome de Lu = `far_atom_alignment` au bit) ; tableau 13 tailles ; figure 13 panneaux (i) au style d'`offset_profiles` (+ (ii) de R9 pour 5…12, npz de R9) ; A.2 : tableau et figure ε − E_D contre 1/N (π, σ ; Lu / plateau ; points seuls) ; A.3 : colonne QE de R9 C.1 avec le plateau, modèles de `c_results.json` tels quels, écart chaîne − QE | 1 tâche, 8 cœurs, 64 Go, 1 h | 10–15 min | — |
| J2 `r10b` | `b0`, `b1`, `b2` | B.0 : porte (R9 B.0 : 9×9 tel quel, 300², pic Γ_T −0,180 (courbe à 5 meV), −Im T̄(K) −0,1775, courbe = `resonance_9x9.npz` à 0,0 ; 240², N_k^int 900 : tel quel E_res −0,2018 (x = 9), pic Γ_T −0,1825 ; plateau −0,1748 (x = 7), −0,1800 ; `b_results(_plateau).json`) ; B.1 : 120², 240², 480², 960² × (tel quel, plateau (i)), N_k^int 900, R_cut 3, η 0,02 : E_res, état, couronne, écarts aux couronnes voisines, pics de la courbe Γ_T (5 et 2,5 meV), \|E_res − pic\| ; tableau ; figure contre 1/N_out ; B.2 : lecture XML + json (porte : E_F et énergies des bandes π, σ relues = `c_results.json` à 0) | 16 cœurs, 128 Go, 3 h | ≈ 1 h en fenêtre complète | afterok J1 (C_N de 10, 11 pour B.2) |

Boucles d'énergie parallélisées du pilote : 1 fil BLAS + fils Python ; g₀ (zgemm) et `eigh` à 16 fils. Figures : `figures/memoire.mplstyle`, `_palette.py`,
français. Rien dans `results/`. Rapport intermédiaire (phase 0, A, B), STOP.

#### GO 2 (C) — après le commit par Greg du code de 0.5 et de la config

| job | contenu | ressources | durée estimée (étalon R6) | dépend de |
|---|---|---|---|---|
| C0 `r10c0` | portes (a) C_N = 0, étiquettes actuelles : 9×9 R_cut 3 médiane Γ·N_cells 3 132,60 meV (0,01 meV), tab:rcut_M 6,679e-2 ; (b) plateau, étiquettes actuelles : 9×9 3 189,01 meV et 2,470e-2 ; 12×12 3 264,51 meV et 3,097e-2 (+ V_loc aligné = R9 au bit) ; (c) Wigner-Seitz : M inchangé à 1e-12 sur la grille MP ; D.1 tel quel = « vraie » de l'audit (valence K 76,956 eV Å²) ; (d) (i) avec C = Lu, étiquettes brutes, contre R9 D.1 exact (81,512 eV Å²) : écart rapporté | 16 cœurs, 128 Go, 2 h | ≈ 45 min | — |
| C1 production (lanceur, sorties dans `results/M2_plateau`, journaux dans `results/M2_plateau/logs`) | specwd ×6 ; nkint ; resigma ×2 ; locality ; analyze (+ `GOLDEN_RESULT` de R6) ; ved + sampling (+ copie md5 de ks_reconstruction) ; resonance 9×9 (+ 6×6, 12×12 : D9) ; c14 (D10) | 16 cœurs, 64–96 Go par tâche | specwd 52 min – 2 h 10 (parallèles) ; res 9×9 6 h 35 ; res 6×6 2 h 39 ; res 12×12 1 h 38 ; autres < 30 min | afterok C0 |
| C1f figures | `make_figures.py`, `level2_families.py`, `nkint_check_post.py`, `make_figures_memoire.py`, `epw_ed_vs_ep.py`, `make_figures_epw.py` (options NOTES_EPW), tous avec `--outdir R10_plateau/fig` quand l'option existe ; pilote `mrr` (m_rcut_resigma.csv) | 16 cœurs, 96 Go, 1 h | ≈ 5 min | afterok C1 |
| C2 `r10c2` | tab:rcut_M final (9×9, 12×12, R_cut 0…6) : tel quel / plateau (étiquettes actuelles) / plateau (Wigner-Seitz) ; porte : la dernière colonne = `m_rcut_convergence.csv` de C1 à 1e-12 ; carte de Kaasbjerg (R9 D.1, D.2) plateau (i) + Wigner-Seitz, disques K, K′, valence, conduction | 16 cœurs, 128 Go, 2 h | ≈ 45 min | afterok C0 |
| C3–C5 | contrôles des sorties touchées par l'audit (P-c2 : anneaux ; P-b2 : abscisses ; tab:tests_M par famille), planches avant/après par figure (`figures/` en lecture contre `R10_plateau/fig/`), table de correspondance v2 tel quel → v2 plateau (format de `table_v1_v2.md`, + ch. 5) | 8 cœurs, 32 Go, 30 min | < 10 min | afterok C1f, C2 |

Total C (arithmétique sur R6) : ≈ 21 h de tâches à 16 cœurs (≈ 340 cœur·h), durée murale ≈ 8 h (critère 9×9). Disque : ≈ 115 Mo dans `results/M2_plateau`,
`cache/` du pilote < 5 Go. Installation dans `figures/` : GO séparé.

#### Décisions attendues (étape G)

1. **D1 — signatures de (a)** : `ws_images`, `ws_phase`, option `ws=` de `Mwr_to_Mwk` et `Mwr_to_Mwk_pairs` (a1) ; `defect_mwr` (a2) ; modules proposés ; tests T1–T7.
2. **D2 — attributions** : Greg écrit (a) (réservé) ; Code écrit (b) après accord sur (a) ; Greg applique (c).
3. **D3 — B.1, chaîne** : « res » de R9 B (un g₀ 900², proposé) ou « sig » du niveau 1 (8 g₀ 900², ≈ 70 min de plus).
4. **D4 — répertoire de résultats** : `results/M2_plateau` (proposé) ; (A) clé `matrices_dir` + changements des chemins de matrices (proposé) ou (B) liens symboliques
   vers les 43 matrices de `results/M2` (aucun changement de script ; risque d'écriture à travers les liens, garde `chmod a-w` sur `results/M2`).
5. **D5 — B.1, fenêtre** : fenêtre complète ±1,5 eV à toutes les grilles (≈ 1 h, 128 Go ; proposé, aucune hypothèse) ou grilles 480² et 960² restreintes (E_res sur
   [−0,5 ; 0], courbe sur les états de [−1,0 ; +0,5] ; justifié par R9 : toutes les E_res dans [−0,284 ; −0,134] eV) avec une porte « restreint = complet » à 240².
6. **D6 — sorties en base de Bloch sur M2 brut** (`M_analysis` : carte, Re M^L à K, max\|M\| ; `lnl_frobenius` ; convention intensive de tab:tests_M) : laissées sur
   M2 brut (pas de M_W ; « M2 restent brutes » ; proposé) ou variante alignée par V(k′)[−C_N·D_N(k − k′)·𝕀₅]V(k)† (approximation (i) ramenée en base de Bloch ;
   arithmétique : ≈ −N_cells·C_N sur la paire π à (K, K), +2,04 eV en 9×9).
7. **D7 — tab:tests_M** : deux lignes par famille seulement (prompt) ou aussi la ligne six tailles avec son libellé exact.
8. **D8 — lecture de B.1 à 960²** : pas de g₀ gelé à 2,5 meV (écarts de couronnes 0,6–6,3 meV) ; aucune variante proposée (un seul changement à la fois), fait rapporté.
9. **D9 — résonance 6×6 et 12×12** dans C.1 (R6 3.4, + 4 h 17 de tâche ; proposé : oui, pour que `results/M2_plateau` contienne les mêmes produits que `results/M2`).
10. **D10 — C14 (`resonance_9x9_shiftL.npz`) et « C14 ancien » (V_sub)** : décalages ±25 meV appliqués au V_loc aligné, V_sub calculé avec C_N = 0 comme en M2
    (proposé) ; ou retirer ces variantes.
11. **D11 — appelants hors chaîne** (`check_onsite_and_NL.py`, `scattering_rate_from_wannier`) : branchés sur `defect_mwr` ou laissés.
12. **D12 — test d'or 5×5** : non rejoué (valeur R6 reprise, proposé) ou rejoué sur la base alignée (3 h 36, 1 nœud).
13. **D13 — critère de résonance 9×9** : mêmes réglages que R6 (16 fils BLAS, 6 h 35 ; proposé, « mêmes paramètres ») ou `OMP_NUM_THREADS=1` par
    l'environnement du lanceur (aucun changement de script ; durée non mesurée).
14. **D14 — plateau sans critère d'arrêt** : C_N = moyenne du plateau (i) pour toutes les tailles quel que soit max\|écart\| (décision du 29 sept. ; le critère
    « sans plateau » de R9 n'est plus appliqué, max\|écart\| rapporté).

**STOP — phase 0 terminée le 2026-09-29.** Rien n'est calculé, rien n'est écrit dans `src/`, `scripts/`, `config/`, `results/`, `figures/` ; attente de l'étape G
puis du GO 1.

## Étape G (2026-09-29) — décisions de Greg

GO 1 : A et B, jobs J1 (`r10a`) et J2 (`r10b`) tels que proposés en 0.7.

- **D1** signatures de (a) acceptées : `ws_images`, `ws_phase`, option `ws=` de `Mwr_to_Mwk` et `Mwr_to_Mwk_pairs` (`wannier_interpolation.py`) ; `defect_mwr`
  (`local_tmatrix.py`) ; tests T1–T7 dans `tests/test_r10_functions.py`.
- **D2** Greg écrit (a) ; Code écrit (b) contre ces signatures après le commit de (a), Greg relit et commit ; Greg applique (c) après A (C_N des 13 tailles lus
  dans `a1_results.json` de R10).
- **D3** chaîne « res » de R9 B (un g₀ 900²). **D4** `results/M2_plateau`, option (A) : clé `matrices_dir` et chemins des matrices ; pas de liens symboliques.
- **D5** fenêtre complète ±1,5 eV à toutes les grilles.
- **D6** variante alignée aussi en base de Bloch : ΔM(k′, k) = −C_N·V(k′) D_N(k − k′) V(k)† (approximation (i) ramenée en base de Bloch), appliquée à la carte
  fig_M_map (π, π*, partie L), à Re M^L à K (tableau des familles) et à tab:L_NL ; brut et aligné côte à côte ; test : grille grossière, bandes de la fenêtre
  gelée, (k, k) → −N_cells·C_N·𝕀 à 1e-12. max|M| et le test intensif restent sur M2 brut (hors k = k′). Changement ajouté à (b) pour `analyze_M.py` et
  `lnl_frobenius_all.py`.
- **D7** tab:tests_M : les deux lignes par famille seulement. **D8** fait rapporté, aucune variante. **D9** résonance 6×6 et 12×12 dans C.1.
- **D10** C14 redéfini : décalages ±rms du plateau autour du V_loc aligné (9×9 : ±9,05 meV, `a1_results.json` de R9) ; « C14 ancien » (V_sub) retiré.
- **D11** appelants hors chaîne laissés tels quels. **D12** test d'or non rejoué, valeur de R6 reprise. **D13** mêmes réglages que R6.
- **D14** confirmé : C_N = moyenne du plateau pour toutes les tailles, rms et max|écart| rapportés, aucun critère d'arrêt.

## GO 1 — soumission (2026-09-29)

- Pilote `r10_driver.py` (sous-commandes a0, a1, a2, a3, atables, b [--parts b0,b1,b2], btables) et lanceur `submit_r10.sh` (tâches a, b) écrits ; diffs
  archivés dans `submitted/<jobid>/`, md5 du pilote dans chaque journal. HEAD da2bc42 (commits EM depuis la phase 0, aucun dans `src/`, `scripts/`, `config/`,
  `results/`).
- Contrôles sur le nœud de connexion avant soumission (journal de test hors de ce répertoire) : P1 5×5 = R9 au bit (C_i −0,0575626062015234 eV ; un atome vrai
  n° 2, +176,1406 meV = audit ; Lu = fonction = publié, écart 0,0) ; couronnes 240² (forme entière x = i² + j² − ij) = `crowns_pi_240` de R9 B (14 couronnes,
  rangs, effectifs, énergies à 0,0) ; lecture des XML de B.2 = R5 C (E_F, gap, ε − E_F des états π et σ) à 0,0.
- Jobs : 22041340 `r10a` (a0 → a1 → a2 → a3 ; 8 cœurs, 64 Go, 1 h) ; 22041341 `r10b` (b0 → b1 → b2 ; 16 cœurs, 128 Go, 3 h ; afterany:22041340, pour qu'un
  échec de A n'arrête pas B ; sans A.1, B.2 laisse C_N vide pour 10 et 11 et B.1 prend C_9 dans R9).

## Rapport intermédiaire — A et B (GO 1, 2026-09-29)

| job | tâche | état | durée | MaxRSS |
|---|---|---|---|---|
| 22041340 `r10a` | a0 → a1 → a2 → a3 | COMPLETED | 4 min 05 | 13 Go |
| 22041341 `r10b` | b0 → b1 → b2 | COMPLETED | 21 min 05 | 19 Go |

Diffs archivés : `submitted/22041340/`, `submitted/22041341/` (pilote md5 0e37fa44c430). Après la soumission, retouches de présentation seulement (vrai signe
moins et virgule dans les tables, cadrage de `offset_profiles_13`, marqueurs de `eres_vs_grid`) ; tables et figures régénérées sans calcul (`atables`, `btables`) ;
ces changements du pilote apparaîtront dans le diff de la prochaine soumission. Chiffres bruts ; tables complètes dans `a/A1_tables.md`, `a/A2_tables.md`,
`a/A3_tables.md`, `b/B_tables.md` ; json `a/a1_results.json`, `a/a2_results.json`, `a/a3_results.json`, `b/b_results.json`.

### A.0 — porte (OK)

P1 sur 5, 6, 7, 8, 9, 12 : moyenne du plateau (i) = `C_i_eV` de R9 **au bit** pour les six tailles ; atomes vraiment les plus loin = audit (2 ; 2, 57 ; 81 ; 20 ;
21, 140 ; 28, 235), décalages à 0,0 meV de l'audit, moyennes +176,14 / −93,56 / −37,04 / +44,78 / −45,88 / −30,61 meV ; P1 à l'atome de Lu = `far_atom_alignment`
au bit et = Lu publié (R5 C) à ≤ 7e-12 meV.

### A.1 — C_N plateau (i) pour les 13 tailles (`a/A1_tables.md`, `fig/offset_profiles_13`)

meV ; plateau = atomes à distance vraie ≥ 0,75 r_max ; « un atome vrai » = atomes à r_max en vraie image (moyenne en cas d'égalité) ; Lu publié = atome de
`far_atom` (axe par axe).

| N | r_max (Å) | C_N = moyenne ; rms ; max\|écart\| ; n | 0,5 Å | un atome vrai : indices ; moyenne | Lu publié : atome ; valeur | C_N − Lu | ⟨ΔV⟩_3D |
|---|---|---|---|---|---|---|---|
| 5×5 | 7,12 | −57,56 ; 94,14 ; 233,70 ; 13 | −59,98 | 2 ; +176,14 | 49 ; −122,33 | +64,77 | +28,62 |
| 6×6 | 8,54 | −50,51 ; 24,56 ; 50,04 ; 26 | −52,33 | 2, 57 ; −93,56 | 71 ; −98,62 | +48,11 | +20,83 |
| 7×7 | 9,97 | −26,91 ; 18,53 ; 35,58 ; 31 | −27,81 | 81 ; −37,04 | 97 ; −38,54 | +11,62 | +14,55 |
| 8×8 | 11,39 | −14,41 ; 18,87 ; 59,19 ; 49 | −14,58 | 20 ; +44,78 | 1 ; −22,97 | +8,56 | +11,23 |
| 9×9 | 12,81 | −25,14 ; 9,05 ; 24,10 ; 53 | −25,29 | 21, 140 ; −45,88 | 161 ; −24,65 | −0,49 | +9,09 |
| **10×10** | 14,24 | −9,34 ; 6,13 ; 12,54 ; 70 | −9,46 | 23 ; −14,29 | 198 ; −8,22 | −1,12 | +7,17 |
| **11×11** | 15,66 | −9,55 ; 5,76 ; 17,57 ; 73 | −9,56 | 26 ; +8,02 | 241 ; −12,16 | +2,61 | +5,95 |
| 12×12 | 17,08 | −18,69 ; 5,42 ; 13,88 ; 95 | −19,08 | 28, 235 ; −30,61 | 287 ; −29,24 | +10,55 | +5,06 |
| **15×15** | 21,35 | −16,22 ; 2,76 ; 8,30 ; 149 | −16,22 | 65, 384 ; −23,34 | 449 ; −15,64 | −0,58 | +1,42 |
| **18×18** | 25,63 | −15,06 ; 1,74 ; 4,99 ; 212 | −15,18 | 115, 570 ; −19,32 | 1 ; −18,18 | +3,12 | −8,48 |
| **21×21** | 29,90 | −14,23 ; 1,03 ; 3,13 ; 275 | −14,33 | 133, 748 ; −16,87 | 881 ; −13,77 | −0,46 | +1,03 |
| **24×24** | 34,17 | −13,48 ; 0,72 ; 2,25 ; 383 | −13,56 | 201, 1000 ; −15,36 | 2 ; −13,12 | −0,36 | −62,74 |
| **27×27** | 38,44 | −13,01 ; 0,53 ; 1,62 ; 461 | −13,06 | 225, 1232 ; −14,36 | 1457 ; −12,77 | −0,24 | −99,91 |

Valeurs « un atome vrai » par atome (nouvelles tailles) : 15×15 −22,16 / −24,52 ; 18×18 −18,58 / −20,06 ; 21×21 −16,38 / −17,37 ; 24×24 −14,99 / −15,73 ;
27×27 −14,09 / −14,63 meV. Contrôles : P1 à l'atome de Lu = `far_atom_alignment` au bit pour les 13 tailles, écart au Lu publié ≤ 1,4e-11 meV. Lecture des Vks :
4 s (10×10) à 32 s (27×27) ; P1 complet ≤ 57 s par taille.

C_N (eV, pleine précision, `a/a1_results.json` → `C_N_eV`) pour le bloc `alignment` de la config (c) : 5×5 −0,0575626062015234 ; 6×6 −0,05050914517629934 ;
7×7 −0,02691429425434127 ; 8×8 −0,014410274533323343 ; 9×9 −0,02514371002267805 ; 10×10 −0,009341716029942201 ; 11×11 −0,009550651557225422 ;
12×12 −0,01868959310083182 ; 15×15 −0,01621854502926427 ; 18×18 −0,015062070057554946 ; 21×21 −0,014234270413329362 ; 24×24 −0,013481540317473475 ;
27×27 −0,0130057717522539. md5 de la copie versionnée `article/R10_plateau/a/a1_results.json` : 3d98bbd244eba8e1c1c645dcf675432b.

### A.2 — niveaux QE réalignés (`a/A2_tables.md`, `fig/levels_vs_invN`)

x = ε − C − E_D (eV) ; x_plateau = x_Lu + Lu − C_N ; w₂ inchangés. R5 C et R7 D1 identiques pour 6, 9, 12 (écart 0,0).

| N | Lu − C_N (meV) | π : x Lu → x plateau | σ : x Lu → x plateau (paire) |
|---|---|---|---|
| 5×5 | −64,77 | −0,3284 → −0,3931 | −0,1308 → −0,1955 |
| 6×6 | −48,11 | −1,0654 → −1,1135 | +0,1560 → +0,1078 |
| 7×7 | −11,62 | −0,4948 → −0,5064 | −0,2823 → −0,2939 |
| 8×8 | −8,56 | −0,3698 → −0,3784 | −0,1954 / −0,1953 → −0,2039 |
| 9×9 | +0,49 | −0,7373 → −0,7368 | +0,1014 → +0,1019 |
| 10×10 | +1,12 | −0,4224 → −0,4212 | −0,2383 → −0,2372 |
| 11×11 | −2,61 | −0,3319 → −0,3345 | −0,1641 → −0,1668 / −0,1667 |
| 12×12 | −10,55 | −0,5508 → −0,5614 | +0,1140 → +0,1034 |
| 15×15 | +0,58 | −0,4593 → −0,4588 | +0,1033 → +0,1039 |
| 18×18 | −3,12 | −0,3881 → −0,3912 | +0,1070 → +0,1038 |
| 21×21 | +0,46 | −0,3437 → −0,3432 | +0,1030 → +0,1035 |
| 24×24 | +0,36 | −0,3074 → −0,3071 | +0,1026 → +0,1030 |
| 27×27 | +0,24 | −0,2785 → −0,2783 | +0,1024 → +0,1027 |

### A.3 — colonne QE de R9 C.1 avec le plateau (`a/A3_tables.md`)

Modèles repliés de R9 non recalculés (« aligné » : C_9 = Lu de R9). Écart chaîne − QE (π), meV, contre QE Lu / QE plateau :

| N | QE π : Lu → plateau | tel quel : π ; écart Lu / plateau | aligné : π ; écart Lu / plateau |
|---|---|---|---|
| 6×6 | −1,0654 → −1,1135 | −1,0438 ; +21,6 / +69,7 | −1,0202 ; +45,2 / +93,3 |
| 9×9 | −0,7373 → −0,7368 | −0,6774 ; +59,8 / +59,3 | −0,6607 ; +76,5 / +76,0 |
| 12×12 | −0,5508 → −0,5614 | −0,5168 ; +34,0 / +44,5 | −0,5022 ; +48,6 / +59,2 |
| 15×15 | −0,4593 → −0,4588 | −0,4277 ; +31,6 / +31,1 | −0,4140 ; +45,4 / +44,8 |
| 18×18 | −0,3881 → −0,3912 | −0,3711 ; +16,9 / +20,0 | −0,3580 ; +30,1 / +33,2 |
| 21×21 | −0,3437 → −0,3432 | −0,3320 ; +11,6 / +11,2 | −0,3193 ; +24,4 / +23,9 |
| 24×24 | −0,3074 → −0,3071 | −0,3033 ; +4,1 / +3,8 | −0,2909 ; +16,5 / +16,2 |
| 27×27 | −0,2785 → −0,2783 | −0,2812 ; −2,7 / −2,9 | −0,2692 ; +9,4 / +9,1 |

σ : QE (paire) +0,1560 → +0,1078 (6), +0,1014 → +0,1019 (9), +0,1140 → +0,1034 (12), +0,1033 → +0,1039 (15), +0,1070 → +0,1038 (18), +0,1030 → +0,1035 (21),
+0,1026 → +0,1030 (24), +0,1024 → +0,1027 (27) ; état pair de la fenêtre des modèles −0,8308 / −0,8132 / −0,8125 (tel quel, 6 / 9 / 12…27) et −0,8061 / −0,7886 / −0,7879
(aligné) : écart chaîne − QE plateau −938,6 à −890,5 meV.

### B.0 — porte (OK) et porte 240² de B.1 (OK)

- B.0 (9×9 tel quel, N_k^int 300, 240²) : pic de la courbe Γ_T −0,1800 eV, pic de −Im T̄(K) −0,1775 eV ; courbe et T̄ = `resonance_9x9.npz` à 0,0 ; médiane des
  états 3 132,8723 meV (R6 3 132,8723).
- 240², N_k^int 900 contre R9 B (`b_results.json`, `b_results_plateau.json`, `b_curves_9x9*.npz`) : tel quel et plateau, même état E_res, écarts 0,0 sur E_res,
  les deux pics, la médiane et la courbe ; couronnes x = 9 (tel quel) et 7 (plateau), comme R9.
- g₀ « res » : 300² en 1,3 min, 900² en 8,3 min (16 fils).

### B.1 — E_res contre la grille de sortie (`b/B_tables.md`, `fig/eres_vs_grid`)

9×9, R_cut 3, η 0,02 eV, N_k^int 900, chaîne « res » (pas de g₀ 2,5 meV), états ±3 eV, E_res = argmax |Γ| des états à ±1,5 eV ; plateau : C_9 = −25,1437 meV.
Énergies en eV relatives à E_D (Wannier). « Voisines » : couronnes adjacentes en énergie moyenne (bande de l'état E_res, bande 3 partout).

| variante | grille | états | E_res | couronne x (rang) ; [min, max] | voisine en dessous : x ; E_res − moyenne | voisine au-dessus : x ; E_res − moyenne | pic Γ_T 5 ; 2,5 meV | \|E_res − pic\| 5 ; 2,5 (meV) |
|---|---|---|---|---|---|---|---|---|
| tel quel | 120² | 10 336 | −0,22689 | 3 (2) ; [−0,2393 ; −0,2269] | 4 ; +42,5 meV | 1 ; −92,5 meV | −0,2300 ; −0,2300 | 3,1 ; 3,1 |
| tel quel | 240² | 41 266 | −0,20184 | 9 (5) ; [−0,2019 ; −0,2018] | 12 ; +31,3 meV | 7 ; −23,9 meV | −0,1800 ; −0,1825 | 21,8 ; 19,3 |
| tel quel | 480² | 165 234 | −0,19993 | 37 (16) ; [−0,2092 ; −0,1999] | 39 ; +10,2 meV | 36 ; +1,9 meV | −0,2000 ; −0,2025 | 0,1 ; 2,6 |
| tel quel | 960² | 660 782 | −0,19050 | 133 (47) ; [−0,1974 ; −0,1905] | 139 ; +7,8 meV | 129 ; +0,5 meV | −0,2000 ; −0,1975 | 9,5 ; 7,0 |
| plateau | 120² | 10 336 | −0,13441 | 1 (1) ; [−0,1344 ; −0,1344] | 3 ; +98,7 meV | — | −0,1350 ; −0,1350 | 0,6 ; 0,6 |
| plateau | 240² | 41 266 | −0,17485 | 7 (4) ; [−0,1810 ; −0,1748] | 9 ; +27,0 meV | 4 ; −40,4 meV | −0,1800 ; −0,1800 | 5,2 ; 5,2 |
| plateau | 480² | 165 234 | −0,17485 | 28 (13) ; [−0,1810 ; −0,1748] | 31 ; +12,4 meV | 27 ; −0,2 meV | −0,1800 ; −0,1800 | 5,2 ; 5,2 |
| plateau | 960² | 660 782 | −0,17210 | 109 (39) ; [−0,1789 ; −0,1721] | 111 ; +5,0 meV | 108 ; +2,6 meV | −0,1800 ; −0,1800 | 7,9 ; 7,9 |

Médianes des états (meV) : tel quel 3 148,67 / 3 136,61 / 3 138,81 / 3 137,74 ; plateau 3 182,76 / 3 191,51 / 3 189,40 / 3 190,31 (120² → 960²). Durée par grille
(deux variantes) : 1 s (120²), 6 s (240²), 49 s (480²), 125 s (960²). Rappel (D8, arithmétique de 0.3) : à 480² et 960², des couronnes voisines sont séparées de
moins que le pas de g₀ (2,5 meV).

### B.2 — familles 5…12 (`b/B_tables.md`)

Contrôle contre R5 C (E_F, gap, ε − E_F des états π et σ) : 0,0 → OK. ε − E_F de la cellule avec lacune (non aligné) ; occupations XML par état (Γ, poids 2),
MV 0,01 Ry.

| N | famille | gap parfaite à Γ (eV) | ΔE_F (meV) | π : ε − E_F ; occupation | σ : ε − E_F ; occupation (chacun) | C_N (meV) | max\|M\| dense (eV) | pic Γ_T plateau 2,5 ; 5 meV (eV) |
|---|---|---|---|---|---|---|---|---|
| 5×5 | non-3m | 2,8093 | −235,8 | −0,0912 ; 0,8778 | +0,1064 ; 0,0611 | −57,56 | 23,850 | −0,2425 ; −0,2400 |
| 6×6 | 3m | 0,0000 | −20,5 | −1,1625 ; 1,0000 | +0,0589 ; 0,1622 | −50,51 | 25,148 | −0,1800 ; −0,1800 |
| 7×7 | non-3m | 3,0030 | −866,1 | −0,0969 ; 0,9016 | +0,1156 ; 0,0492 | −26,91 | 23,818 | −0,2450 ; −0,2450 |
| 8×8 | non-3m | 1,9761 | −583,3 | −0,0817 ; 0,8346 | +0,0928 ; 0,0827 | −14,41 | 24,311 | −0,2425 ; −0,2450 |
| 9×9 | 3m | 0,0000 | −2,4 | −0,7784 ; 1,0000 | +0,0602 ; 0,1584 | −25,14 | 25,348 | −0,1775 ; −0,1800 |
| 10×10 | non-3m | 2,0389 | −508,7 | −0,0858 ; 0,8534 | +0,0983 ; 0,0733 | −9,34 | — | — |
| 11×11 | non-3m | 1,5093 | −209,3 | −0,0788 ; 0,8205 | +0,0890 ; 0,0897 | −9,55 | — | — |
| 12×12 | 3m | 0,0000 | +2,5 | −0,6015 ; 1,0000 | +0,0633 ; 0,1497 | −18,69 | 25,721 | −0,1775 ; −0,1750 |

### Fichiers

- Répertoire de travail : `r10_driver.py`, `submit_r10.sh`, `submitted/`, `JOBID`, `r10_log.txt`, `slurm-r10-*` ; `a/` (json, tables, `profiles_<S>.npz` ≤ 60 ko) ;
  `b/` (json, table, `b1_curves.npz`) ; `fig/` (`offset_profiles_13`, `levels_vs_invN`, `eres_vs_grid`, pdf et png) ; `cache/` (g₀ 300² et 900² « res », 2 × 813 Mo,
  220 Mo sur le disque ; TEST, reconstructible en 10 min ; aucune suppression, manifeste en fin de campagne).
- Copie `article/R10_plateau/` : rapport, README, pilote, lanceur, tables, json, npz de `a/` et `b/`, figures (ni slurm, ni `cache/`).
- Rien dans `src/`, `scripts/`, `config/`, `results/`, `figures/`.

**STOP — GO 1 terminé le 2026-09-29.** Suite (étape G) : Greg écrit (a) et commit ; Code écrit (b) contre ces signatures (y compris D6 pour `analyze_M.py` et
`lnl_frobenius_all.py`), Greg relit et commit ; Greg applique (c) avec les C_N ci-dessus ; puis GO 2.

## (b) — changements de `scripts/` écrits par Code (2026-09-29 ; non commités)

Base : (a) commité par Greg (`c61c118`, HEAD `c1d9792`). Signatures relues dans le diff : `ws_images` n'accepte que `A_cols` (3, 3) ; `Mwr_to_Mwk(_pairs)(ws=)`
exige `len(ws["dist"]) == len(R)` (étiquettes brutes avec `R_center = R_d`, ou recentrées avec `R_center = 0`) ; `defect_mwr` soustrait C_N sur M_W(R, R) pour tous
les w et toutes les mailles de la boîte (étiquettes brutes), puis `recenter_mwr`, sans conversion d'unités. (c) n'est pas appliqué : les scripts importent
`matrices_dir` et `alignment_C` de `config.py` tels que proposés en 0.5 (c).

`git diff --stat` : 11 fichiers de `scripts/`, +182 −99. Rien dans `src/`, `config/`, `tests/`, `results/`, `figures/`.

| fichier | changement |
|---|---|
| `compute_spectral_wannier.py` | `lt.defect_mwr(M, U, U_dis, k_coarse, MP, n_box=N, C_N=alignment_C(cfg, S))` (dense et grossier) à la place de Mbk_to_Mwk → Mwk_to_Mwr → recenter_mwr ; M_ed grossier lu dans `matrices_dir` ; ligne `[align]` dans le journal |
| `rcut_resigma.py`, `resonance_criteria.py` | idem (`defect_mwr`, C_N de la config, ligne de journal) |
| `resonance_metrics.py` | `vloc_from` par `defect_mwr` ; D10 : V_sub retiré (ML_diag_mean, lecture de M_L, matrice identité de la taille de M, t_S, G_S, GS_e, ligne « alignment ») ; clés npz retirées `Gamma_T_noshift`, `G_S`, `ML_diag_mean` (aucun lecteur, `grep`) ; clé ajoutée `C_N_eV` ; `--shift-L-meV` inchangé dans son principe (constantes en meV autour du V_loc aligné), aide : « ±rms du plateau (9x9 : 9.05,-9.05) » |
| `m_rcut_convergence.py` | `defect_mwr` ; `ws = ws_images(Rn, 0, MP, A_uc)` (cellule du .save dense, bohr) ; Pm, Pp = `ws_phase(kf, ws, nR, ∓1)` ; masques R_cut sur Rn inchangés ; ligne `[ws]` (étiquettes à égalité, nombre d'images) |
| `mwr_locality_coarse_vs_dense.py` | `defect_mwr` (dense et grossier) ; abscisse `ws_images(Rn, 0, MP, A_cols)["dist"]` (A de la l. 19 en colonnes, a₃ sans effet, R_z = 0) ; M_ed dans `matrices_dir` |
| `analyze_Ved.py` | r = `al.true_min_image_dist` sur (i/n₁, j/n₂, z_vac), centre s_vac, cellule en Å (deux boucles) ; X, Y, S1, S2 devenus inutiles retirés (la carte utilise IC, JC, inchangés) |
| `analyze_M.py` | chemins des matrices par `MAT = matrices_dir(cfg)` (l. 21, M_ed, nb16, coarsecheck, M_L) ; `corners` et `map_corners` supprimés ; fermeture par `defect_mwr(…, C_N=0.0)` ; D7 : deux lignes par famille (tailles lues dans `cfg["families"]`, libellé construit sur les tailles présentes) ; D6 : carte alignée (`map_Vpi_aligned`, `map_Vpistar_aligned`, `map_Lpar_aligned`, `map_Npar_aligned`, `map_Labs_aligned`, `map_C_N_eV`), Re M^L à K aligné (`lnl_<S>_ReL_aligned`, `lnl_<S>_C_N_eV`), portes D6 (ci-dessous, clés `align_gate_*`) ; clés brutes inchangées |
| `lnl_frobenius_all.py` | D6 : colonnes ajoutées en fin de ligne `C_N_eV`, `ratio_full_aligned`, `mean_fL_aligned_eV`, `ratio_diag_aligned`, `ratio_min_aligned`, `ratio_max_aligned`, `ratio_median_pairs_aligned` ; colonnes brutes inchangées ; contrôle refusant de D_N contre la somme explicite sur la boîte (trois kets, 1e-12) |
| `level2_families.py` | colonne `ReML_K_aligned_eV` ajoutée en fin de ligne (lue dans `M_analysis.npz`) — hors des deux fichiers nommés par D6, voir point 3 ci-dessous |
| `make_figures_memoire.py` | fig_Ved (c) : `axhline(alignment_C(cfg, S))` pointillé, couleur de la taille, pour chaque taille tracée |
| lignes `RES = results_dir(…)` | commentaire « R6 : results/M2 » remplacé par « R10 : produits (results/M2_plateau) » dans les fichiers touchés |

Variante alignée en base de Bloch (D6), forme fermée utilisée dans `analyze_M.py` et `lnl_frobenius_all.py` :
ΔM[m, k′, n, k] = −C_N·D_N(k − k′)·Σ_w V(k′)[m, w]·V(k)[n, w]*, D_N(q) = Σ_{R ∈ [0, N)²} e^{2πi q·R} calculé sur les indices MP entiers ; ΔM s'ajoute à M et à M^L,
pas à M^NL (contrôle : M_dense = M_L_dense + M_NL_dense en v2, 8,8e-17 sur une colonne 9×9).

Portes D6 dans `analyze_M.py` (taille de référence, refus au-delà du seuil) :
(i) forme fermée contre la chaîne de Wannier (`defect_mwr` sur M = 0 en jauge identité, puis `Mwr_to_Mwk_pairs`), 729 k′ contre K, seuil 1e-12 ;
(ii) sous-grille grossière 9×9 de la grille dense 27×27, jauge de Wannier, toutes les paires (k′, k) : ΔM_W = −N_cells·C_N·δ_kk′·𝕀₅, seuil 1e-12 ;
(iii) bandes de la fenêtre gelée ([−25 ; −1,74] eV, 4 ou 5 par k), (k, k) de la sous-grille, base de Bloch : ΔM = −N_cells·C_N·𝕀, seuil 1e-9 (point 1).

### Essais (nœud de connexion, 16 Go ; `config.py` de substitution conforme à 0.5 (c) dans le scratchpad, `results_dir` → scratchpad ; rien écrit dans `results/`)

- pyflakes sur les 11 fichiers, avant et après : aucun message nouveau (les messages préexistants, imports inutilisés, restent identiques).
- `analyze_M.py`, en-tête + section 1 + bloc D6 exécutés tels quels (9×9) : portes (i) 2,9e-15, (ii) 1,3e-15, (iii) 1,4e-10 (max|V V† − 𝕀| sur les bandes gelées
  1,37e-10) ; clés brutes de la carte (`map_Vpi`, `map_Vpistar`, M^L_∥, M^NL_∥, |M^L|, |M^NL|, kx) = `M_analysis.npz` de M2 à 0,0 ; `map_corners` absent.
  Carte alignée (eV Å²) : Ṽ_π à K 107,22 → 114,89 ; M^L_∥ à K 84,39 → 92,09 ; ⟨M^L_∥⟩ 89,53 → 89,63 ; max|ΔM| (ligne π, paire à K) 10,72 = A_cell·N_cells·|C_9|.
- Re M^L à K, lignes de la section 3 (eV) : 9×9 +11,3051 → +13,3417 ; 5×5 +10,5343 → +11,9733 ; 12×12 +11,2220 → +13,9133 ; décalages +2,0366 / +1,4391 / +2,6913
  = −N_cells·C_N (paire π à K dans la fenêtre gelée pour les trois tailles ; 5×5 : point de grille le plus proche de K à 0,0555 Å⁻¹).
- D7 sur les `scale_<S>_dense` de M2 : 3m (N = 6, 9, 12) 2,3e-2 ; non-3m (N = 5, 7, 8) 2,1e-2 (R9 : 2,26e-2 et 2,05e-2).
- `lnl_frobenius_all.py 5x5,9x9` (28 s) : les 13 colonnes brutes = `results/M2/lnl_frobenius.csv` (chaînes identiques) ; contrôle D_N passé ; aligné :
  ⟨‖M^NL‖_F⟩/⟨‖M^L‖_F⟩ 5×5 0,258 → 0,257 (diagonale k′ = k 0,287 → 0,270) ; 9×9 0,252 → 0,252 (0,278 → 0,255).
- `mwr_locality_coarse_vs_dense.py` restreint à 5×5 : p_z–p_z sur site de la lacune +31,0023 → +31,0599 eV (dense), +1,1235 → +1,1811 eV (grossier), décalage
  +57,563 meV = −C_5 ; poids hors site identiques à 0,0 ; abscisse max 20,785 → 13,892 a (dense, D = 25, borne D/√3 = 14,43 a) ; 3,464 → 2,646 a (grossier).
- `m_rcut_convergence.py --size 5x5 --nf 24 --rcuts 0,1,2,3` (15 s) : 24 étiquettes à égalité, 649 images pour 625 étiquettes ; max|ΔM|/max|M| 2,898e-1 / 1,583e-1 /
  6,161e-2 / 2,553e-2 (aucune référence pour cette taille : exécution seulement ; la porte chiffrée est celle de C.2 sur 9×9 et 12×12).
- Lignes de `rcut_resigma.py` jusqu'au recentrage, 5×5 et 9×9 : V_loc (R_cut 3) aligné = V_loc(C = 0) − C_N·diag(in_box) **au bit** ; R_d et étiquettes recentrées
  identiques à la chaîne sans alignement ; 5×5 : 25 des 29 mailles du R_cut 3 dans la boîte 5×5 (les 4 autres, dans la zone de remplissage de la grille 25×25,
  ne reçoivent pas C_N : approximation (i)) ; 9×9 : 29 sur 29.
- `analyze_Ved.py 9x9` : anneau 9,625 Å 534 points, 11,075 Å 518 (audit : 534, 518) ; anneaux r < 0,433 a_sc = 9,609 Å (192) : 1,2e-14 (profil), 4,3e-14 (profil
  masqué) de M2, pas 0,0 ; les 30 anneaux au-delà changent ; les 24 autres clés 9×9 = M2 à 0,0. Cause de l'écart de 1e-14 : r nouveau − ancien ≤ 5,3e-15, les
  49 543 points restent dans le même anneau, mais `np.histogram` pondéré somme par différences de sommes cumulées dans l'ordre trié des r ; une somme par anneau
  dans l'ordre du tableau (`bincount`) donne 0,0.
- `make_figures_memoire.py --outdir <scratchpad>` sur des copies des npz de M2 : six figures écrites ; fig_Ved (c) avec les six traits C_N (7×7 −26,91 et 9×9 −25,14 meV
  presque superposés).
- Non exécutés ici (mémoire du nœud de connexion ou durée) : `compute_spectral_wannier.py`, `resonance_metrics.py`, `resonance_criteria.py`, `analyze_M.py` complet
  (sections 2, 4, 5), `mwr_locality` sur les six tailles, `m_rcut_convergence.py` 9×9 et 12×12, `level2_families.py` sur des produits alignés. Ils passent en C.0/C.1
  derrière leurs portes.

### Écarts au plan et points à trancher

1. **Test D6 « (k, k) → −N_cells·C_N·𝕀 à 1e-12 »** : impossible à la lettre en base de Bloch ; V V† = 𝕀 sur les bandes gelées à 1,37e-10 seulement (`wannier_u.mat`
   et `wannier_u_dis.mat` écrits avec 10 décimales ; U†U − 𝕀 1,7e-10, U_dis†U_dis − 𝕀 1,4e-10). Écrit : portes (i) et (ii) à 1e-12 (exactes, jauge de Wannier), porte
   (iii) en base de Bloch à **1e-9** (seuil de Code, à confirmer ou à remplacer).
2. **Grille grossière du test** : pas de wannierisation grossière de la 9×9 (`wannier/9x9` absent ; `wannier/` : 5×5, 7×7, 8×8 grossiers) ; le test utilise la
   sous-grille 9×9 de la grille dense 27×27 (même algèbre : D_9(k − k′) = 81·δ_kk′ sur cette sous-grille).
3. **`level2_families.py`** : colonne alignée ajoutée pour que « Re M^L à K (tableau des familles) » ait brut et aligné côte à côte ; fichier hors des deux nommés par
   D6 ; à garder ou à retirer.
4. **fig_M_map** : `make_figures.py` et `make_figures_memoire.py` inchangés (carte brute) ; la variante alignée est dans `M_analysis.npz` ; la planche côte à côte
   brut/aligné est prévue au pilote (C.4). Des panneaux alignés dans la figure du mémoire seraient un changement séparé.
5. **fig_Ved (c)** : les traits C_N n'ont pas d'entrée de légende (légende à six colonnes inchangée) ; à dire dans la légende de la figure.
6. **Contrôle C.3 des anneaux de `analyze_Ved`** : « anneaux < 0,433 a_sc redonnés à 0,0 » donne 1,2e-14 et 4,3e-14 (cause ci-dessus) ; proposé pour C.3 : appartenance
   aux anneaux identique (indices) et valeurs à 1e-13, ou recalcul par `bincount` à 0,0.
7. **Ordre des commits** : les scripts de (b) importent `matrices_dir` et `alignment_C` ; entre un commit de (b) seul et (c), ils échouent à l'import. Commit de (b)
   et (c) ensemble, ou (c) d'abord.
8. **Scripts hors chaîne (D11, laissés)** : 14 scripts lisent des matrices par `f"{RES}/M_…"` et chercheront dans `results/M2_plateau` après (c) : `check_ML_coarse_kernel`,
   `check_onsite_and_NL`, `check_M_dense_vs_coarse`, `check_M_dense_nb20_vs_nb16`, `compute_tmatrix`, `compute_spectral`, `migrate_M_norm`, `_old_vs_new_7x7`,
   `_eta_scan`, `_normtest`, `tag_vacancy_sublattice`, `_mcheck`, `test_local_tmatrix_real`, `validate_ML_grid_7x7`. Aucun n'est utilisé en C.

Ajouts aux contraintes du lanceur de C (0.6) : C14 appelé avec `--shift-L-meV 9.05,-9.05` (D10) ; les lignes de la porte A.2 que R6 ajoutait à
`M_tests_summary.csv` (`article/R6_production_corrigee/etape3/r6_tests_gate_row.py`, `results_dir` de la config) seront reprises dans le pilote pour
`results/M2_plateau`.

Fichiers de l'essai (scratchpad, temporaires) : config de substitution, bancs d'essai, sorties ; hors dépôt et hors répertoire de campagne.

**STOP — (b) écrit le 2026-09-29, non commité.** Suite : relecture et commit par Greg ; (c) par Greg (points 1, 3 et 7 ci-dessus) ; puis GO 2.

## (c) — config appliquée par Code à la demande de Greg (2026-09-29 ; non commitée)

(b) commité par Greg (`ea91d61`, « R10 b »). Greg demande à Code d'appliquer (c) et de vérifier trois points ; Greg commit et pousse.

Fichiers : `config/production.json` (`results_dir` → `results/M2_plateau` ; `matrices_dir` = `results/M2` ; `results_dir_frozen` = [`results/M`, `results/M2`] ;
bloc `alignment` : version « plateau (i), R10, 2026-09-29 », méthode, `C_N_eV` des 13 tailles, `source` = `article/R10_plateau/a/a1_results.json`, `source_md5`,
`labels_offgrid` = `wigner_seitz` ; reste du fichier inchangé à l'octet, même mise en forme) ; `src/electron_defect_interaction/config.py` (clés obligatoires
`matrices_dir`, `alignment` ; `matrices_dir()` ; `alignment_C()`, `KeyError` si la taille manque ; `dense_paths(…)["mfile"]` → `matrices_dir` ; docstring de
`results_dir` ; ligne `verbose` avec `matrices=` et `C_N=`) ; `.gitignore` (mêmes exceptions pour `results/M2_plateau/` que pour `results/M2/`).

Vérifications (chargeur réel `load_production`) :
1. C_N de 5×5, 6×6, 7×7, 8×8, 9×9, 12×12 = `C_i_eV` de R9 (`article/R9_controles/a/a1_results.json`) au bit ; arrondis −57,56 / −50,51 / −26,91 / −14,41 / −25,14 /
   −18,69 meV ; écrits en eV (|C_N| < 0,1). Le `C_retenu_eV` de R9 (autre valeur pour 5×5, 6×6, 7×7, 8×8, 12×12) n'est pas utilisé.
2. `source_md5` 3d98bbd244eba8e1c1c645dcf675432b = md5 de `article/R10_plateau/a/a1_results.json` = md5 de `R10_plateau/a/a1_results.json` ; les 13 C_N de la config =
   `C_N_eV` du fichier au bit.
3. `matrices_dir` → `results/M2`, `results_dir` → `results/M2_plateau` ; `dense_paths("9x9")["mfile"]` = `results/M2/M_dense_9x9.npy` (existe) ; taille absente
   (13×13) → `KeyError`.

Contrôles : import de `config` des 11 scripts de (b) OK ; `git check-ignore` : `M_analysis.npz`, `specwd_*_prod.npz`, csv, README de `results/M2_plateau/` versionnés,
`logs/` et autres fichiers ignorés ; `pytest tests` : 192 passed (7 min 53 s). `results/M2_plateau/` n'existe pas encore (créé au GO 2).

Points de (b) toujours ouverts (section « (b) ») : seuil 1e-9 de la porte D6 (iii), colonne de `level2_families.py`, contrôle C.3 des anneaux (1e-14).

**STOP — (c) appliquée le 2026-09-29, non commitée.** Suite : commit et push par Greg ; GO 2.
