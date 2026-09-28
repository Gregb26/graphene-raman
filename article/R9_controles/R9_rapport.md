# R9 — Contrôles avant le chapitre 4 — rapport de campagne

Statut **TEST** (post-traitement seul). Prompt R9 (Greg, 2026-09-27). Ordre : phase 0 → STOP → étape G (Greg) → (GO) A → B → C → D →
rapport → STOP. Dépôt `graphene-raman`, HEAD e4122ab (« R6 clos », arbre propre) ; racine `$GRAPHENE_RAMAN`. Répertoire de travail
`graphene/qe/defects/R9_controles/`, copie versionnée `article/R9_controles/`. Chiffres bruts, sans interprétation.

## Phase 0 (2026-09-27 ; rien n'est calculé)

Opérations faites : lectures (rapports R4, R5, R6, R7 ; `src/`, `scripts/`, pilotes `r4_driver.py`, `r5_driver.py`, `r6_d4.py`, `r7_driver.py` ;
`config/production.json` ; journaux `results/M2/logs/specwd_*` ; clés de `results/M2/resonance_9x9_shiftL.npz` et de
`R7_tailles_3m/d1_8pts/window_*.npz`) et `sacct` des jobs M^L denses de septembre. Aucun job, aucune écriture hors de ce répertoire et de
sa copie. Les lignes marquées « arithmétique » sont des produits de chiffres déjà publiés (R4–R7) ou du cône de Dirac (ħv_F = 5,459 eV Å,
R5 C ; a = 2,4659 Å ; |b| = 4π/(√3 a) = 2,9422 Å⁻¹ ; A_cell = 5,266 Å²), écrites pour dimensionner le plan ; ce ne sont pas des mesures.

### 0.1 Alignement : fonction, définition, signe, valeurs par taille

**Fonction** : `far_atom_alignment(V_d, V_p, x_red_d, x_red_p, A_cols, radii)` de `src/electron_defect_interaction/defects/alignment.py`
(R4, commit 75c8656), appelée par `r4_driver.py` (d1, d5), `r5_driver.py` (c) et `r7_driver.py` (d1) avec `radii = (0,5 ; 1,0)` Å ;
valeur de travail : 1,0 Å.

Définition exacte :
1. Potentiels : `graphene/qe/defects/super_cell/NxN/{defective,pristine}/Vks_NxN_{d,p}` (pp.x `plot_num = 1`, potentiel local total avec XC,
   filplot, Ry), lus par `qe_io.get_pot(…, subtract_mean=False, to_hartree=True)`, × HA2EV, `transpose(2, 1, 0)` → V[ix, iy, iz] en eV
   (`load_pot_eV` de r4/r5). Positions `qe_io.get_x_red`, cellule `get_A_volume` × Bohr → Å.
2. Site de la lacune (`vacancy_site`) : l'atome de la parfaite dont la distance (image minimale) au plus proche atome de la cellule avec
   lacune est maximale.
3. « Loin du défaut » = **un seul atome** (`far_atom`) : l'atome de la cellule avec lacune le plus éloigné (image minimale) du site de la
   lacune (`argmax`, premier indice en cas d'égalité) ; dans la parfaite, l'atome le plus proche de cette position (même position pour une
   géométrie non relaxée). 9×9 : atome 161 (1-based), 18,51 Å.
4. `sphere_average` : moyenne arithmétique de V sur les points de grille à distance < r du centre (image minimale, sphère 3D, tous les plans z) ;
   9×9 : 8 611 points à 1,0 Å, 1 039 à 0,5 Å.
5. Décalage « Lu » = ⟨V_d⟩_sphère − ⟨V_p⟩_sphère. **Signe** : c'est la valeur de ΔV = V_d − V_p autour de l'atome loin ; R4/R5/R7 le
   retranchent aux valeurs propres de la cellule avec lacune (ε_QE − Lu). Identité R5 A.2 / R6 1.2 (30 états QE) : constante s = +24,66 meV = −Lu.
6. À part : ⟨V_d⟩ − ⟨V_p⟩ 3D = composante G = 0 de ΔV sur toute la cellule (vide compris) (`ved_metrics`, `mean3d`).

Aucun profil de ΔV autour de **tous** les atomes n'a été calculé jusqu'ici (un atome ; ligne a₁ et plans frontière en R4 D5). Aucune
soustraction dans la production de M (`compute_M_dense_stages.py`, `subtract_mean=False`) ; le seul retrait existant est la variante
« C14 ancien » de `resonance_metrics.py` (moyenne diagonale de M^L, 5 427,3 meV en v2 9×9).

**Convention retenue pour R9** : C_N ≡ valeur de ΔV loin du défaut (même signe que Lu) ; ΔV_aligné = ΔV − C_N. 9×9 avec Lu :
C_9 = −24,65 meV → ΔV_aligné = ΔV + 24,65 meV. La variante C14 « C = +25 meV » de R6 3.5 (V_loc + 25 meV·𝕀 sur les 29 mailles de R_cut 3)
est donc, au signe près et à 0,35 meV près, la variante alignée 9×9 dans l'approximation (i) de 0.2 : médiane 3 187,96 meV (tel quel 3 132,87),
E_res −0,175, pic de la courbe Γ_T −0,180, écart relatif maximal de Γ_T 19,2 % (R6 3.5). La valeur de C_N utilisée en A.3–D sera celle de
A.1 (décision 3 ci-dessous).

Valeurs déjà calculées (meV ; « bord » et « plans » : R4 D5, ΔV le long de a₁ de la lacune à la demi-boîte et max|ΔV − Lu| sur les plans
frontière dans le plan des atomes) :

| N | source | atome loin (1-based) / distance (Å) | Lu 1,0 Å | Lu 0,5 Å | ⟨ΔV⟩_3D | bord ligne a₁, aligné | plans frontière, plan des atomes, aligné |
|---|---|---|---|---|---|---|---|
| 5 | R4 D5 = R5 C | 49 / 9,97 | −122,33 | −131,64 | +28,62 | −19,3 | 135,0 |
| 6 | idem | 71 / 12,81 | −98,62 | −125,83 | +20,83 | −42,5 | 78,3 |
| 7 | idem | 97 / 14,24 | −38,54 | −37,93 | +14,55 | −2,2 | 77,5 |
| 8 | idem | 1 / 17,08 | −22,97 | −24,62 | +11,23 | −2,9 | 26,5 |
| 9 | R4 D1, D5 | 161 / 18,51 | −24,65 | −21,59 | +9,09 | +0,7 | 17,0 |
| 10 | R4 D5 = R5 C | 198 / 19,93 | −8,22 | −7,80 | +7,17 | −3,5 | 10,4 |
| 11 | idem | 241 / 22,78 | −12,16 | −12,09 | +5,95 | +0,7 | 21,2 |
| 12 | idem | 287 / 25,63 | −29,24 | −36,07 | +5,06 | −10,5 | 16,8 |
| 15 | R7 D1 | — | −15,64 | −14,54 | +1,42 | — | — |
| 18 | R7 D1 | — | −18,18 | −20,34 | −8,48 | — | — |
| 21 | R7 D1 (8 points) | — | −13,77 | −13,44 | +1,03 | — | — |
| 24 | idem | — | −13,12 | −12,94 | −62,74 | — | — |
| 27 | idem | — | −12,77 | −12,69 | −99,91 | — | — |

Tailles de M2 (denses et grossiers) : 5, 6, 7, 8, 9, 12 (10×10 grossier seulement, 11×11 exclu en R6).

### 0.2 Le terme constant dans M2

**M grossiers (k commensurables).** Noyau `compute_ML_R_mpi` : M^L_{n′k′,nk} = dV Σ_{r∈SC} ψ*_{n′k′} ΔV ψ_{nk} avec ψ = u e^{2πi k·r}/√Ω_uc,
dV = Ω_uc/N_uc. Pour ΔV = C : M^L[C] = C·(1/N_uc) Σ_{R∈SC} e^{2πi(k−k′)·R} Σ_{p∈maille} e^{2πi(k−k′)·p} u*_{n′k′}(p) u_{nk}(p). Pour k, k′ sur la
grille N×N, Σ_R = N_cells δ_kk′ et (1/N_uc) Σ_p u*_{n′k} u_{nk} = Σ_G C*C = δ_nn′ (Parseval discret) : **M^L[C] = N_cells·C·𝕀 exactement**
(7×7 compris : le rééchantillonnage de Fourier sur 217 conserve une constante). M^NL n'est pas touché. Donc M_aligné = M2 − N_cells·C_N·𝕀
(C_N en Hartree dans les fichiers). En base de Wannier grossière (repliée, période N) : C_N·𝕀₅ sur **toutes** les mailles ; dans (a1)
(H = diag ε + M/N_cells) : décalage rigide C_N (R5 A.5 l'a vérifié pour 76,58 meV).
Arithmétique, N_cells·|C_N| avec Lu 1,0 Å : 3,06 / 3,55 / 1,89 / 1,47 / 2,00 / 4,21 eV pour N = 5 / 6 / 7 / 8 / 9 / 12, à comparer à
max|M| (bandes 1–16) 23–27 eV (R6 3.2).

**M denses (ΔV plongé dans la cellule zéro-paddée p·N).** Le noyau `compute_ML_R_mpi_shared` intègre sur la seule grille de la super-cellule
N×N avec k, k′ sur la grille D×D (D = pN) : ΔV est nul hors de la **boîte** des mailles R = (i, j), 0 ≤ i, j < N, de la grille pp.x. Pour
ΔV = C·1_boîte :

  M^L[C·1_boîte]_{n′k′,nk} = C · D_N(k − k′) · S_{n′k′,nk}(k − k′),
  D_N(q) = Σ_{i,j=0}^{N−1} e^{2πi(q₁i + q₂j)} = e^{iπ(N−1)(q₁+q₂)} Π_a sin(πNq_a)/sin(πq_a),
  S_{n′k′,nk}(q) = (1/N_uc) Σ_{p∈maille} e^{2πi q·p} u*_{n′k′}(p) u_{nk}(p)   (p en coordonnées réduites de la grille de maille).

Sur la grille D×D : q = 0 → N_cells·δ_nn′ (comme le grossier) ; q ∈ (1/N)ℤ² non nul (paires de k coïncidents modulo la grille N×N) → 0 ;
tout autre q → non nul, enveloppe en sinus cardinal de largeur ~1/N autour de k′ = k : c'est le **facteur de forme de la boîte N×N**. En base
de Wannier : C·F(R, R′), F(R, R′) = ⟨w_R|1_boîte|w_R′⟩ (exact dans les 20 bandes, qui contiennent les w) : ≈ δ_RR′ pour une maille au cœur de la
boîte, → 0 hors de la boîte, fractionnaire pour les mailles du bord (fonctions de Wannier à cheval sur le bord).

Boîte en étiquettes recentrées Rn = R − R_d (R_d de `recenter_mwr`, journaux `specwd_*` de R6) et sites de l'amas V_loc hors de la boîte :

| N | R_d | boîte (Rn₁, Rn₂) | sites hors boîte, R_cut 3 (29 sites) | idem R_cut 4 (49 sites) | mailles du bord dans l'amas R_cut 3 |
|---|---|---|---|---|---|
| 5 | [2, 2] | [−2, 2]² | 4 | 24 | couronne \|Rn\| = 2 (bord des deux côtés) |
| 6 | [2, 2] | [−2, 3]² | 2 | 14 | Rn₁ ou Rn₂ = −2 ou 3 |
| 7 | [3, 3] | [−3, 3]² | 0 | 4 | (±3, 0), (0, ±3) |
| 8 | [4, 4] | [−4, 3]² | 0 | 2 | (3, 0), (0, 3) |
| 9 | [4, 4] | [−4, 4]² | 0 | 0 | aucune (bord à \|Rn_a\| = 4, dans R_cut 4) |
| 12 | [5, 5] | [−5, 6]² | 0 | 0 | aucune |

**Proposition (i)**, en base de Wannier, coût nul : M_W(R, R) − C_N·𝕀₅ pour les mailles R de la boîte, rien hors de la boîte ni hors diagonale ;
sur l'amas : V_loc,aligné = V_loc − C_N·P_boîte (P_boîte = 𝕀 sur les 5 orbitales des sites de l'amas situés dans la boîte).

**Proposition (ii)**, exacte : M^L[1_boîte] par le noyau partagé corrigé (`compute_ML_R_mpi_shared`, 09280bb) avec ΔV = 1 sur la grille N×N,
**sans modifier `src/`** : le pilote passe au noyau un module `io` enveloppe de `qe_io` dont `get_pot` renvoie des 1 (chemin « défectueux ») et
des 0 (chemin « parfait ») sur la grille du vrai fichier (ΔV = 1 − 0 ; l'argument `io=` du noyau est utilisé tel quel) ; autre voie : un
argument `Ved=None` ajouté au noyau (P5, changement de `src/`). Puis M_aligné,exact = M2 − C_N·M^L[1_boîte] (20 bandes, tous les k) et
F_W = V†M^L[1_boîte]V en base de Wannier (V_loc,exact = V_loc − C_N·P F_W P). Coût par taille = étape `ml` de production (1 nœud exclusif,
32 rangs × 6 fils, mesures `sacct` de septembre) :

| N | D | étape ml (python) | job complet | MaxRSS étape | fichier (20 bandes) |
|---|---|---|---|---|---|
| 5 | 25 | 14 min | 23 min | 42 Go | 2,50 Go |
| 6 | 24 | 16 min | 28 min | 39 Go | 2,12 Go |
| 7 | 28 | 41 min | 1 h 03 | 57 Go | 3,93 Go |
| 8 | 32 | 1 h 30 | 2 h 06 | 74 Go | 6,71 Go |
| 9 | 27 | 53 min | 1 h 25 | 50 Go | 3,40 Go |
| 12 | 24 | 59 min | 1 h 41 | 39 Go | 2,12 Go |

Six tailles : ≈ 4 h 30 de nœud exclusif (étapes ml), 20,8 Go. Proposé : 9×9 seulement (validation de (i), A.2), les autres sur décision.

Arithmétique (ordre de grandeur du seul terme constant dans tab:rcut_M, hypothèses : C_N = Lu 1,0 Å sur chaque diagonale de M_W dans la boîte,
maximum de |ΔM| en k′ = k, où le terme constant de la paire π reconstruite vaut C_N × nombre de mailles gardées) : R_cut 3, 9×9 :
(81 − 29) × 24,65 meV / 25,44 eV = 5,0e-2 ; 12×12 : (144 − 29) × 29,24 meV / 25,94 eV = 1,30e-1. Valeurs mesurées en R6 : 6,68e-2 / 9,02e-2.

### 0.3 Chaîne repliée

**Existant** :
- `supercell_fold.fold_hwr_to_supercell(Hwr, Rw, ndegen, N)` : H à Γ de la N×N pour tout N (classes R mod N, ndegen divisé) ; porte 1
  `wannier_gate` (R6 : 4,9e-14 eV à N = 9).
- `supercell_fold.fold_mwr_to_supercell(Mwr, R_mwr, N, R_local=None)` : somme sur **toutes les paires d'images** (R ≡ c, R′ ≡ c′ mod N) =
  élément de matrice à Γ de la somme périodique des copies ; correct aussi quand les images se recouvrent (l'amas R_cut 3 s'étend sur 7 mailles :
  recouvrement pour N ≤ 6). V_loc (disposition L·nw + w) s'y passe directement : `V.reshape(nL, nw, nL, nw).transpose(1, 0, 3, 2)` avec
  `R_mwr = R_loc` : **aucune fonction nouvelle** pour replier H(R) (`wannier/27x27`) + V_loc(R_cut 3) dans une N×N quelconque.
- `r6_d4.py cmd_d4`, variante (c-3) : H_SC = fold(H(R)) + fold(M_W, R_local = R_loc) à N = 9 ; `r4.spectrum_blocks` (parité : 3 sp² pairs,
  2 p_z impairs) ; E_D = Wannier à K (= quadruplet de la parfaite repliée pour N = 3m) ; w₂ par `w_from_r` (coefficients C du `.save` dense et
  `folded_density_2d`, possible seulement si la grille N×N ⊂ 27×27 : N = 9, 27) ; `wsite` (site + voisins + sp²) ; `alignment.rigid_shift_fit`
  (états QE < E_D − 4 eV).
- QE : `R7_tailles_3m/d1_8pts/window_NxN.npz`, N = 6…27 : spectre complet `e_d_all`, états de la fenêtre [−3, +1] (ε alignée Lu, ε − E_D,
  parité, w₂, w₁), E_D, Lu.
- Voisins de la lacune (V_loc 9×9, lacune A en Rn = 0) : p_z(B) (WF 4) des mailles (0,0,0), (−1,0,0), (0,−1,0) (`NN_CELLS` de R4) ;
  indice `cell_of(R mod N)·5 + 4` dans le modèle replié, `L·5 + 4` dans l'amas.

**Manquant** :
- LDOS sur les p_z des trois voisins : **n'existe nulle part** (ni dans `src/`, ni dans `scripts/`, ni dans les pilotes R4–R7 ; `projwfc_io`
  lit seulement les pdos de projwfc, R1 nspin 2). À ajouter : LDOS du modèle replié à partir des paires propres (lorentziennes η) (P4) et LDOS de
  la limite diluée G = g₀ + g₀ t g₀ sur des orbitales de l'amas (P3).
- Côté QE, aucune projection p_z pour R7 (projwfc = calcul QE, exclu). Substitut possible : LDOS pondérée par w₂ (disque de 2 Å, qui contient
  les trois voisins) à partir des npz de R7 — ce n'est pas la même grandeur (décision 7).
- Critère de localisation du modèle sans w₂ (N ≠ 9, 27) : proposé = poids sur les p_z des trois voisins, seuil 3 × moyenne des états de la
  fenêtre de la parfaite repliée (analogue du seuil 3⟨w₂⟩_P de R4) ; « π quasi-lié » = état impair de plus grand poids (comme R7 avec w₂).

**Coût** : dimension 5N², bloc π 2N² (V_loc hors bloc σ–π = 0,0 exactement, R6 2.0 ; H(R) σ–π ≤ 5,5e-10 eV, R4 ; couplage pair–impair de
(c-3) 1,7e-7 eV, R6 2.1). Le bloc π suffit pour l'état π et la LDOS p_z ; le couplage σ–π abandonné (≤ 1,7e-7 eV) sera mesuré à N = 27
(complet contre bloc π). Au-delà de N = 27, chaque maille est couplée aux 729 mailles de H(R) (27×27) : ≈ 1 458 non-nuls par ligne du bloc π,
pas d'avantage au creux sans tronquer H(R) (approximation nouvelle, non proposée) ; le shift-invert près de −0,2 eV passe par la même
factorisation. Première couronne de la grille repliée autour de K (arithmétique, ħv_F|b|/N = 16,06 eV / N) :

| N | 5N² | Mo (complexe, complet) | bloc π 2N² | Mo (bloc π) | ħv_F\|b\|/N (eV) |
|---|---|---|---|---|---|
| 6–27 | 180–3 645 | ≤ 213 | 72–1 458 | ≤ 34 | 2,68–0,595 |
| 36 | 6 480 | 672 | 2 592 | 107 | 0,446 |
| 45 | 10 125 | 1 640 | 4 050 | 262 | 0,357 |
| 54 | 14 580 | 3 401 | 5 832 | 544 | 0,297 |
| 63 | 19 845 | 6 301 | 7 938 | 1 008 | 0,255 |
| 72 | 25 920 | 10 750 | 10 368 | 1 720 | 0,223 |
| 81 | 32 805 | 17 219 | 13 122 | 2 755 | 0,198 |
| 108 | 58 320 | 54 420 | 23 328 | 8 707 | 0,149 |

`eigh` complet ≈ 3 × la matrice en mémoire, temps ∝ n³ ; pas de mesure disponible ici : calibration dans le job sur N = 27, 36, 45 (bloc π),
extrapolation n³, arrêt si un N dépasse 1 h. **N_max proposé : 81 en dense sur le bloc π** (dimension 13 122, 2,8 Go ; 3m : 30, 36, 45, 54,
63, 72, 81) ; 108 en option si la calibration le permet. Pour N ≤ 27 (C.1), matrice complète (σ et π).

### 0.4 Kaasbjerg Fig. 3 : les deux valeurs déjà citées

**(a) R6 3.2, « |M|·A_cell à K = 107,22 eV Å² »** : `scripts/analyze_M.py` §1 (l. 39–56 ; job post_analyze 21862371). M2 dense 9×9
(`results/M2/M_dense_9x9.npy`, 27×27, 20 bandes, eV). k = K exactement (indice 495, (2/3, 1/3)). pairK = les deux bandes de plus grand poids p_z
|V_{n,pz(A)}|² + |V_{n,pz(B)}|² (V = U_dis·U) à K = bandes (3, 4) (0-based), dégénérées. Pour chaque k′ de la grille : (a, b) = paire π à k′ par
poids p_z, triée en énergie (a = π, b = π*). Ṽ_π(k′) = A_cell·[Σ_{n∈pairK} |M_{a k′, n K}|²]^{1/2}, A_cell = |a₁ × a₂| = 5,266 Å². « À K » :
k′ = K, a = un des deux membres de la paire dégénérée (une **ligne** du bloc 2 × 2, norme sur la paire du ket). Maximum sur les 729 k′ :
116,46 (π) ; π* : 146,02, à K 109,79. v1 : 23,87 à K.

**(b) « M̄ intravallée 76 eV Å² » (R5)** : la valeur n'apparaît telle quelle ni dans `R5_rapport.md`, ni dans les json/tables de R5, ni dans
NOTES_TGAMMA. La seule définition des dépôts qui la redonne : R4 D2 (`R4_rapport.md`, table « Base de Bloch ») — M̄ = ½ Tr du bloc 2 × 2 de la
paire π/π* dégénérée à (K, K) du M dense 27×27, × A_cell. v1 : 3,2426 eV = 17,07 eV Å² (L 0,1396 ; NL 3,1030). Avec M^L × 81 (convention M2) :
81 × 0,1396 + 3,1030 = 14,41 eV → 75,9 eV Å² ; c'est la grandeur d'`analyze_M.py` §3 (½ Re Tr à K ; R6 niveau 2, 9×9 : Re M^L(K) + Re M^NL(K)
= 11,305 + 3,103 = 14,408 eV → 75,87 eV Å²). À confirmer par Greg.

**Différence** : même M, même bloc (K, K), contraction différente : ½ Re Tr (invariant de jauge) contre norme d'une ligne (le bra est un membre de la
paire dégénérée : dépend de la jauge du `.save` à K). Rapport 107,22 / 75,87 = 1,413 (√2 = 1,414). Arithmétique sur v1 (R4 D2) : ‖bloc‖_F = 6,489 eV
= 2,001 × M̄ ; un bloc hermitien 2 × 2 avec ‖·‖_F = 2|½ Tr| a pour valeurs propres (0, Tr) (rang un). D.2 imprimera, pour le bloc M2, les valeurs
propres, ½ Tr, les deux normes de ligne et ‖·‖_F.

**(c) Grandeur de Kaasbjerg** (définitions du prompt) : Ṽ = A_cell·|M^{nn}_{k′k}| (norme maille), intrabande (valence = π, conduction = π*),
k = K + δx̂ fixe (δ ≪ 2π/a), carte en k′ ; « ~70 eV Å² » = valeurs près de K et K′ ; V₀ ≈ 27 eV ; super-cellule 11×11. Faits : (a) et (b) sont
à k = K exactement (paire dégénérée), aucun des deux n'est un élément intrabande à k ≠ K ; (a) somme sur les deux membres de la paire du ket,
(b) trace la paire. Pour D.1 : δ = 0,01|b₁| = 0,0294 Å⁻¹ (0,27 pas de la grille 27×27, dont le pas est |b₁|/27 = 0,109 Å⁻¹ ; arithmétique :
ħv_F δ = 0,161 eV), k hors grille → interpolation de Wannier de production (5 WF, M_W dense 27×27, boîte MP simple comme `wannier_interpolate`),
paire π par poids p_z ; direction x̂ = x cartésien de la cellule QE (angle avec ΓK imprimé) ; carte k′ sur une grille 240² repliée dans la
première zone ; k′ = K et K′ exacts (dégénérés) exclus et comptés ; disques de rayon 0,05|b₁| = 0,147 Å⁻¹ (≈ 520 points de la grille 240²).

### 0.5 Fonctions à ajouter, plan, jobs, coûts

#### Fonctions proposées pour `src/` (non écrites ; Greg les écrit ou autorise Code)

| # | module | signature | équation | formes |
|---|---|---|---|---|
| P1 | `defects/alignment.py` | `atom_sphere_shifts(V_d, V_p, x_red_d, x_red_p, A_cols, radius) -> dict` | pour chaque atome a de la cellule avec lacune : Δ_a = ⟨V_d⟩_{S(τ_a, ρ)} − ⟨V_p⟩_{S(τ_p(a), ρ)}, S(τ, ρ) = {r de la grille : \|r − τ\|_image < ρ} (même ensemble de points que `sphere_average`, parcouru sur une boîte d'indices locale), p(a) = atome parfait le plus proche ; la valeur à `far_atom` redonne `far_atom_alignment` exactement (contrôle intégré) | V (n1, n2, n3) eV ; x_red (n, 3) ; A_cols (3, 3) Å ; sortie : `s_vac` (3,), `dist` (n_d,) Å (image minimale au site de la lacune), `shift`, `mean_d`, `mean_p` (n_d,) eV, `npts_d`, `npts_p`, `i_p` (n_d,) int |
| P2 | `wannier/wannier_interpolation.py` | `Mwr_to_Mwk_pairs(Mwr, R, k_bra, k_ket) -> Mwk` | Mwk[w, k′, W, k] = Σ_{R,R′} e^{−2πi k′·R} Mwr[w, R, W, R′] e^{+2πi k·R′} (convention de `Mwr_to_Mwk`, qu'elle redonne pour k_bra = k_ket) | Mwr (nw, nR, nw, nR) ; R (nR, 3) int ; k (n′, 3), (n, 3) ; sortie (nw, n′, nw, n) |
| P3 | `defects/many_body/local_tmatrix.py` | `cluster_ldos(g0, V_loc, idx) -> (rho, rho0)` | G = g₀ + g₀ t g₀ = (1 − g₀V)⁻¹ g₀ sur l'amas (exact pour des orbitales de l'amas) ; ρ_i = −(1/π) Im G_ii, ρ⁰_i = −(1/π) Im g₀,ii | g0 (nE, n, n) de `local_green_batch` ; V_loc (n, n) ; idx liste d'indices L·nw + w ; sorties (nE, len(idx)), états/eV/orbitale/spin |
| P4 | `wannier/supercell_fold.py` | `ldos_from_eigenpairs(e, v, idx, egrid, eta) -> rho` | ρ_i(E) = Σ_j \|v_ij\|² (η/π)/((E − e_j)² + η²) (même lorentzienne que `resonance_metrics`) | e (n,) ; v (n, n) vecteurs en colonnes ; idx ; egrid (nE,) ; sortie (nE, len(idx)) |
| P5 (option) | `defects/local_R.py` | `compute_ML_R_mpi_shared(…, Ved=None)` et `prep_realspace_inputs(…, Ved=None)` | ΔV fourni remplace V_d − V_p (même grille) | — ; **non nécessaire** si l'enveloppe `io` de 0.2 (ii) est acceptée (recommandé : noyau de production inchangé) |

Reste dans le pilote (pas destiné à `src/`) : reprises de logique inline des scripts de production — niveau 1 par `lt.scattering_rate_fast` ;
Σ complexe comme `rcut_resigma.py` ; courbe Γ_T, T̄(K), Born comme `resonance_metrics.py` ; `fine_pi` de `m_rcut_convergence.py` ; `pi_pair`
d'`analyze_M.py` ; `w_from_r` de `r6_d4.py` ; M_W périodique « M2 dense ⊂ N² k » comme la variante (b) de R6 2.1 ; enveloppe `io`. Chaque reprise
est vérifiée par la porte de sa partie contre le chiffre de R6 avant usage. Le pilote importe `r4_driver`, `r5_driver` et `r6_d4` comme l'a fait R6.

#### Plan par partie

- **A.0** : niveau 1 9×9, M2 tel quel, R_cut 3, 240², η 0,02, N_k^int 300 → 3 132,60 meV, −0,175 eV (écart toléré 0,01 meV, E_res identique).
- **A.1** : (i) P1 pour N = 5, 6, 7, 8, 9, 12 à 1,0 Å (0,5 Å rapporté), profil Δ_a(dist) ; (ii) diagonales de M_W(R, R) par orbitale (p_z(A), p_z(B),
  sp²) en fonction de la distance de l'atome de l'orbitale au site de la lacune, deux versions : **périodique** (M2 dense restreint aux N² k coïncidents,
  20 bandes, rotation U aux k coïncidents, double TF sur la grille N×N : c'est le M auquel s'applique l'identité exacte) et **boîte** (M_W dense de
  production, mailles de la boîte, distance au bord de la boîte indiquée) ; (iii) ⟨ΔV⟩_3D. Plateau et C_N : décision 2. Tableau N × (i, ii, iii,
  Lu), figure `offset_profiles`. Arrêt selon la décision 2.
- **A.2** : grossier 9×9, bandes 0–3, `prep_realspace_inputs` puis `prep["Ved"] = 1`, `compute_ML_R_mpi` (16 rangs) contre 81·𝕀 (1e-12 relatif) ;
  dense 9×9 (ii) → F_W : max|F_W − 𝕀₅| sur l'amas R_cut 3, sur R_cut 4, sur toute la boîte, par distance au bord ; hors boîte ; ‖C_9(F_W − P_boîte)‖
  sur V_loc.
- **A.3** : variantes alignées sans recalcul de M :
  convention intensive — grossiers 5, 7, 8, 9, bandes 1–16, M2 − N_cells C_N 𝕀 (exact) ; **écart relevé** : la ligne de tab:tests_M (7,7e-2) est
  calculée par `analyze_M.py` l. 170 sur les M **denses** des **six** tailles, alors que son libellé (l. 171) dit « N = 5, 7, 8, 9 » ;
  arithmétique sur les maxima de R6 3.2 : denses six tailles (25,72 − 23,82)/24,70 = 7,7e-2 ; denses 5, 7, 8, 9 : 6,3e-2 ; grossiers 5, 7, 8, 9 :
  8,0e-2 (décision 5) ;
  tab:rcut_M 9×9 et 12×12, R_cut 0…6, (i) (et exact pour 9×9) ; niveau 1 six tailles × R_cut 3, 4 × (tel quel, aligné (i)) : médiane |Γ|·N_cells,
  E_res, pic de la courbe Γ_T, Re Σ médian, écarts des familles 3m / non-3m ; résonance 9×9 (pic Γ_T, pic −Im T̄(K), Re T̄(E_D), Born/T médian,
  bloc π : min|det|, min|λ|, 300²) ; C14 de R6 : énergie de l'écart relatif maximal de Γ_T (`Gamma_T_shiftp25`, `Gamma_T_shiftm25` contre `Gamma_T`
  de `resonance_9x9_shiftL.npz`, lecture seule).
- **B** : B.0 (−0,180 avec la courbe production `eg = egrid[::2]`, pas 5 meV ; −0,177 pour −Im T̄(K), 300²) ; B.1 9×9 et 12×12 × N_k^int 300, 450,
  600, 900 × (tel quel, aligné) : courbe Γ_T et −Im T̄(K) au pas de 2,5 meV (η/8, sans sous-échantillonnage) sur [−1,5, +1,5] eV, Γ_T(E_D) (états à K,
  définition de `nkint_check_post.py`), E_res (argmax des états 240² à ±1,5 eV) avec l'espacement mesuré des niveaux distincts de la grille 240²
  près de −0,2 eV ; B.2 tableau et figure `resonance_vs_nkint` (positions contre 1/N_k^int). Arithmétique (cône de Dirac, couronnes x = i² + ij + j²
  autour de K, énergie ħv_F|b|√x/N) : grille 240² → 0,134 (x = 4), 0,177 (7), 0,201 (9), 0,232 (12), 0,241 (13), 0,268 (16) eV ; les E_res de R6
  (−0,175 ; −0,202 ; −0,227 ; −0,238) sont à ≤ 5 meV des couronnes 7, 9, 12, 13. Grilles internes, écart entre couronnes près de 0,2 eV : 300² 8–21 meV,
  450² 3–15, 600² 2–9, 900² 1,6–3,3 meV (η = 20 meV).
- **C** : C.0 (c-3) 9×9 non aligné −0,677 eV (w₂ 0,227) à 1 meV ; C.1 N = 6…27 (3m), V_loc 9×9 de production (cache R6 `Vloc_M2_9x9.npz`),
  non aligné et aligné (C_9 sur les 29 sites) : E_D (quadruplet replié), états localisés (parité, poids p_z des trois voisins, ε − E_D), LDOS
  (P4, η 20 meV) et son maximum sur [−1, 0] eV, décalage rigide (QE R7 < E_D − 4) ; tableau N × (QE R7, non aligné, aligné) ; C.2 bloc π
  N = 30…81, ajustements 1/N, 1/N² sur 6…27 (comme R7) et valeurs extrapolées contre calculées ; C.3 P3 avec g₀ 600² sur [−1,1, +0,1] eV (pas
  2,5 meV), maximum, contre B (−Im T̄(K)) et C.2 ; figure `folded_vs_R7`.
- **D** : D.0 (proposée) : 107,22 eV Å² (définition (a), M2 dense sur grille) ; D.1 P2 + rotation H(k), carte Ṽ_π et Ṽ_π* (240²), disques K et
  K′, tel quel et aligné (exact 9×9, F_W de A.2) ; D.2 (a) et (b) recalculées, bloc 2 × 2 à (K, K) complet ; figure `kaasbjerg_fig3_map`.

Portes : échec → STOP sur la partie, rapporté. Paramètres de production (config v2) sauf N_k^int en B et les grilles de C/D. Boucles d'énergie
parallélisées à 1 fil BLAS (+ fils Python, `R4_EIG_WORKERS`, leçon R4/R6) ; g₀ (zgemm) et `eigh` des grandes matrices à 16 fils BLAS
(`threadpoolctl` si présent dans l'environnement, sinon étapes séparées). Figures : `figures/memoire.mplstyle`, `_palette.py`, texte en français.

#### Jobs (rrg-cotemich-ac ; lanceur : décision 10)

| job | sous-commandes | contenu | ressources, durée estimée | dépend de |
|---|---|---|---|---|
| J1 `r9a` | `a0`, `a1`, `a2c` | porte A.0 ; A.1 (i), (ii) périodique et boîte, (iii) ; décision d'arrêt ; A.2 grossier (srun 16 rangs) | 16 cœurs, 128 Go, 2 h | — |
| J2 `r9box` | `a2d` | M^L[1_boîte] 9×9 dense (noyau partagé, enveloppe `io`) → `cache/ML_box_9x9.npy` (3,4 Go) ; F_W | nœud exclusif 32 × 6, 2 h (étape ml 53 min en septembre) | — |
| J3 `r9a3` | `a3` | A.3 complet (+ variantes 9×9 exactes) | 16 cœurs, 128 Go, 4 h | afterok J1, J2 |
| J4 `r9b` | `b` | B.0–B.2 (g₀ 900² : ≈ 9 × le coût à 300²) | 16 cœurs, 96 Go, 3 h | afterok J1 |
| J5 `r9c` | `c` | C.0–C.3 (calibration n³ avant N ≥ 54) | 16 cœurs, 128 Go, 4 h | afterok J1 |
| J6 `r9d` | `d` | D.0–D.2 | 8 cœurs, 64 Go, 1 h | afterok J2 |

Total (bornes hautes) ≈ 14 h de tâches à 8–16 cœurs + 2 h de nœud exclusif (+ ≈ 4 h 30 si M^L[1_boîte] pour les six tailles). Disque : `cache/`
≈ 3,4 Go (+ 17 Go pour les six tailles), npz de résultats < 1 Go ; rien dans `results/`. Aucune suppression.

#### Décisions attendues (étape G)

1. **P1–P4** (et P5 si l'enveloppe `io` est refusée) : écrites par Greg ou autorisées pour Code.
2. **Plateau de A.1 et portée du STOP.** Proposé : plateau (i) = atomes de la couronne extérieure (dist ≥ d_max − a) ; critère : max − min ≤ 5 meV et
   |moyenne(couronne extérieure) − moyenne(couronne précédente)| ≤ 5 meV ; (ii) idem sur les diagonales p_z (A et B) du M_W périodique ; STOP
   **par taille** (tailles sans plateau exclues des variantes alignées et rapportées) plutôt que sur toute la partie A. Faits de R4 D5 : max|ΔV − Lu|
   sur les plans frontière (plan des atomes) 135,0 / 78,3 / 77,5 / 26,5 / 17,0 / 16,8 meV pour N = 5 / 6 / 7 / 8 / 9 / 12.
3. **C_N retenu** : valeur de plateau de (i) (proposé), Lu (un atome) rapporté à côté ; (ii) comme contrôle.
4. **M^L[1_boîte]** : 9×9 seulement (proposé) ou six tailles (+ 4 h 30 de nœud exclusif).
5. **Convention intensive** : ensemble à rejouer (grossiers 5, 7, 8, 9 comme le prompt ; denses six tailles comme la valeur de tab:tests_M ; les deux).
6. **C** : V_loc 9×9 de production pour tous les N (aligné avec C_9 sur ses 29 sites) ; C.2 sur le bloc π ; N_max 81.
7. **Colonne QE de la LDOS en C.1** : substitut pondéré par w₂ (npz R7) ou colonne vide.
8. **Si A s'arrête** : B, C, D en « tel quel » seulement (proposé) ou attente.
9. **D.1** : grille 240², exclusion de K et K′ exacts, variante alignée exacte (F_W de J2) ou (i).
10. **Soumission** : `submit_r9.sh` (lanceur seul) ou `sbatch --wrap`.

**STOP — phase 0 terminée le 2026-09-27.** Rien n'est calculé, rien n'est écrit dans `src/`, `scripts/` ni `results/` ; attente de l'étape G
puis du GO.

## Étape G (2026-09-28) — décisions de Greg, fonctions écrites, premiers jobs

### Décisions (Greg, 2026-09-28)

- Fonctions : P1–P4 écrites par Code (lisibles, commentées, docstring avec l'équation, test minimal) ; **P5 refusée** : méthode exacte par le
  module `io` substitut dans le pilote seulement, noyau de production intouché.
- 1. STOP par taille ; une taille sans plateau est rapportée, sa variante alignée est calculée avec Lu et marquée « sans plateau » dans toutes les tables.
- 2. C_N = plateau de (i), Lu à côté ; plateau = moyenne des décalages de sphère (P1) sur les atomes à distance ≥ 0,75·r_max de la lacune ;
  moyenne, rms, max|écart|, nombre d'atomes ; critère (i)/(ii) de 5 meV inchangé.
- 3. M^L[1_boîte] exact : 9×9 et 5×5.
- 4. Convention intensive : grossiers 5, 7, 8, 9 et denses six tailles, aligné et non aligné, chaque ligne étiquetée par son ensemble.
- 5. LDOS QE en C.1 : colonne vide ; comparaison à QE sur les états (ε − E_D, parité, poids p_z des voisins, w₂), référence d'énergie écrite par colonne.
- 6. Lanceur unique `submit_r9.sh` (tâche par argument), diff avant chaque soumission.
- Amendements : B — observables principaux = pic de la courbe Γ_T et pic de −Im T̄(K) (pas ≤ 2,5 meV) ; E_res rapporté avec l'indice et l'énergie de sa
  couronne 240², de même pour les E_res de R6. C.2 — famille 3m 36, 45, 54, 63, 72, 81 (bloc π), première couronne à côté de ε et du maximum de LDOS.
  Écart 1 (convention intensive) consigné, rien corrigé dans le dépôt. 76 eV Å² = A_cell·½ Tr du bloc π à (K, K).
- GO immédiat : jobs M^L[1_boîte] 9×9 et 5×5, B non aligné.

### Fonctions écrites dans `src/` (non commitées ; Greg commit)

| # | fichier | fonction | contenu |
|---|---|---|---|
| P1 | `defects/alignment.py` | `atom_sphere_shifts(V_d, V_p, x_red_d, x_red_p, A_cols, radius)` + aides `true_min_image_dist`, `_sphere_points` | Δ_a = ⟨V_d⟩_S(τ_a) − ⟨V_p⟩_S(τ_p(a)) autour de chaque atome ; même ensemble de points que `sphere_average` (critère recopié, parcours sur une boîte d'indices locale) ; renvoie `dist` (vraie image minimale), `dist_axis` (convention de `far_atom`), `shift`, moyennes, nombres de points, partenaire, `i_far_axis` |
| P2 | `wannier/wannier_interpolation.py` | `Mwr_to_Mwk_pairs(Mwr, R, k_bra, k_ket)` | M[w, k′, W, k] = Σ e^{−2πi k′·R} Mwr e^{+2πi k·R′}, listes k′ et k indépendantes |
| P3 | `defects/many_body/local_tmatrix.py` | `cluster_ldos(g0, V_loc, idx)` | ρ_i = −(1/π) Im G_ii, G = g₀ + g₀ T g₀ = [1 − g₀V]⁻¹ g₀ ; ρ⁰_i = −(1/π) Im g₀,ii |
| P4 | `wannier/supercell_fold.py` | `ldos_from_eigenpairs(e, v, idx, egrid, eta)` | ρ_i(E) = Σ_n \|⟨i\|n⟩\|² (η/π)/((E − ε_n)² + η²) |

Tests : `tests/test_r9_functions.py` (4 tests : P1 redonne `far_atom_alignment` à 1e-12 et un ΔV constant partout ; P2 = `Mwr_to_Mwk` et sous-bloc
rectangulaire ; P3 = bloc de l'amas de (z − H₀ − V)⁻¹ à 1e-12 ; P4 = −(1/π) Im de la résolvante à 1e-12). `PYTHONPATH=src pytest tests/` : 9 passés
(5 anciens + 4). Aucune autre ligne de `src/` modifiée (en-têtes de module complétés).

**Écart relevé (définition, rien corrigé)** : `alignment.min_image_dist`, donc `far_atom` (atome de Lu, R4/R5/R7), réduit chaque coordonnée réduite
séparément ; dans une cellule à 60° ce n'est pas l'image minimale. Distances à la lacune de l'atome de Lu, axe par axe → vraie : 5×5 9,97 → 6,21 Å ;
6×6 12,81 → 7,40 ; 9×9 18,51 → 11,12 ; 12×12 25,63 → 14,80. Vrai r_max (atome le plus éloigné) : 7,12 / 8,54 / 12,81 / 17,08 Å (= a_sc/√3) ; atomes à
≥ 0,75 r_max : 13 / 26 / 53 / 95. La valeur de Lu n'en dépend pas (sphère de 1 Å) ; R9 utilise la vraie image minimale pour r_max et le plateau et
rapporte les deux distances.

### Pilote et lanceur

`r9_driver.py` (sous-commandes a2c, a2d, a2dpost, a0, a1, a3, a3pole, b, btables, c, d) ; `submit_r9.sh` (tâches box, b, a, a3, c, d ; hors job :
diff de `r9_driver.py` et de `submit_r9.sh` depuis la soumission précédente, puis sbatch ; copie `submitted/<jobid>/` avec `diff.txt`, ligne dans `JOBID`).
Reprises de logique de production dans le pilote, chacune vérifiée par une porte : niveau 1 (`scattering_rate_fast` / `rcut_resigma.py`, grille « sig »),
courbe Γ_T, T̄(K), Born (`resonance_metrics.py`, grille « res »), `fine_pi` (`m_rcut_convergence.py`), `pi_pair` (`analyze_M.py`), `w_from_r` (`r6_d4.py`).
Contrôle sur le nœud de connexion (cache V_loc de R6, H(R) 27×27) : le repliement H(R) + V_loc de C redonne (c-3) de R6 2.1 à la sixième décimale
(π −0,677446, σ −0,813248, −1,765234, −1,762682, +0,302370) ; bloc π seul = bloc impair complet (écart 0,0).

### Jobs du GO immédiat (2026-09-28 08:22, `JOBID`, `submitted/`)

| job | tâche | contenu |
|---|---|---|
| 21954787 `r9box9x9` | box 9x9 | a2c (M^L d'une constante, 9×9 grossière, bandes 0–3, 32 rangs) ; a2d (M^L[1_boîte] dense 9×9) ; a2dpost |
| 21954788 `r9box5x5` | box 5x5 | a2d (M^L[1_boîte] dense 5×5) ; a2dpost |
| 21954789 `r9b` | b --variants brut | porte B.0 puis B.1 non aligné (9×9, 12×12 ; 300, 450, 600, 900) |

Premiers résultats bruts :
- **A.2 grossier (a2c)** : max|M^L[1] − 81·𝕀|/81 = **1,56e-13** (seuil 1e-12 → OK) ; diagonale 81,000000000000 (min = max) ; hors diagonale max 1,26e-11 ;
  hermiticité 9,8e-17 ; 20 s.
- **Porte B.0** : pic de la courbe Γ_T −0,180 eV (R6 −0,180), pic de −Im T̄(K) −0,1775 (R6 −0,177), courbe et T̄ identiques à `resonance_9x9.npz`
  (écart relatif 0,0), médiane des états 3 132,87 meV (R6 3 132,87) → **OK**.

## Résultats (GO du 2026-09-28)

Préférence de Greg (2026-09-28) : pas de GO à attendre à chaque diff ; les diffs restent archivés (`submitted/<jobid>/diff.txt`) et chaque job écrit
le md5 du pilote dans son journal. Réponses de Greg aux points ouverts : « c'est correct comme ça » ; appliqués tels qu'annoncés — « sans plateau » si
max|écart| au plateau > 5 meV ; une taille dont (i) et (ii) diffèrent de plus de 5 meV est traitée comme « sans plateau » ; convention intensive,
denses 6, 7, 8, 12 alignés sur les seuls blocs k = k′ (exact : M^L[1_boîte]_{n′k,nk} = N_cells δ_nn′), 5×5 et 9×9 exacts.

Jobs (`JOBID`) : 21954787 `r9box9x9` COMPLETED 54 min ; 21954788 `r9box5x5` COMPLETED 14 min ; 21954789 `r9b` (B non aligné) COMPLETED 42 min ;
21955512 `r9a` FAILED (erreur d'exécution dans a1 : masque p_z de la version « boîte » construit sur les N² mailles périodiques au lieu des D² ; corrigé,
diff archivé) et ses quatre dépendants annulés par SLURM ; resoumis : 21955644 `r9a` COMPLETED 2 min 40, 21955645 `r9a3`, 21955646 `r9b` (aligné, exact),
21955647 `r9c`, 21955648 `r9d` COMPLETED 25 s.

### A.0 — porte (job 21955644)

Niveau 1, 9×9, M2 tel quel, R_cut 3, 240², η 0,02, N_k^int 300 (chaîne « sig » de `scattering_rate_fast` / `rcut_resigma.py`) : médiane Γ·N_cells
**3 132,60 meV** (specwd 3 132,5968, écart < 1e-9 relatif), E_res **−0,1748 eV** (R6 −0,175), Re Σ médian −26,65 meV (R6 −26,65), |Re Σ|/Γ médian 0,125,
min Γ +1,75 meV, 41 266 états → **OK**.

### A.1 — C_N par taille (job 21955644 ; `a/a1_results.json`, `a/A1_tables.md`, `a/a1_profiles_*.npz`, `fig/offset_profiles`)

Plateau = atomes à distance vraie ≥ 0,75 r_max de la lacune (décision 2). (i) sphères de 1,0 Å (P1) ; (ii) diagonale p_z de M_W « périodique » (M2 dense
restreint aux N² k coïncidents, 20 bandes, U aux mêmes k, TF sur la grille N×N ; position de chaque orbitale = son atome, écart aux atomes de la parfaite
≤ 2,1e-9 Å) et « boîte » (M_W dense de production, mailles de la boîte à ≥ 1 maille du bord) ; (iii) ⟨V_d⟩ − ⟨V_p⟩ 3D. meV.

| N | r_max (Å) | (i) 1,0 Å : moyenne ; rms ; max\|écart\| ; n | (i) 0,5 Å | (ii) p_z périodique : moyenne ; rms ; max\|écart\| ; n | (ii) sp² périodique | (ii) p_z boîte ; n | (iii) | Lu 1,0 Å | \|(i) − (ii)\| | C_N retenu |
|---|---|---|---|---|---|---|---|---|---|---|
| 5×5 | 7,12 | −57,56 ; 94,14 ; 233,70 ; 13 | −59,98 | −51,22 ; 90,36 ; 235,74 ; 13 | −32,13 | −111,54 ; 1 | +28,62 | −122,33 | 6,34 | −122,33 (sans plateau, Lu) |
| 6×6 | 8,54 | −50,51 ; 24,56 ; 50,04 ; 26 | −52,33 | −48,55 ; 24,86 ; 52,65 ; 26 | −50,51 | −44,04 ; 4 | +20,83 | −98,62 | 1,96 | −98,62 (sans plateau, Lu) |
| 7×7 | 9,97 | −26,91 ; 18,53 ; 35,58 ; 31 | −27,81 | −25,96 ; 17,84 ; 33,36 ; 31 | −22,49 | −31,16 ; 5 | +14,55 | −38,54 | 0,96 | −38,54 (sans plateau, Lu) |
| 8×8 | 11,39 | −14,41 ; 18,87 ; 59,19 ; 49 | −14,58 | −13,63 ; 18,37 ; 59,52 ; 48 | −14,07 | +7,22 ; 10 | +11,23 | −22,97 | 0,78 | −22,97 (sans plateau, Lu) |
| 9×9 | 12,81 | −25,14 ; 9,05 ; 24,10 ; 53 | −25,29 | −24,71 ; 9,38 ; 25,35 ; 53 | −24,76 | −25,99 ; 13 | +9,09 | −24,65 | 0,43 | −24,65 (sans plateau, Lu) |
| 12×12 | 17,08 | −18,69 ; 5,42 ; 13,88 ; 95 | −19,08 | −18,38 ; 5,49 ; 14,65 ; 94 | −18,71 | −18,02 ; 36 | +5,06 | −29,24 | 0,31 | −29,24 (sans plateau, Lu) |

- Aucune taille n'a max|écart| ≤ 5 meV sur son plateau : les six tailles sont « sans plateau » ; leurs variantes alignées (A.3, B, C, D) utilisent Lu (1,0 Å),
  étiqueté. La 5×5 dépasse aussi le critère (i)/(ii) (6,34 meV) ; les cinq autres sont à ≤ 1,96 meV.
- P1 à l'atome de Lu redonne `far_atom_alignment` (même valeur au bit) et les valeurs publiées (R4 D5, R5 C) pour les six tailles.
- Profils (figure) : (i) et (ii) se suivent atome par atome à toutes les distances ; la moyenne du plateau et Lu diffèrent de 64,8 (5×5), 48,1 (6×6),
  11,6 (7×7), 8,6 (8×8), 0,5 (9×9) et 10,6 meV (12×12).

### A.2 — identité et facteur de forme de la boîte (jobs 21954787, 21954788 ; `a2/`, `cache/ML_box_{9x9,5x5}.npy`, `cache/Fw_{9x9,5x5}.npz`)

- Grossier 9×9 (a2c, bandes 0–3, `compute_ML_R_mpi`, ΔV = 1) : max|M^L[1] − 81·𝕀|/81 = **1,56e-13** (seuil 1e-12) → OK ; diagonale 81,000000000000 ;
  hors diagonale ≤ 1,26e-11 ; 20 s.
- Dense, M^L[1_boîte] par `compute_ML_R_mpi_shared` inchangé (module `io` substitut) : 9×9 en 53,1 min (même durée que l'étape ml de septembre),
  5×5 en 13,4 min. Contrôles : hermiticité 3,5e-16 / 2,9e-16 ; diagonale = N_cells à 3,4e-13 / 1,1e-13 ; bloc des k coïncidents = N_cells·𝕀 à 8,1e-8 (9×9,
  bruit des k du XML) / 3,0e-15 (5×5) ; règle de somme Σ_R F_ww(R, R) = N_cells pour les cinq orbitales à 3e-10 près (81 / 25).
- F_W(R, R) par distance au bord de la boîte (mailles ; identiques à 1e-5 près pour 9×9 et 5×5) — diagonale min par orbitale (sp² 1, sp² 2, sp² 3, p_z(A), p_z(B)) :
  dernière maille (bord 0) 0,521 / 0,521 / 0,975 / 0,925 / 0,900 ; bord 1 : ≥ 0,9972 ; bord 2 : ≥ 0,9999 ; bord ≥ 3 : ≥ 0,99998 ; première maille extérieure :
  jusqu'à 0,469 (sp²) et 0,050 (p_z(B)) ; ≥ 2 mailles à l'extérieur : ≤ 2,6e-3. Hors site : max|F(R, R′)| = 4,76e-2 (paires à cheval sur le bord).
- Sur l'amas de V_loc : 9×9, R_cut 3 (29 sites, aucun hors boîte, bord ≥ 1) : max|F − P_boîte| = **2,58e-3** (hors site 2,5e-4 ; deux sites à |F − 1| > 1e-3,
  Rn = (−3, 0) et (0, −3), sp² 0,9974) → écart (i)/exact ≤ 24,65 meV × 2,6e-3 = 0,06 meV ; R_cut 4 : 4,69e-1 (mailles du bord ±4). 5×5, R_cut 3 (4 sites hors
  boîte) : 4,79e-1 (hors site 2,3e-2 ; 26 sites à |F − P_boîte| > 1e-3).

### D — grandeur de Kaasbjerg (job 21955648 ; `d/d_results.json`, `d/D_tables.md`, `d/d1_maps.npz`, `fig/kaasbjerg_fig3_map`)

- Porte D.0 : Ṽ_π(K) (définition (a), `analyze_M.py` §1) = **107,22 eV Å²**, max **116,46** (R6 107,22 / 116,46) → OK ; A_cell = 5,2658 Å².
- D.2, bloc 2 × 2 de la paire π/π* dégénérée à (K, K), × A_cell (eV Å²) :

| partie | variante | valeurs propres | ½ Tr | normes des deux lignes | ‖·‖_F |
|---|---|---|---|---|---|
| total | M2 tel quel | −1,71 ; 153,45 | **75,87** | **107,22** ; 109,79 | 153,46 |
| total | aligné exact | 8,81 ; 163,96 | 86,38 | 114,73 ; 117,47 | 164,20 |
| L | tel quel | −1,71 ; 120,77 | 59,53 | 84,39 ; 86,41 | 120,78 |
| L | aligné exact | 8,81 ; 131,28 | 70,04 | 91,94 ; 94,13 | 131,58 |
| NL | — | −0,00 ; 32,68 | 16,34 | 22,83 ; 23,38 | 32,68 |

  (b) = ½ Tr = 75,87 ; (a) = norme de la première ligne = 107,22. Valeurs propres −1,71 et 153,45 : rapport 1,1 %. Alignement exact : ½ Tr +10,51 eV Å²
  = 81 × 24,65 meV × A_cell.
- D.1, Ṽ = A_cell·|M^{nn}_{k′k}|, k = K + δx̂ (δ = 0,01|b₁| = 0,0294 Å⁻¹ ; x̂ à 30,0° de ΓK dans la cellule QE ; ε(k) − E_D = −0,1614 (π) et +0,1604 eV (π*)),
  interpolation de Wannier de production (M_W dense 27×27, boîte MP de `Mwk_to_Mwr`, P2), carte 240² ; K et K′ exacts exclus (2 points) ; disques de
  rayon 0,05|b₁| = 0,147 Å⁻¹ (512 points chacun) :

| variante | bande | disque K : moyenne [min, max] | disque K′ : moyenne [min, max] |
|---|---|---|---|
| M2 tel quel | valence (π) | 76,96 [75,51 ; 78,67] | 81,18 [80,49 ; 82,31] |
| M2 tel quel | conduction (π*) | 78,16 [76,14 ; 80,37] | 82,30 [81,50 ; 83,70] |
| aligné exact | valence | 81,51 [77,36 ; 86,15] | 81,17 [80,53 ; 82,28] |
| aligné exact | conduction | 82,71 [77,79 ; 87,32] | 82,29 [81,50 ; 83,68] |

  Carte entière (240², hors K et K′) : valence 72,07–85,66 (aligné 72,07–86,15), conduction 73,75–104,88 (aligné 73,72–104,81) eV Å². Kaasbjerg : ~70 eV Å² près de K et K′, V₀ ≈ 27 eV, super-cellule 11×11.
- Exécution : la première figure dessinait le contour de zone avec la liste de coins d'`analyze_M.py` (§1, `map_corners`), qui commence par (1/3, 1/3) :
  dans la cellule QE (a₁, a₂ à 60°) ce ne sont pas les points K (K = (2/3, 1/3), `production.json`) ; **écart relevé dans le script de production**
  (contour de `fig_M_map` à vérifier ; rien corrigé dans le dépôt). Figure R9 refaite avec les coins K/K′ (`r9_driver.py dfig`, sans recalcul) ; les
  énergies à k du json étaient absolues sous le nom `eps_k_minus_ED` (corrigé : E_D Wannier de la même H(R), gardé sous `eps_k_abs`).

### A.3 — variantes alignées sans recalcul de M (job 21955645, 4 h 42 ; `a/a3_results.json`, `a/A3_tables.md`)

C_N utilisé partout : Lu (1,0 Å), étiqueté « sans plateau » (A.1). « Aligné (i) » = V_loc − C_N·P_boîte ; « exact » = V_loc − C_N·P F_W P (9×9, 5×5).
Toutes les lignes « M2 tel quel » redonnent `specwd_<S>_prod.npz` de R6 (médiane à 1e-6 relatif, E_res identique) : régression OK pour les 12 couples
(taille, R_cut).

**C14 de R6** (lecture de `resonance_9x9_shiftL.npz`) : C = +25 meV, écart relatif maximal de Γ_T 19,2 % à −0,115 eV (médian 3,40 %) ; C = −25 meV, 16,3 %
à −0,120 eV (médian 4,02 %).

**Convention intensive** (max|M|, bandes 1–16, eV ; C_N = Lu ; dans toutes les lignes alignées le maximum est sur un élément k = k′, dans aucune ligne brute) :

| ensemble | 5×5 | 6×6 | 7×7 | 8×8 | 9×9 | 12×12 | (max − min)/moyenne |
|---|---|---|---|---|---|---|---|
| grossiers 5, 7, 8, 9 — tel quel | 25,198 | — | 23,250 | 24,311 | 24,801 | — | 7,98e-2 |
| grossiers 5, 7, 8, 9 — aligné (exact, N_cells C_N 𝕀) | 26,625 | — | 24,498 | 24,444 | 25,561 | — | 8,63e-2 |
| denses six tailles — tel quel | 23,850 | 25,148 | 23,818 | 24,311 | 25,348 | 25,721 | 7,71e-2 (= tab:tests_M) |
| denses six tailles — aligné (5, 9 exacts ; 6, 7, 8, 12 blocs k = k′) | 25,717 | 27,454 | 24,498 | 24,444 | 25,677 | 27,939 | 1,347e-1 |
| denses 5, 7, 8, 9 — tel quel / aligné | | | | | | | 6,29e-2 / 5,07e-2 |

**tab:rcut_M** (max|ΔM|/max|M| de la paire π, grille fine 60²) :

| R_cut | sites | 9×9 tel quel | 9×9 aligné (i) | 9×9 aligné exact | 12×12 tel quel | 12×12 aligné (i) |
|---|---|---|---|---|---|---|
| 0 | 1 | 3,157e-1 | 3,160e-1 | 3,160e-1 | 8,673e-1 | 8,747e-1 |
| 1 | 5 | 2,129e-1 | 1,632e-1 | 1,624e-1 | 2,172e-1 | 2,318e-1 |
| 2 | 13 | 1,139e-1 | 5,583e-2 | 5,576e-2 | 1,249e-1 | 1,681e-1 |
| **3** | 29 | **6,679e-2** | **2,476e-2** | 2,496e-2 | **9,021e-2** | **7,308e-2** |
| 4 | 49 | 4,419e-2 | 1,713e-2 | 1,722e-2 | 7,570e-2 | 5,198e-2 |
| 5 | 81 | 1,749e-2 | 1,367e-2 | 1,398e-2 | 5,356e-2 | 3,260e-2 |
| 6 | 113 | 1,207e-2 | 1,201e-2 | 1,210e-2 | 3,288e-2 | 2,168e-2 |

max|M_π| (référence, toutes R) : 9×9 25,440 (tel quel), 25,561 eV (aligné) ; v1 9×9 R_cut 3 : 2,27e-2.

**Niveau 1** (240², η 0,02, N_k^int 300) et courbe Γ_T (grille de `resonance_metrics.py`) :

| taille | R_cut | variante | médiane Γ·N_cells (meV) | E_res (eV) | Re Σ médian (meV) | pic Γ_T 5 ; 2,5 meV | pic −Im T̄(K) | Re T̄(E_D) | sites hors boîte |
|---|---|---|---|---|---|---|---|---|---|
| 5×5 | 3 | tel quel | 3 288,90 | −0,284 | −1 429,90 | −0,305 ; −0,3025 | −0,2700 | +6,519 | 4 |
| 5×5 | 3 | aligné (i) | 3 140,34 | −0,175 | +1 174,43 | −0,180 ; −0,1800 | −0,1775 | +12,025 | 4 |
| 5×5 | 3 | aligné exact | 3 130,22 | −0,175 | +1 133,42 | −0,180 ; −0,1800 | −0,1775 | +11,948 | 4 |
| 5×5 | 4 | tel quel / aligné (i) / exact | 3 302,74 / 3 143,55 / 3 143,25 | −0,284 / −0,181 / −0,181 | −1 437,73 / +1 169,89 / +1 172,10 | −0,305 / −0,180 / −0,180 | −0,2700 / −0,1775 / −0,1775 | +6,452 / +11,921 / +11,889 | 24 |
| 6×6 | 3 | tel quel / aligné (i) | 3 155,98 / 3 439,97 | −0,227 / −0,134 | −695,92 / +1 560,98 | −0,240 / −0,175 | −0,2200 / −0,1675 | +9,218 / +14,805 | 2 |
| 6×6 | 4 | tel quel / aligné (i) | 3 206,64 / 3 556,00 | −0,227 / −0,134 | −996,33 / +1 999,29 | −0,240 / −0,175 | −0,2225 / −0,1675 | +8,764 / +15,260 | 14 |
| 7×7 | 3 | tel quel / aligné (i) | 3 052,03 / 2 927,30 | −0,269 / −0,227 | −1 033,37 / −62,97 | −0,300 / −0,240 | −0,2675 / −0,2625 | +7,165 / +8,961 | 0 |
| 7×7 | 4 | tel quel / aligné (i) | 3 146,18 / 2 944,65 | −0,269 / −0,227 | −1 398,98 / +141,59 | −0,300 / −0,240 | −0,2675 / −0,2625 | +6,696 / +9,185 | 4 |
| 8×8 | 3 | tel quel / aligné (i) | 3 020,42 / 2 958,71 | −0,269 / −0,227 | −796,51 / −216,25 | −0,245 / −0,240 | −0,2650 / −0,2250 | +7,930 / +9,036 | 0 |
| 8×8 | 4 | tel quel / aligné (i) | 3 084,43 / 2 970,15 | −0,269 / −0,227 | −1 141,13 / −172,31 | −0,275 / −0,240 | −0,2650 / −0,2250 | +7,486 / +9,057 | 2 |
| **9×9** | 3 | tel quel / aligné (i) / exact | **3 132,60** / 3 185,03 / 3 185,02 | −0,175 / −0,175 / −0,175 | −26,65 / +592,62 / +592,55 | −0,180 / −0,180 / −0,180 | −0,1775 / −0,1725 / −0,1725 | +11,122 / +12,594 / +12,594 | 0 |
| 9×9 | 4 | tel quel / aligné (i) / exact | 3 147,47 / 3 196,41 / 3 196,02 | −0,175 ×3 | −410,29 / +674,74 / +671,57 | −0,180 ×3 | −0,1775 / −0,1750 / −0,1750 | +10,556 / +12,606 / +12,603 | 0 |
| 12×12 | 3 | tel quel / aligné (i) | 3 162,78 / 3 332,74 | −0,175 / −0,134 | +448,29 / +1 166,10 | −0,180 / −0,175 | −0,1750 / −0,1675 | +12,429 / +14,323 | 0 |
| 12×12 | 4 | tel quel / aligné (i) | 3 161,49 / 3 395,10 | −0,175 / −0,134 | +249,19 / +1 514,67 | −0,180 / −0,175 | −0,1750 / −0,1675 | +12,104 / +14,722 | 0 |

Familles (médianes Γ·N_cells ; moyenne, (max − min)/moyenne) : R_cut 3 tel quel 3m 3 150,5 meV, 0,96 % ; non-3m 3 120,4, 8,60 % — aligné (i) 3m 3 319,2, 7,68 % ;
non-3m 3 008,8, 7,08 %. R_cut 4 tel quel 3m 3 171,9, 1,87 % ; non-3m 3 177,8, 6,87 % — aligné 3m 3 382,5, 10,63 % ; non-3m 3 019,4, 6,59 %.
Γ_T(E_D) (états à E_D ± 1e-6 eV) : non défini pour 5×5, 7×7, 8×8 (aucun état de la grille 240² à moins de 1e-6 eV du E_D de leurs wannierisations) ;
9×9 3 592,7 / 4 224,8 meV, 6×6 2 788,4 / 5 202,9, 12×12 4 206,5 / 5 149,8 (tel quel / aligné, R_cut 3).

**Résonance 9×9** (R_cut 3, 300²) : Born/T médian 46,49 (min 1,11, max 211,7 ; Γ_Born médian 177 204,7 meV = R6 3.4) tel quel ; 45,84 (1,05 ; 203,1) aligné (i)
et exact ; max de la courbe Γ_T 37,84 → 40,56 eV ; max de −Im T̄(K) 24,45 → 26,87 eV ; |Re Σ|/Γ médian 0,125 → 0,215. Critère de pôle (fenêtre ±3 eV) :

| variante | bloc π : min\|det\|/max (ε − E_D) | bloc π : min\|λ\| (ε − E_D) ; λ | complet : min\|det\|/max | complet : min\|λ\| ; λ |
|---|---|---|---|---|
| tel quel | 2,330e-2 (−0,172) | 0,3967 (−0,170) ; +0,2186 + 0,3310 i | 1,318e-4 (−0,812) | 0,0019 (−0,812) ; +0,0019 i |
| aligné (i) = exact | 1,717e-2 (−0,170) | 0,3450 (−0,127) ; +0,2507 + 0,2370 i | 1,260e-4 (−0,787) | 0,0019 (−0,787) ; +0,0019 i |

(tel quel = R6 2.3 / 3.4 aux chiffres imprimés.)

### B — position de la résonance contre N_k^int (jobs 21954789 et 21955646 ; `b/b_results.json`, `b/B_tables.md`, `b/b_curves_*.npz`, `fig/resonance_vs_nkint`)

Porte B.0 : OK (voir Étape G). Courbes au pas de 2,5 meV (grille de `resonance_metrics.py` sans sous-échantillonnage) ; C_N = Lu (« sans plateau ») ; énergies
en eV relatives à E_D (Wannier) ; E_res avec sa couronne x = |i b₁ + j b₂|²/|b₁|² de la grille 240² (rang parmi les x réalisés, intervalle [min, max] des ε de la couronne).

| taille, variante | N_k^int | pic de la courbe Γ_T | pic de −Im T̄(K) | Γ_T(E_D) (meV) | E_res ; couronne |
|---|---|---|---|---|---|
| 9×9 tel quel | 300 / 450 / 600 / 900 | −0,1800 / −0,1825 / −0,1825 / −0,1825 | −0,1775 / −0,1825 / −0,1900 / −0,1925 | 3 592,7 / 3 446,4 / 3 419,2 / 3 412,0 | −0,1748 (x = 7) / −0,2018 (x = 9) / −0,2018 (x = 9) / −0,2018 (x = 9) |
| 9×9 aligné (i) | idem | −0,1775 / −0,1800 / −0,1800 / −0,1800 | −0,1725 / −0,1750 / −0,1750 / −0,1750 | 4 224,8 / 4 052,1 / 4 019,9 / 4 011,4 | −0,1748 (x = 7) ×4 |
| 9×9 aligné exact | idem | identique à (i) aux décimales imprimées | identique | 4 224,8 / 4 052,1 / 4 019,9 / 4 011,4 | −0,1748 (x = 7) ×4 |
| 12×12 tel quel | idem | −0,1775 / −0,1800 / −0,1800 / −0,1800 | −0,1750 / −0,1750 / −0,1775 / −0,1750 | 4 206,5 / 4 034,6 / 4 002,5 / 3 994,0 | −0,1749 (x = 7) ×4 |
| 12×12 aligné (i) | idem | −0,1775 ×4 | −0,1675 / −0,1500 / −0,1550 / −0,1550 | 5 149,8 / 4 939,4 / 4 900,1 / 4 889,7 | −0,1749 (x = 7) ×4 |

- Couronnes π de la grille 240² près de −0,2 eV (9×9 ; 12×12 identique à 0,1 meV) : x = 4 : −0,1344 ; x = 7 : [−0,1810 ; −0,1748] (24 états) ; x = 9 : −0,2018
  (12 états) ; x = 12 : [−0,2393 ; −0,2269] ; x = 13 : [−0,2472 ; −0,2382]. Écarts entre moyennes de couronnes voisines : 43,5 (4 → 7), 23,9 (7 → 9), 31,3 meV (9 → 12).
- E_res de R6 en regard (9×9) : N_k^int 150 (−0,238) → couronnes x = 13 ou 12 (intervalles qui se chevauchent) ; 300 (−0,175) → x = 7 ; 450 et 600 (−0,202) →
  x = 9 ; R_cut 1 (−0,227) → x = 12. 12×12 (−0,175) → x = 7.
- Γ_T(E_D) de B (grille « res », point le plus proche par arrondi) diffère de celui de R6 3.3 (grille de `rcut_resigma.py`, argmin) : 3 592,7 contre 3 630,7 meV à 300².
- Médiane des états (grille « res ») 9×9 : 3 132,87 / 3 135,92 / 3 136,20 / 3 136,61 meV (tel quel) ; 3 185,66 / 3 195,06 / 3 189,45 / 3 189,49 (aligné).

### C — chaîne repliée contre R7 (job 21955647, 2 h 08 ; `c/c_results.json`, `c/C_tables.md`, `c/c*_ldos_*.npz`, `fig/folded_vs_R7`)

V_loc = 9×9 de production (R_cut 3) ; aligné = V_loc − C_9·P_boîte avec C_9 = Lu = −24,65 meV (« sans plateau ») ; H(R) = `wannier/27x27`. Références d'énergie :
QE (R7) = ε − E_D(quadruplet de la parfaite QE) après alignement Lu 1,0 Å ; modèle = ε − E_D(quadruplet de la parfaite repliée) (= E_D Wannier à K à 1,4e-14 eV ;
quadruplet dégénéré à 6,1e-8 eV ; porte 1 de R4 ≤ 1,3e-12 eV pour N = 6 … 27). Couplage σ–π résiduel 1,7e-7 eV.

- **Porte C.0** : N = 9 non aligné, π localisé **−0,6774 eV, w₂ 0,2268** (R6 (c-3) −0,677, 0,227) → OK ; les cinq états localisés de (c-3) redonnés
  (−1,7652 ; −1,7627 ; −0,6774 ; +0,3024 impairs, −0,8132 pair).

**C.1** (matrice complète ; poids = Σ|c|² sur les p_z(B) des trois voisins ; « localisés » = impairs de la fenêtre [−3, +1] à poids > 3 × moyenne des impairs de la parfaite repliée) :

| N | 1ʳᵉ couronne π | QE : π quasi-lié ; w₂ | QE : paire σ ; w₂ | modèle tel quel : π ; poids ; n loc. | tel quel : pair de la fenêtre ; poids site | tel quel : max LDOS | aligné : π ; poids | aligné : pair | aligné : max LDOS | décalage rigide (meV) tel quel ; aligné |
|---|---|---|---|---|---|---|---|---|---|---|
| 6 | −2,658 | −1,0654 ; 0,284 | +0,156 ×2 ; 0,711 | −1,0438 ; 0,281 ; 4 | −0,831 ; 0,80 | −0,0000 (bord) | −1,0202 ; 0,280 | −0,806 | −0,9975 (bord) | +60,9 ; +41,7 |
| 9 | −1,803 | −0,7373 ; 0,246 | +0,101 ×2 ; 0,713 | −0,6774 ; 0,246 ; 4 | −0,813 ; 0,79 | −0,6775 | −0,6607 ; 0,244 | −0,789 | −0,6600 | +6,8 ; −0,9 |
| 12 | −1,357 | −0,5508 ; 0,223 | +0,114 ×2 ; 0,713 | −0,5168 ; 0,227 ; 9 | −0,813 ; 0,79 | −0,5175 | −0,5022 ; 0,225 | −0,788 | −0,5025 | +18,4 ; +13,9 |
| 15 | −1,085 | −0,4593 ; 0,207 | +0,103 ×2 ; 0,713 | −0,4277 ; 0,215 ; 11 | −0,813 ; 0,79 | −0,4275 | −0,4140 ; 0,213 | −0,788 | −0,4150 | +8,2 ; +5,4 |
| 18 | −0,904 | −0,3881 ; 0,194 | +0,107 ×2 ; 0,712 | −0,3711 ; 0,205 ; 26 | −0,813 ; 0,79 | −0,3700 | −0,3580 ; 0,203 | −0,788 | −0,3575 | +12,9 ; +10,8 |
| 21 | −0,774 | −0,3437 ; 0,183 | +0,103 ×2 ; 0,711 | −0,3320 ; 0,197 ; 31 | −0,813 ; 0,79 | −0,3325 | −0,3193 ; 0,196 | −0,788 | −0,3200 | +9,8 ; +8,2 |
| 24 | −0,677 | −0,3074 ; 0,174 | +0,103 ×2 ; 0,712 | −0,3033 ; 0,189 ; 38 | −0,813 ; 0,79 | −0,3025 | −0,2909 ; 0,189 | −0,788 | −0,2900 | +9,9 ; +8,7 |
| 27 | −0,601 | −0,2785 ; 0,166 | +0,102 ×2 ; 0,712 | −0,2812 ; 0,182 ; 48 | −0,813 ; 0,79 | −0,2800 | −0,2692 ; 0,182 | −0,788 | −0,2700 | +10,2 ; +9,3 |

Colonne LDOS de QE vide (décision 5). N = 6 : le π quasi-lié est hors de [−1, 0] eV, le maximum de LDOS est au bord de la fenêtre.

**C.2** (bloc π ; N = 27 : bloc π = matrice complète à 1,4e-14 eV près) :

| N | dim | 1ʳᵉ couronne π | tel quel : π quasi-lié ; poids | tel quel : max LDOS | aligné : π ; poids | aligné : max LDOS | eigh (s), tel quel ; aligné |
|---|---|---|---|---|---|---|---|
| 36 | 2 592 | −0,450 | −0,2370 ; 0,163 | −0,2375 | −0,2260 ; 0,165 | −0,2250 | 16 ; 16 |
| 45 | 4 050 | −0,360 | −0,2094 ; 0,144 | −0,2100 | −0,1995 ; 0,148 | −0,2000 | 64 ; 64 |
| 54 | 5 832 | −0,299 | −0,1896 ; 0,126 | −0,1900 | −0,1808 ; 0,132 | −0,1800 | 182 ; 183 |
| 63 | 7 938 | −0,257 | −0,1742 ; 0,109 | −0,1750 | −0,1664 ; 0,116 | −0,1675 | 451 ; 446 |
| 72 | 10 368 | −0,224 | −0,1613 ; 0,093 | −0,1625 | −0,1546 ; 0,101 | −0,1550 | 990 ; 979 |
| 81 | 13 122 | −0,199 | −0,1503 ; 0,078 | −0,1500 | −0,1445 ; 0,088 | −0,1450 | 1 992 ; 2 064 |

Ajustements ε(N) = ε∞ + a/N^p (π quasi-lié), comme R7 :

| série | points | 1/N : ε∞ ; a ; rms | 1/N² : ε∞ ; rms | 1/N extrapolé à N = 54 / 81 | calculé à 54 / 81 |
|---|---|---|---|---|---|
| modèle tel quel | 6 … 27 | −0,0486 ; −5,848 ; 0,0148 | −0,2805 ; 0,0298 | −0,1568 / −0,1207 | −0,1896 / −0,1503 |
| modèle aligné | 6 … 27 | −0,0402 ; −5,762 ; 0,0141 | −0,2688 ; 0,0298 | −0,1469 / −0,1113 | −0,1808 / −0,1445 |
| QE (R7) | 6 … 27 | −0,0523 ; −6,089 ; 0,0045 | −0,2963 ; 0,0465 | −0,1651 / −0,1275 | — |
| modèle tel quel (supplément) | 9 … 27 | −0,0781 ; −5,334 ; 0,0048 | −0,2514 ; 0,0141 | −0,1769 / −0,1440 | −0,1896 / −0,1503 |
| modèle aligné (supplément) | 9 … 27 | −0,0683 ; −5,273 ; 0,0048 | −0,2397 ; 0,0140 | −0,1660 / −0,1334 | −0,1808 / −0,1445 |
| QE (supplément) | 9 … 27 | −0,0499 ; −6,131 ; 0,0047 | −0,2495 ; 0,0187 | −0,1634 / −0,1256 | — |

Ajustements du maximum de LDOS sur 6 … 27 : pollués par le point N = 6 (maximum au bord de la fenêtre ; rms 0,178 tel quel) ; sur 9 … 27 (supplément) :
tel quel 1/N ε∞ −0,0768, rms 0,0047 ; aligné −0,0690, rms 0,0045 (`c/C_tables.md`).

**C.3** (limite diluée, P3, même V_loc, g₀ des caches de B) : maximum de la LDOS des trois p_z(B) voisins sur [−1, 0] eV : tel quel **−0,1925 eV** (0,747 états/eV,
600² ; 0,748, 900²) ; aligné **−0,1775** (600², 0,791) et −0,1750 (900², 0,790). En regard : pic de −Im T̄(K) (B) tel quel −0,1900 (600²) / −0,1925 (900²),
aligné −0,1750 / −0,1750 ; C.2 à N = 81 : max LDOS −0,1500 (tel quel), −0,1450 (aligné), première couronne −0,199.

### Fichiers, état du dépôt

- Répertoire de travail : `r9_driver.py`, `submit_r9.sh`, `submitted/<jobid>/` (pilote, lanceur, diff), `JOBID`, `r9_log.txt`, `slurm-r9-*`, `a/`, `a2/`, `b/`, `c/`, `d/`
  (json, tables, npz < 2 Mo), `fig/` (offset_profiles, resonance_vs_nkint, folded_vs_R7, kaasbjerg_fig3_map ; pdf et png), `cache/` (12 Go : ML_box 9×9 3,4 Go et
  5×5 2,5 Go, Fw, Mwr des six tailles, g₀) — TEST, supprimable après consignation (CLAUDE.md, règle 2), aucune suppression faite.
- Copie `article/R9_controles/` : rapport, README, pilote, lanceur, tables md, json, figures (pas de npz ni de slurm, pas de `cache/`).
- Dépôt (non commité, Greg) : `src/…/defects/alignment.py` (P1), `src/…/wannier/wannier_interpolation.py` (P2), `src/…/defects/many_body/local_tmatrix.py` (P3),
  `src/…/wannier/supercell_fold.py` (P4), `tests/test_r9_functions.py`, `article/R9_controles/`. Rien d'écrit dans `results/`, `scripts/`, `config/`.
- Écarts relevés, rien corrigé : (1) libellé de la ligne « convention intensive » d'`analyze_M.py` (l. 171 : « N = 5, 7, 8, 9 », calcul sur les denses six tailles) ;
  (2) `alignment.min_image_dist` / `far_atom` : image minimale axe par axe (distances de l'atome de Lu) ; (3) `analyze_M.py` §1 `map_corners` : coins de zone
  (1/3, 1/3), … qui ne sont pas les points K de la cellule QE (contour de `fig_M_map` à vérifier).

**STOP — R9 terminé le 2026-09-28.** Chiffres bruts ; le choix « sans plateau » (seuil 5 meV sur max|écart|) fait utiliser Lu pour les six tailles ; un rejeu avec la
moyenne du plateau comme C_N ne demande que les caches (A.3, B, C, D).

## Clôture — rejeu C_N = plateau (2026-09-28)

C_N = moyenne du plateau (i) : estimateur choisi après lecture de A.1, avant tout résultat aligné avec cet estimateur.

Prompt « R9 — Clôture » (Greg, 2026-09-28). Statut TEST (inchangé). Ordre : phase 0 → STOP → (GO) R → rapport → README → manifeste → STOP. Aucune suppression.
Cadre de R9 ; lecture seule de `results/`, `config/`, `scripts/`, `src/` (P1–P4 compris) ; caches de `cache/` réutilisés, aucun recalcul de M ni de M^L[1_boîte].

### Phase 0 (rien n'est calculé ; le pilote n'est pas modifié)

**0.1 Valeurs de C_N** (`a/a1_results.json` : `C_i_eV` = moyenne (i), 1,0 Å, atomes à ≥ 0,75 r_max ; `C_retenu_eV` = Lu, utilisé en R9). meV ; arithmétique :
N_cells·(plateau − Lu), terme extensif ajouté à la diagonale de M.

| N | plateau (i) | Lu (R9) | plateau − Lu | N_cells·(plateau − Lu) |
|---|---|---|---|---|
| 5×5 | −57,56 | −122,33 | +64,77 | +1,619 eV |
| 6×6 | −50,51 | −98,62 | +48,11 | +1,732 eV |
| 7×7 | −26,91 | −38,54 | +11,62 | +0,570 eV |
| 8×8 | −14,41 | −22,97 | +8,56 | +0,548 eV |
| 9×9 | −25,14 | −24,65 | −0,49 | −0,040 eV |
| 12×12 | −18,69 | −29,24 | +10,55 | +1,519 eV |

Valeurs attendues du prompt retrouvées (−57,56 / −50,51 / −26,91 / −14,41 / −25,14 / −18,69).

**0.2 Changement du pilote** (proposé, non appliqué : `cloture/phase0_r9_driver.diff`, 492 lignes de diff dont 269 changées ; version complète `cloture/r9_driver_propose.py` ;
lanceur `cloture/phase0_submit_r9.diff`, `cloture/submit_r9_propose.sh`) :
- option `--cn lu | plateau` sur `a3`, `a3pole`, `b`, `btables` ; `C_of(S)` renvoie `C_i_eV` (« plateau (i) ») ou `C_retenu_eV` (« Lu (1,0 Å) », avec assertion que
  c'est bien Lu) ; sans `--cn`, comportement de R9 inchangé ;
- étiquette dans toutes les sorties : fichiers suffixés (`a3_results_plateau.json`, `A3_tables_plateau.md`, `b_results_plateau.json`, `B_tables_plateau.md`,
  `b_curves_<S>_plateau.npz`, figure intermédiaire `resonance_vs_nkint_plateau`) et clés de variantes `aligne_plateau`, `exact_plateau` (`aligne_lu`, … pour `--cn lu`) ;
  `V_variant` lit la base de la clé (`aligne_plateau` → `aligne`) ; les R9 existants (`a3_results.json`, `b_results.json`) ne sont ni lus en écriture ni réécrits ;
- avec `--cn`, `a3` et `a3pole` ne recalculent pas le « tel quel » (déjà dans R9) ; E_res de niveau 1 porte désormais l'état et sa couronne 240² (`level1_stats` renvoie
  aussi l'indice de l'état ; couronnes π sur [−1,5 ; 0] eV stockées par taille) ; familles calculées pour toutes les variantes présentes ;
- sous-commande `r0` (porte R.0) et sous-commande `synth` (R.4, aucune donnée recalculée) ; lanceur : tâche `r0`, `--cn` transmis à `a3pole`.
- Aucune autre modification de logique (formules, grilles, caches, seuils inchangés).

Note d'exécution : le `r9_driver.py` actuel (md5 f2ba6adb9a44) diffère déjà de la dernière soumission (`submitted/last`, afce56ad675b) par les post-traitements
sans calcul ajoutés après les jobs de R9 (`dfig`, `cfig`, coins de zone, légende) ; le diff de la prochaine soumission les montrera.

**0.3 Test intensif hors k = k′** (nouveau) : fonction `offdiag_intensive(cfg)` appelée par `a3` quand `--cn` est donné ; elle lit les deux C_N dans `a1_results.json`
(indépendamment de `--cn`) et calcule max|M| sur les éléments k ≠ k′ (bandes 1–16, eV) pour tel quel, Lu, plateau :
- grossiers 5, 7, 8, 9 : paires k ≠ k′ ; l'alignement exact (N_cells C_N 𝕀) ne touche que k = k′ : trois valeurs égales par construction ;
- denses six tailles : hors blocs k = k′ ; 5×5 et 9×9 alignés exactement (M2 − C_N M^L[1_boîte], qui a des éléments k ≠ k′ par le facteur de forme de la boîte) ;
  6×6, 7×7, 8×8, 12×12 : alignement des seuls blocs k = k′ → trois valeurs égales par construction (colonne « mode » dans le tableau).
Sortie : `offdiag_intensive` dans `a3_results_plateau.json`, tableau N × variante et (max − min)/moyenne par ensemble.

**0.4 Jobs** (rrg-cotemich-ac ; mêmes ressources que R9) :

| job | tâche | contenu | ressources | durée estimée | précédent |
|---|---|---|---|---|---|
| J1 `r9r0` | `r0` | R.0 : --cn lu, 9×9 et 12×12, R_cut 3 : médiane alignée (i) et tab:rcut_M contre `a3_results.json` (1e-9 relatif) | 16 cœurs, 96 Go, 1 h | ≈ 15 min | — |
| J2 `r9a3` | `a3 --cn plateau` puis `a3pole --cn plateau` | R.1, R.2 | 16 cœurs, 160 Go, 5 h | ≈ 2 h (caches g₀ et M_W : niveau 1 R_cut 3 ≈ 2,5 min et R_cut 4 ≈ 11 min par variante, 8 + 8 variantes) | A.3 de R9 : 4 h 42 dont ≈ 1 h 20 de g₀ |
| J3 `r9b` | `b --cn plateau --variants aligne,exact` | R.3 (+ porte B.0 rejouée) | 16 cœurs, 96 Go, 2 h | ≈ 15 min | 21955646 : 15 min 07 |
| — | `synth` (nœud de connexion) | R.4, figures | lecture de json | < 1 min | — |

J2 et J3 en afterok de J1.

**0.5 Parties non rejouées (C, D)** : lecture du pilote — `cmd_c` (l. 1083) et `cmd_d` (l. 1350) prennent `C9, tag9 = C_of("9x9")` et n'utilisent que C_9
(`c/c_results.json` et `d/d_results.json` : C9 = −24,6514 meV, Lu). Plateau − Lu = −0,49 meV à 9×9, soit 2,0 % de C_9. Ordre de grandeur attendu (arithmétique : les
écarts aligné − tel quel de R9 × 0,0200 ; pas de calcul) :
- C : V_loc,aligné(plateau) − V_loc,aligné(Lu) = +0,49 meV sur les 145 éléments diagonaux ; π quasi-lié : N = 9 : 16,7 meV × 0,02 = 0,33 meV ; N = 27 : 0,24 ; N = 54 : 0,18 ;
  N = 81 : 0,12 meV ; pair σ (−0,813 → −0,788) : 0,5 meV ; maximum de LDOS de C.3 (15 meV de décalage aligné − tel quel) : 0,3 meV, sous le pas de 2,5 meV ;
- D : ½ Tr du bloc (K, K) aligné exact, linéaire en C : 81 × A_cell × 0,49 meV = +0,21 eV Å² (86,38 → 86,59) ; moyennes des disques autour de K (écart exact − tel quel
  +4,55 eV Å²) : ≈ +0,09 eV Å² ; autour de K′ : < 0,01.

**STOP — phase 0 de la clôture terminée le 2026-09-28.** Rien n'est calculé ni appliqué ; attente du GO.

### GO (2026-09-28) et exécution

GO de Greg : appliquer le diff proposé, r0 → (afterok) `a3 --cn plateau` + `a3pole` et `b --cn plateau --variants aligne,exact` → `synth`. Précision R.2 appliquée
avant soumission (étiquettes (a)/(b) des lignes « identique par construction », note en tête du tableau). Diff de soumission archivé (`submitted/21976369/diff.txt`).

| job | tâche | état | durée |
|---|---|---|---|
| 21976369 `r9r0` | porte R.0 | COMPLETED | 5 min 11 |
| 21976370 `r9a3` | `a3 --cn plateau` puis `a3pole --cn plateau` | COMPLETED | 2 h 49 min 51 (niveau 1 à R_cut 4 ≈ 17 min par variante) |
| 21976371 `r9b` | `b --cn plateau --variants aligne,exact` | FAILED à la toute fin (voir ci-dessous) ; résultats complets écrits | 22 min 41 |

Erreur d'exécution (corrigée) : la figure intermédiaire de `b` ne retenait que les tailles ayant une variante « brut » ; le fichier `b_results_plateau.json` n'en contient
pas → `ValueError` de matplotlib après l'écriture du json et des tables. Sélection des tailles corrigée (`b_figure`), figure régénérée par `btables --cn plateau`
(sans calcul). Présentation de `resonance_vs_nkint` corrigée ensuite (légende commune sous les panneaux, « C_N » en LaTeX). `synth` exécuté sur le nœud de connexion
(lecture des json ; `cloture/synthese.md`, 220 lignes ; version R9 de la figure gardée sous `fig/resonance_vs_nkint_R9`).

### R.0 — porte (job 21976369)

`--cn lu` redonne R9 à toutes les décimales : 9×9, R_cut 3, C_N −24,6514 meV, médiane alignée (i) **3 185,0251 meV** (R9 3 185,0251), tab:rcut_M R_cut 3
**2,476216e-2** (R9 2,476216e-2) ; 12×12 : C_N −29,2351 meV, **3 332,7408 meV** (R9 3 332,7408), **7,308052e-2** (R9 7,308052e-2) → **OK**.

### R.1 — A.3 avec C_N = plateau (job 21976370 ; `a/a3_results_plateau.json`, `a/A3_tables_plateau.md`, `cloture/synthese.md`)

Niveau 1 (240², η 0,02, N_k^int 300), R_cut 3 ; tel quel / aligné (i) Lu / aligné (i) plateau (E_res avec sa couronne x de la grille 240² ; toutes les grandeurs et
R_cut 4 : `cloture/synthese.md`) :

| taille | médiane Γ·N_cells (meV) | E_res (eV ; couronne) | pic Γ_T (2,5 meV) | pic −Im T̄(K) | Re Σ médian (meV) | Re T̄(E_D) (eV) |
|---|---|---|---|---|---|---|
| 5×5 | 3 288,90 / 3 140,34 / 2 997,14 | −0,2843 (19) / −0,1748 (7) / −0,2268 (12) | −0,3025 / −0,1800 / −0,2425 | −0,2700 / −0,1775 / −0,2250 | −1 429,90 / +1 174,43 / −226,16 | +6,519 / +12,025 / +8,940 |
| 6×6 | 3 155,98 / 3 439,97 / 3 149,15 | −0,2270 (12) / −0,1345 (4) / −0,1749 (7) | −0,2375 / −0,1775 / −0,1800 | −0,2200 / −0,1675 / −0,1775 | −695,92 / +1 560,98 / +474,13 | +9,218 / +14,805 / +11,884 |
| 7×7 | 3 052,03 / 2 927,30 / 2 942,85 | −0,2694 (16) / −0,2270 (12) / −0,2694 (16) | −0,3000 / −0,2425 / −0,2450 | −0,2675 / −0,2625 / −0,2625 | −1 033,37 / −62,97 / −356,03 | +7,165 / +8,961 / +8,407 |
| 8×8 | 3 020,42 / 2 958,71 / 2 974,18 | −0,2695 (16) / −0,2270 (12) / −0,2270 (12) | −0,2450 / −0,2425 / −0,2425 | −0,2650 / −0,2250 / −0,2625 | −796,51 / −216,25 / −435,27 | +7,930 / +9,036 / +8,618 |
| 9×9 | 3 132,60 / 3 185,03 / 3 189,01 | −0,1748 (7) ×3 | −0,1800 / −0,1775 / −0,1775 | −0,1775 / −0,1725 / −0,1725 | −26,65 / +592,62 / +605,08 | +11,122 / +12,594 / +12,624 |
| 12×12 | 3 162,78 / 3 332,74 / 3 264,51 | −0,1749 (7) / −0,1345 (4) / −0,1749 (7) | −0,1775 / −0,1775 / −0,1775 | −0,1750 / −0,1675 / −0,1700 | +448,29 / +1 166,10 / +913,84 | +12,429 / +14,323 / +13,619 |

Exact (F_W), R_cut 3 : 5×5 plateau 2 999,20 meV, E_res −0,2268 (12) ; 9×9 plateau 3 189,00, −0,1748 (7). Médianes R_cut 4 (tel quel / Lu / plateau) : 5×5 3 302,74 / 3 143,55 /
2 997,07 ; 6×6 3 206,64 / 3 556,00 / 3 161,15 ; 7×7 3 146,18 / 2 944,65 / 2 953,28 ; 8×8 3 084,43 / 2 970,15 / 2 991,24 ; 9×9 3 147,47 / 3 196,41 / 3 200,00 ;
12×12 3 161,49 / 3 395,10 / 3 283,17.

Familles (médianes Γ·N_cells ; moyenne, (max − min)/moyenne) :

| R_cut | famille | tel quel | Lu | plateau |
|---|---|---|---|---|
| 3 | 3m (6, 9, 12) | 3 150,5 ; 0,96 % | 3 319,2 ; 7,68 % | 3 200,9 ; 3,60 % |
| 3 | non-3m (5, 7, 8) | 3 120,4 ; 8,60 % | 3 008,8 ; 7,08 % | 2 971,4 ; 1,83 % |
| 4 | 3m | 3 171,9 ; 1,87 % | 3 382,5 ; 10,63 % | 3 214,8 ; 3,80 % |
| 4 | non-3m | 3 177,8 ; 6,87 % | 3 019,4 ; 6,59 % | 2 980,5 ; 1,47 % |

tab:rcut_M (max|ΔM|/max|M| de la paire π, grille fine 60² ; `fig/rcut_aligned`) :

| R_cut | 9×9 tel quel | 9×9 Lu | 9×9 plateau | 9×9 exact Lu | 9×9 exact plateau | 12×12 tel quel | 12×12 Lu | 12×12 plateau |
|---|---|---|---|---|---|---|---|---|
| 0 | 3,157e-1 | 3,160e-1 | 3,171e-1 | 3,160e-1 | 3,171e-1 | 8,673e-1 | 8,747e-1 | 8,678e-1 |
| 1 | 2,129e-1 | 1,632e-1 | 1,630e-1 | 1,624e-1 | 1,621e-1 | 2,172e-1 | 2,318e-1 | 1,892e-1 |
| 2 | 1,139e-1 | 5,583e-2 | 5,568e-2 | 5,576e-2 | 5,562e-2 | 1,249e-1 | 1,681e-1 | 1,250e-1 |
| **3** | 6,679e-2 | 2,476e-2 | **2,470e-2** | 2,496e-2 | 2,495e-2 | 9,021e-2 | 7,308e-2 | **3,097e-2** |
| 4 | 4,419e-2 | 1,713e-2 | 1,718e-2 | 1,722e-2 | 1,726e-2 | 7,570e-2 | 5,198e-2 | 1,840e-2 |
| 5 | 1,749e-2 | 1,367e-2 | 1,358e-2 | 1,398e-2 | 1,389e-2 | 5,356e-2 | 3,260e-2 | 1,286e-2 |
| 6 | 1,207e-2 | 1,201e-2 | 1,199e-2 | 1,210e-2 | 1,209e-2 | 3,288e-2 | 2,168e-2 | 1,272e-2 |

Convention intensive (max|M|, bandes 1–16, eV ; tel quel / Lu / plateau) : grossiers 5 25,198 / 26,625 / 25,198 ; 7 23,250 / 24,498 / 23,928 ; 8 24,311 / 24,444 / 24,311 ;
9 24,801 / 25,561 / 25,601 ; denses (5, 9 exacts ; 6, 7, 8, 12 blocs k = k′) : 5 23,850 / 25,717 / 24,097 ; 6 25,148 / 27,454 / 25,722 ; 7 23,818 / 24,498 / 23,928 ;
8 24,311 / 24,444 / 24,311 ; 9 25,348 / 25,677 / 25,717 ; 12 25,721 / 27,939 / 26,420. (max − min)/moyenne : grossiers 5, 7, 8, 9 **7,98e-2 / 8,63e-2 / 6,76e-2** ;
denses six tailles **7,71e-2 / 1,347e-1 / 9,96e-2** ; denses 5, 7, 8, 9 6,29e-2 / 5,07e-2 / 7,30e-2.

Résonance 9×9 (R_cut 3, 300² ; tel quel / Lu / plateau ; exact = (i) aux décimales imprimées, sauf max de la courbe Γ_T exact Lu 40,557) : Born/T médian 46,49 / 45,84 / 45,85 ; max de la courbe Γ_T
37,835 / 40,558 / 40,598 eV ; max de −Im T̄(K) 24,45 / 26,87 / 26,91 eV ; |Re Σ|/Γ médian 0,125 / 0,215 / 0,217. Critère de pôle : bloc π min|det|/max 2,330e-2 (−0,172) /
1,717e-2 (−0,170) / 1,706e-2 (−0,170) ; min|λ| 0,3967 (−0,170) / 0,3450 (−0,127) / 0,3440 (−0,127) ; complet 1,318e-4 (−0,812) / 1,260e-4 (−0,787) / 1,261e-4 (−0,785),
|λ| 0,0019 (−0,812 / −0,787 / −0,787).

### R.2 — test intensif hors k = k′ (job 21976370 ; `offdiag_intensive`)

Lignes « identique par construction » : (a) grossiers, identité exacte (M^L[C] ∝ δ_kk′) ; (b) denses 6, 7, 8, 12, alignement approché sur les seuls blocs k = k′.
Seuls les denses 5×5 et 9×9 (alignement exact) mesurent un effet de l'alignement hors diagonale.

| ensemble | taille | tel quel | Lu | plateau | mode |
|---|---|---|---|---|---|
| grossiers | 5×5 / 7×7 / 8×8 / 9×9 | 25,198 / 23,250 / 24,311 / 24,801 | idem | idem | identique par construction (a) |
| denses | 5×5 | 23,850 | 25,546 | 24,038 | exact : effet mesuré |
| denses | 6×6 / 7×7 / 8×8 / 12×12 | 25,148 / 23,818 / 24,311 / 25,721 | idem | idem | identique par construction (b) |
| denses | 9×9 | 25,348 | 25,506 | 25,534 | exact : effet mesuré |

(max − min)/moyenne : grossiers 7,98e-2 (trois variantes) ; denses six tailles 7,71e-2 / 7,61e-2 / 7,69e-2.

### R.3 — B avec C_N = plateau (job 21976371 ; `b/b_results_plateau.json`, `b/B_tables_plateau.md`, `fig/resonance_vs_nkint`)

Porte B.0 rejouée : OK (écart 0,0 à `resonance_9x9.npz`). N_k^int 300 / 450 / 600 / 900 ; tel quel / Lu / plateau (exact plateau 9×9 = aligné plateau aux décimales imprimées,
sauf Γ_T(E_D) 450² : 4 065,5 contre 4 065,6) :

| taille | grandeur | tel quel | Lu | plateau |
|---|---|---|---|---|
| 9×9 | pic Γ_T | −0,1800 / −0,1825 / −0,1825 / −0,1825 | −0,1775 / −0,1800 / −0,1800 / −0,1800 | −0,1775 / −0,1800 / −0,1800 / −0,1800 |
| 9×9 | pic −Im T̄(K) | −0,1775 / −0,1825 / −0,1900 / −0,1925 | −0,1725 / −0,1750 / −0,1750 / −0,1750 | −0,1725 / −0,1750 / −0,1750 / −0,1725 |
| 9×9 | Γ_T(E_D) (meV) | 3 592,7 / 3 446,4 / 3 419,2 / 3 412,0 | 4 224,8 / 4 052,1 / 4 019,9 / 4 011,4 | 4 238,8 / 4 065,6 / 4 033,2 / 4 024,7 |
| 9×9 | E_res (couronne) | −0,1748 (7) / −0,2018 (9) / −0,2018 (9) / −0,2018 (9) | −0,1748 (7) ×4 | −0,1748 (7) ×4 |
| 12×12 | pic Γ_T | −0,1775 / −0,1800 / −0,1800 / −0,1800 | −0,1775 ×4 | −0,1775 ×4 |
| 12×12 | pic −Im T̄(K) | −0,1750 / −0,1750 / −0,1775 / −0,1750 | −0,1675 / −0,1500 / −0,1550 / −0,1550 | −0,1700 / −0,1525 / −0,1600 / −0,1625 |
| 12×12 | Γ_T(E_D) (meV) | 4 206,5 / 4 034,6 / 4 002,5 / 3 994,0 | 5 149,8 / 4 939,4 / 4 900,1 / 4 889,7 | 4 780,6 / 4 585,1 / 4 548,6 / 4 539,0 |
| 12×12 | E_res (couronne) | −0,1749 (7) ×4 | −0,1749 (7) ×4 | −0,1749 (7) ×4 |

### R.4 — synthèse

`cloture/synthese.md` (220 lignes) : chaque chiffre de R.1–R.3 en trois colonnes (tel quel, Lu, plateau) plus les deux colonnes exactes (9×9, 5×5) ; porte R.0 ; C_N par taille.
Figures : `fig/resonance_vs_nkint` (trois variantes ; version R9 gardée sous `fig/resonance_vs_nkint_R9`), `fig/rcut_aligned` (tab:rcut_M, 9×9 et 12×12, R_cut 0…6, trois
variantes + exactes pour la 9×9) ; `fig/offset_profiles` inchangée. C et D non rejoués (0.5).

Écart d'exécution relevé après coup (corrigé) : `a3_tables` écrivait toujours `a/A3_tables.md`, sans suffixe, contrairement à la phase 0 de la clôture ; le rejeu
`--cn plateau` avait remplacé la table de R9 par celle du plateau (les json n'étaient pas touchés). Nom de sortie corrigé (`A3_tables{suffixe}.md`), sous-commande
`a3tables [--cn]` ajoutée (tables depuis le json, sans calcul) ; `a/A3_tables.md` régénéré depuis `a3_results.json` est identique au bit à la version commitée
(5a4bc94) ; la table du plateau est `a/A3_tables_plateau.md`. Les autres sorties du rejeu étaient suffixées (json, courbes, tables et figure de B).

### README, manifeste, git

- `README.md` (répertoire de travail et copie) : format du CLAUDE.md, ligne « Lecture : » laissée vide.
- `manifeste_R9.md` (répertoire de travail et copie) : versionné (copie `article/R9_controles/`, P1–P4, tests : état git par fichier), gardé dans le répertoire de travail
  (npz, slurm, `submitted/`, `JOBID`, `r9_log.txt` ; ≈ 9,7 Mo), supprimable (`cache/`, 72 fichiers, 50,0 Go en taille apparente, 12 Go sur le disque ; rôle, commande et
  durée de reconstruction par fichier). Aucun fichier n'est promu en production. Aucune suppression faite.
- git (lecture) : `src/` (P1–P4) et `tests/test_r9_functions.py` sont dans le commit 5a4bc94 (« R9 checkpoint », poussé : `main...origin/main`), inchangés depuis ;
  restent à commiter par Greg les fichiers de `article/R9_controles/` listés à la fin du manifeste.

**STOP — clôture de R9 terminée le 2026-09-28.** La suppression de `cache/` attend un GO séparé sur `manifeste_R9.md`.
