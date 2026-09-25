# R5 — Base ou M ? Diagnostic de la marche (a1) de R4 ; vérité DFT selon la taille : rapport de campagne

Statut TEST. Prompt d'origine : R5 (Greg, 2026-09-25). Répertoire de travail : `graphene/qe/defects/R5_base_vs_M/` (hors dépôt).
Copie versionnée : `graphene-raman/article/R5_base_vs_M/` (au rapport final). Chemins relatifs à `/project/6004866/gregb26/`
(`$PROJECTS`) ; le dépôt est `graphene-raman/` ; le scratch est `/scratch/gregb26/qe_tmp/` (= `links/scratch/qe_tmp/`). Unités :
énergies en eV sauf mention ; les XML de QE et les fichiers M sont en Hartree ; les potentiels pp.x en Rydberg.

## Phase 0 — Localisation et interfaces (2026-09-25)

Rien de la campagne n'est calculé. Seules lectures faites sur le nœud de connexion : XML et attributs HDF5 des `.save`
(k, nbnd, igwx, gamma_only), en-têtes des `.out`, fichiers json/npz de R4, `sacct` des jobs de référence. Aucun `.save` ouvert en
écriture ; aucun job soumis ; git non touché.

### 0.1 Chemins

**Super-cellules N = 5…12** (scf à Γ, `K_POINTS gamma`, ecutwfc 100 Ry, `assume_isolated = '2D'`, c = 30 bohr, conv_thr 1e-10,
degauss 0,01 Ry mv ; inputs `graphene/qe/defects/super_cell/NxN/{defective,pristine}/scf.in`, `submit.scf` : 64 rangs MPI, 4 G/cœur).
`.save` sur le scratch `qe_tmp/defect_NxN_{d,p}/defect_NxN_{d,p}.save/`, miroir `graphene/qe/qe_tmp_backup/defect_NxN_{d,p}/`
(md5 2026-09-17, `MD5SUMS_2026-09-17.txt`), liens `graphene-raman/data/graphene/supercell/qe/defect_NxN_{d,p}.save`. Tous **gamma_only**
(demi-sphère), nks 1, npol 1, un seul `wfc1.hdf5` ; wfc lisibles par tranche de bandes avec `qe_gamma_io.read_wfc_gamma` (R4).

| N | atomes d / p | nbnd d / p | ondes planes (demi-sphère) | grille FFT | `wfc1.hdf5` d (Go) | E_F d / p (eV) | scf 64 rangs d / p (s ; itérations) |
|---|---|---|---|---|---|---|---|
| 5 | 49 / 50 | 130 / 130 | 119 010 | 150 × 150 × 192 | 0,25 | −4,59085 / −4,35505 | 12,9 / 12,9 (22 / 16) |
| 6 | 71 / 72 | 172 / 174 | 171 510 | 180 × 180 × 192 | 0,47 | −4,23901 / −4,21850 | 29,4 / 32,9 (23 / 26) |
| 7 | 97 / 98 | 224 / 226 | 233 409 | 216 × 216 × 192 | 0,84 | −4,67548 / −3,80934 | 50,4 / 41,0 (20 / 16) |
| 8 | 127 / 128 | 284 / 286 | 304 854 | 240 × 240 × 192 | 1,39 | −4,54876 / −3,96547 | 118 / 98,9 (24 / 17) |
| 9 | 161 / 162 | 352 / 354 | 385 710 | 270 × 270 × 192 | 2,18 | −4,22195 / −4,21952 | 200 / 161 (24 / 19) |
| 10 | 199 / 200 | 428 / 430 | 476 289 | 300 × 300 × 192 | 3,27 | −4,58372 / −4,07501 | 301 / 332 (23 / 18) |
| 11 | 241 / 242 | 512 / 514 | 576 360 | 360 × 360 × 192 | 4,73 | −4,50379 / −4,29446 | 721 / 473 (26 / 18) |
| 12 | 287 / 288 | 604 / 606 | 685 905 | 360 × 360 × 192 | 6,64 | −4,21722 / −4,21973 | 907 / 673 (27 / 17) |

Potentiels pp.x (plot_num 1, Ry) : `super_cell/NxN/{defective,pristine}/Vks_NxN_{d,p}` pour les 8 tailles (9×9 : 241 Mo ; 12×12 : 428 Mo) ;
`unit_cell/NxN/Vks_uc_NxN` aussi présents (3 Mo). Les grilles des Vks sont celles du tableau ; pour N = 7 la grille 216 n'est pas un
multiple de 7 (la chaîne M rééchantillonne à 217 ; sans effet ici : le protocole D1 travaille sur la grille du `.save`).

**Mailles unitaires N×N** (2 atomes, nbnd 16, ecutwfc 100 Ry, grille FFT 30 × 30 × 192, 9 629 ondes planes à Γ, nosym + noinv,
gamma_only = .FALSE., une `wfc<ik>.hdf5` par k) — scratch `qe_tmp/<répertoire>/<prefix>.save`, miroir `qe_tmp_backup/<répertoire>/`
(tous présents, `_wann` compris), liens `data/graphene/unit_cell/qe/defect_NxN.save` :

| N | répertoire / prefix | calcul | k | convention des k | remarque |
|---|---|---|---|---|---|
| 5 | `defect_unit_cell_5x5` (lien `defect_5x5.save`) ; `defect_unit_cell_5x5_wann` | scf ; nscf | 25 ; 25 | [−0,5, 0,5) ; liste | `unit_cell/5x5/bands.dat` (chemin Γ–K–M–Γ, 3 décimales, 2026-06-18) |
| 6 | `defect_unit_cell_6x6` | scf | 36 | [−0,5, 0,5) | K ∈ grille |
| 7 | `defect_unit_cell_7x7/defect_7x7.save` (prefix `defect_7x7`, lien `defect_7x7.save`) ; `defect_unit_cell_7x7_wann` | scf ; nscf | 49 ; 49 | [−0,5, 0,5) ; liste | — |
| 8 | `defect_unit_cell_8x8` ; `defect_unit_cell_8x8_wann` | scf ; nscf | 64 ; 64 | [−0,5, 0,5) ; liste | — |
| 9 | `defect_unit_cell_9x9` (production de `M_ed_9x9`, `M_L_9x9`) | scf | 81 | [−4/9, 4/9] | K ∈ grille ; **pas de nscf.in** |
| 10 | `defect_unit_cell_10x10` | nscf | 100 | [0, 1) | — |
| 11 | `defect_unit_cell_11x11` | **bands** | 181 | chemin Γ–K–M–Γ, 60 points par segment | la grille 121 k (nscf.out du 2026-06-19 11:32, valeurs propres non imprimées) a été **écrasée** par le run bands (12:40, même prefix/outdir) ; miroir identique ; K exact à l'indice 61 (k saisi à 6 décimales : 0,666667, 0,333333) |
| 12 | `defect_unit_cell_12x12` | scf | 144 | [−0,5, 0,5) | K ∈ grille |
| dense 27 | `defect_uc_dense_27` (production de `M_dense_9x9`, Wannier 27×27) | nscf | 729 | [0, 26/27] | nbnd 20, diago_full_acc, 7 × « c_bands: 1 eigenvalues not converged » |

**Input nscf de production de la maille 9×9 : il n'y en a pas.** `unit_cell/9x9/` ne contient que `scf.in` (calculation scf,
`K_POINTS automatic 9 9 1 0 0 0`, nbnd 16, conv_thr 1e-14, mixing_beta 0,3, electron_maxstep 400, nosym, noinv, assume_isolated 2D,
pseudo_dir `abinit_processing/pseudo`), `submit.scf` (64 rangs, 1 G/cœur, 2 h ; 13,1 s de mur, 3 nœuds), `pp.in`, `Vks_uc_9x9`.
Le nscf de production le plus proche est `unit_cell/27x27/nscf.in` (K_POINTS crystal 729, nbnd 20, `diago_full_acc = .true.`, pas de
`diago_thr_init`, conv_thr 1e-14 ; `submit_dense_nb20.sh` : 64 rangs, 2 G/cœur, 4 h ; 207 s de mur). L'input de B est donc dérivé de
`scf.in` (§0.5), avec un changement de plus que les quatre prévus : `calculation` scf → nscf.

**Espace disque d'un `.save` à 128 bandes × 81 k** : Σ_k igwx = 771 419 (9 456 à 9 629 par k) → 128 × 771 419 × 16 o = **1,58 Go** de
wfc (16 bandes : 0,197 Go calculé, 207 Mo mesuré), + XML (~0,4 Mo) + charge-density (2,1 Mo) ; pic pendant le run ≤ 2 × (fichiers wfc
temporaires) ≈ 3,2 Go. Scratch : 113 Go / 20 To utilisés. Remarque : 771 419 = 2 × 385 710 − 1 = sphère complète de la 9×9 à Γ :
la réunion des 81 bases d'ondes planes de maille est exactement la base de la super-cellule (voir 0.2).

**Données R4 réutilisées** (`R4_quasi_lie/`) : `cache/Mwr_9x9.npz`, `cache/Vloc_9x9.npz` ; `d1/d1_states.npz` (E_D −4,238471, seuil
0,0812, s_vac, états 297–326 de la 9×9 lacune : x, parité, w₂, w₁) et `d1/d1_results.json` (décalage Lu D 1,0 Å : −24,65 meV) ;
`d4/d4_results.json` (variantes a1_tot/L/NL, a2_16, a3_20, b_5wf, c_all, c_3, QE_D1 ; portes ; marches), `d4/d4_spectra.npz` (spectres
complets) ; `prep/M_coarse_{tot,L,NL}.npy` (Ha, 16 × 81 × 16 × 81), `prep/M_NL_coarse_new_9x9.npy` + sidecar, `prep/prep_results.json`
(D2, moyenne diagonale M^L grossier 76,58 meV). Les vecteurs propres de H_(a1) ne sont pas stockés (recalcul : 1 296 × 1 296, secondes).

Valeurs de référence du prompt, retrouvées : `d4_results.json` → QE_D1 pairs +0,101 ×2 (w₂ 0,713), impairs −1,759 ×2 (0,114),
−0,737 (0,246), +0,269 (0,102) ; a1_tot impairs −1,3496 (0,2260), −1,8030 ×2 (0,1009), pairs −2,812 ×2 (0,287) ; a2_16 impair −1,3496
(**0,0836**) ; a3_20 −1,359 ; b_5wf = c_all −1,315 ; `prep_results.json` → p_z–p_z 6,6166 = 0,3113 + 6,3053 eV, ⟨M^L⟩_diag grossier
76,58 meV ; `d1_results.json` → ⟨ΔV⟩_3D +9,09 meV, décalage Lu (1,0 Å) −24,65 meV, E_D(SC) −4,23847 eV.

### 0.2 Produit scalaire ⟨nk|Ψ_QE⟩

Ψ_QE est stocké dans `defect_9x9_d.save/wfc1.hdf5` en demi-sphère gamma_only : coefficients D(g) sur 385 710 indices de Miller g de
la super-cellule, D(−g) = D*(g), norme |D(0)|² + 2 Σ_{g≠0} |D(g)|² = 1. Un état de maille |nk⟩ (C_nk(G), Σ_G |C|² = 1 sur la maille,
k sur la grille 9×9) s'écrit sur la super-cellule avec les **mêmes coefficients** aux indices g_sc = 9(k + G) : ψ_nk = Ω_uc^{−1/2} Σ_G C e^{i(k+G)·r}
= √81 · Ω_sc^{−1/2} Σ_G C e^{i g_sc·r}. Normalisation : l'état normalisé sur la super-cellule est |nk⟩_sc = ψ_nk/√N_cells (N_cells = 81),
de coefficients C_nk(G) en g_sc ; donc

- c_nk = ⟨nk|Ψ⟩_sc = Σ_G C*_nk(G) D(9(k + G)) (D complété par conjugaison hors de la demi-sphère) ; Σ_nk |c_nk|² ≤ 1, et
  δ = 1 − Σ|c|² est **exactement** le poids hors des bandes retenues : la réunion des 81 sphères de maille (Σ igwx = 771 419) coïncide
  avec la sphère complète de la super-cellule (2 × 385 710 − 1 = 771 419), même ecut ; aucune onde plane n'est perdue.
- H_in = diag(ε_nk) + M/81 est bien l'Hamiltonien dans la base |nk⟩_sc (M en norme unit_cell, divisé une fois par N_cells, comme (a1)
  de R4) ; Ψ_in = Σ c_nk |nk⟩_sc, Ψ_out = Ψ − Ψ_in.

La routine de R4 (`supercell_fold.sc_planewave_index` + `folded_density_2d`) sert **en partie** : `sc_planewave_index` donne l'index
plat de g_sc = 9(k + G) sur la grille FFT 270 × 270 × 192 (utilisable tel quel pour le produit scalaire, par lecture de D placé sur
la grille, et pour Ψ_in), à condition de passer le k du `.save` qui fournit C (voir 0.3) ; `folded_density_2d` ne rend que |Ψ|² sommé sur
z (pas de Ψ(r) 3D ni de coefficients). Manquent : la mise sur grille des coefficients gamma_only complétés (faite en interne par
`qe_gamma_io.density_2d`, non exposée), le produit scalaire, et Ψ(r) 3D normalisé (∫|Ψ|² dV = 1 avec dV = Ω_sc/N_grid, Ω_sc = 81 Ω_uc).

Jauge : les ψ_nk doivent venir du `.save` qui a produit le M utilisé, sinon les phases (et les rotations dans les sous-espaces
dégénérés, bandes 15–17) rendent c†Mc faux. Vérifié par les dates et par R4 :

| base | M | `.save` de M | ψ_nk lus dans | constat |
|---|---|---|---|---|
| 16 bandes | `M_ed_9x9.npy` (= `M_L_9x9` 2026-06-18 17:29 + M^NL juin ; copies `prep/M_coarse_*.npy`) | `defect_unit_cell_9x9.save` (wfc du 2026-06-18 13:54, jamais régénéré) | le même | R4 J1 : M^NL recalculé sur ce `.save` = M_ed − M_L à 1,1e-16 relatif |
| 20 bandes | `M_dense_9x9.npy` restreint aux 81 k (indices 0, 3, 6, … de la 27×27) ; `M_L_dense` 2026-09-04 13:15, `M_NL_dense`/`M_dense` 13:47 | `defect_uc_dense_27.save` (wfc du 2026-09-04 12:09 ; `U`, `U_dis` de pw2wannier 12:19, manifeste `baa17b88b1b69e51`) | le même | dates cohérentes ; ε(16 bandes) − ε(dense, 81 k) ≤ 0,92 meV (R4) |

La jauge de Ψ_QE elle-même est indifférente (|c|², quotients de Rayleigh, résidus).

### 0.3 Cause du w₂ faux des variantes (a2)–(c) de R4

**Identifiée** (lecture seule, XML des deux `.save`) : les 81 k de la maille grossière sont écrits dans [−4/9, 4/9] (grille `automatic`,
QE ramène k − nint(k)) et les 729 k du `.save` dense dans [0, 26/27] (liste explicite du `nscf.in`). Le pilote R4 sélectionne les 81 k
du dense par égalité modulo 1 (`idx81`), puis construit l'index d'ondes planes de la super-cellule avec `sc_planewave_index(k81, G20, …)`,
c'est-à-dire g_sc = 9(k81 + G) alors que les coefficients C20 lus dans le dense sont référés à k27 = k81 + ΔG :

- ΔG ∈ {(0,0), (0,1), (1,0), (1,1)} ; **56 des 81 k** ont ΔG ≠ 0 (écart fractionnaire résiduel 4,4e-9) ; les ensembles {k + G} des deux
  `.save` sont identiques (vérifié au k n° 44 : 146 vecteurs en plan, nG 9 504 des deux côtés) ; attribut `xk` des `wfc<ik>.hdf5` du dense
  conforme au XML (cartésien).
- Effet : chaque ψ_nk du dense est multiplié par e^{2πi·9ΔG·x} (phase périodique de module 1) : |ψ_nk|² est inchangé (états de Bloch purs
  à w₂ = 0,027–0,030 dans toutes les variantes, ce qui a masqué l'erreur), mais les interférences entre k de ΔG différents sont fausses
  pour toute superposition → w₂ 0,084 au lieu de 0,226 à −1,350 eV pour (a2), même énergie, même vecteur propre à une jauge près.
- (a1) utilise C16 et k81 du même `.save` (grossier) : correct. (b) et (c) passent par `flat20` = le même index fautif : mêmes w₂ faux
  (0,083 / 0,106 / 0,088) ; les « poids WF site + voisins » (0,247, 0,448), les énergies, les portes 1 et 2 et les marches en énergie de
  R4 ne dépendent pas de cet index et restent valables. La correction R4 de J6 (étiquettes recentrées R_d) portait sur un autre bug, réel.

Correction (**erreur d'exécution**, dans le pilote R5, R4 non modifié) : passer à `sc_planewave_index` le k du `.save` qui fournit C
(`qe_io.get_k_red(dense)[idx81]`) et ajouter le garde-fou `supercell_fold.check_k_reference(k_used, k_save)` (refus si k_used − k_save
n'est pas nul à 1e-6 près — un multiple entier non nul est précisément l'erreur). Vérification en A.4 : w₂(a2) = w₂(a1) à 1e-3 pour les
états localisés (attendu : −1,350 → 0,226 et la paire σ à −2,812 → ≈ 0,287, alors classée « sous le seuil » dans R4) ; le tableau D4 est
refait avec les w₂ corrigés pour a2_16, a3_20, b_5wf, c_all, c_3 (les lignes « la paire paire passe sous le seuil » et « (a2) → (a3) … 
toujours sous le seuil » de R4 sont à considérer comme caduques ; une ligne d'erratum dans `R4_rapport.md` est à la discrétion de Greg).

### 0.4 Application de ΔV à un état de super-cellule en ondes planes

Existant (production, en lecture) :

- `qe_io.get_pot(path, subtract_mean, to_hartree)` : lit un filplot pp.x (Ry → Ha) ; `V.transpose(2,1,0)` donne [ix, iy, iz] ;
  `local_R.prep_realspace_inputs` fait ΔV = V_d − V_p (subtract_mean = False en production dense) sur la grille du `.save` (270 × 270 × 192
  pour la 9×9, commensurable : 9 × 30) ; ΔV^L par produit sur la grille = Σ_r Φ*(r) ΔV(r) Φ′(r) dV, dV = Ω_sc/N_grid, exactement le
  noyau de `compute_ML_R` (`M += dV (Ψ* ΔV) Ψᵀ`).
- `qe_gamma_io.read_wfc_gamma(save, bands=(b0, b1))` (demi-sphère), `density_2d` (complète la sphère en interne : A[idx] = C,
  A[−idx] = C*, puis ifftn), `mirror_parity_z`, `inplane_disc_mask`.
- `supercell_fold.sc_planewave_index`, `folded_density_2d` (Σ_nk d_nk C_nk(G) placé sur la grille, ifftn, |·|² sommé sur z).
- Partie non locale (`defects/non_local.py`) : `build_K_vectors(k_red, G_red, keep, B)` (K = k + G cartésiens, normes, directions),
  `compute_phase(K, τ)` (e^{−iK·τ}), `compute_angular_part(K̂, lmax)` (Y_lm), `pseudo_io.read_upf` (D_ii Ry → Ha, r·β_i, grille r ;
  C.upf : l = 0, 1, deux canaux par l), `fq_from_fr` (Hankel), préfacteur 4π/√Ω, B_lim = pref Σ_g C*(g) F_li(|K|) Y_lm(K̂) e^{−iK·τ},
  M_NL[n, j] = Σ_li D_li B[n] B*[j], et M^NL = M_d − M_p (mêmes atomes sauf l'atome retiré ⇒ ΔV^NL = − V^NL de l'atome retiré).
  Toutes acceptent nkpt = 1 (Γ de la super-cellule) et une liste de g arbitraire ; l'atome retiré est donné par `alignment.vacancy_site`
  (indice i_vac_p dans la parfaite, position s = (13/27, 13/27, 0)).

Manquant (fonctions nouvelles, §0.6) : (i) mise sur grille des coefficients gamma_only complétés et Ψ(r) 3D normalisé ; (ii) produit
scalaire c_nk ; (iii) Ψ_in(r) 3D à partir de d_nk (variante de `folded_density_2d` sans la somme sur z) ; (iv) ⟨Φ|ΔV^L|Φ′⟩ sur la grille ;
(v) projections B_lim(Φ) de l'atome retiré sur la sphère complète de la super-cellule (771 419 g, nkpt = 1, préfacteur 4π/√Ω_sc ; avec
Ω_sc = 81 Ω_uc, M^NL_sc = M^NL_uc/81, cohérent avec M/N_cells) et ⟨Φ|ΔV^NL|Φ′⟩ = − Σ_li D_li B_li(Φ) B*_li(Φ′). Vérification du chemin
(A.2) : pour des états de Bloch purs (c = δ), ⟨nk|ΔV^L|n′k′⟩ direct contre M^L/81 et ⟨nk|ΔV^NL|n′k′⟩ contre M^NL/81 sur quelques paires
(attendu ≤ 1e-6 eV ; même grille, mêmes coefficients), puis l'identité c†H_in c = ε_QE c†c − ⟨Ψ_in|ΔV|Ψ_out⟩ + constante d'alignement.
Mémoire : grille 14,0 M points → 224 Mo par tableau complexe (Ψ, Ψ_in), 112 Mo pour ΔV ; sphère complète 771 419 g : K̂ 18 Mo, Y 74 Mo ;
30 états × 2 bases : minutes.

### 0.5 Input nscf pour B et coûts

Fichier préparé (**non lancé**) : `b/nscf_nb128.in`, dérivé de `unit_cell/9x9/scf.in` ; diff complet (`b/nscf_nb128.diff`) :

```diff
--- unit_cell/9x9/scf.in	2026-06-18 13:53:53
+++ b/nscf_nb128.in	2026-09-25
@@ -1,10 +1,10 @@
-! ground state SCF of pristine graphene unit cell
-! 9x9 k grid
+! R5 (B) : nscf 128 bandes de la maille unitaire, grille 9x9 (81 k, nosym, noinv), dérivé de unit_cell/9x9/scf.in
+! outdir NEUF : copie de charge-density.hdf5 + data-file-schema.xml du .save de production (jamais écrit)
 
 &CONTROL
-  calculation = 'scf'       ! SCF ground state calculation
+  calculation = 'nscf'      ! R5 : non self-consistent sur la densité copiée
   prefix      = 'defect_unit_cell_9x9'        ! prefix of output files
-  outdir      = '/home/gregb26/links/scratch/qe_tmp/defect_unit_cell_9x9'   ! directory of tempory files
+  outdir      = '/home/gregb26/links/scratch/qe_tmp/R5_uc9x9_nb128'   ! R5 : outdir neuf (copie du .save)
   pseudo_dir  = '/home/gregb26/links/projects/rrg-cotemich-ac/gregb26/abinit_processing/pseudo' ! directory of pseudo potentials
 /
 &SYSTEM
@@ -15,7 +15,7 @@
   occupations = 'smearing' 
   smearing    = 'mv'     ! Marzari-Vanderbilt cold smearing
   degauss     = 0.01    ! value of gaussian smearing in Ry
-  nbnd        = 16       ! number of bands
+  nbnd        = 128      ! R5 : 128 bandes
   assume_isolated = '2D' ! helps for 2d systems
   nosym = .true.
   noinv = .true. 
@@ -24,6 +24,8 @@
   electron_maxstep =  400   ! number of iterations in SCF step
   conv_thr         =  1.d-14 ! convergence threshold for self-consistency
   mixing_beta      =  0.3
+  diago_full_acc   = .true.   ! R5 : bandes vides convergées comme les occupées
+  diago_thr_init   = 1.0d-12  ! R5 : seuil explicite (défaut nscf = conv_thr/nelec/10 = 1.25e-16, inatteignable)
 /
 
 ATOMIC_SPECIES
```

Cinq changements (les quatre du prompt + `calculation`, inévitable puisque la production est un scf) ; prefix, pseudo_dir, cellule,
positions, ecutwfc, smearing, nosym/noinv, assume_isolated, conv_thr et `K_POINTS automatic 9 9 1 0 0 0` inchangés (la grille automatique
en nscf avec nosym + noinv regénère les 81 k dans l'ordre du scf ; B.1 le vérifie explicitement). `diago_thr_init = 1e-12 Ry` : le défaut
nscf (conv_thr/nelec/10 = 1,25e-16) est inatteignable et a produit les « eigenvalues not converged » du 27×27 ; 1e-12 Ry ≪ 1e-6 eV demandé.
`b/submit_nscf_nb128.sh` (non lancé, GO « nscf » requis) : 1 nœud, 64 rangs MPI, 2 G/cœur, 1 h ; refuse si l'outdir contient déjà des wfc ;
copie `charge-density.hdf5` + `data-file-schema.xml` du `.save` de production vers `qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.save/`
(md5 imprimés), puis `srun pw.x < nscf_nb128.in > nscf_nb128.out`. Le `.save` source n'est jamais ouvert en écriture.

Coûts estimés :

| étape | référence mesurée | loi | estimation |
|---|---|---|---|
| nscf 128 bandes × 81 k | 27×27 : 729 k × 20 bandes, 207 s de mur, 64 rangs (0,0142 s par k·bande) | ∝ nbnd … ∝ nbnd² | 2,5 min … 16 min (64 rangs) ; demande 1 h |
| M^L grossier 128 bandes (`compute_ML_R_mpi_shared`, 32 rangs × 6 fils, 1 nœud) | coarse check 16 bandes (Bk = 1 296) : 71 s, MaxRSS 4,9 G (job 20203627) ; dense 20 bandes × 729 k (Bk = 14 580) : 53 min, MaxRSS 50 G (job 20212562) | ∝ Bk² (Bk = 10 368) | 64 × 71 s = 76 min (borne haute) ; (10 368/14 580)² × 53 min = 27 min ; mémoire : u_nk partagé 128 × 81 × 172 800 × 16 o = 28,7 G + 32 × M_local 1,72 G = 55 G + blocs Ψ ≈ 90 G |
| M^NL grossier 128 bandes (`compute_M_NL`, processus frais) | R4 J1 : 16 bandes, 172 s, 16 fils | ∝ nbnd (B), ∝ nbnd² (produit final, négligeable) | ≈ 25 min ; mémoire ≤ 8 × celle de R4 |
| fichiers | (128 × 81)² × 16 o | — | 1,72 Go par M (L, NL, tot), dans le répertoire de travail, hors dépôt |

Remarque : `scripts/compute_M.py --stage ml` (juin, `compute_ML_R_mpi`) a subi des OOM/corruptions à nk ≥ 81 (logs `M_ed_4_*`) ; la
chaîne à utiliser est le noyau R partagé de septembre (`compute_ML_R_mpi_shared`, validé = juin à 1,3e-15 par R4 J1) puis `compute_M_NL`
en processus frais, puis somme ; les scripts de production codent les liens `data/` en dur, le pilote appelle donc les routines de
bibliothèque avec le chemin du `.save` neuf (aucun script de production édité).

### 0.6 Plan d'exécution

Problèmes relevés en phase 0 (règle « deux sortes ») :

| n° | problème | sorte | traitement |
|---|---|---|---|
| 1 | w₂ de (a2)–(c) de R4 faux : index d'ondes planes construit avec k81 au lieu de k27 (0.3) | exécution (R4) | corrigé dans le pilote R5 + garde-fou ; A.4 refait le tableau D4 ; R4 non modifié |
| 2 | pas d'input nscf de production pour la 9×9 | écart | dérivation depuis `scf.in`, un changement de plus (`calculation`), diff rapporté |
| 3 | `.save` de la maille 11×11 = run bands (grille 121 k écrasée) | écart | E_D(11) par la route « maille » prise sur le chemin (K exact à l'indice 61) ; la route « repliement 121 k » n'existe pas pour N = 11 ; rapporté |
| 4 | `PYTHONPATH=src` écrase le chemin `EBPYTHONPREFIXES` (h5py invisible) | exécution | `sys.path.insert` dans le pilote (comme R4) |

Définitions retenues (à confirmer au GO) :

- Fenêtre ε − E_D ∈ [−3, +1] eV ; parité ⟨σ_h⟩ (z₀ = 0) ; w₂ / w₁ = fraction de |Ψ|² (sommée sur z) dans le disque de 2,0 / 1,0 Å ;
  seuil de localisation : 3 × ⟨w₂⟩ de la parfaite de **même taille** (R4 : 0,0812 pour la 9×9 ; la fraction d'aire du disque varie en 1/N²).
- Alignement : Lu 1,0 Å (`alignment.far_atom_alignment`) avec `Vks_NxN_{d,p}` ; 9×9 : −24,65 meV (R4).
- E_D : 9×9 = quadruplet de la parfaite (−4,23847 eV). Pour C, N = 3m (6, 9, 12) : quadruplet de la parfaite N×N ; N ≠ 3m (5, 7, 8, 10, 11),
  K n'est pas sur la grille N×N : E_D(N) := E_K(maille 9×9, −4,238470) + Δ_Γ + Δ_fold, avec Δ_Γ = médiane sur les 8 états les plus bas à Γ
  de [ε(maille N×N) − ε(maille 9×9)] (Γ est sur toutes les grilles ; deux SCF de la même maille ne diffèrent que par la grille k) et
  Δ_fold = médiane sur les états sous E_D − 4 eV de [ε(parfaite N×N) − ε(maille N×N repliée)] (R4 : 1,3e-6 eV pour la 9×9 dans la
  fenêtre). Pour N = 11, Δ_fold est remplacé par l'alignement du seul état le plus bas (Γ, bande 1) et E_K est lu directement sur le
  chemin du run bands ; pour N = 6, 9, 12 les deux routes (quadruplet, maille) sont rapportées. c = 1/(2N²).
- A.1 : c_nk = ⟨nk|Ψ⟩_sc (0.2) ; δ = 1 − Σ|c|² ; R = c†H_in c / c†c ; répartition sur les 5 états propres de H_in de plus grand
  |⟨φ_j|c⟩|²/c†c ; résidu ‖(H_in − ε_QE) c‖/‖c‖ ; par parité (H_in diagonalisé par bloc, comme R4) ; ε_QE aligné Lu, rapporté à
  E_D de la base (16 bandes : −4,238470 ; 20 bandes : −4,238895 ; R4 §D4).
- A.2 : parties L et NL séparées (`M_coarse_L`, `M_coarse_NL` ; `M_L_dense`, `M_NL_dense` restreints) ; identité c†H_in c = ε_QE c†c
  − ⟨Ψ_in|ΔV|Ψ_out⟩ + s c†c, où s (constante d'alignement : E_D de la base − E_D(SC) + décalage Lu, et moyenne G = 0 de ΔV non
  soustraite) est calculé et rapporté, pas ajusté.
- A.3 : comptages sous E_D − 3 eV par parité (QE défaut : 352 états ; (a1) : 1 296 ; (a3) : 1 620) ; spectres impairs triés sur [−4, +1]
  côte à côte ; min_k [λ_k(a1) − λ_{k−s}(QE)] pour s = 0, 1, 2 ; idem pairs.
- A.5 : M^L − 76,58 meV·𝕀 sur la diagonale (a1) ; attendu : décalage rigide −76,58 meV ; rapporté max|Δε| des états localisés après
  retrait du décalage ajusté (`alignment.rigid_shift_fit`).
- B.3 : sous-blocs n ∈ {16, 24, 32, 48, 64, 96, 128} du M à 128 bandes ; ψ_nk et parité du `.save` neuf ; A.1 et A.3 à chaque n ;
  w₂ avec les ψ_nk du `.save` neuf (jauge cohérente avec le M à 128 bandes ; le garde-fou de 0.3 s'applique : k du `.save` neuf).
- Coût 15×15 / 18×18 (C) : série scf 5…12 à 64 rangs, ajustement t ∝ N^4,96 (défaut) / N^4,59 (parfaite) → 15×15 ≈ 44 / 31 min,
  18×18 ≈ 110 / 71 min ; à partir de la 9×9 seule (200 s) avec l'exposant 5 : 42 / 104 min, avec l'exposant 3 (N_at^1,5) : 15 / 27 min ;
  `wfc1.hdf5` ∝ N⁴ : 16 Go (15×15), 34 Go (18×18) par super-cellule. À finaliser dans C.

Jobs (rrg-cotemich-ac, 1 nœud, `module restore qe; module load mpi4py/4.0.3 scipy-stack`, `.venv/bin/python`, sorties ici) :

| job | sous-commande du pilote | contenu | ressources | dépendance | écart à `config/production.json` |
|---|---|---|---|---|---|
| J1 `r5a` | `a` | A.4 (vérification 0.3 : w₂(a2) = w₂(a1) à 1e-3, tableau D4 refait), A.1, A.2 (vérification états purs ≤ 1e-6 eV puis les 30 états, L / NL), A.3, A.5 ; figures | 16 cœurs, 64 G, 1 h | — | aucun paramètre de production utilisé (pas de matrice T) |
| J2 `r5c` | `c` | C : protocole D1 pour N = 5…12 (valeurs propres, parité, w₂, w₁, Lu 1,0 Å, E_D deux routes, comptages), ajustements 1/N et 1/N², figure size, coût 15×15 / 18×18 | 16 cœurs, 64 G, 2 h (12×12 : grille 25 M points, wfc 6,6 Go lus par tranche) | — (parallèle à J1) | fenêtre [−3, +1] |
| J3 `r5nscf` | `b/submit_nscf_nb128.sh` (pw.x) | B.1 : nscf 128 bandes dans l'outdir neuf | 64 rangs MPI, 2 G/cœur, 1 h | **GO nscf** | calcul QE unique autorisé |
| J4 `r5m128` | `b2ml` (srun, MPI) → `b2nl` (processus frais) → `b2sum` | B.2 : M^L, M^NL, M à 128 bandes + sidecars (unit_cell, hartree, nbnd 128, save = outdir neuf) ; vérifications (bloc 16 × 16 contre `M_ed_9x9` par valeurs singulières, hermiticité) | 1 nœud exclusif (`--mem=0`), 32 rangs × 6 fils, 3 h | afterok J3 | — |
| J5 `r5b` | `b1 b3` | B.1 (16 premières bandes contre le `.save` grossier ≤ 1e-6 eV, ordre des k, bande 128 à Γ/K/M, parité 128 × 81), B.3 (n ∈ {16…128}), B.4 ; figures delta_vs_n, ladder_bands | 16 cœurs, 128 G (H à 128 bandes : 10 368², 1,7 Go ; eigh par bloc de parité ≈ 5 184²) | afterok J4 | — |

Ordre : J1 et J2 après le GO (indépendants) ; J3 après le GO nscf ; J4, J5 enchaînés par afterok. Aucun job n'est relancé
automatiquement ; un échec de J3 laisse l'outdir neuf en l'état (jamais supprimé sans GO).

Pilote unique `r5_driver.py` (sous-commandes `a`, `c`, `b1`, `b2ml`, `b2nl`, `b2sum`, `b3`, `tables`, `figs` ; `sys.path.insert(0, src)` ;
chargement de `config/production.json` pour E_D Wannier et K ; journal `r5_log.txt` ; résultats `a/`, `c/`, `b/` (json + npz), `fig/`),
`submit_r5.sh` (16 cœurs / 64–128 G, une soumission par sous-commande) et `submit_r5_m128.sh` (nœud exclusif, srun pour `b2ml`).
Aucun script jetable dans le dépôt ; aucun fichier de production modifié.

Fonctions à ajouter (modules dans `src/electron_defect_interaction/` ; les modules R4 sont non commités et peuvent recevoir des ajouts) :

| module | signature | rôle |
|---|---|---|
| `wavefunctions/sc_projection.py` (nouveau) | `sc_state_grid(C, mill, ngfft, gamma_only) -> A (n1, n2, n3)` | coefficients d'un état de super-cellule sur la grille FFT, sphère complétée (A[−g] = A*[g]) |
| idem | `grid_to_real(A, Omega_sc) -> psi (n1, n2, n3)` | Ψ(r) = N_grid · ifftn(A)/√Ω_sc, ∫|Ψ|² dV = 1 |
| idem | `bloch_overlaps(A, C_nkg, nG, flat_idx) -> c (nb, nk)` | c_nk = Σ_G C*_nk(G) A[flat(k, G)] |
| idem | `bloch_superposition_grid(d_nk, C_nkg, nG, flat_idx, ngfft) -> A` | Σ_nk d_nk C_nk(G) sur la grille (Ψ_in), sans somme sur z |
| idem | `sc_full_sphere(C, mill) -> (D, mill_full)` | liste complète (771 419 g) pour la partie non locale |
| `defects/deltav_pw.py` (nouveau) | `expect_local(A_bra, A_ket, dV_grid, Omega_sc) -> complex` | ⟨Φ|ΔV^L|Φ′⟩ = Σ_r Φ* ΔV Φ′ dV (noyau de `compute_ML_R`) |
| idem | `removed_atom_projections(D, mill_full, B_sc, tau_vac, upf, Omega_sc) -> B (l, i, m)` | B_lim(Φ) avec `build_K_vectors`, `compute_phase`, `compute_angular_part`, `fq_from_fr` (nkpt = 1, préfacteur 4π/√Ω_sc) |
| idem | `expect_nonlocal(B_bra, B_ket, ekb_li) -> complex` | ⟨Φ|ΔV^NL|Φ′⟩ = − Σ_li D_li B_li(Φ) B*_li(Φ′) |
| idem | `check_pure_bloch(M_L, M_NL, N_cells, pairs, …) -> dict` | vérification états purs contre M/N_cells (≤ 1e-6 eV) |
| `wannier/supercell_fold.py` (R4) | `check_k_reference(k_used, k_save, tol=1e-6)` | refuse un k décalé d'un vecteur entier (bug 0.3) |
| `defects/basis_diagnostics.py` (nouveau) | `basis_diagnostics(H_blocks, c, e_ref) -> dict(delta, rayleigh, residual, top5)` | A.1 par bloc de parité |
| idem | `interlacing(eps_ref, eps_model, shifts=(0, 1, 2), lo, hi) -> dict` | A.3 : comptages, min_k [λ_k − λ_{k−s}] |
| idem | `band_subblock(M, eps, C, n) -> (M_n, eps_n, C_n)` | B.3 : sous-blocs n bandes du M à 128 |
| `defects/alignment.py` (R4) | `dirac_quadruplet(eps, ef, n=4) -> (E_D, idx, spread)` | E_D des parfaites N = 3m |
| idem | `gamma_shift(eps_a, eps_b, nlow=8) -> (shift, resid)` ; `lowest_state_shift(eps_a, eps_b)` | Δ_Γ et Δ_fold (route « maille ») |
| idem | `size_fits(N, eps, families, powers=(1, 2)) -> dict(eps_inf, a, resid)` | ajustements ε_∞ + a/N^p par famille |
| `io/qe_gamma_io.py` (R4) | `read_wfc_gamma_ecut_check(save) -> igwx` (optionnel) | contrôle 771 419 = Σ igwx (0.2) |

Tests unitaires prévus au GO (nœud de connexion, secondes) : `bloch_overlaps` d'un état de Bloch pur = δ_nk et Σ|c|² = 1 ;
`check_k_reference` refuse les 56 k décalés ; `grid_to_real` normé à 1e-12 ; `expect_nonlocal` d'un état pur contre `compute_M_NL` sur une
maille 1×1 fictive (ou report de la vérification à J1 avec la vraie 9×9) ; `size_fits` sur une loi exacte.

**STOP — phase 0 terminée le 2026-09-25 ; attente du GO (couvre J1, J2 ; J3 = nscf demande le GO « nscf » séparé ; J4, J5 suivent J3).**

## Phase GO — Méthodes (GO du 2026-09-25 ; J1–J5)

Précisions du GO appliquées : erratum daté ajouté à la fin de `article/R4_quasi_lie/R4_rapport.md` (et de la copie de travail,
identiques) ; A.2 : états purs avant les 30 états, STOP A.2 si échec ; C : gap à Γ de la parfaite (comptage d'électrons, n_occ = nelec/2),
ħv_F·min|k − K| de la grille N×N, π dans le gap pour N ≠ 3m, ajustements par famille (3m en tête, N ≠ 3m à part), coût 15×15 / 18×18 sans
rien lancer ; B : |Δε_π| entre n successifs, premier n d'un état pair localisé à moins de 0,3 eV de +0,101 ; `.save` 128 bandes gardé,
manifeste en fin de rapport.

### Code (répertoire de travail, hors dépôt — Greg a refusé l'écriture dans `src/`)

Modules de campagne : `r5_sc_projection.py` (sphère complète gamma_only, mise sur grille, Ψ(r) normé, recouvrements c_nk, superposition
Ψ_in, contrôle de l'union des ondes planes), `r5_deltav_pw.py` (⟨Φ|ΔV^L|Φ′⟩ sur la grille = noyau de `compute_ML_R` ; projecteur KB de
l'atome retiré sur la sphère complète de la super-cellule = noyau de `compute_M_NL` avec Ω_sc, k = 0, un atome ; vérification sur états
purs), `r5_basis_diagnostics.py` (δ, quotient de Rayleigh, résidu, poids sur les états propres ; entrelacement), `r5_alignment_ext.py`
(`check_k_reference`, quadruplet, décalages à Γ et par l'état le plus bas, ajustements 1/N^p). Pilote `r5_driver.py` (sous-commandes `a`,
`c`, `b1`, `b2ml`, `b2nl`, `b2sum`, `b3`, `atables`, `ctables`, `btables`) ; il importe les aides du pilote R4 (`r4_driver` : `setup`,
`read_wfc_subset`, `uc_parity`, `spectrum_blocks`, caches M_W / V_loc) en redirigeant leur journal vers `r5_log.txt` (R4 non modifié) ;
`submit_r5.sh` (16 cœurs, 96 G), `submit_r5_m128.sh` (nœud exclusif, 32 rangs × 6 fils), `b/submit_nscf_nb128.sh` (64 rangs).
Routines de production réutilisées telles quelles : `qe_io`, `matrix_io`, `config`, `Mbk_to_Mwk`, `Hwr_to_Hwk`, `compute_ML_R_mpi_shared`,
`compute_M_NL`, `read_upf`, `build_K_vectors`, `compute_phase`, `compute_angular_part`, `fq_from_fr`. Aucun fichier de production modifié.

Tests sur le nœud de connexion (2026-09-25) : union des 81 sphères de maille = sphère complète de la 9×9 (771 419 indices distincts) ;
état de Bloch pur → c = δ à 1e-15 et norme 1 sur la grille ; ∫|Ψ|² dV = 1 à 1e-15 ; `check_k_reference` refuse les 56 k décalés du dense
et accepte les k du `.save` neuf (max |Δk| 9e-15 contre la maille grossière) ; ⟨nk|ΔV^NL|n′k′⟩ direct = M^NL/81 à 1e-15 eV sur six paires ;
`size_fits` exact sur une loi 1/N ; `basis_diagnostics` : état propre exact → δ = 6e-16, résidu 4e-15, poids 1. (Le FFT multi-fils est
refusé sur le nœud de connexion, « Resource temporarily unavailable » ; sans effet dans les jobs.)

### Jobs (rrg-cotemich-ac, 1 nœud, `module restore qe; module load mpi4py/4.0.3 scipy-stack`, `.venv/bin/python`)

| job | contenu | ressources | identifiant | durée / MaxRSS |
|---|---|---|---|---|
| J3 | nscf pw.x 128 bandes, grille 9×9, outdir neuf `qe_tmp/R5_uc9x9_nb128/` (copie de charge-density + XML, md5 imprimés) | 64 rangs, 2 G/cœur | 21808943 | 117 s de mur (2 min 21 de job), 0,89 G par rang ; ethr 1e-12, 34,3 itérations Davidson en moyenne, 0 « not converged » ; `.save` 1,5 Go + 64 fichiers `.wfc*` temporaires (1,4 Go) |
| J4 | `b2ml` (srun MPI) → `b2nl` → `b2sum` : M^L, M^NL, M à 128 bandes | nœud exclusif, 32 × 6 | 21810702 | voir §B |
| J1 | `a` : A.1–A.5 | 16 cœurs, 96 G | 21810853 (plantage dans les tables), 21811069 (calcul complet, 127 s ; tables et figures refaites par `atables` depuis le json) | ≈ 3 min |
| J2 | `c` : C | 16 cœurs, 96 G | 21810854, 21811070 (plantages dans les tables), 21811149 (complet, 2 min 10), puis resoumis pour E_D(11) (§C) | ≈ 3 min |
| J5 | `b1 b3` : B.1, B.3, B.4 | 16 cœurs, 96 G | 21810855 (afterok J4) | voir §B |

Problèmes rencontrés (règle « deux sortes ») :

| n° | problème | sorte | traitement |
|---|---|---|---|
| 1 | test unitaire : ⟨nk|ΔV^NL|n′k′⟩ direct 27,2 fois trop petit (D_li de `read_upf` en Hartree, M en eV) | exécution | facteur `energy_scale = HA2EV` dans le projecteur ; accord 1e-15 ensuite |
| 2 | J1, J2 : plantages dans l'écriture des tables (complexes / clés entières relues du json) | exécution | corrigés ; sous-commandes `atables`, `ctables`, `btables` (tables depuis les json) |
| 3 | C : gap à Γ défini par E_F comptait le quadruplet entier comme occupé (N = 3m : gap 1,7 eV au lieu de 0) | exécution | n_occ = nelec/2 |
| 4 | C, N = 11 : le seul `.save` de maille est le run bands, dont la référence d'énergie diffère de +2,37 eV (Δ_Γ +2 373,7 meV, Δ_fold −2 384,6 meV) ; la route E_K9 + Δ_Γ + Δ_fold laisse un résidu de −10,9 meV | écart de données | E_D(11) pris sur le point K du même run (indice 61) + Δ_fold du même run (les références s'annulent) ; les deux valeurs sont rapportées |
| 5 | **A.2 : ⟨nk|ΔV^L|n′k′⟩ direct = 81,0 × M^L/81 sur toutes les paires, NL exact** (§A.2) | **écart avec la production** | STOP A.2 (convention de production) ; variante diagnostique « M^L × 81 » ajoutée à A.1, A.2, A.4 et B.3, étiquetée hors prompt ; production non touchée |

### Définitions

Comme en phase 0 (§0.6), plus : ε_QE sur l'échelle d'une base = e_QE − décalage Lu − E_D(SC) + E_D(base) ; constante d'alignement s de A.2
définie par c†H_in c = (e_QE − Lu) c†c − ⟨Ψ_in|ΔV|Ψ_out⟩ + s c†c (Ψ_in = Σ c_nk |nk⟩_sc, Ψ_out = Ψ − Ψ_in, ΔV = V_d − V_p de pp.x sans
soustraction, ΔV^NL = − projecteurs KB de l'atome retiré) ; « M^L × 81 » = M^L multiplié par N_cells avant la somme avec M^NL.

## A — Base ou M, avec les données existantes (J1, 21811069 ; 127 s)

Tables complètes : `a/A_tables.md` ; résultats `a/a_results.json`, `a/a_overlaps.npz` (c_nk des 30 états sur les bases 16 et 20),
`a/a_qe_all.npz` (parité des 352 états QE), `a/a_spectra.npz` ; figures `fig/a_delta`, `fig/a_ladder_corrected`.

Bases : 16 bandes = maille grossière (`defect_unit_cell_9x9.save`, M_ed de juin, E_D(K) = −4,238470) ; 20 bandes = dense restreint aux 81 k
(`defect_uc_dense_27.save`, M dense restreint, E_D(K) = −4,238895). Union des ondes planes repliées = sphère de la super-cellule (771 419)
pour les deux bases ; M_tot − M_L − M_NL = 1,3e-15 ; couplage pair–impair de H_in : 9,8e-4 (16), 8,0e-6 (20) eV. Les 30 états QE
(bandes 297–326, gamma_only) ont une norme 1 à 1e-12 et la même parité qu'en R4 (dev 0). Bug 0.3 confirmé : 56 des 81 index d'ondes planes
diffèrent entre k81 et k27[idx81] ; `check_k_reference` refuse (« 56 of 81 k differ … by a non-zero integer vector »).

### A.2 — vérification sur états purs : **ÉCHEC pour la convention de production** (STOP A.2), facteur 81,0 sur M^L

⟨nk|ΔV|n′k′⟩_sc calculé directement (états normés sur la super-cellule) contre M/81, six paires (16 bandes) et sept (20 bandes) dont
intra-k (3, K, 3, K), (3, K, 4, K) et inter-k (3, K, 3, K′), (2, 5, 7, 40), (10, 17, 1, 60) :

| partie | base 16 : max\|Δ\| | base 20 : max\|Δ\| | rapport direct / (M/81) |
|---|---|---|---|
| non locale (M^NL) | 3,4e-15 eV | 1,5e-9 eV (bruit 1,3e-7 des k du XML dense) | 1,0000 |
| locale (M^L) | 1,7e-1 eV | 1,5e-1 eV | **81,0000** sur toutes les paires non nulles (ex. (3, K, 3, K) : direct +0,155642, M^L/81 = +0,001922, M^L du fichier = +0,155642) |

Le même test avec M^L × 81 : max|Δ| = 6,8e-15 (16) et 3,5e-9 eV (20), L et NL. **Constat** : dans `M_ed_9x9` (juin) comme dans
`M_dense_9x9` (septembre), M^L est l'élément de matrice entre états normés sur la super-cellule (le noyau `compute_ML_R*` construit
Ψ = u_nk e^{ik·r}/√Ω_sc, de norme 1 sur la super-cellule) tandis que M^NL est l'élément entre états normés sur la maille (préfacteur
4π/√Ω_uc de `compute_M_NL`) ; les sidecars déclarent `unit_cell` pour les deux. Les deux parties de M diffèrent donc d'un facteur
N_cells = 81 l'une par rapport à l'autre. Conséquence directe (pas d'interprétation) : dans H_(a1) = diag ε + M/81, la partie locale
pèse 1/81 de ce que la partie non locale suppose ; (a1) total ≈ (a1) M^NL seul (R4 : −1,350 contre −1,370 ; −2,812 contre −2,847).
Aucun fichier de production n'est touché ; la suite de A rapporte la convention de production **et** la variante diagnostique « M^L × 81 »
(M^L × N_cells + M^NL, les deux parties alors en norme `unit_cell`), étiquetée hors prompt.

Les 30 états (convention diagnostique seulement, la production étant arrêtée à A.2) : (i) c†(M/81)c et (ii) application directe de ΔV à
Ψ_in coïncident à 2,7e-15 (L) et 3,6e-14 eV (NL) sur la base 16, 6,7e-9 / 8,6e-9 sur la base 20 ; l'identité
c†H_in c = (e_QE − Lu) c†c − ⟨Ψ_in|ΔV|Ψ_out⟩ + s c†c donne **s = +24,66 meV** pour les 30 états (min +24,66, max +24,67 ; base 20 :
+25,06 à +25,17), soit s = −(décalage Lu) à 0,01 meV près (Lu = −24,65 meV) : avec e_QE non aligné, la constante est nulle
(⟨ΔV⟩_3D = +9,09 meV n'y entre pas). États en évidence (base 16 ; eV) :

| bande | parité | ε − E_D QE | Σ\|c\|² | c†(M/81)c tot (L ; NL) | ⟨Ψ_in\|ΔV^L\|Ψ_out⟩ | ⟨Ψ_in\|ΔV^NL\|Ψ_out⟩ | c†H_in c | (e_QE − Lu) c†c − ⟨in\|ΔV\|out⟩ |
|---|---|---|---|---|---|---|---|---|
| 318, 319 | π | −1,759 | 0,9998 | +0,0284 (0,0284 ; 0,0000) | −0,0091 | +0,0000 | −6,0124 | −5,9878 |
| 320 | π | −0,737 | 0,9991 | +0,3666 (0,3383 ; 0,0283) | −0,0382 | −0,0151 | −4,9426 | −4,9180 |
| 324, 325 | σ | +0,101 | 0,9638 | +5,8490 (5,0523 ; 0,7967) | −1,5566 | −0,5255 | −1,9290 | −1,9052 |
| 326 | π | +0,269 | 0,9995 | +0,1480 (0,1338 ; 0,0143) | −0,0195 | −0,0077 | −3,9647 | −3,9401 |

### A.1 — recouvrements, quotient de Rayleigh, résidu, répartition

δ = 1 − Σ|c_nk|² ; fuite vers l'autre parité ≤ 1e-9 ; R − ε_QE et résidu ‖(H_in − ε_QE) c‖/‖c‖ sur l'échelle de la base ; « top » = poids
|⟨φ_j|c⟩|²/c†c sur les états propres de H_in (ε_j − E_D). Convention de production (« prod ») et diagnostique (« M^L × 81 ») :

| bande | parité | ε − E_D QE | base | δ | prod : R − ε (meV) ; résidu (eV) ; premiers poids | M^L × 81 : R − ε (meV) ; résidu (eV) ; premiers poids |
|---|---|---|---|---|---|---|
| 324, 325 | σ | +0,101 | 16 | 0,0362 | −3 042 ; 5,04 ; (−2,812 ; 0,31), (−2,812 ; 0,21), (−3,294 ; 0,11) | +2 136 ; 5,65 ; (+1,481 ; 0,63), (+1,480 ; 0,34), (+37,7 ; 0,01) |
| 324, 325 | σ | +0,101 | 20 | 0,0241 | −2 986 ; 5,49 ; (−2,848 ; 0,44), (−3,314 ; 0,19), (−3,967 ; 0,08) | +1 410 ; 4,61 ; (+1,068 ; 0,97), (+1,068 ; 0,02), (+45,0 ; 0,01) |
| 320 | π | −0,737 | 16 | 0,0009 | −306 ; 1,82 ; (−1,350 ; 0,81), (+0,108 ; 0,09), (−2,654 ; 0,02) | +29 ; 0,86 ; (**−0,727 ; 0,999**), (+0,263 ; 0,001) |
| 320 | π | −0,737 | 20 | 0,0007 | −305 ; 1,87 ; (−1,359 ; 0,80), (+0,106 ; 0,09), (−2,656 ; 0,02) | +16 ; 0,76 ; (**−0,734 ; 0,999**), (+0,259 ; 0,000) |
| 318, 319 | π | −1,759 | 16 | 0,0002 | −44 ; 0,20 ; (−1,803 ; 0,86), (−1,803 ; 0,13) | −16 ; 0,05 ; (−1,775 ; 0,85), (−1,775 ; 0,15) |
| 318, 319 | π | −1,759 | 20 | 0,0002 | −43 ; 0,22 ; (−1,803 ; 0,95), (−1,803 ; 0,04) | −16 ; 0,05 ; (−1,776 ; 0,81), (−1,776 ; 0,19) |
| 326 | π | +0,269 | 16 | 0,0005 | −130 ; 1,25 ; (+0,108 ; 0,90), (−1,350 ; 0,05) | +3 ; 0,62 ; (**+0,263 ; 0,999**) |
| 326 | π | +0,269 | 20 | 0,0003 | −129 ; 1,28 ; (+0,106 ; 0,89), (−1,359 ; 0,06) | −4 ; 0,55 ; (**+0,259 ; 0,999**) |

Les 30 états : δ ≤ 0,0009 pour les 28 états π (base 16) et ≤ 0,0007 (base 20) ; δ = 0,036 / 0,024 pour la paire σ (`fig/a_delta`).
Avec M^L × 81, le résidu des états π reste de 0,05 à 0,9 eV bien que le poids sur un seul état propre soit 0,999 : la partie hors base
(δ ~ 1e-3) est couplée par ΔV (voir ⟨Ψ_in|ΔV|Ψ_out⟩ ci-dessus). Tableau complet des 30 états (deux conventions, deux bases) dans `a/A_tables.md`.

### A.3 — entrelacement (énergies ε − E_D ; QE aligné Lu ; modèles sur leur E_D)

QE défaut : 352 états, 243 pairs, 109 impairs (parité exacte à 7e-7). Sous E_D − 3 eV : QE 241 pairs et 55 impairs ; (a1) 16 bandes
241 / 55 ; (a3) 20 bandes 241 / 55 (mêmes comptages). Fenêtre [−3, +1] : 2 pairs + 28 impairs dans QE, (a1) et (a3).
min_k [λ_k(modèle) − λ_{k−s}(QE)] (fenêtre ; tout le spectre commun) :

| modèle | parité | s = 0 | s = 1 | s = 2 |
|---|---|---|---|---|
| (a1) 16 | σ | −2,913 ; −2,913 | −2,913 ; −2,913 | +0,321 ; −0,226 |
| (a1) 16 | π | −0,612 ; −0,612 | −0,044 ; −0,057 | −0,005 ; −0,007 |
| (a3) 20 | σ | −2,950 ; −2,950 | −2,950 ; −2,950 | +0,285 ; −0,239 |
| (a3) 20 | π | −0,621 ; −0,621 | −0,044 ; −0,057 | −0,005 ; −0,007 |

Spectres impairs triés côte à côte sur [−4, +1] (28 états chacun dans [−3, +1] ; listes complètes dans `a/A_tables.md`) : QE
−2,793, −2,793, −2,790, −2,777, −2,777, −2,535, −2,105, −2,105, −2,103, −2,086, −2,086, −2,001, −1,815, −1,811, −1,811, −1,808, −1,807,
−1,807, −1,798, −1,798, −1,798, −1,759, −1,759, −0,737, +0,003, +0,015, +0,015, +0,269 ; (a1) −2,790 ×5, −2,654, −2,100 ×5, −2,021,
−1,804 ×9, −1,803 ×2, −1,350, −0,000 ×3, +0,108.

### A.4 — w₂ corrigés (bug 0.3) : tableau D4 refait

Vérification : w₂(a2) = w₂(a1) sur les cinq états localisés de (a1) : max|Δw₂| = **2,7e-4** (0,2871 → 0,2874 ×2, 0,1009 → 0,1009 ×2,
0,2260 → 0,2260) → OK. Tableau (ε − E_D ; w₂ corrigé ; w₂ de R4 entre parenthèses) :

| variante | pairs (σ) localisés | impairs (π) localisés |
|---|---|---|
| QE (D1) | +0,101 ×2 (0,713) | −1,759 ×2 (0,114) ; −0,737 (0,246) ; +0,269 (0,102) |
| (a1) tot | −2,812 ×2 (0,287) | −1,803 ×2 (0,101) ; −1,350 (0,226) |
| (a1) M^L seul | aucun | −1,803 ×2 (0,101) ; −1,780 (0,103) |
| (a1) M^NL seul | −2,847 ×2 (0,255) | −1,370 (0,222) |
| (a2) 16 b., M dense | −2,812 ×2 (0,287 ; R4 0,058 / 0,074) | −1,803 ×2 (0,101 ; 0,091) ; −1,350 (0,226 ; 0,084) |
| (a3) 20 b. | −2,848 ×2 (0,254 ; 0,067 / 0,055) | −1,803 ×2 (0,101 ; 0,091) ; −1,359 (0,225 ; 0,084) |
| (b) 5 WF | −2,514 ×2 (0,489 ; 0,106 / 0,088) | −1,803 ×2 (0,101 ; 0,091) ; −1,315 (0,232 ; 0,083) |
| (c-all) | −2,514 ×2 (0,489) | −1,803 ×2 (0,101) ; −1,315 (0,232) |
| (c-3) | −2,516 / −2,515 (0,490) | −1,803 ×2 (0,100 / 0,101) ; −1,312 (0,232) |
| **(a1) M^L × 81** (diagnostique) | aucun dans [−3, +1] (paire σ à +1,481 / +1,480, hors fenêtre ; poids 0,63 + 0,34 de l'état QE +0,101) | −1,775 ×2 (0,116) ; **−0,727 (0,236)** ; **+0,263 (0,107)** |
| **(a3) M^L × 81** | aucun (paire σ à +1,068 ; poids 0,97 + 0,02) | −1,776 ×2 (0,116) ; **−0,734 (0,238)** ; **+0,259 (0,106)** |

Écarts (a1) M^L × 81 contre QE, états π localisés : −0,727 contre −0,737 (+10 meV ; w₂ 0,236 contre 0,246), +0,263 contre +0,269 (−6 meV),
−1,775 contre −1,759 (−16 meV) ; base 20 : +3, −10, −17 meV. Paire σ : +1,48 (16 b.) / +1,07 (20 b.) contre +0,101, δ = 0,036 / 0,024.
Erratum ajouté à `R4_rapport.md` (deux copies) avec ce tableau (sans les lignes diagnostiques). Figure `fig/a_ladder_corrected`.

### A.5 — moyenne diagonale de M^L soustraite

Moyenne diagonale de M^L grossier recalculée : 76,58 meV (= R4 J1). Soustraite de la diagonale de M (norme `unit_cell`, base 16) :
décalage rigide observé −0,945 meV = −76,58/81 exactement (médian, min, max identiques), max|Δε − médiane| = 3,3e-13 eV, Δx des états
localisés après retrait du décalage ≤ 1e-7 meV. Lecture « ΔV − 76,58 meV » (constante × 81 sur la diagonale) : décalage rigide −76,579 meV,
max|Δε − médiane| = 2,1e-13 eV. Dans les deux cas : décalage rigide et rien d'autre, à la précision machine (la constante est proportionnelle
à l'identité dans la base repliée).

## C — Vérité DFT selon la taille, N = 5…12 (J2, 21811223 ; 2 min 20 ; sans QE)

Tables complètes : `c/C_tables.md` ; résultats `c/c_results.json` ; figure `fig/size`. Protocole D1 de R4 par taille : valeurs propres
et parité ⟨σ_h⟩ (z₀ = 0 ; déviation ≤ 1e-6 partout) des états gamma_only de la fenêtre ε − E_D ∈ [−3, +1], w₂ (disque de 2,0 Å) et w₁ (1,0 Å)
autour du site de la lacune, alignement Lu (1,0 Å, valeur de travail ; 0,5 Å rapporté) avec `Vks_NxN_{d,p}`, seuil de localisation
3 × ⟨w₂⟩ des états de la fenêtre de la parfaite N×N (fraction d'aire du disque de 2 Å : 9,6 % (5×5) … 1,7 % (12×12) ; seuils 0,262 … 0,046).
Références : E_K de la maille 9×9 = −4,238470 eV (bandes 3–4 à K dégénérées à 1,5e-8) ; ħv_F = 5,459 eV Å = (ε₅ − ε₄)/(2|δk|) au point
K + b₁/27 de la maille dense (|δk| = 0,1090 Å⁻¹).

E_D par taille : N = 3m → quadruplet de la parfaite (étalement ≤ 3,2e-7 eV) ; N ≠ 3m → route « maille » E_K9 + Δ_Γ + Δ_fold
(Δ_Γ = médiane sur les 8 états les plus bas à Γ de ε(maille N×N) − ε(maille 9×9) ; Δ_fold = médiane sur les états sous E_D − 4 eV de
ε(parfaite N×N) − ε(maille N×N repliée), résidu ≤ 1,2e-5 eV) ; N = 11 (seul `.save` de maille = run bands, référence d'énergie décalée de
+2,37 eV) → K du même run (indice 61) + Δ_fold par l'état le plus bas du même run ; la route E_K9 y donne −4,24936 (résidu de −10,9 meV entre
Δ_Γ +2 373,67 et Δ_fold −2 384,56 meV), la route K −4,23851 (retenue). Pour N = 6, 9, 12 les deux routes coïncident à 0,64, 0,00 et
0,12 meV du quadruplet.

| N | c = 1/(2N²) | E_D (eV) [source] | route maille : Δ_Γ ; Δ_fold (meV) | Lu 1,0 / 0,5 Å (meV) | ⟨V_d⟩ − ⟨V_p⟩ (meV) | seuil w₂ | état π localisé (ε − E_D ; w₂ ; w₁ ; bande) | doublet σ (ε − E_D ; w₂) | localisés σ / π | fenêtre P ; D (σ/π) |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0,02000 | −4,23138 [maille] | +7,08 ; +0,00 | −122,3 / −131,6 | +28,6 | 0,262 | −0,328 ; 0,267 ; 0,055 ; 98 | −0,131 ×2 ; 0,702 | 2 / 1 | 6 (0/6) ; 9 (2/7) |
| 6 | 0,01389 | −4,23745 [quadruplet] | +0,37 ; +0,00 (maille −4,23809) | −98,6 / −125,8 | +20,8 | 0,183 | −1,065 ; 0,284 ; 0,053 ; 140 | +0,156 ×2 ; 0,711 | 2 / 1 | 19 (0/19) ; 21 (2/19) |
| 7 | 0,01020 | −4,23904 [maille] | −0,59 ; +0,02 | −38,5 / −37,9 | +14,6 | 0,133 | −0,495 ; 0,301 ; 0,061 ; 194 | −0,282 ×2 ; 0,703 | 2 / 1 | 12 (0/12) ; 15 (2/13) |
| 8 | 0,00781 | −4,23768 [maille] | +0,79 ; −0,01 | −23,0 / −24,6 | +11,2 | 0,103 | −0,370 ; 0,223 ; 0,046 ; 254 | −0,195 ×2 ; 0,708 | 2 / 2 | 33 (0/33) ; 34 (2/32) |
| 9 | 0,00617 | −4,23847 [quadruplet] | +0,00 ; −0,00 (maille −4,23847) | −24,7 / −21,6 | +9,1 | 0,081 | −0,737 ; 0,246 ; 0,047 ; 320 | +0,101 ×2 ; 0,713 | 2 / 4 | 28 (0/28) ; 30 (2/28) |
| 10 | 0,00500 | −4,23889 [maille] | −0,42 ; −0,01 | −8,2 / −7,8 | +7,2 | 0,066 | −0,422 ; 0,241 ; 0,050 ; 398 | −0,238 ×2 ; 0,708 | 2 / 2 | 39 (0/39) ; 42 (2/40) |
| 11 | 0,00413 | −4,23851 [K du run bands] | (+2 373,67 ; −2 384,56 : −4,24936) | −12,2 / −12,1 | +5,9 | 0,054 | −0,332 ; 0,191 ; 0,039 ; 482 | −0,164 ×2 ; 0,707 / 0,708 | 2 / 2 | 42 (0/42) ; 45 (2/43) |
| 12 | 0,00347 | −4,23869 [quadruplet] | −0,09 ; −0,01 (maille −4,23857) | −29,2 / −36,1 | +5,1 | 0,046 | −0,551 ; 0,223 ; 0,044 ; 572 | +0,114 ×2 ; 0,713 | 2 / 9 | 55 (0/55) ; 57 (2/55) |

Gap à Γ de la parfaite N×N (n_occ = nelec/2 ; HOMO = état n_occ, LUMO = état n_occ + 1) et distance de la grille N×N à K :

| N | HOMO − E_D (eV) | LUMO − E_D (eV) | gap Γ (eV) | n_occ (états ≤ E_F) | min\|k − K\| (Å⁻¹) | ħv_F·min\|k − K\| (eV) | π localisé dans le gap ? | E_F p / d (eV) |
|---|---|---|---|---|---|---|---|---|
| 5 | −1,500 | +1,309 | 2,809 | 100 (100) | 0,340 | 1,855 | oui | −4,35505 / −4,59085 |
| 6 | 0 | 0 | 1,7e-7 | 144 (146) | 0 | 0 | — | −4,21850 / −4,23901 |
| 7 | −1,508 | +1,495 | 3,003 | 196 (196) | 0,243 | 1,325 | oui | −3,80934 / −4,67548 |
| 8 | −1,025 | +0,951 | 1,976 | 256 (256) | 0,212 | 1,159 | oui | −3,96547 / −4,54876 |
| 9 | 0 | 0 | 7,7e-8 | 324 (326) | 0 | 0 | — | −4,21952 / −4,22195 |
| 10 | −1,027 | +1,012 | 2,039 | 400 (400) | 0,170 | 0,927 | oui | −4,07501 / −4,58372 |
| 11 | −0,773 | +0,736 | 1,509 | 484 (484) | 0,154 | 0,843 | oui | −4,29446 / −4,50379 |
| 12 | 0 | 0 | 4,0e-8 | 576 (578) | 0 | 0 | — | −4,21973 / −4,21722 |

Pour N ≠ 3m, l'état π localisé et le doublet σ sont tous les deux dans le gap à Γ de la parfaite (σ : −0,131, −0,282, −0,195, −0,238,
−0,164 contre HOMO −1,50, −1,51, −1,02, −1,03, −0,77). Le gap à Γ vaut 1,5 à 1,6 fois ħv_F·min|k − K| (2,809 / 1,855 = 1,51 ; 3,003 / 1,325 = 2,27 ;
1,976 / 1,159 = 1,70 ; 2,039 / 0,927 = 2,20 ; 1,509 / 0,843 = 1,79).

Tous les états localisés (w₂ > seuil) : N = 5 : π −0,328 (0,267), σ −0,131 ×2 (0,702) ; 6 : π −1,065 (0,284), σ +0,156 ×2 (0,711) ;
7 : π −0,495 (0,301), σ −0,282 ×2 (0,703) ; 8 : π −1,598 (0,151), −0,370 (0,223), σ −0,195 ×2 (0,708) ; 9 : π −1,759 ×2 (0,114), −0,737 (0,246),
+0,269 (0,102), σ +0,101 ×2 (0,713) ; 10 : π −1,856 (0,091), −0,422 (0,241), σ −0,238 ×2 (0,708) ; 11 : π −1,166 (0,137), −0,332 (0,191),
σ −0,164 ×2 (0,707 / 0,708) ; 12 : π −2,621 ×2 (0,048), −2,438 ×2 (0,060), −1,963 (0,054), −1,314 ×2 (0,064), −0,551 (0,223), +0,186 (0,079),
σ +0,114 ×2 (0,713). w₂ du doublet σ : 0,702 à 0,713 pour les huit tailles ; w₁ 0,266 à 0,281.

Ajustements ε(N) = ε_∞ + a/N^p (eV ; état π = impair localisé de plus grand w₂ ; σ = moyenne du doublet) :

| état | famille | N | ε − E_D | p = 1 : ε_∞ ; a ; rms ; max résidu | p = 2 : ε_∞ ; a ; rms ; max résidu |
|---|---|---|---|---|---|
| π | **3m** (extrapolation de tête) | 6, 9, 12 | −1,065, −0,737, −0,551 | **−0,046 ; −6,137 ; 0,007 ; 0,010** | **−0,409 ; −23,94 ; 0,024 ; 0,033** |
| σ | 3m | 6, 9, 12 | +0,156, +0,101, +0,114 | +0,055 ; +0,572 ; 0,012 ; 0,017 | +0,087 ; +2,352 ; 0,011 ; 0,015 |
| π | N ≠ 3m (à part) | 5, 7, 8, 10, 11 | −0,328, −0,495, −0,370, −0,422, −0,332 | −0,417 ; +0,208 ; 0,062 ; 0,108 | −0,411 ; +1,147 ; 0,061 ; 0,107 |
| σ | N ≠ 3m | 5, 7, 8, 10, 11 | −0,131, −0,282, −0,195, −0,238, −0,164 | −0,264 ; +0,471 ; 0,050 ; 0,085 | −0,239 ; +1,969 ; 0,049 ; 0,083 |

Famille 3m : trois points, deux paramètres ; le rms distingue 1/N (0,007) de 1/N² (0,024) pour π, pas pour σ (0,012 / 0,011). Famille
N ≠ 3m : rms 0,05–0,06 pour les deux lois, sans tendance monotone en N (5 → 11 : −0,328, −0,495, −0,370, −0,422, −0,332).
Figure `fig/size` : ε − E_D des états localisés contre 1/N (ronds N = 3m, carrés N ≠ 3m ; bleu σ, orange π ; droites d'ajustement π).

Coût de 15×15 (449 atomes) et 18×18 (647 atomes), rien lancé, à partir des murs scf 5…12 (64 rangs MPI ; défaut 13, 29, 50, 118, 200,
301, 721, 907 s ; parfaite 13, 33, 41, 99, 161, 332, 473, 673 s) : t ∝ N^4,96 (défaut) / N^4,59 (parfaite) → 15×15 : 44 / 31 min, 18×18 :
110 / 71 min ; depuis la 9×9 seule avec ces exposants : 42 / 104 (défaut), 28 / 65 min (parfaite) ; avec l'exposant 3 (N_at^1,5) : 15 / 27 et
12 / 21 min. `wfc1.hdf5` ∝ N⁴ : 16,2 Go (15×15) et 33,6 Go (18×18) par super-cellule (12×12 : 6,64 Go) ; pp.x et lecture des wfc en
proportion. Les deux tailles sont 3m.

## B — Convergence en bandes (J3 nscf 21808943 ; J4 21810702, 50 min ; J5 21810855, 21 min)

Tables complètes : `b/B_tables.md` ; résultats `b/b1_results.json` (+ `b1_parity.npz`), `b/b2_results.json`, `b/b3_results.json`
(+ `b3_overlaps.npz`) ; figures `fig/delta_vs_n`, `fig/ladder_bands` (colonnes « n× » = M^L × 81).

### B.1 — le `.save` à 128 bandes

- Liste des 81 k identique à la maille grossière, même ordre (max|Δk| 8,8e-15) ; `check_k_reference` accepte.
- 16 premières bandes contre `defect_unit_cell_9x9.save` : bandes 1–4 (occupées) ≤ 5e-7 eV ; bandes 5–14 ≤ 7 µeV ; bande 15 : 0,17 meV ;
  bande 16 : **0,81 meV** (attendu ≤ 1e-6 eV : non retrouvé pour les deux bandes du haut du scf grossier ; la bande 16 y est dégénérée
  avec la 17ᵉ à 4e-9 eV dans le nscf ; même ordre de grandeur que l'écart scf grossier / nscf dense de R4, 0,92 meV). E_D(K) = −4,238471 eV.
- Bande 128 (ε − E_D) : Γ +74,16, K +80,48, M +75,63 eV ; la base à 128 bandes couvre [−24,0, +80,5] eV autour de E_D.
- Parité des 128 × 81 états : max|1 − |⟨σ_h⟩|| = 6,2e-12 ; 5 547 pairs, 4 821 impairs, 0 indéterminé (comptes par bande dans `B_tables.md`).
- nscf : 117 s de mur, ethr 1e-12, 0 « not converged » ; `.save` 1,59 Go (81 wfc) + 64 fichiers `.wfc*` temporaires (1,4 Go) dans l'outdir.

### B.2 — M à 128 bandes (chaîne grossière de production sur le `.save` neuf)

`compute_ML_R_mpi_shared` (32 rangs × 6 fils, u_nk partagé 28,7 G, Bk = 10 368) : 1 646 s, MaxRSS 34,4 G ; `compute_M_NL` (processus frais) :
1 324 s ; fichiers `b/M_L_9x9_nb128.npy`, `b/M_NL_9x9_nb128.npy`, `b/M_ed_9x9_nb128.npy` (1,72 Go chacun, Ha, sidecars `unit_cell`, nbnd 128,
`.save` source ; hors dépôt). max|M^L| 1,03e-2 Ha, max|M^NL| 0,378 Ha (les deux parties portent les conventions de A.2 : M^L norme
super-cellule, M^NL norme maille) ; hermiticité 2,6e-18 (L), 9,6e-15 (NL). Bloc 16 × 16 contre le M grossier de juin :

| partie | valeurs singulières, écart relatif max (bloc 16 bandes) | idem, bloc bandes 1–8 (non dégénérées) | valeurs propres de H_(a1), écart max (eV) | diagonale, écart max (Ha) |
|---|---|---|---|---|
| L | 8,5e-4 | 7,2e-7 | 4,4e-4 | 2,2e-3 |
| NL | 1,29e-3 | 1,19e-6 | 1,06e-2 | 5,2e-2 |
| tot | 1,31e-3 | 1,27e-6 | 1,15e-2 | 5,3e-2 |

L'attendu (≤ 1e-6) est retrouvé sur le bloc des bandes 1–8 ; le bloc 16 bandes s'en écarte comme en R4 (bandes 16–17 dégénérées : le
sous-espace tronqué n'est pas le même) et par les 0,8 meV de la bande 16 ; les diagonales diffèrent dans les paires dégénérées (jauge).

### B.3 — (a1) à n ∈ {16, 24, 32, 48, 64, 96, 128} bandes (ε − E_D ; w₂ ; seuil 0,0812 ; ψ_nk et parité du `.save` neuf)

Convention de production (M = M^L + M^NL tel que produit par la chaîne) :

| n | dim | état π localisé (ε − E_D ; w₂) | Δε_π vs n précédent (meV) | doublet σ (ε − E_D ; w₂) | σ à moins de 0,3 eV de +0,101 | médiane \|Δε\| vs QE σ / π (eV) | médiane \|Δε\| vs n précédent σ / π |
|---|---|---|---|---|---|---|---|
| 16 | 1 296 | −1,350 ; 0,226 | — | −2,812 ×2 ; 0,287 | aucun | 2,913 / 0,007 | — |
| 24 | 1 944 | −1,368 ; 0,225 | −18,0 | −2,867 ×2 ; 0,236 | aucun | 2,968 / 0,007 | 0,055 / 0,000 |
| 32 | 2 592 | −1,379 ; 0,224 | −11,0 | −2,901 ×2 ; 0,201 | aucun | 3,003 / 0,007 | 0,034 / 0,000 |
| 48 | 3 888 | −1,397 ; 0,222 | −18,4 | −2,909 ×2 ; 0,193 | aucun | 3,010 / 0,007 | 0,007 / 0,000 |
| 64 | 5 184 | −1,408 ; 0,220 | −11,1 | −2,922 ×2 ; 0,180 | aucun | 3,023 / 0,007 | 0,013 / 0,000 |
| 96 | 7 776 | −1,423 ; 0,218 | −14,5 | −2,933 ×2 ; 0,168 | aucun | 3,035 / 0,007 | 0,011 / 0,000 |
| 128 | 10 368 | −1,431 ; 0,216 | −8,8 | −2,940 ×2 ; 0,160 | aucun | 3,041 / 0,007 | 0,007 / 0,000 |

Premier n avec un état pair localisé à moins de 0,3 eV de +0,101 : **aucun** jusqu'à 128. QE : π −0,737 (0,246), σ +0,101 ×2 (0,713).
Autres états π localisés : −1,803 ×2 (0,101) à tous les n (QE −1,759 ×2). La fenêtre compte 2 σ + 28 π à tous les n.

Variante diagnostique M^L × 81 + M^NL (hors prompt), mêmes sous-blocs :

| n | état π localisé (ε − E_D ; w₂) | doublet σ (ε − E_D ; w₂) | autres π localisés | σ à moins de 0,3 eV de +0,101 |
|---|---|---|---|---|
| 16 | −0,727 ; 0,236 | aucun dans [−3, +1] (paire à +1,48) | −1,775 ×2 (0,116) ; +0,263 (0,107) | non |
| 24 | −0,741 ; 0,240 | +0,843 ×2 ; 0,674 | −1,777 ×2 (0,115) ; +0,256 (0,105) | non |
| 32 | −0,747 ; 0,242 | +0,428 ×2 ; 0,697 | −1,779 ×2 (0,115) ; +0,252 (0,105) | non |
| 48 | −0,754 ; 0,244 | +0,345 ×2 ; 0,700 | −1,780 ×2 (0,114) ; +0,248 (0,104) | **oui (premier n)** |
| 64 | −0,757 ; 0,245 | +0,211 ×2 ; 0,708 | −1,782 ×2 (0,114) ; +0,247 (0,103) | oui |
| 96 | −0,759 ; 0,245 | +0,137 ×2 ; 0,713 | −1,783 ×2 (0,114) ; +0,246 (0,103) | oui |
| 128 | −0,760 ; 0,245 | +0,116 ×2 ; 0,714 | −1,783 ×2 (0,114) ; +0,245 (0,103) | oui |

Avec M^L × 81 : Δε_π entre n successifs −14, −6, −7, −3, −2, −1 meV (−0,727 → −0,760 ; QE −0,737) ; doublet σ +0,843 → +0,116 (QE +0,101),
w₂ 0,674 → 0,714 (QE 0,713) ; +0,263 → +0,245 (QE +0,269) ; −1,775 → −1,783 (QE −1,759).

A.1 à chaque n (états en évidence ; convention de production ; δ indépendant de M) :

| état | n | δ | R − ε_QE (meV) | résidu (eV) | 3 plus grands poids (ε_j − E_D ; poids) |
|---|---|---|---|---|---|
| σ 324 (+0,101) | 16 / 32 / 64 / 128 | 0,0363 / 0,0057 / 0,0017 / 0,0003 | −3 042 / −2 820 / −2 748 / −2 717 | 5,04 / 6,53 / 7,10 / 7,49 | 128 : (−3,390 ; 0,19), (−2,940 ; 0,15), (−2,940 ; 0,14) |
| π 320 (−0,737) | 16 / 32 / 64 / 128 | 0,0009 / 0,0003 / 0,0001 / 0,0000 | −306 / −302 / −301 / −301 | 1,82 / 1,96 / 2,08 / 2,14 | 128 : (−1,431 ; 0,76), (+0,090 ; 0,11), (−2,672 ; 0,03) |
| π 326 (+0,269) | 16 / 32 / 64 / 128 | 0,0005 / 0,0001 / 0,0000 / 0,0000 | −130 / −128 / −128 / −128 | 1,25 / 1,35 / 1,43 / 1,47 | 128 : (+0,090 ; 0,87), (−1,431 ; 0,06) |

(δ(128) : σ 0,0003, π ≤ 5e-5 ; `fig/delta_vs_n` ; valeurs à tous les n dans `B_tables.md`.) Le quotient de Rayleigh R − ε_QE et le résidu ne
diminuent pas avec n (convention de production) : la partie hors base ne pèse plus rien (δ → 0) et l'écart reste dans H_in.

A.3 à chaque n (convention de production) : comptages sous E_D − 3 identiques à QE (σ 241, π 55) pour tous les n ; fenêtre 2 σ + 28 π ;
min_k[λ_k(n) − λ_{k−s}(QE)] pour π : s = 0 : −0,612 (16), −0,630, −0,641, −0,660, −0,671, −0,685, −0,694 (128) ; s = 1 : −0,044 ; s = 2 : −0,005
(tous n) ; pour σ : s = 0 et 1 : −2,913 → −3,041 ; s = 2 : +0,322 → +0,193 (fenêtre).

### B.4 — n = 128 ne reproduit pas QE à 0,1 eV pour π (convention de production : −1,431 contre −0,737)

Rapporté sans conclure : δ(128) de l'état π = 3e-5 (Σ|c|² = 0,99997) ; ⟨Ψ_in|ΔV|Ψ_out⟩ restant = −0,0008 (L) − 0,0008 (NL) eV
(n = 16 : −0,0382 − 0,0151) ; ⟨Ψ_in|ΔV^L|Ψ_in⟩ direct = 0,2816 eV contre c†M^Lc/81 = 0,0035 (rapport 81,0 ; NL : 0,0083 = 0,0083 à 1e-16).
Pour la paire σ : δ(128) = 3e-4, ⟨Ψ_in|ΔV|Ψ_out⟩ = −0,026 (L) − 0,016 (NL) eV (n = 16 : −1,56 − 0,53), ⟨Ψ_in|ΔV^L|Ψ_in⟩ = 2,769 contre 0,034
(c†M^Lc/81). Avec M^L × 81 (diagnostique) : π à −0,760 (écart QE −23 meV), σ à +0,116 (+15 meV) à n = 128.

## Fichiers de la campagne et manifestes

Répertoire de travail `graphene/qe/defects/R5_base_vs_M/` : `README.md`, `R5_rapport.md`, modules `r5_*.py` (4), pilote `r5_driver.py`,
`submit_r5.sh`, `submit_r5_m128.sh`, `b/nscf_nb128.in` + `.diff` + `submit_nscf_nb128.sh` + `nscf_nb128.out`, `JOBID`, `r5_log.txt`,
`a/` (json, npz, `A_tables.md`), `c/` (json, `C_tables.md`), `b/` (json, npz, `B_tables.md`, **`M_{L,NL,ed}_9x9_nb128.npy` 3 × 1,72 Go**),
`fig/` (pdf + png : `a_delta`, `a_ladder_corrected`, `size`, `delta_vs_n`, `ladder_bands`), `slurm-r5-*`. Copie versionnée
`graphene-raman/article/R5_base_vs_M/` : tout sauf les `.npy`, `slurm-*`, `JOBID`, `.seen*`, `__pycache__` et les fichiers > 5 Mo.
Fichiers de production modifiés : aucun ; `article/R4_quasi_lie/R4_rapport.md` et sa copie de travail : section « Erratum (R5, 2026-09-25) »
ajoutée en fin de fichier. Git : rien ajouté, rien commité (Greg fait add, commit, push).

Données nouvelles volumineuses (aucune suppression ; GO séparé pour l'un ou l'autre manifeste) :

| fichier | taille | manifeste proposé |
|---|---|---|
| `qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.save/` (81 wfc + XML + charge-density copiée + C.upf) | 1,59 Go | **à miroiter avec md5** vers `qe_tmp_backup/R5_uc9x9_nb128/` si le `.save` à 128 bandes doit servir (B.3 : avec la convention de production, n = 128 ne rapproche pas le π de QE ; avec M^L × 81, la paire σ n'entre à 0,3 eV de +0,101 qu'à partir de n = 48) ; sinon manifeste de nettoyage |
| `qe_tmp/R5_uc9x9_nb128/defect_unit_cell_9x9.wfc1…64` (temporaires pw.x) | 1,4 Go | nettoyage (redondants avec le `.save`) |
| `R5_base_vs_M/b/M_{L,NL,ed}_9x9_nb128.npy` | 3 × 1,72 Go | même sort que le `.save` (recalcul 50 min) |
| `R5_base_vs_M/cache` | aucun | — |

**STOP — rapport R5 terminé le 2026-09-25.** Décisions attendues de Greg : (1) la convention relative de M^L et M^NL dans les fichiers M de
production (A.2 : facteur N_cells = 81 entre les deux parties ; touche `M_ed_*`, `M_dense_*`, M_W, V_loc et tout ce qui en découle) ;
(2) le sort du `.save` à 128 bandes et des M à 128 bandes ; (3) commit de l'erratum R4 et de la copie `article/R5_base_vs_M/`.
