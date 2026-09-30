# EM2 — Couplage électron-photon : références indépendantes (DFT directe, postw90 `kubo`)

Date : 2026-09-29. Série EM (§2.5 du mémoire), après EM1 et M4. Prompt : `memoire/EM/EM2_prompt.md`.
Statut : **PRODUCTION**. Répertoire de travail : `$P/graphene/qe/electron_photon/EM2/` (hors dépôt) ; copie versionnée `memoire/EM/EM2/`.
Référence en LECTURE SEULE (rien écrit) : `$P/graphene/qe/defects/unit_cell/27x27/` et `$S/qe_tmp/defect_uc_dense_27/`.
Ni commit ni push ; `src/`, `tests/` et `memoire/EM/EM.md` non modifiés. GO par étape levé par Greg en début de séance.

Chemins abrégés : `$P` = `/home/gregb26/links/projects/rrg-cotemich-ac/gregb26`, `$S` = `/home/gregb26/links/scratch`.
Dépôt à `7a64154` (après le pull de Greg ; `memoire/EM/M4_sigma/` présent). Binaires : module `quantumespresso/7.5`
(`module restore qe`) : `pw.x`, `bands.x` (QE 7.5), `postw90.x` (Wannier90 3.1.0, `postw90.x -v`).

## Résumé

- **A (DFT directe)** : 84 k (72 sur les anneaux, 12 sur un cercle de Dirac), pw.x `bands` + `bands.x lp` (job 22034896, 1 min 33).
  ε_π, ε_π* QE − Wannier ≤ **0,09 meV**. |ħv_cv| dans le plan, k par k : **complet −0,042 % à +0,002 %** ; centres seuls −1,8 à −3,3 % ;
  sans Berry −47 à +29 %. Cercle de Dirac : ħv_F des pentes QE 5,46977 = ⟨|ħv_cv|⟩ de p 5,46978 eV·Å (**10⁻⁶**) → p de `bands.x` = m_e v,
  commutateur [V_NL, r] **inclus** (confirmé dans le source qe-7.5) ; Wannier 5,4692 (+1×10⁻⁴).
- **B (postw90 3.1.0 `kubo`, 1201², réglages du prompt)** : σ_xx/σ₀ = 1,00235 à 0,2 eV (après **×2 de spin, absent de postw90**) ;
  écart à M4 −1,4 / −2,2 / −2,7×10⁻³ aux lasers (yy : −0,4 / −0,7 / −0,9×10⁻³), 2,9×10⁻² au max près de van Hove.
  Décomposé exactement : **diagonale de r** (postw90 par défaut : éq. 44 sans log ; `_tb.dat` : −Im ln M) −1,2 à −2,2×10⁻³ ;
  **préfacteur** 1/(ε_c − ε_v) contre 1/ħω de `kubo.py` −2,6 à −3,5×10⁻⁴ ; `use_ws_distance` −0,2 à −1,3×10⁻⁴ ; grille 10⁻¹⁵ ;
  **résidu à conventions égales ≤ 5×10⁻⁵ aux lasers** (2,1×10⁻⁴ au max).
- Largeur postw90 = √2 η = 0,0565685425 eV (vérifié) ; `kubo_eigval_max` par défaut = E_D + 3,17 eV (fixé à 1000) ; `postw90.x` du module = série.

## Étape 0 — Inventaire (lecture seule)

- F17 `ring_kpoints_crystal(tb, K, mu, npoints)` présent (`ring.py:233`), signature de l'appel du prompt ; `memoire/EM/M4_sigma/` présent après le pull.
- Référence : `nscf.in` (calculation `nscf`, prefix `defect_uc_dense_27`, outdir `$S/qe_tmp/defect_uc_dense_27`, ecutwfc 100 Ry,
  smearing mv 0,01 Ry, nbnd 20, `assume_isolated = '2D'`, `nosym`, `noinv`, 8 électrons, pseudo `$P/abinit_processing/pseudo/C.upf`,
  NC, MD5 `34a24e64…`) ; `wannier.win` (num_wann 5, num_bands 20, dis_froz_max −1,74 eV, mp_grid 27 27 1) ; `wannier.mmn` (87 Mo),
  `.chk`, `.eig` présents (postw90 relit le `.mmn` pour r(R), `get_AA_R`).
- **Cellule** : `CELL_PARAMETERS bohr` du nscf = `unit_cell_cart bohr` du `.win`, mêmes trois vecteurs, même ordre
  (4.0354919061, ∓2.3298923383, 0 ; 0, 0, 30). Revérifié par script (`em2_kpoints.py`, égalité exacte des flottants).
- `.save` de référence : 2,1 Go, 729 k ; `data-file-schema.xml` = celui du nscf, `fermi_energy` = −0,15508 Ha = **−4,21994 eV** (E_D + 19,0 meV).

## A. DFT directe aux anneaux

### A.1 Liste de k (`em2_kpoints.py` → `em2_kpoints_crystal.txt`, `em2_kpoints_table.txt`, `em2_kpoints.log`)

- `ring_kpoints_crystal(tb, K, E_D, {2.33: 48, 1.96: 12, 2.54: 12})` sur `wannier/27x27/wannier_tb.dat` (sha256 `baa17b88…`),
  `K = make_grid_tb(tb, 3).K` = (2/3, 1/3, 0) à 10⁻¹² ; premier point à 2,33 eV : (0.73959906, 0.40626573, 0) ✔.
- + 12 k sur un cercle de rayon q = 0,005 Å⁻¹ autour de K (θ = 0, 30, …, 330°), k_red = k·a_i/2π.
- 84 k au total, ordre : 2,33 (48), 1,96 (12), 2,54 (12), Dirac (12) ; 14 décimales, poids 1. La table donne ħω, θ, k cristallin et cartésien, |k − K|.
- Anneaux centrés sur K à 3×10⁻⁷ – 10⁻⁶ Å⁻¹ ; |k − K| : 0,1638–0,2102 (1,96), 0,1922–0,2617 (2,33), 0,2081–0,2940 Å⁻¹ (2,54).
- `tb.lattice` = cellule QE × 0,529177210903 à 6,9×10⁻⁸ Å : Wannier90 3.1.0 est compilé avec les constantes CODATA 2006
  (bohr = 0,52917720859 Å) ; sans effet sur les coordonnées cristallines (sans dimension, calculées dans le même repère).

### A.2 Run QE `bands` (job **22034896**, 32 tâches, `-nk 4`, 1 min 33 s dont pw.x 15,6 s et bands.x 43 s)

- `bands.in` = `nscf.in` de référence avec trois changements seulement (`bands.in.diff`) : `calculation = 'bands'`,
  `outdir = '$S/qe_tmp/em2_bands_27'`, liste de k. Mêmes ecutwfc, smearing, nbnd 20, `nosym`, `noinv`, conv_thr, `diago_full_acc`.
- outdir NEUF : `charge-density.hdf5` et `data-file-schema.xml` copiés (`cp -p`) du `.save` de référence, md5 identiques
  (`e5704084…`, `b0dfc11c…`) ; rien écrit dans la référence.
- Le potentiel est celui de la référence : même référence d'énergie que le `.eig` de Wannier90 (vérifié : ε QE − ε Wannier ≤ 0,09 meV, A.4).
- `.wfcN` en vrac (24 fichiers) supprimés après le run ; `.save` (84 `wfc*.hdf5`, 255 Mo) miroité :
  `$P/graphene/qe/qe_tmp_backup/em2_bands_27/` (`rsync -a -r --no-o --no-g --open-noatime --chmod=Dg+s`), **md5 OK, 88 fichiers**,
  `qe_tmp_backup/MD5SUMS_EM2_2026-09-29.txt` ; groupe `rrg-cotemich-ac` partout.

### A.3 `bands.x` avec `lp = .true.` (`bands_pp.in`, `filp = 'em2_p_avg.dat'`) — ce que calcule QE 7.5

Lu dans les sources au tag `qe-7.5` (`$P/codes/qe`, `git show qe-7.5:…` ; `write_p_avg.f90` et `compute_ppsi.f90` identiques au checkout) :

- **Commutateur non local : INCLUS.** `compute_ppsi` appelle `commutator_Hx_psi` (PW/src) puis multiplie par i/2 :
  `ppsi = (i/2)[H, x_α]|ψ_v⟩`. `commutator_Hx_psi` = terme cinétique −2i(k+G)_α·tpiba + terme non local (dérivées des
  projecteurs, `gen_us_dj`, `gen_us_dy`, coefficients D de `compute_deff`) ; le `goto 111` qui l'enlèverait est commenté.
  En unités Rydberg (ħ = 1, m = 1/2) : ppsi = m v|ψ⟩, soit « p » = m_e v̂ avec v̂ = (i/ħ)[Ĥ, r̂] complet (Ismail-Beigi 2001).
  Pseudo NC (`pseudo_type="NC"`) : la branche ultradouce (`okvan`) ne sert pas.
- **Unités** : p en ħ/bohr ; ħv = (ħ²/m_e)·p = 2 Ry·bohr × p = **14,39964548 eV·Å × p** (constantes CODATA 2018 de QE 7.5).
- **Format** : bloc `&p_mat nbnd=20, nks=84 /`, puis pour chaque k : `xk` (cartésien, 2π/alat, `3f10.6`) et `nbnd_occ` ;
  pour α = 1, 2, 3 (x, y, z cartésiens de la cellule QE = ceux de Wannier90) : une ligne par bande vide c = nbnd_occ+1…nbnd,
  contenant |⟨ψ_c|p_α|ψ_v⟩|² pour v = 1…nbnd_occ (`5f15.8`). Seuls les **modules au carré** sont écrits, pas les phases.
- **Pas de moyenne sur les états dégénérés** (malgré le nom `p_avg`) : chaque paire (c, v) est écrite séparément.
  Sans conséquence ici : π et π* ne sont dégénérés avec aucune autre bande sur les 84 k.
- **Occupation** : `nbnd_occ` = dernière bande avec `et ≤ ef`, ef lu dans le XML (celui du nscf, E_D + 19 meV ; un run `bands` ne le recalcule pas).
  `lp` n'a **pas** refusé le smearing (aucun test d'occupation dans `bands.f90` ni `write_p_avg`) : pas de relance en `occupations = 'fixed'`.
  **nbnd_occ = 4 aux 84 k**, y compris sur le cercle de Dirac (π* à E_D + 27 meV > ef = E_D + 19 meV).
- `bands_pp.out` (138 Ko) : l'essentiel est l'analyse de symétrie de `lsym` (défaut `.true.`), sans objet ici.

### A.4 Comparaison k par k (`em2_A_compare.py` → `em2_A_summary.txt`, `em2_A_table.txt`, `em2_A.npz`, `em2_A_compare.log`)

- ε QE : XML du run `bands` (Hartree, pleine précision, × 27,211386245988). Ordre des k vérifié contre la liste (XML à 10⁻⁹, fichier p à 10⁻⁶).
- π, π* choisis **par l'énergie** autour de E_D, dans QE et dans Wannier : ce sont les bandes **4 et 5 de QE aux 84 k** (assert), la dernière
  colonne de valence et la première ligne vide du fichier p (nbnd_occ = 4).
- |ħv_cv| dans le plan = 14,39965 × √(|p^x_cv|² + |p^y_cv|²) pour QE ; √(|ħv^x_cv|² + |ħv^y_cv|²) de `compute_velocity` pour les trois
  variantes (complet `tb, 'berry'` ; `centres_only(tb), 'berry'` ; `tb, 'no_berry'`). |p^z_cv| = 0 exactement (miroir σ_h).

**Énergies (QE − Wannier, meV)** :

| Ensemble | n | ε_π min / moy / max | ε_π* min / moy / max | \|ħv_cv\|_QE (eV·Å) |
|---|---|---|---|---|
| anneau 1,96 eV | 12 | −0,02 / +0,01 / +0,03 | −0,01 / +0,02 / +0,05 | 4,3772–6,9223 |
| anneau 2,33 eV | 48 | −0,01 / +0,00 / +0,03 | −0,03 / +0,03 / +0,09 | 4,2030–7,2421 |
| anneau 2,54 eV | 12 | −0,03 / −0,01 / +0,02 | −0,03 / −0,01 / +0,01 | 4,1080–7,4297 |
| cercle de Dirac | 12 | −0,00 / −0,00 / −0,00 | +0,00 / +0,00 / +0,01 | 5,4348–5,5049 |

**|ħv_cv|_Wannier / |ħv_cv|_QE − 1 (%), min / moy / max, k par k** :

| Ensemble | n | complet | centres seuls | sans Berry |
|---|---|---|---|---|
| anneau 1,96 eV | 12 | −0,042 / −0,013 / +0,002 | −1,956 / −1,830 / −1,752 | −33,8 / −0,09 / +22,5 |
| anneau 2,33 eV | 48 | −0,040 / −0,017 / −0,000 | −2,747 / −2,551 / −2,394 | −42,2 / −0,02 / +28,1 |
| anneau 2,54 eV | 12 | −0,035 / −0,017 / −0,005 | −3,252 / −3,002 / −2,782 | −47,3 / +0,05 / +29,4 |
| cercle de Dirac | 12 | −0,016 / −0,011 / −0,005 | −0,018 / −0,013 / −0,006 | −0,73 / −0,01 / +0,70 |

**Composantes** (|ħv^x_cv|² et |ħv^y_cv|² sont chacune invariantes de jauge ; écart max / (ħv_F)², x ; y) : complet ≤ 1,1×10⁻³ ; 7,4×10⁻⁴ ;
centres seuls ≤ 0,10 ; 0,076 ; sans Berry ≤ **0,88** ; 0,076 (sans Berry, la composante y est celle des centres seuls, comme en F15).

- Contrôle d'unités près de K : |ħv_cv|_QE sur le cercle de Dirac = 5,43–5,50 eV·Å (gauchissement trigonal ±0,64 % à q = 0,005),
  moyenne 5,46978 → ħv_F = 5,469 ✔.
- **C₃** : |ħv_cv|_QE est invariant par θ → θ + 120° à ≤ 4×10⁻⁶ près ; Wannier complet à 0,6–1,5×10⁻⁴ près (même ordre que
  |σ_yy − σ_xx| = 1,9×10⁻⁴ de M4). Le motif en θ de l'écart complet − QE (figure, panneau b) vient donc de Wannier.

### A.5 Contrôle d'unités indépendant de Wannier (cercle de Dirac, q = 0,005 Å⁻¹, 12 k)

| Grandeur | eV·Å |
|---|---|
| ħv_F^DFT = ⟨(ε_π* − ε_π)/(2q)⟩_θ (valeurs propres QE) | **5,46977** (dispersion en θ 3,5×10⁻², gauchissement trigonal) |
| ⟨\|ħv_cv\|⟩_θ de `bands.x lp` | **5,46978** (dispersion 7,0×10⁻²) |
| rapport p / pentes − 1 | **+1,1×10⁻⁶** |
| Wannier, `fermi_velocity(q = 0,005)` inter / intra ; ⟨\|ħv_cv\|⟩ aux 12 k | 5,46916 / 5,46921 ; 5,46916 |
| ħv_F^DFT / 5,469188 (M4, q = 10⁻³) − 1 | +1,06×10⁻⁴ |

- Accord à 10⁻⁶ ≪ 10⁻³ : **p de `bands.x` = m_e v, commutateur [V_NL, r] compris**, conversion 14,39965 eV·Å juste.
  Ce que le commutateur pèse n'est pas chiffré (il faudrait un `bands.x` recompilé avec le `goto 111`, hors cadre).
- La moyenne sur 12 angles annule les termes en cos 3θ ; le biais isotrope en q² est ≈ (q a)² ~ 5×10⁻⁵.

Figure de contrôle : `em2_A_anneaux.{pdf,png}` (`em2_A_figure.py`, style du mémoire) — (a) |ħv_cv|(θ) sur l'anneau de 2,33 eV,
DFT (48 points) et les trois variantes Wannier (720 angles) ; (b) écart relatif complet − DFT sur les trois anneaux laser.

## B. postw90 `kubo`

### B.1 Mise en place (`postw90/`)

- `wannier.chk` et `wannier.eig` copiés (`cp -p`) d'EM1 (md5 = référence : `42311f82…`, `9cb7e96a…`), `wannier.mmn` en lien symbolique
  vers la référence (87 Mo, lecture seule), `wannier.win.ref` = `.win` de référence (md5 `86d73765…`).
- Un sous-répertoire par run (postw90 écrit des noms fixes), `.win` produit par `mkwin.sh DIR N [lignes]` (bloc kubo inséré après
  `write_hr`, liens, `wannier.win.diff`). Lanceur `submit_postw90.sh DIR…` (1 tâche, `srun -n 1 postw90.x wannier` par répertoire).

### B.2 Conventions vérifiées dans le source de Wannier90 **v3.1.0**

Source : archive du tag v3.1.0 (GitHub, sha256 `40651a98…`, dans le scratchpad de la séance). `$P/codes/wannier90` est un checkout
*develop* restructuré (`postw90_readwrite.F90`), différent de 3.1.0 : il n'a pas servi.

| Point | Ce que fait postw90 3.1.0 | Réglage / conséquence |
|---|---|---|
| Mots-clés | `berry`, `berry_task = kubo`, `berry_kmesh`, `fermi_energy`, `kubo_freq_{min,max,step}`, `kubo_adpt_smr`, `kubo_smr_type`, `kubo_smr_fixed_en_width`, `kubo_eigval_max` (tous confirmés dans `parameters.F90`) | 0,2 → 6,0 eV par 0,01 : 581 fréquences, même grille que M4 (à la précision `E16.8` du fichier) |
| **Largeur** | `delta = utility_w0gauss((ε_m−ε_n−ħω)/w, 0)/w` avec `utility_w0gauss(x, 0) = exp(−x²)/√π` (`utility.F90:1049`, `berry.F90`) : gaussienne d'écart-type w/√2 | **w = √2·η = 0,0565685425 eV** (ta supposition était juste) |
| Largeur fixe | `kubo_adpt_smr = false` ; la ligne « adaptive smearing » de `berry_main` n'est écrite que si le lissage est adaptatif : absente du `.wpout` ✔ | Le `.wpout` affiche « Using global smearing parameters » : **bug de précédence** dans `param_write` (`kubo_adpt_smr .eqv. adpt_smr .and. …` se lit `F .eqv. (… .and. w == 0)` = vrai) ; sans effet sur le calcul |
| **`kubo_eigval_max`** | défaut = `dis_froz_max + 0,6667` = −1,073 eV = **E_D + 3,17 eV** s'il y a une fenêtre gelée (`parameters.F90:1991`) : couperait les transitions dont le π* est au-dessus | fixé à **1000 eV**. (Le défaut de `kubo_freq_max` est aussi `dis_froz_max − E_F + 0,6667` = 3,17 eV ; fixé à 6,0) |
| **Spin** | Kubo sur les états de Wannier sans spin, occupations 0/1 (`pw90common_get_occ`), aucun facteur 2 (`berry.tex`, éq. sig-H ; « Number of electrons per state : 2 » est seulement affiché) | **× g_s = 2** appliqué à postw90 ; contrôle : **0,50118 brut → 1,00235** à 0,2 eV |
| Unités | `fac = 1e8 e²/(ħ V_c)` → S/cm pour la cellule 3D ; colonne 2 de `-kubo_S_ab.dat` = Re σ^H (partie absorptive, symétrisée) | σ_2D/σ₀ = σ_3D × 100 × c × 2/σ₀, c = 30 bohr = 15,8753163 Å, σ₀ = e²/4ħ = 6,085337×10⁻⁵ S avec les constantes de la compilation (bloc CODATA 2006 de `w90_constants`, bohr = 0,52917720859 Å) : les constantes s'annulent |
| `use_ws_distance` | défaut `.true.` (`parameters.F90:1509`) | gardé pour le run du prompt ; mesuré séparément |
| Hermitisation de A | partie hermitienne de A(q) sur la grille grossière, avant la TF (`get_oper.F90:468`) | équivalente à la nôtre en k (décision 3 d'EM.md) |
| **Diagonale de r** | `get_AA_R` : A_nm(q) = i Σ_b w_b b S_nm pour **tous** n, m, diagonale comprise (éq. 44 WYSV06), sauf si `transl_inv = true` : diagonale = −Σ_b w_b b Im ln S_nn (MV97 éq. 31). Le `_tb.dat` de wannier90.x (`hamiltonian_write_tb`) prend **toujours** −Σ w b Im ln M_nn sur la diagonale | **postw90 par défaut n'utilise pas le même r(R) que notre `_tb.dat`** ; `transl_inv = true` les aligne (testé, ci-dessous) |
| **Préfacteur** | chaque transition pondérée par (ε_m − ε_n)\|A_nm\|² = \|ħv_nm\|²/(ε_c − ε_v) | `kubo.py` (M4) pondère par \|ħv_cv\|²/ħω : même limite η → 0, écart O(η²) ≈ (η/ħω)²(σ − ω dσ/dω) |
| Parallélisme | **`postw90.x` du module = exécutable SÉRIE** (« Running in serial (with serial executable) », compilé sans `-DMPI` bien que lié à OpenMPI) | 1 tâche ; avec 8 tâches, 8 copies identiques écrivent les mêmes fichiers |

### B.3 Runs (tous : binaire du module, 1 tâche, ~105 Mo de mémoire)

| Répertoire | berry_kmesh | Options en plus | Job | Durée postw90 |
|---|---|---|---|---|
| `k301` | 301 301 1 | — (réglages du prompt) | 22035645 | 158 s |
| `k301_nows` | 301 301 1 | `use_ws_distance = false` | 22035645 | 61 s |
| `k301_ti` | 301 301 1 | `transl_inv = true` | 22036408 | 148 s |
| `k301_ti_nows` | 301 301 1 | `transl_inv = true`, `use_ws_distance = false` | 22036408 | 61 s |
| **`k1201`** | **1201 1201 1** | — (**réglages du prompt, production**) | **22035981** | **2318 s** |
| `k1201_ti` | 1201 1201 1 | `transl_inv = true` | 22036409 | 2339 s |
| `k1201_ti_nows` | 1201 1201 1 | `transl_inv = true`, `use_ws_distance = false` | 22036410 | 1001 s |

- Ressources du vrai calcul, déduites de l'essai à 301² (156 s) : ×16 en k → ~40 min en série ; demandé 1 tâche, 4 Go, 1 h 30 ; mesuré 38 min 40, 105 Mo.
- Jobs ratés, tous de cette campagne : **22035375** (8 tâches = 8 copies série ; annulé après 1 min 38, son `.wpout` gardé en
  `k301/wannier.wpout.22035375_annule`) ; **22035532** (`/usr/bin/time` absent des nœuds, échec en 2 s) ; **22036356–58**
  (`berry_kmesh` corrompu par mon premier générateur de `.win`, annulés avant de démarrer ; générateur réécrit en `mkwin.sh`, qui contrôle N).
- Décomposition côté dépôt : `em2_B_ours.py` (job **22035650**, 8 cœurs, 9 min 19, 3,9 Go ; `submit_B_ours.sh`) recalcule notre σ
  (mode complet, `compute_velocity`, `gaussian_eta`, `kubo_normalize` du paquet ; `src/` non modifié) : N = 1200 décalé en 1/ħω
  (**redonne le npz de M4 à 2,8×10⁻¹⁴**), en 1/Δ, puis grille centrée en Γ N = 1201 et 301 en 1/Δ → `em2_B_ours.npz`.

### B.4 postw90 (réglages du prompt, 1201²) contre M4 (`em2_B_compare.py` → `em2_B_compare_k1201_k1201_ti_k1201_ti_nows.txt`, `em2_B_compare.log`)

| ħω (eV) | σ_xx postw90 | σ_xx M4 | écart | σ_yy postw90 | σ_yy M4 | écart |
|---|---|---|---|---|---|---|
| 0,20 | 1,00235 | 1,00265 | −2,9×10⁻⁴ | 1,00236 | 1,00246 | −1,0×10⁻⁴ |
| 1,96 | 1,26227 | 1,26367 | −1,40×10⁻³ | 1,26317 | 1,26360 | −4,3×10⁻⁴ |
| 2,33 | 1,40336 | 1,40555 | −2,19×10⁻³ | 1,40469 | 1,40536 | −6,7×10⁻⁴ |
| 2,54 | 1,50880 | 1,51148 | −2,68×10⁻³ | 1,51044 | 1,51135 | −9,1×10⁻⁴ |
| pic de van Hove | 4,06 eV, 6,8778 | 4,05 eV, 6,8992 | | | | |

- Toute la courbe : max \|écart\| xx **2,9×10⁻²** à 4,01 eV (moyenne 4,2×10⁻³) ; yy 1,5×10⁻² à 4,00 eV (moyenne 1,5×10⁻³).
- σ_xy : \|σ_xy\| ≤ 2,8×10⁻⁸ (M4 : 9×10⁻¹²) ; σ_zz ~ 6×10⁻²² (miroir σ_h).
- Contrôle σ/σ₀ → 1 à 0,2 eV : **1,00235** (0,50118 avant le facteur de spin).

### B.5 Décomposition de l'écart (`em2_B_decompose.py` → `em2_B_decomposition.txt` ; la somme des termes = l'écart total à 10⁻¹²)

σ_xx/σ₀ :

| Terme | 0,20 eV | 1,96 | 2,33 | 2,54 | 4,05 | max \|terme\| (eV) |
|---|---|---|---|---|---|---|
| préfacteur 1/Δ au lieu de 1/ħω (notre code) | −1,8×10⁻⁴ | −2,6×10⁻⁴ | −3,1×10⁻⁴ | −3,5×10⁻⁴ | −1,7×10⁻³ | 1,1×10⁻² (4,11) |
| grille centrée en Γ 1201 au lieu de décalée 1200 | +8×10⁻⁷ | 10⁻¹⁵ | 10⁻¹⁵ | 10⁻¹⁵ | 10⁻¹⁴ | 8×10⁻⁷ (0,20) |
| résidu : postw90 (`transl_inv`, sans ws) − nous (même grille, même préfacteur) | −2,5×10⁻⁵ | +3,5×10⁻⁵ | +5×10⁻⁶ | −5,3×10⁻⁵ | −4×10⁻⁵ | 2,1×10⁻⁴ (3,89) |
| `use_ws_distance` | −7×10⁻⁵ | −1,5×10⁻⁵ | −1,3×10⁻⁴ | −9,4×10⁻⁵ | −1,2×10⁻³ | 1,7×10⁻³ (4,02) |
| diagonale de r (éq. 44 sans log, défaut de postw90) | −1×10⁻⁵ | **−1,2×10⁻³** | **−1,7×10⁻³** | **−2,2×10⁻³** | −1,9×10⁻² | 1,9×10⁻² (4,06) |
| **total postw90 − M4** | −2,9×10⁻⁴ | −1,4×10⁻³ | −2,2×10⁻³ | −2,7×10⁻³ | −2,1×10⁻² | 2,9×10⁻² (4,01) |

σ_yy : même décomposition dans `em2_B_decomposition.txt` ; diagonale de r −2,6, −4,2, −5,4×10⁻⁴ aux lasers (3 à 4 fois moins qu'en xx),
`use_ws_distance` +5×10⁻⁵ à +1×10⁻⁴ (signe opposé à xx), préfacteur identique à xx.

**Isotropie** \|σ_yy − σ_xx\| (max sur la courbe ; aux trois lasers) :

| Calcul | max | 1,96 / 2,33 / 2,54 eV |
|---|---|---|
| M4 (et notre σ en 1/Δ, grille Γ) | 2,7×10⁻³ | 7×10⁻⁵ / 1,9×10⁻⁴ / 1,3×10⁻⁴ |
| postw90, réglages du prompt | 1,2×10⁻² | 9×10⁻⁴ / 1,3×10⁻³ / 1,6×10⁻³ |
| postw90, `transl_inv` | **2,6×10⁻⁵** | 6×10⁻⁷ / 2×10⁻⁶ / 2×10⁻⁶ |
| postw90, `transl_inv`, sans ws | 2,7×10⁻³ | 6×10⁻⁵ / 2,3×10⁻⁴ / 1,5×10⁻⁴ |

- La petite anisotropie de M4 (déjà notée dans EM.md, 1,9×10⁻⁴ à 2,33 eV) est celle de postw90 **sans** `use_ws_distance` ;
  avec, elle tombe à 2×10⁻⁶ aux lasers. La boîte de R simple (décision 4 d'EM.md) n'est pas symétrique par C₃, les images de
  Wigner-Seitz le sont. Même ordre que la brisure de C₃ de |ħv_cv| Wannier vue en A (0,6–1,5×10⁻⁴) : lien plausible, non testé.
- La diagonale « éq. 44 sans logarithme » de postw90 par défaut brise davantage l'isotropie (1,2×10⁻²) : elle n'est pas invariante
  par translation (c'est le sens de `transl_inv`).
- Le **résidu** (2,1×10⁻⁴ max sur le flanc de van Hove, ≤ 5×10⁻⁵ aux lasers) est le même à 301² et à 1201² : ce n'est pas du bruit
  de grille. Sur le flanc (dσ/dħω ≈ 10 eV⁻¹), il équivaut à ~2×10⁻⁵ eV. Origine non cherchée (r(R) reconstruit par postw90 depuis le
  `.mmn` et le `.chk` au lieu du `m_matrix` de wannier90.x, H(R) de `get_HH_R`…).
- L'essai à 301² n'est pas convergé (notre σ, grille Γ : 301 − 1201 jusqu'à 5,5×10⁻² à 3,34 eV) ; il sert au temps de calcul et au
  test de l'hypothèse `transl_inv` (`em2_B_compare_k301_k301_nows_k301_ti_k301_ti_nows.txt`).

Figure de contrôle : `em2_B_sigma.{pdf,png}` (`em2_B_figure.py`) — (a) σ_xx/σ₀ de M4 et de postw90 (réglages du prompt et `transl_inv`) ;
(b) \|écart\| à M4 en échelle log : préfacteur seul, postw90, postw90 `transl_inv`, résidu.

`em2_postw90_sigma.npz` (hw, sigma (581, 3, 3) en σ/σ₀, g_s = 2 appliqué, format de `scripts/make_figures_em.py`) est écrit depuis
**`k1201` (réglages du prompt)**. `make_figures_em.py` n'a pas été relancé (EM3).

## C. Ce qui reste à décider (Greg)

1. **Quelle courbe postw90 montrer** (EM3) : réglages du prompt (`k1201`, −2,2×10⁻³ à 2,33 eV, surtout la diagonale de r) ou
   `transl_inv = true` (`k1201_ti`, −4,4×10⁻⁴ à 2,33 eV, dont −3,1×10⁻⁴ de préfacteur), qui utilise le même r(R) que notre `_tb.dat`.
   Le npz écrit est celui des réglages du prompt ; celui de `k1201_ti` se produit en passant son répertoire en premier dans le code (une ligne).
2. **Critère « concordant avec postw90 »** de M4 (EM.md §8) : pas de seuil chiffré dans EM.md. Aux lasers, l'écart est dominé par deux
   conventions (diagonale de r, préfacteur), et le résidu à conventions égales est ≤ 5×10⁻⁵.
3. **Préfacteur 1/ħω ou 1/Δ** dans `kubo.py` : les deux sont exacts quand η → 0 ; l'écart (−3×10⁻⁴ aux lasers, 1,1×10⁻² près de van Hove)
   fait partie du biais en η du pilote M4 (σ(η) ≈ σ(0) + 0,42 η²). Le biais de la forme 1/Δ n'a pas été mesuré.
4. **`use_ws_distance`** : EM.md §10 le donne « négligeable (≤ 10⁻⁵) » ; mesuré ici −1,3×10⁻⁴ à 2,33 eV (xx), jusqu'à 1,7×10⁻³ près de
   van Hove, et c'est lui qui rend σ isotrope à 2×10⁻⁶ près. Décision 4 d'EM.md (pas de Wigner-Seitz) à réexaminer ou à documenter.
5. **EM.md à réaligner** (session locale) : §10 EM2 (facteur de spin ×2 absent de postw90, largeur √2 η, conversion σ_3D × 100 × c,
   `kubo_eigval_max` et `kubo_freq_max` par défaut liés à `dis_froz_max`, binaire série, `transl_inv`) ; §2/§5 : p de `bands.x`
   = m_e v avec [V_NL, r] (répond à la question ouverte de §0 « Potentiel non local ») ; statut EM2 dans §11.

## Fichiers, commandes, md5

- **Commandes** (depuis `$P/graphene/qe/electron_photon/EM2/`, environnement `source /etc/profile; module restore qe; module load mpi4py/4.0.3 scipy-stack;`
  `export PYTHONPATH=$GRAPHENE_RAMAN/src:$PYTHONPATH`, Python = `$GRAPHENE_RAMAN/.venv/bin/python`) :
  `python em2_kpoints.py` → `sbatch submit_bands.sh` (22034896) → `python em2_A_compare.py` → `python em2_A_figure.py` ;
  `cd postw90; ./mkwin.sh DIR N [lignes]; sbatch [--time=…] submit_postw90.sh DIR…` (22035645, 22035981, 22036408–10) ;
  `sbatch submit_B_ours.sh` (22035650) ; `python em2_B_compare.py k1201 k1201_ti k1201_ti_nows` (et la série k301) ;
  `python em2_B_decompose.py` ; `python em2_B_figure.py` ; `./em2_sync_copy.sh` (copie versionnée).
- **md5** : `MD5SUMS_EM2_2026-09-29.txt` (liste de k, `bands.in`, `em2_p_avg.dat`, npz, `.chk`/`.eig`/`.win` de postw90, `kubo_S` de
  production, XML du run bands) ; miroir du `.save` : `qe_tmp_backup/MD5SUMS_EM2_2026-09-29.txt` (88 fichiers, `md5sum -c` OK).
- **Scratch** : `$S/qe_tmp/em2_bands_27/` (`.save` 255 Mo + `defect_uc_dense_27.xml`), miroité ; rien d'unique sur le scratch.
- **Copie versionnée** `memoire/EM/EM2/` (102 fichiers, 3,2 Mo, `em2_sync_copy.sh`) : tout sauf les `.chk`/`.eig` (= référence),
  le lien `.mmn`, les sorties SLURM (`em2_*_<job>.out/.err` ; la sortie d'`em2_B_ours.py` est gardée en `em2_B_ours.log`), les `JOBID_*`,
  `bands_pp.out` et `em2_bands.dat*` (sorties annexes de `bands.x`), les `jdos`/`kubo_A` de postw90 et le `.wpout` du run annulé.
  Aucun fichier > 5 Mo. Rien de commité.
- Rien n'a été écrit dans `unit_cell/27x27/`, dans `$S/qe_tmp/defect_uc_dense_27/`, dans `src/`, `tests/` ni `memoire/EM/EM.md` ;
  aucun job antérieur touché (seul `R7c_relax_27` 21852177 était en file, en attente, non touché).

## Résumé en dix lignes (pour la session locale)

1. EM2 fait le 2026-09-29 sur rorqual, PRODUCTION ; travail `graphene/qe/electron_photon/EM2/`, copie `memoire/EM/EM2/` (non commitée).
2. A : 84 k (anneaux 1,96/2,33/2,54 eV : 12/48/12 ; cercle de Dirac q = 0,005 Å⁻¹ : 12), QE 7.5 `bands` + `bands.x lp`, 1 min 33.
3. ε_π, ε_π* QE − Wannier ≤ 0,09 meV ; |ħv_cv| complet : −0,042 à +0,002 % k par k ; centres seuls −1,8 à −3,3 % ; sans Berry −47 à +29 %.
4. `bands.x lp` : p = (i/2)[H, x] en Ry = m_e v, [V_NL, r] inclus (source qe-7.5) ; ħv = 14,39965 eV·Å × p ; |p|² seulement, pas de moyenne dégénérée.
5. Contrôle sans Wannier : ħv_F des pentes QE 5,46977 = ⟨|ħv_cv|⟩ de p 5,46978 eV·Å (10⁻⁶) ; Wannier 5,4692 (+1×10⁻⁴).
6. B : postw90 3.1.0 (série), 1201², w = √2·0,04, `kubo_eigval_max` = 1000 (défaut = E_D + 3,17 eV), ×2 de spin (absent de postw90) : 1,00235 à 0,2 eV.
7. postw90 − M4 (xx) : −1,4 / −2,2 / −2,7×10⁻³ aux lasers, 2,9×10⁻² max près de van Hove ; yy environ trois fois moins.
8. Décomposé exactement : diagonale de r (postw90 éq. 44 sans log vs `_tb.dat` −Im ln M ; `transl_inv = true` aligne) −1,2 à −2,2×10⁻³ ;
   préfacteur 1/(ε_c−ε_v) vs 1/ħω −3×10⁻⁴ ; `use_ws_distance` ≤ 1,3×10⁻⁴ ; grille 10⁻¹⁵ ; résidu à conventions égales ≤ 5×10⁻⁵.
9. `use_ws_distance` rend σ isotrope (|σ_yy − σ_xx| 2×10⁻⁶ aux lasers, contre 1,9×10⁻⁴ pour M4) : l'anisotropie de M4 vient de la boîte de R simple.
10. À décider : courbe postw90 de la figure (réglages du prompt écrits dans `em2_postw90_sigma.npz`, ou `transl_inv`), seuil de concordance,
    préfacteur, `use_ws_distance` ; EM.md §10/§11 à réaligner (spin, largeur, conversion ×100·c, défauts liés à `dis_froz_max`, binaire série).
