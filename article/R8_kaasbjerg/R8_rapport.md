# R8 — Reproduction des Fig. 13 et 14 de Kaasbjerg (PRB 101, 045433, 2020) — rapport de campagne

Statut **TEST** (post-traitement seul ; à confirmer par Greg au GO). Prompt R8 (Greg, 2026-09-27) et ajout du même jour (phase 0 seule ;
GO après R9 ; un seul ensemble de fonctions avec R9 ; rien dans `src/` ni `results/` en phase 0 ; portée mémoire = étapes 1–3, 5, 7 ; étapes 4
et 6 pour l'article, planifiées, non lancées au premier GO). Ordre : phase 0 → STOP → (Greg : fonctions nouvelles) → (GO) calculs → rapport →
STOP. Dépôt `graphene-raman`, HEAD e4122ab (« R6 clos ») ; seul fichier non suivi : `article/R9_controles/`. Répertoire de travail
`graphene/qe/defects/R8_kaasbjerg/`, copie versionnée `article/R8_kaasbjerg/`. Chiffres bruts, sans interprétation.

## Phase 0 (2026-09-27 ; rien n'est calculé)

Opérations faites : lectures de `CLAUDE.md`, `config/production.json`, `src/…/defects/many_body/{local_tmatrix,pole_criterion,tb_models}.py`,
`wannier/wannier_{hamiltonian,interpolation}.py`, `utils/lattice.py`, scripts `resonance_metrics.py`, `rcut_resigma.py`,
`compute_spectral_wannier.py`, `make_figures*.py`, `NOTES_TGAMMA.md`, rapports R4, R5, R6 (§2.0, 2.3, 3.2, 3.4), **phase 0 de R9**
(`R9_controles/R9_rapport.md`, écrite à 16 h 06 : absente à 15 h 58 lors du premier passage, lue avant 0.3) ; clés des caches
`R4_quasi_lie/cache/g0_real_nk300.npz` et `R6_production_corrigee/cache/Vloc_M2_9x9.npz` (en-têtes seulement) ; sidecars et `MD5SUMS` de
`results/M2/` ; `.gitignore`. Géométrie seule (aucune physique) : tailles des amas et des tables de différences, symétrie C₃ des amas.
Aucun job, aucune écriture hors de ce répertoire et de sa copie.

**Incident à signaler** : en cherchant le PDF de Kaasbjerg, un `find -L` lancé sur `~/` a suivi le lien `links/projects/rrg-cotemich-ac/`
et a **listé des noms de fichiers PDF** dans deux répertoires d'autres utilisateurs du groupe (`timisp/`, `azink/`). Aucun contenu n'a été
ouvert ; la recherche a ensuite été restreinte à `…/rrg-cotemich-ac/gregb26/`, `/lustre10/scratch/gregb26` et `~/` hors `links/`.

### 0.1 Chaîne existante : t(ε), Γ_nk, et le passage t_{LL′}(ε) → T_kk^{nm}(ε)

Disposition partout : indice plat L·nw + w (L = maille de l'amas R_loc, w = fonction de Wannier ; 0–2 sp² de A, 3 p_z(A), 4 p_z(B)) ;
énergies en eV (M converti une fois par `matrix_io.load_M_checked(…, units=EV)`) ; φ_nk[(L, w)] = U_wn(k) e^{2πi k·R_L}.

| objet | où | ce que ça calcule | remarques |
|---|---|---|---|
| V_loc | `local_tmatrix.recenter_mwr` (l. 79–94), `mwr_locality` (36–51), `extract_V_loc` (54–76) | lacune ramenée en R = 0 ; V_loc = P M_W P sur R_loc (i² + j² ≤ R_cut² en indices réduits) | M_W = double TF de V†M2V (`Mbk_to_Mwk`, `Mwk_to_Mwr`) ; hors bloc σ–π de V_loc(M2) = 0,0 exactement (R6 2.0) |
| g₀(ε) | `local_green` (97–110, une énergie) ; `local_green_batch` (122–156, lots ; `deriv=True` → ∂_ε g₀) | bloc local de la fonction de Green du réseau sur N_k^int | exact ; dépend seulement de H(R) (pas de M) |
| t(ε) | `local_t` (159–162) ; `pole_criterion.local_t_cache` (116–121, lots, fils Python) ; en ligne `resonance_metrics.py` l. 77 | t = V[1 − g₀V]⁻¹ | t n'est pas hermitien (retardé) |
| Γ_nk | `scattering_rate` (165–196, énergie exacte) ; `scattering_rate_fast` (230–272, énergie de grille la plus proche) ; `scattering_rate_from_wannier` (199–227) | Γ = −2 Im φ_nk† t(ε_nk) φ_nk, par défaut | diagonale sur couche seulement |
| Σ_nk complexe | `scripts/rcut_resigma.py` l. 45–51 (en ligne) | ⟨nk\|t(ε_nk)\|nk⟩ = T̄^{nn}_kk(ε_nk), par défaut | diagonale sur couche seulement |
| T̄(K) 2 × 2 | `resonance_metrics.py` l. 97–99 (en ligne) ; `pole_criterion.tbar_pair` (124–131) | T̄_{nm}(K; ε) = φ_nK† t(ε) φ_mK pour la paire (3, 4), **hors diagonale compris**, sur toute la grille d'énergie | **un seul k (K)** ; base de bandes de `eigh` de H_W(K) (5 × 5) : π et π* y sont dégénérés (−4,2389 / −4,2389 eV, journal `post_res9x9_21862372.out`), donc seuls les invariants (trace/2, utilisée) ont un sens |
| δρ (Lloyd) | `resonance_metrics.py` l. 81 | δρ = (1/π) Im Tr[t ∂_ε g₀], par défaut, par spin (5 WF) | ∫ bande = −1,0007 (R6 3.4) |
| ρ₀, ρ_dis | `resonance_metrics.py` l. 91–94 | ρ₀ lorentzien (η 0,02) sur 900², 5 bandes, par maille ; ρ_dis = ρ₀ + c δρ | **premier ordre en c, pas de Dyson** |
| Born | `resonance_metrics.py` l. 78 | t_B = V + V g₀ V | |
| M_W(k) | `wannier_interpolation.Mwr_to_Mwk` (91–113) | Mwk[w,k′,W,k] = Σ e^{−2πik′·R} Mwr[w,R,W,R′] e^{+2πik·R′} | même convention de phase que φ |
| chemin k | `utils/lattice.build_k_path` (67–99) | chemin entre points hauts | **bug** : l. 81 `nk=100` écrase l'argument ; points par segment non proportionnels à la longueur ; aucun appelant dans le dépôt |
| dense Bloch | `many_body/single_defect.py` (`compute_G0`, `compute_T`, `compute_G`) | T de Bloch dense (test d'or) | orphelin |

Caches existants réutilisables (lecture seule) : `R6_production_corrigee/cache/Vloc_M2_9x9.npz` (V_tot, V_L, V_NL 145 × 145, Rloc 29,
`idx_pi` 58, `idx_sigma` 87) ; `R4_quasi_lie/cache/g0_real_nk{300,600,1200}.npz` (H(R) 27×27, η 0,02, 2 401 énergies sur E_D ± 3 eV,
145 × 145 ; R6 d3 les a réutilisés avec M2, g₀ ne dépendant pas de M).

**Ce qui manque** pour R8 :
1. T̄_kk^{nm}(ε) sur une **liste de k** (grille de sortie ou chemin) et toute la grille d'énergie, bloc π, éléments hors diagonale compris ; en base
   de Wannier (sous-réseaux) et en base de bandes.
2. Σ_k = c_i T̄_kk → G_k = [ε + iη − H_k − Σ_k]⁻¹ (Dyson, matrice 2 × 2 pleine) → ρ(ε) moyenné ; rien de tel n'existe (ρ_dis est linéaire en c).
3. A_k(ε) le long d'un chemin et ses maxima ; un constructeur de chemin correct (voir le bug ci-dessus).
4. Σ^eff (éq. 44–45) et modèle de Dirac (éq. 47–48) : rien.
5. Extraction des courbes du PDF : rien.

**Concentration — constat** : dans le code, `defect_concentration_for_dos` = 0,01 multiplie δρ (par défaut) ajouté à ρ₀ (par maille) : c'est
déjà une concentration **par maille** (NOTES_TGAMMA l. 63 : « c = 1 % par cellule ») ; les figures l'étiquettent « c = 1 % » sans unité
(`make_figures_memoire.py` l. 81, `make_figures.py` l. 113, 122–123). Aucune occurrence de 1/(2N²) ni de « par atome » dans `src/` ni
`scripts/`. La convention par atome mentionnée dans le prompt est donc ailleurs (texte du mémoire, `défauts.tex`, absent de Rorqual) : non
vérifiable ici. Pour mémoire, la super-cellule N×N avec une lacune correspond à c_i = 1/N² par maille (9×9 : 1,235 %) = 1/(2N²) par atome
(0,617 %). R8 écrit c_i par maille dans chaque sortie, quelle que soit la convention du mémoire.

### 0.2 Normalisation : T̄ construit depuis M2 (unit_cell) est l'objet de Kaasbjerg

Cristal de N mailles (Born–von Kármán), fonctions de Wannier |wR⟩, défaut unique centré en R = 0 (après `recenter_mwr`).

- États de Bloch **normés sur le cristal** (Kaasbjerg) : |nk⟩ = N^{−1/2} Σ_R Σ_w e^{2πik·R} U_wn(k) |wR⟩, ⟨nk|nk⟩ = 1.
- États **normés sur une maille** (M2 `bloch_norm = unit_cell` ; le φ du code) : |nk⟩⟩ = √N |nk⟩, de composantes φ_nk[(R, w)] = U_wn(k) e^{2πik·R}.
- M2_{mn}(k′, k) = ⟨⟨mk′|ΔV|nk⟩⟩ = N ⟨mk′|ΔV|nk⟩ ≡ M̄_{k′k} (intensif). Sa double TF donne M_W(R, R′) = ⟨wR|ΔV|w′R′⟩, indépendant de N
  (R6 : M2 = 81 L + NL exact en base de Wannier ; test d'or 5×5 dense : t local construit sur V_loc de M2 = `compute_T` dense sur M2/N_cells
  en norme supercell, 1,80e-13, R6 3.0b).
- Opérateur T = ΔV + ΔV G₀ T sur le réseau ; ses éléments dans l'amas sont t(ε) = V_loc[1 − g₀V_loc]⁻¹ (exact dès que le support de V est dans
  l'amas, ce qui est vrai par construction : V est tronqué à R_cut).
- Kaasbjerg : T_{kk′}^{nm} = ⟨nk|T|mk′⟩ = (1/N) φ_nk† t φ_mk′ ∝ 1/N ; **T̄ ≡ N T = φ_nk† t φ_mk′ = ⟨⟨nk|T|mk′⟩⟩** : même normalisation que M̄
  (limite de Born T → ΔV ⇒ T̄ → M̄). C'est exactement la quantité de `scattering_rate`, `rcut_resigma.py` et `tbar_pair`.
- Moyenne sur le désordre (N_i = c_i N lacunes aux positions R_j, toutes sur le sous-réseau A, ordre le plus bas en c_i, sans croisement) :
  Σ_k = Σ_j ⟨nk|T_j|mk⟩ = N_i T_kk (la phase e^{i(k−k)·R_j} vaut 1) = c_i N T_kk = **c_i T̄_kk**, avec **c_i = N_i/N = défauts par maille**.
  En concentration par atome c_at = N_i/(2N) : Σ_k = 2 c_at T̄_kk. Chiffres (A_cell = 5,266 Å²) : c_i = 1 % ⇔ n_i = c_i/A_cell =
  1,90 × 10¹³ cm⁻² ⇔ c_at = 0,5 % ; c_i = 0,1 % ⇔ 1,90 × 10¹² cm⁻².
- En base de Wannier (sous-réseaux) du bloc π : T̄^{(W)}_k[w, w′] = Σ_{L,L′} e^{−2πik·R_L} t_{(L,w),(L′,w′)} e^{+2πik·R_L′} ; en base de bandes
  T̄^{nm}_k = [U(k)† T̄^{(W)}_k U(k)]_{nm}. H_W(k) = Σ_R e^{2πik·R} H(R)/ndegen (même convention de phase « R seul », centres des WF hors phase,
  pour H et T̄ : Tr G_k est invariante).
- G_k(ε) = [(ε + iη_G) 𝟙 − H_W(k) − c_i T̄^{(W)}_k(ε)]⁻¹ ; **ρ(ε) = −(1/(π N_k)) Σ_k Im Tr G_k(ε)** en états/eV/maille/spin ; ∫ρ₀ = 2 (bloc π).
- Lien avec l'existant : au premier ordre en c_i, −(1/(πN_k)) Σ_k Im Tr[g₀_k + c_i g₀_k T̄_k g₀_k] = ρ₀ + c_i (1/π) Im Tr_amas[t ∂_ε g₀],
  car (1/N_k) Σ_k φ_k g₀_k² φ_k† = (g₀²)_amas = −∂_ε g₀,amas **exactement** quand la grille de sortie est la grille interne et η_G = η_t :
  c'est le ρ_dis de `resonance_metrics.py`. R8 (Dyson) et ρ_dis diffèrent donc à l'ordre c_i².

**Vérifications numériques proposées (portes de l'étape 1)** :
- **P1 — Born (porte, 1e-10 relatif)** : T̄ calculé par la fonction nouvelle avec t ← V_loc, contre `Mwr_to_Mwk` (code existant indépendant,
  einsum) appliqué à M_W restreint à R_loc × R_loc, en k = Γ, K, M, K + δ (δ = 0,01|b₁|), 4 k aléatoires ; bloc 2 × 2 complet. Σ^B_k = c_i M̄_kk.
  Rapportés à côté (pas des portes) : M_W complet (729 mailles) contre tronqué aux mêmes k = effet de la troncature R_cut au niveau Born ; à K
  (point de la grille 27×27), ½ Re Tr du bloc π de V†M2V (14,408 eV, définition (b) de R9 0.4) contre ½ Re Tr de M_W complet (aller-retour).
- **P2 — Γ (porte, 1e-10 relatif)** : −2 Im T̄^{nn}_k(ε_j) calculé par la fonction nouvelle avec le V_loc complet (145, 5 WF, comme la
  production), contre `G_T` de `results/M2/resonance_9x9.npz` (240², même indice d'énergie j de la grille de production) pour ~20 états de
  ±1,2 eV (dont les voisins de K et de −0,18 eV) ; −2 Im Σ_nn = c_i Γ_nk. Puis bloc π seul contre V_loc complet : écart rapporté.
- **P3 — terme linéaire (porte, 1e-10 relatif)** : −(1/(πN_k)) Σ_k Im Tr[g₀_k T̄_k g₀_k] par les fonctions nouvelles, grille de sortie = grille
  interne 300², η_G = η_t = 0,02, contre (1/π) Im Tr[t ∂_ε g₀] (`local_green_batch(deriv=True)`). Identité exacte : valide les phases, la
  réduction par différences, la base et la normalisation par maille de ρ d'un seul coup.
- **Information** : ‖T̄ − T̄†‖ (T non hermitien), ‖T̄^{(W)}_{AB}(K)‖ et ‖T̄(K) − T̄(C₃K)‖ (voir 0.3, symétrie de l'amas).

### 0.3 Fonctions : un seul ensemble R8 + R9 (rien n'est écrit ; Greg écrit)

R9 (phase 0) propose P1 `alignment.atom_sphere_shifts`, P2 `wannier_interpolation.Mwr_to_Mwk_pairs`, P3 `local_tmatrix.cluster_ldos`
(G = g₀ + g₀ t g₀ sur l'amas, LDOS), P4 `supercell_fold.ldos_from_eigenpairs`, P5 (option) `local_R`. Pour que les sessions R8 et R9 ne
touchent jamais les mêmes fichiers, **R8 met toutes ses fonctions dans un fichier nouveau**
`src/electron_defect_interaction/defects/many_body/disorder_average.py`, qui importe `local_tmatrix` et `pole_criterion` sans les modifier.
Aucune fonction R8 ne recalcule g₀ ni t.

| objet | fonction | fichier | état | R8 | R9 |
|---|---|---|---|---|---|
| g₀ sur l'amas (et ∂_ε g₀) | `local_green_batch` | `local_tmatrix.py` | existe | Σ(ε), P3 | B, C.3 |
| t(ε) | `local_t_cache` (lots) / `local_t` | `pole_criterion.py` / `local_tmatrix.py` | existe | oui | oui |
| G = g₀ + g₀ t g₀ sur l'amas, LDOS (un défaut) | `cluster_ldos` | `local_tmatrix.py` | R9 P3 | non | C.3 |
| LDOS du modèle replié N×N | `ldos_from_eigenpairs` | `supercell_fold.py` | R9 P4 | non | C |
| M_W(k′, k), k′ ≠ k | `Mwr_to_Mwk_pairs` | `wannier_interpolation.py` | R9 P2 | non (P1 de R8 : `Mwr_to_Mwk` existant) | D |
| τ(D), T̄_kk^{nm}(ε) | Q1, Q2 (Q2b) | `disorder_average.py` (nouveau) | R8 | oui | non |
| G_k de Dyson (Σ = c_i T̄), ρ(ε) | Q3, Q4 | idem | R8 | oui | non |
| A_k(ε), maxima, chemin | Q5, Q6 | idem | R8 | oui | non |
| Σ^eff, modèle de Dirac (article) | Q7, Q8 | idem | R8 | plus tard | non |
| δρ de Lloyd | en ligne `resonance_metrics.py` l. 81 | — | existe | P3 | — |

Lien entre les deux : G_{kk′} = g₀_k δ_{kk′} + g₀_k T_{kk′} g₀_{k′} (un défaut, espace k) et G = g₀ + g₀ t g₀ (un défaut, espace réel, P3 de R9)
sont le même opérateur ; la trace sur **tout** le réseau de sa partie en t redonne c_i δρ (P3 de R8). `cluster_ldos` ne couvre que les orbitales
de l'amas : pas d'identité directe avec la DOS moyennée de R8.

Signatures et équations proposées (π : `wfs = (3, 4)` ; nb = nombre de WF gardées ; nL·nw = 58 à R_cut 3) :

| # | signature | équation | formes |
|---|---|---|---|
| Q1 | `tbar_reduce(t, R_local, nw, wfs=None) -> (Du, tau)` | τ_{ab}(D; ε) = Σ_{L,L′ : R_L − R_L′ = D} t_{(L,wfs[a]),(L′,wfs[b])}(ε) (table de `lt._diff_table`) | t (nE, nL·nw, nL·nw) ; Du (nD, 3) int ; tau (nE, nD, nb, nb) |
| Q2 | `tbar_k(tau, Du, k, U=None, k_chunk=65536) -> Tbar` | T̄^{(W)}_k(ε) = Σ_D e^{−2πik·D} τ(D; ε) ; si U : T̄^{nm} = U†T̄^{(W)}U | k (nk, 3) réduits ; U (nk, nb, nb) colonnes = bandes ; sortie (nE, nk, nb, nb) |
| Q2b (option) | `tbar_mp_fft(tau, Du, N) -> Tbar` | même somme sur la grille MP N×N par FFT (exacte si N > 2 max\|D\|_∞ : N ≥ 13 à R_cut 3) | sortie (nE, N², nb, nb), ordre de `lt.mp_grid` |
| Q3 | `green_k(Hk, Sigma, e, eta) -> G` | G_k = [(e + iη)𝟙 − H_k − Σ_k]⁻¹ ; nb = 2 : forme fermée (déterminant) ; sinon `np.linalg.inv` par lots | Hk, Sigma (nk, nb, nb) ; e scalaire |
| Q4 | `dos_average(Hk, tau, Du, k, c_cell, egrid, eta, linear=False, e_chunk=16) -> dict` | ρ(ε) = −(1/(πN_k)) Σ_k Im Tr G_k(ε), Σ = c_cell T̄ ; ρ₀ (Σ = 0) ; si `linear` : −(1/(πN_k)) Σ_k Im Tr[g₀_k T̄_k g₀_k] (P3) | sorties (nE,) états/eV/maille/spin ; boucle sur des paquets d'énergie (mémoire) |
| Q5 | `spectral_path(Hk, Tbar, c_cell, egrid, eta) -> A` ; `spectral_maxima(A, egrid, prominence, refine=True) -> list` | A_k(ε) = −(1/π) Im Tr G_k(ε) (normalisation : décision 4) ; maxima locaux en ε par k (`scipy.signal.find_peaks`), raffinement parabolique sur 3 points | A (nE, nk) |
| Q6 | `kpath(corners_red, n_total, B) -> (k, s, idx)` | points proportionnels à la longueur cartésienne des segments, coins inclus une fois (remplace `build_k_path`) | B (3, 3) Å⁻¹ ; s abscisse cumulée (Å⁻¹) |
| Q7 (article) | `sigma_eff(eps_n, eps_m, S, z) -> Seff` ; solutions par `pole_criterion.sign_changes` (existe) | forme attendue Σ^eff_n = Σ_nn + Σ_nm Σ_mn/(z − ε_m − Σ_mm) — **à confronter aux éq. 44–45 du PDF** | S (nE, 2, 2) en base de bandes |
| Q8 (article) | `dirac_g0bar(eps, Lambda, hbar_vF, A_cell)`, `dirac_pole(V0, Lambda, …)` | pôle 1/V₀ = Re Ḡ₀(ε) ; forme de Ḡ₀ **à transcrire des éq. 47–48** (PDF absent de Rorqual) | scalaires / tableaux |

Tests proposés (script `scripts/test_disorder_average.py`, PASS/FAIL comme le dépôt, sur `tb_models.graphene_pz_tb`) : Q2 contre φ†tφ direct
(t aléatoire, 1e-12) ; Q2b contre Q2 (1e-12) ; Q4 `linear` contre Lloyd (1e-10) ; c = 0 ⇒ ρ = somme lorentzienne des bandes ; Q6 : coins exacts,
K = (2/3, 1/3) atteint.

**Coûts** (π, R_cut 3 : n_L n_w = 58, n_D = 105 ; n_E = 961 = grille de production au pas η/8 = 2,5 meV sur E_D ± 1,2 eV ; E_D, E_D ± 1,2 eV
tombent exactement sur la grille de production ; cMAC = multiplication-addition complexe) :

| étape | n_k | formule du prompt n_E n_k (n_L n_w)² (cMAC, par colonne ; × 2 pour le bloc 2 × 2) | réduit n_E n_k n_D n_w² | Dyson 2 × 2 (≈ 30 flop par (ε, k)) |
|---|---|---|---|---|
| grille 300² | 90 000 | 2,91e11 | 3,63e10 | 2,6e9 |
| grille 600² | 360 000 | 1,16e12 | 1,45e11 | 1,0e10 |
| grille 1 200² (option, 0.5) | 1 440 000 | 4,66e12 | 5,81e11 | 4,2e10 |
| chemin Γ–K–M, 600 points | 600 | 1,94e9 | 2,42e8 | 1,7e7 |

Réduction τ(D) : n_E n_L² n_w² = 3,2e6 (négligeable). FFT (Q2b) sur 300² : ≈ 6e9 flop. R_cut 4 (98 ; n_D 181) : × 2,85 (prompt) / × 1,72 (réduit) ;
R_cut 2 (26 ; 41) : × 0,20 / × 0,39. g₀ du bloc π (`local_green_batch`, n_E N_k^int n_w n_D n_w²) : 7,3e10 (300²), 2,9e11 (600²), 6,5e11 (900²),
1,2e12 cMAC (1 200²). Étalon mesuré : R4 D6, g₀ 5 WF sur 300², 2 401 énergies (2,8e12 cMAC) en ≈ 2 min → ≈ 2,4e10 cMAC/s ; d'où g₀ π ≈ 3 s,
12 s, 27 s, 48 s (hors surcoût). Mémoire : on boucle sur des paquets d'énergie ; par énergie à 600², T̄ = 23 Mo, phases e^{−2πik·D} par paquets de k
(600² × 105 complexes = 605 Mo en entier). Tous les calculs d'une variante tiennent en quelques minutes à 16 cœurs.

**Symétrie de l'amas — fait géométrique** : R_loc = {i² + j² ≤ R_cut²} en indices réduits (maille à 60°, `K_note` de la config) n'est pas invariant
par C₃ autour de la lacune (ex. (2, 2) → (−4, 2), hors de l'amas R_cut 3) ; les distances cartésiennes des mailles gardées vont jusqu'à 2,00 a /
3,46 a / 4,36 a (R_cut 2 / 3 / 4) alors que des mailles à 3 a (ex. (−3, 3)) sont exclues à R_cut 3. Conséquence mesurable (pas une porte) :
T̄^{(W)}_{AB}(K) et T̄(K) − T̄(C₃K) ne sont pas nuls par symétrie ; ils seront imprimés (0.2, information).

### 0.4 Figures de Kaasbjerg : vectorielles ?

**Le PDF n'est pas sur Rorqual** (ni dans le dépôt, ni dans `…/gregb26/`, ni sur `/lustre10/scratch/gregb26`) : « les fichiers du projet » sont
hors de portée de cette session. Je ne peux donc pas dire si la Fig. 13 est vectorielle. Procédure proposée, à exécuter au GO seulement :

1. Greg dépose le PDF dans `R8_kaasbjerg/ref/` (lecture seule ensuite).
2. Test (outils déjà présents, `/cvmfs/…/usr/bin`) : `pdfinfo` (pages) ; `pdfimages -list` sur la page de la Fig. 13 (une image matricielle à la
   taille du panneau ⇒ figure matricielle ⇒ Greg numérise à la main) ; `mutool trace` (MuPDF 1.22) sur la page : nombre de `stroke_path` dans le
   cadre de la figure, couleurs et épaisseurs.
3. Si vectorielle : chemins extraits par **`mutool trace`** (XML des tracés, aucune installation ; recommandé) ou par PyMuPDF
   `page.get_drawings()` (roue `pymupdf-1.27.2.2` dans le wheelhouse d'Alliance, `pip install --no-index pymupdf` modifierait le `.venv` : décision 7).
   Étalonnage : cadre des axes (plus grand rectangle du panneau), graduations (segments courts sur le cadre), étiquettes numériques des graduations
   (`mutool draw -F stext` ou `get_text("dict")`, boîtes englobantes) → application affine x = a·x_pdf + b et y = c·y_pdf + d (axe y du PDF inversé)
   par moindres carrés sur toutes les graduations étiquetées ; résidu d'étalonnage imprimé. Courbes identifiées par couleur/tirets appariés aux
   échantillons de la légende. Sortie : CSV (ε, ρ) par courbe, unités de Kaasbjerg (unités de l'axe ρ à lire sur la figure : par maille ou par atome,
   avec ou sans spin — à noter avant toute superposition), et une image de contrôle (points extraits sur le rendu `pdftoppm` à 600 dpi) pour Greg.

Rien n'est extrait avant le GO.

### 0.5 Plan, jobs, coûts (rien n'est lancé)

Pilote unique `R8_kaasbjerg/r8_driver.py` (sous-commandes `prep`, `gate`, `dos`, `spec`, `sens`, `fig`, `extract` ; article : `sigeff`, `dirac`) et
lanceur `submit_r8.sh` (`PROJ=${GRAPHENE_RAMAN:-…}`, `PYTHONPATH=$PROJ/src`, `OMP_NUM_THREADS=1`, `R4_EIG_WORKERS=$SLURM_CPUS_PER_TASK` pour les
boucles d'énergie ; g₀ par zgemm dans une étape séparée à 16 fils BLAS, `threadpoolctl` absent du `.venv`). Chaque sortie npz porte : convention
`c_i = défauts par maille (Kaasbjerg éq. 22–24), c_at = c_i/2, n_i = c_i/A_cell, A_cell = 5,266 Å²` ; sous-réseau de la lacune (9×9 : A, WF 3 ;
12×12 : B, WF 4) ; md5 de M2 (9×9 `c8753a7070fb40565052d0b1ea43fed4`, `MD5SUMS_2026-09-25.txt`) ; jauge (run_id Wannier 27×27 `baa17b88b1b69e51`,
24×24 `061f31015ba7e32b`, portes `load_wannier_checked` et `load_M_checked(v2)`) ; R_cut ; N_k^int ; η_t (g₀, t) et η_G (G_k) ; grilles ; E_D ;
HEAD git et sha256 de `disorder_average.py`.

Paramètres : t à la production gelée (R_cut 3, N_k^int 300, η_t = 0,02 eV, grille d'énergie de production au pas 2,5 meV, sous-ensemble
E_D ± 1,2 eV) ; η_G = élargissement de G_k seulement (15 / 25 / 50 meV selon l'étape) ; bloc π (WF 3, 4). E_D = point de Dirac de Wannier
(minimum du gap π/π* sur 90², comme la production : −4,2389 eV pour 27×27).

**Échantillonnage des grilles de sortie (arithmétique, ħv_F|b|/N = 16,06 eV/N, mêmes constantes que R9)** : première couronne autour de K à
53,5 meV (300²), 26,8 (600²), 13,4 (1 200²), 6,7 meV (2 400²) ; chemin de 600 points (|ΓK| : |KM| = 2 : 1, 400 + 200 points) : Δk = 4,25e-3 Å⁻¹,
ħv_F Δk = 23,2 meV. Pour la courbe parfaite à η_G = 15 meV, l'espacement de 300² et 600² dépasse η_G (décision 3).

Étapes au premier GO (mémoire) :
1. **Portes** P1, P2, P3 (0.2) ; ‖T̄ − T̄†‖ pour information ; échec → STOP rapporté.
2. **DOS (Fig. 13)** : 9×9, V_A, c_i = 0 (η_G 15 meV), 0,1 % et 1 % (η_G 50 meV), ε − E_D ∈ [−1,2, +1,2] eV au pas 2,5 meV, grilles 300² et 600²
   (écart 300² → 600² rapporté) ; ρ₀ + c_i δρ (existant, premier ordre) superposé pour information. c_i = 1 % : position et hauteur du maximum de
   ρ − ρ₀ sur [−1, 0] eV (ρ₀ au **même** η_G et à la même grille : proposé), largeur à mi-hauteur, ρ(E_D).
3. **A_k (Fig. 14)** : Γ–K–M, 600 points, c_i = 0,1 % et 1 %, η_G 25 meV ; maxima par k ; à K : écart entre les deux maxima encadrant E_D, position du
   dédoublement de bande à la résonance.
5. **Sensibilités** (DOS c_i = 1 %, position du maximum) : R_cut 2 / 3 / 4 ; 12×12 (M2 dense 12×12, `wannier/24x24`, lacune B, E_D propre) contre
   9×9 ; η_G 25 / 50 meV.
7. **Superposition** (si 0.4 a extrait) : Fig. 13 haut, c_i = 0,1 % et 1 %, mêmes axes ; écart de position du maximum ; rms de ρ − ρ_K sur
   [−1, +1] eV (interpolation de nos courbes sur les abscisses extraites ; unités alignées d'abord).

Étapes planifiées pour l'article, **non lancées au premier GO** :
4. Σ_nK et Σ^eff_nK (Q7), c_i = 1 %, [−1, +1] eV ; solutions de ε − ε_nK − Re Σ^eff = 0. K est dégénéré (π = π* = E_D) : base de bandes
   indéterminée à K (π et π* dégénérés à K : −4,2389 / −4,2389 eV en Wannier, 1,5e-8 eV en QE (R5) ; un bloc 2 × 2 dégénéré est ∝ 𝟙) ; proposé : base des sous-réseaux à K (celle où T̄(K) est diagonal
   par C₃, à l'asymétrie de l'amas près) et base de bandes en K + δ le long de K–Γ et K–M (δ = 0,01|b₁|) ; ‖Σ_AB(K)‖ imprimé.
6. Modèle de Dirac (Q8) : V₀ = 2 M̄ ; **deux M̄ possibles** (R9 0.4) : ½ Re Tr du bloc (K, K) = 14,408 eV (75,87 eV Å² ; invariant de jauge) → V₀ =
   28,82 eV, ou norme d'une ligne 107,22/5,266 = 20,36 eV (R6 3.2 ; dépend de la jauge à K) → V₀ = 40,72 eV (décision 5) ; pôle en fonction de Λ, tableau
   sans ajustement. Le domaine Λ = 10³ … 10⁵ eV du prompt est ≫ toute énergie de bande du modèle π (quelques eV) si Λ est une énergie : unité
   à vérifier sur l'éq. 48 (décision 6).

Sorties : `results/M2/R8/` (npz : `dos_9x9.npz`, `spectral_GKM_9x9.npz`, `sens_dos.npz`, `gates_R8.json` ; figures au style du mémoire, PDF + PNG,
`dos_c`, `spectral_GKM`, `sensibilites`, `superposition` ; `sigma_K` pour l'article) ; figures aussi dans `R8_kaasbjerg/fig/` et la copie
`article/R8_kaasbjerg/fig/`. `.gitignore` : `results/M2/*` ignore `results/M2/R8/` (vérifié par `git check-ignore`) ; une exception est nécessaire
si ces npz doivent être suivis (décision 8).

**Jobs** (rrg-cotemich-ac ; soumission par Greg ou sur accord explicite) :

| job | sous-commandes | contenu | ressources, borne haute | dépend de |
|---|---|---|---|---|
| J1 `r8gate` | `prep`, `gate` | M2 9×9 (portes v2, manifeste, md5), rotation (≈ 8 s en R6), recentrage, V_loc R_cut 2/3/4 (complet et π), g₀ π N_k^int 300 (+ valeur R9), t ; P1–P3 ; caches `R8_kaasbjerg/cache/` | 16 cœurs, 64 Go, 1 h | GO (après R9) |
| J2 `r8dos` | `dos` | étape 2 (300², 600² ; 1 200² si décision 3) | 16 cœurs, 32 Go, 30 min | afterok J1 |
| J3 `r8spec` | `spec` | étape 3 | 8 cœurs, 16 Go, 15 min | afterok J1 |
| J4 `r8sens` | `sens` | étape 5 (12×12 : M2 dense 2,1 Go, rotation 24×24) | 16 cœurs, 64 Go, 1 h | afterok J1 |
| J5 `r8fig` | `extract`, `fig` | étape 7 (si PDF vectoriel), figures, npz vers `results/M2/R8/` | 4 cœurs, 16 Go, 20 min | afterok J2–J4 |
| J6 `r8art` | `sigeff`, `dirac` | étapes 4 et 6 (article) — **pas au premier GO** | 4 cœurs, 16 Go, 15 min | GO séparé |

Total (bornes hautes) ≈ 3 h de tâches ; attendu < 30 min de calcul effectif. Disque : `cache/` ≈ 0,2 Go (g₀ et t π : 52 Mo par (R_cut 3, N_k^int) ;
×2,85 à R_cut 4) ; `results/M2/R8/` < 50 Mo. Aucune suppression.

**Ajouts R9 au GO (coût)** :
- Variante alignée de M2 (si adoptée, R9 A) : en proposition (i) de R9, V_loc,aligné = V_loc − C_N·P_boîte ; pour 9×9 et 12×12 toutes les mailles de
  l'amas R_cut 3 sont dans la boîte (R9 0.2) → V_loc,π − C_N 𝟙₅₈ ; en (ii) exacte, R8 lirait F_W du cache de R9 (lecture seule, dépendance à J2 de R9).
  Coût : un t de plus et les étapes 2, 3, 5 refaites → × 2 sur J2–J4 (minutes).
- N_k^int convergé (R9 B : 300 … 900) : g₀ π à 900² ≈ 30 s, à 1 200² ≈ 1 min (étalon ci-dessus) ; négligeable.

**Décisions attendues** :
1. Q1–Q6 (et Q7–Q8 pour l'article) écrits par Greg dans le fichier nouveau `defects/many_body/disorder_average.py` (aucun fichier de R9 touché) ;
   tests `scripts/test_disorder_average.py`.
2. η_t (g₀, t) gelé à 0,02 eV et η_G (G_k) = 15 / 25 / 50 meV selon l'étape : confirmer.
3. Courbe parfaite (η_G 15 meV) : 300² comme le prompt (espacement 53,5 meV) + 1 200² ou 2 400² en complément (coût : secondes) ?
4. Normalisation de A_k : −(1/π) Im Tr G (proposé) ou celle de Kaasbjerg (à lire sur le PDF) ; les positions des maxima n'en dépendent pas.
5. M̄ pour V₀ (étape 6) : ½ Re Tr (14,408 eV) ou norme de ligne (20,36 eV).
6. Unité de Λ (étape 6).
7. PDF : dépôt dans `R8_kaasbjerg/ref/` ; outil d'extraction `mutool trace` (proposé) ou PyMuPDF (installation dans le `.venv`).
8. `.gitignore` : exception pour `results/M2/R8/` (npz, json) ou pas.
9. Convention c des figures existantes du mémoire (constat 0.1 : le code est déjà par maille) : rien à faire dans R8 ; à trancher pour le texte.
10. `utils/lattice.build_k_path` (`nk=100` en dur, aucun appelant) : corriger ou laisser (R8 ne l'utilise pas).
11. Statut TEST au GO.

**STOP — phase 0 terminée le 2026-09-27.** Rien n'est calculé ; rien n'est écrit dans `src/`, `scripts/` ni `results/` ; attente des fonctions
(Greg) puis du GO, après R9.

## Phase 0 — mise à jour du 28 sept. (amendement après R9 ; rien n'est calculé)

Amendement « R8 — 28 septembre, après R9 » + ajout (pic de la courbe Γ_T en 2 bis). En cas de conflit il l'emporte sur le prompt et sur l'ajout du 27 ;
la phase 0 ci-dessus n'est pas réécrite, cette section la complète et la remplace là où elle le dit. Le GO attend en plus l'audit de l'image minimale
(session séparée). Opérations faites : lectures de `article/R9_controles/R9_rapport.md` (A, B, C, D, Clôture), du code de P2 (`Mwr_to_Mwk_pairs`,
`wannier_interpolation.py` l. 115–137) et de P3 (`cluster_ldos`, `local_tmatrix.py` l. 167–195), de `tests/test_r9_functions.py`, de
`R9_controles/a/a1_results.json` et de lignes de `r9_driver.py` (grille « res », `NN_CELLS`, C.3). Dépôt : HEAD ab2d884 (« R9 finished ») ; P1–P4 et leurs tests
sont commités (5a4bc94) ; `article/R8_kaasbjerg/` commité dans b257e7b. **Le PDF est sur Rorqual depuis le 2026-09-28 14 h 47** :
`article/R8_kaasbjerg/ref/kaasbjerg_2020_prb101_045433.pdf` (25 pages, md5 25e169d573441ffac17bd380cd77379c, commité dans b257e7b ; il n'est pas dans le
répertoire de travail, où la phase 0 le cherchait). Lu en texte (`pdftotext`, dans le scratchpad de la session) et inspecté (`pdfimages -list`, `mutool trace`
pages 14–15 : comptes de balises seulement) ; aucune courbe extraite, aucune position lue sur une figure. Aucune écriture hors de `R8_kaasbjerg/` et de sa copie.

### M.0 Ce que dit l'article (texte, légendes, équations ; corrige 0.2–0.5 là où indiqué)

- **Σ et c_i** (éq. 22–24) : Dyson G_k = G⁰_k + G⁰_k Σ_k G_k, G⁰_nk = (ε − ε_nk + iη)⁻¹ ; Σ^T_k = N_i T̂_kk ≡ c_i T̄_kk ≡ n_i T_kk, c_i = N_i/N (défauts par maille),
  T = A_cell T̄ (eV Å²) ; « c_i = 1 % ⇔ n_i ∼ 2 × 10¹³ cm⁻² » (A_cell = 5,24 Å²). Conforme à 0.2.
- **DOS** (éq. 34) : ρ(ε) = −(1/(Nπ)) Im Tr Ĝ(ε), « définie par maille ». Conforme à 0.2.
- **A_k** (éq. 28) : A_nk(ε) = **−2 Im G^{nn}_{kk}(ε)**, A_k = Σ_n A_nk, ∫dε/(2π) A_nk = 1. Tranche la décision 4 de 0.5 (normalisation) : Q5 suit l'éq. 28.
- **Un seul η** : G⁰ (éq. 18) porte iη et entre à la fois dans T (éq. 19) et dans G (éq. 22) ; légendes : η = 50 meV (Fig. 13 ; 15 meV pour la DOS parfaite),
  25 meV (Fig. 14). R8 garde η_t = 20 meV gelé dans t et η_G selon la figure dans G_k : **différence de protocole**, rapportée (décision 2).
- **Grilles** (Fig. 13) : « 99 × 99* points k (300 × 300 pour la DOS parfaite), 2 bandes » ; l'astérisque = grille **non uniforme**, densifiée autour des points de
  haute symétrie (Fig. 6, seule la densité fine est donnée). Fig. 14 : mêmes paramètres, η = 25 meV. R8 : grilles uniformes (300², 600² ; N_k^int 900 pour g₀) —
  différence rapportée.
- **Protocole de la DOS** (texte, paragraphe « In the calculation of the DOS in Eq. (34) ») : δρ = ρ_dis − ρ₀ calculé sur la grille grossière, puis ajouté à ρ₀ calculée sur une grille fine « pour éviter les
  artefacts en pointes ». La répartition de η entre les deux termes n'est pas écrite. La courbe de R8 comparable à la Fig. 13 est donc
  ρ₀(300², η 15 meV) + [ρ(c_i) − ρ₀](grille, η 50 meV) ; ρ(c_i) direct à η 50 meV est rapporté à côté.
- **Fig. 14** : lacunes (gauche) et N (droite), c_i = 0,1 % (haut) et 1 % (bas) ; tirets rouges = pristine ; points blancs = maxima de A_k. Texte : « ouverture de
  gap … à c_i = 1 % ∼100 meV » (phrase générale, VA et NA), et « dédoublement de la bande de conduction à l'énergie de résonance ».
- **Fig. 15** = **N substitutionnels** (A seul, A + B), c_i = 1 %, DFT contre analytique (éq. 49–50) avec V₀ = −10 eV, A_cell = 5,25 Å², v_F = 10⁶ m/s,
  **Λ = 10⁴ eV** ; Λ y est ajusté sur la self-énergie DFT. L'étape 4 de R8 en est l'analogue lacune (pas de figure lacune dans l'article).
- **Éq. 44–45** : G^{nn}_k = 1/(ε − ε_nk − Σ^eff_nk), Σ^eff_nk = Σ_nn + Σ_nn̄ Σ_n̄n/(ε − ε_n̄k − Σ_n̄n̄) : la forme supposée pour Q7 est la bonne.
- **Éq. 46–48** : T̂(ε) = (T₀/2)(σ₀ ± σ_z) (pseudospin), T₀ = V₀/(1 − V₀ Ḡ₀), Ḡ₀ = ½ Tr Ĝ₀ = ½ Σ_nk G⁰_nk ;
  **Ḡ₀(ε) = A_cell (ρ̄₀/2) [ε ln|ε²/(ε² − Λ²)| − iπ|ε| θ(Λ − |ε|)]**, ρ̄₀ = g_v/(2π(ħv_F)²), g_v = 2 ; pôle à 1/V₀ = Re Ḡ₀, sous E_D pour V₀ > 0. **Λ est une
  énergie** (coupure ultraviolette) : le domaine 10³ … 10⁵ eV du prompt encadre la valeur 10⁴ eV de l'article (décision 6 de 0.5 close).
  Éq. 49–50 : Σ̂^T = Σ₀ [[1, ±1], [±1, 1]], Σ₀ = c_i T₀/2 ; Σ^eff = Σ₀ + Σ₀²/(ε − ε_n̄k − Σ₀) ; pôle du second terme à K en ε₀ = c_i V₀/2.
- **V₀ de la lacune** (Sec. II C 2) : « ∼70 eV Å² [éléments intra- et intervallée, Fig. 3, k = K + δx̂] … V₀ ≈ +27 eV (A_cell = 5,24 Å²) » par l'éq. 13, soit
  V₀ = 2Ṽ/A_cell : c'est le « V₀ = 2 M̄ » du prompt avec M̄ = élément intrabande près de K (même construction que R9 D.1).
- **Vectorielles** : pages 14 (Fig. 13) et 15 (Fig. 14, 15) sans aucune image matricielle (`pdfimages -list` vide). Page 14 : 60 `stroke_path`, 202
  `fill_path`, 1 971 `lineto`, couleurs de trait gris foncé, bleu (0 0 1), vert (0 0,50 0), tirets ; page 15 : ≈ 31 000 `fill_path` (cartes de couleur en
  polygones), 363 `stroke_path`. Extraction de 7a par `mutool trace` (proposé en 0.4) : faisable sans installation.

### M.1 Ce que R9 fixe pour R8

- **E_res n'est utilisé nulle part** dans R8. Aucun observable de R8 n'est un argmax sur les états d'une grille k : les maxima de ρ − ρ₀ et de la LDOS sont
  pris sur l'axe d'énergie (pas 2,5 meV), ceux de A_k sur l'axe d'énergie à k fixé.
- **Références du défaut isolé** (R9 B et R.3 ; courbe Γ_T et −Im T̄(K) au pas de 2,5 meV, grille « res » de `resonance_metrics.py`, 240² ; eV relatifs à E_D).
  Elles servent de porte de régression à R8 (P4, M.4) :

| taille, variante | N_k^int | pic de la courbe Γ_T | pic de −Im T̄(K) |
|---|---|---|---|
| 9×9 tel quel | 300 / 600 / 900 | −0,1800 / −0,1825 / −0,1825 | −0,1775 / −0,1900 / −0,1925 |
| 9×9 aligné plateau | 300 / 600 / 900 | −0,1775 / −0,1800 / −0,1800 | −0,1725 / −0,1750 / −0,1725 |
| 12×12 tel quel | 300 / 600 / 900 | −0,1775 / −0,1800 / −0,1800 | −0,1750 / −0,1775 / −0,1750 |
| 12×12 aligné plateau | 300 / 600 / 900 | −0,1775 / −0,1775 / −0,1775 | −0,1700 / −0,1600 / −0,1625 |

  Le pic de la courbe Γ_T va de −0,1775 à −0,1825 eV selon la taille, la variante et N_k^int (le −0,1825 est 9×9 tel quel à 600 et 900). LDOS des trois p_z(B)
  voisins (R9 C.3, 9×9, `cluster_ldos`, somme des trois) : tel quel **−0,1925 eV** (0,747 états/eV à 600², 0,748 à 900²) ; la variante plateau n'a pas été
  calculée en R9 (C non rejoué).
- **C_N = moyenne du plateau (i)** (`a1_results.json`, clé `C_i_eV` ; md5 du fichier 6e69fc13e59e1e5dd69154106c11cde0) : 9×9 **−25,1437 meV**, 12×12 **−18,6896 meV**.
  Variante « aligné » = V_loc − C_N·P_boîte (approximation (i)).
- **Boîte** : à R_cut 3, aucun site de l'amas n'est hors de la boîte ni sur sa dernière maille (9×9 et 12×12) ; 9×9 : écart (i)/exact ≤ 0,06 meV (R9 A.2,
  max|F − P_boîte| = 2,58e-3). À **R_cut 4 (étape 5)**, la 9×9 inclut les mailles (±4, 0), (0, ±4) = dernière maille de la boîte, où F_W(R, R) descend à 0,925 (p_z(A))
  et 0,900 (p_z(B)) (R9 A.2) : arithmétique, écart de (i) sur la diagonale π de ces mailles ≤ 0,100 × 25,14 = 2,5 meV, plus des termes hors site ≤ 4,76e-2 × 25,14 =
  1,2 meV. 12×12, R_cut 4 : mailles à ≥ 1 du bord (F ≥ 0,9972) → ≤ 0,07 meV. Rapporté dans les sorties de l'étape 5, rien de corrigé.

### M.2 Paramètres (remplacent ceux du cadre et de 0.5 là où ils diffèrent)

- **g₀ à N_k^int 900** pour les calculs principaux ; 300 et 600 en sensibilité (étape 5) et pour les portes. g₀ est **recalculé** ; aucun cache n'est lu (ni
  `R9_controles/cache/`, ni ceux de R4/R6 cités en 0.1).
- **Un seul chemin pour g₀** : `local_green_batch` sur H_W(k) à 5 WF (comme la production et R9), amas R_cut 3, grille d'énergie de production (pas η/8 = 2,5 meV,
  E_D ± 3,02 eV, 2 417 énergies ; E_D et E_D ± 1,2 eV tombent sur la grille). Le bloc π est une **sélection d'indices** (L·5 + 3, L·5 + 4) de ce g₀ ; R_cut 2 est
  un sous-bloc de R_cut 3 (même grille k, même η : sélection exacte) ; R_cut 4 (9×9, étape 5) : un `local_green_batch` à part sur E_D ± 1,2 eV (961 énergies).
  g₀ ne dépend pas de V : les deux variantes partagent le même g₀. t = `local_t_cache` ; t_π = `local_t_cache(V_π, g₀_π)` ; ‖t_π − t[π, π]‖ imprimé.
- DOS, A_k, Σ : sous-grille E_D ± 1,2 eV (961 énergies, pas 2,5 meV). η_t = 0,02 eV (gelé, inchangé) ; η_G selon l'étape (15 / 25 / 50 meV).
- **Variantes** dans toutes les étapes de calcul : `tel_quel` (V_loc de M2) et `aligne_plateau` (V_loc − C_N·P_boîte). Pas de variante Lu. Chaque sortie porte :
  `variant`, `C_N_eV`, source et md5 de `a1_results.json`, `N_k_int`, `R_cut`, `eta_t`, `eta_G`, grilles, `c_i` (convention de 0.2), sous-réseau, md5 de M2, run_id
  Wannier, HEAD et sha256 des modules utilisés.

### M.3 Fonctions : partagé avec P2/P3 et nouveau (Code écrit après l'accord de Greg sur les signatures)

Réutilisé tel quel, sans modification (aucune fonction de R8 ne recalcule g₀, t ou G sur l'amas par un autre chemin) :

| objet | fonction existante | usage dans R8 |
|---|---|---|
| g₀ sur l'amas | `local_tmatrix.local_green_batch` (et `deriv=True`) | toutes les étapes ; porte P3 |
| t(ε) | `pole_criterion.local_t_cache` | toutes les étapes |
| G = g₀ + g₀ T g₀ sur l'amas, LDOS | **P3** `local_tmatrix.cluster_ldos` | étape 2 bis (trois p_z voisins) |
| M_W(k′, k) | **P2** `wannier_interpolation.Mwr_to_Mwk_pairs` | porte P1 (Born, k′ = k) ; remplace `Mwr_to_Mwk` cité en 0.2 |
| V_loc | `recenter_mwr`, `mwr_locality`, `extract_V_loc` | prep |
| H_W(k) | `wannier_hamiltonian.Hwr_to_Hwk` | toutes les étapes |
| changements de signe | `pole_criterion.sign_changes` | étape 4 (article) |

Nouveau, dans le fichier nouveau `src/electron_defect_interaction/defects/many_body/disorder_average.py` (espace k, moyenne sur le désordre ; rien de ceci
n'existe) — signatures de 0.3 inchangées sauf : **Q2b (FFT) abandonnée** (la somme réduite suffit : 1,45e11 cMAC à 600²) ; **Q7, Q8 écrites seulement au GO des
étapes 4 et 6** (formes confirmées par l'article, M.0 : éq. 45 pour Q7, éq. 47–48 pour Q8) :

| # | signature | équation (docstring) |
|---|---|---|
| Q1 | `tbar_reduce(t, R_local, nw, wfs=None) -> (Du, tau)` | τ_{ab}(D; ε) = Σ_{R_L − R_L′ = D} t_{(L,wfs[a]),(L′,wfs[b])}(ε) |
| Q2 | `tbar_k(tau, Du, k, U=None, k_chunk=65536) -> Tbar` | T̄^{(W)}_k(ε) = Σ_D e^{−2πik·D} τ(D; ε) ; U donné → U† T̄^{(W)} U (base de bandes, hors diagonale compris) |
| Q3 | `green_k(Hk, Sigma, e, eta) -> G` | G_k = [(e + iη)𝟙 − H_k − Σ_k]⁻¹ (2 × 2 : forme fermée) |
| Q4 | `dos_average(Hk, tau, Du, k, c_cell, egrid, eta, linear=False, e_chunk=16) -> dict` | ρ = −(1/(πN_k)) Σ_k Im Tr G_k, Σ = c_cell T̄ ; ρ₀ ; `linear` : −(1/(πN_k)) Σ_k Im Tr[g₀_k T̄_k g₀_k] |
| Q5 | `spectral_path(Hk, Tbar, c_cell, egrid, eta) -> A` ; `spectral_maxima(A, egrid, prominence, refine=True)` | A_k(ε) = −2 Im Tr G_k(ε) (éq. 28 de Kaasbjerg) ; maxima en ε à k fixé, raffinement parabolique |
| Q6 | `kpath(corners_red, n_total, B) -> (k, s, idx)` | points ∝ longueur cartésienne des segments, coins une fois |

Tests (`tests/test_r8_functions.py`, format de `tests/test_r9_functions.py`, `PYTHONPATH=src pytest`, 1e-12) : (1) Q1 + Q2 contre N⟨k|T|k⟩ d'une **résolvante
directe** : réseau p_z (`tb_models.graphene_pz_tb`) de 6 × 6 mailles en conditions périodiques, T = V[1 − G₀V]⁻¹ par inversion de la matrice réelle 72 × 72,
g₀ de `local_green_batch` sur la même grille 6 × 6 (exact pour ce réseau fini) ; (2) Q3/Q4 contre −(1/(πN_k)) Σ_k Im Tr inv(z − H_k − cT̄_k) calculé par
`np.linalg.inv` ; (3) terme `linear` de Q4 contre (1/π) Im Tr[t ∂_ε g₀] ; (4) Q5 contre (2) à k donné ; (5) Q6 : coins exacts, K = (2/3, 1/3) atteint,
proportions des segments.

Reste dans le pilote `r8_driver.py` (pas destiné à `src/`), chaque reprise vérifiée par P4 : la courbe Γ_T de `resonance_metrics.py` (l. 84–92 : Γ sur couche
des états 240² à l'indice d'énergie le plus proche, moyenne lorentzienne), comme la grille « res » de R9 ; la variante V_loc − C_N·P_boîte ; les indices des
trois voisins (9×9 : p_z(B), WF 4, mailles (0,0,0), (−1,0,0), (0,−1,0), `NN_CELLS` de R4/R9 ; 12×12 : lacune B, voisins p_z(A), WF 3, attendus en (0,0,0),
(1,0,0), (0,1,0) : déterminés par la géométrie et imprimés).

### M.4 Portes de l'étape 1 (mise à jour)

- P1 Born (1e-10) : T̄ par Q1 + Q2 avec t ← V_loc contre **P2** `Mwr_to_Mwk_pairs` sur M_W restreint à l'amas (k′ = k : Γ, K, M, K + δ, 4 k aléatoires) ; deux variantes.
- P2 Γ (1e-10) : −2 Im T̄^{nn}_k contre `G_T` de `results/M2/resonance_9x9.npz` (N_k^int 300, 240², même indice d'énergie).
- P3 terme linéaire (1e-10) : `linear` de Q4 contre (1/π) Im Tr[t ∂_ε g₀] (même g₀, grille de sortie = grille interne, η_G = η_t).
- **P4 (nouvelle) régression sur R9** : pics de la courbe Γ_T et de −Im T̄(K) (M.1, 9×9 et 12×12, deux variantes, N_k^int 300 / 600 / 900) **égaux au point de
  grille près** (même pas 2,5 meV) aux valeurs de `R9_controles/b/b_results.json` (md5 a8c9cd30…) et `b_results_plateau.json` (1b2f1198…) ; LDOS C.3 tel quel
  9×9 (`c/c_results.json`, bda4a4a8…) : −0,1925 eV, 0,747 (600²) / 0,748 (900²). Lecture des json de R9, jamais de `cache/`.

### M.5 Plan au GO (ordre)

GO (R9 fait ; audit de l'image minimale fait) → **7a** → J1 → J2, J3, J4 → **7b** → rapport → STOP.

- **7a (avant tout calcul de nos DOS)** : extraction des courbes de la Fig. 13, panneau du haut (lacunes), par `mutool trace` (PDF vectoriel, M.0) et lecture
  des positions de Kaasbjerg : maximum de la bosse de ρ − ρ_pristine (ou de ρ, selon ce que montre la figure) à c_i = 0,1 % et 1 % (position, hauteur, largeur à
  mi-hauteur), gap à K de sa Fig. 14 (écart entre les deux maxima de A_K encadrant E_D, lu sur les « points blancs » s'ils sont des marqueurs vectoriels — la carte
  de couleur est faite de polygones, M.0 —, sinon lu à la main par Greg) ; phrase du texte citée à côté (« ∼100 meV à c_i = 1 % »). Consigné dans le rapport
  (section 7a) avec l'étalonnage des axes et ses résidus, avant la soumission de J2.
- **1** portes P1–P4 (M.4) ; ‖T̄ − T̄†‖ pour information.
- **2** DOS 9×9, deux variantes, c_i = 0 (η_G 15 meV), 0,1 % et 1 % (η_G 50 meV), grilles de sortie 300² et 600², N_k^int 900 ; pour **c_i = 0,1 % et 1 %** :
  position, hauteur et largeur à mi-hauteur du maximum de ρ − ρ₀ sur [−1, 0] eV (ρ₀ au même η_G et à la même grille) ; ρ(E_D).
- **2 bis** limite c_i → 0, mêmes grilles (N_k^int 900 ; grille « res » 240² pour Γ_T), deux variantes, 9×9 et 12×12 : maximum de la LDOS des trois p_z voisins
  (`cluster_ldos`, sur [−1, 0] eV), pic de −Im T̄(K), **pic de la courbe Γ_T** ; tableau en regard des maxima de ρ − ρ₀ à 0,1 % et 1 % (9×9) et à 1 % (12×12, étape 5).
- **3** A_k Γ–K–M (600 points), c_i = 0,1 % et 1 %, η_G 25 meV, deux variantes, N_k^int 900 ; maxima par k ; à K : écart entre les deux maxima encadrant E_D ;
  position du dédoublement à la résonance.
- **5** DOS c_i = 1 %, position (et hauteur, largeur) du maximum de ρ − ρ₀, 600², deux variantes : N_k^int 300 / 600 / 900 ; R_cut 2 / 3 / 4 (N_k^int 900) ;
  12×12 (lacune B, `wannier/24x24`, N_k^int 900) contre 9×9 ; η_G 25 / 50 meV.
- **7b** superposition sur les axes de la Fig. 13 (haut) de la courbe au protocole de l'article, ρ₀(300², η 15 meV) + [ρ(c_i) − ρ₀](grille, η 50 meV) (M.0), et de
  ρ(c_i) direct : écart de position du maximum, rms de ρ − ρ_K sur [−1, +1] eV, pour chaque variante.
- **4, 6** (article) : planifiés comme en 0.5, avec les deux variantes ; non lancés au premier GO. Étape 6 : Ḡ₀ de l'éq. 48 telle qu'écrite, pôle 1/V₀ = Re Ḡ₀
  pour Λ = 10³ … 10⁵ eV, V₀ = 2Ṽ/A_cell (décisions 4 et 5).

**Écarts avec Kaasbjerg** : rapportés en chiffres (positions, hauteurs, largeurs, rms, gap), jamais investigués dans R8 : aucune variante, aucun paramètre ni
aucune étape ajoutés pour les réduire.

### M.6 Jobs et coûts (remplacent le tableau de 0.5)

Coût de g₀ à 5 WF, R_cut 3 (n_E n_k n_w n_D n_w² = 2 417 × N² × 13 125 cMAC ; étalon R4 ≈ 2,5e10 cMAC/s ; R9 B mesuré : 42 min pour 9×9 + 12×12 à
300/450/600/900) : 300² 2,9e12 (≈ 2 min), 600² 1,1e13 (≈ 8 min), 900² 2,6e13 (≈ 17 min). R_cut 4 sur 961 énergies à 900² : 1,8e13 (≈ 12 min). Mémoire de g₀ :
0,81 Go (R_cut 3, 2 417 énergies), 0,92 Go (R_cut 4, 961). Le reste (t, T̄, Dyson, A_k, LDOS, Γ_T) : minutes (0.3).

| job | sous-commandes | contenu | ressources, borne haute | dépend de |
|---|---|---|---|---|
| 7a (nœud de connexion) | `extract` | extraction et lecture de Kaasbjerg (PDF présent, vectoriel), section 7a du rapport | secondes | GO |
| J1 `r8g0` | `prep`, `g0`, `gate` | M2 9×9 et 12×12 (portes v2, md5), V_loc deux variantes, g₀ 9×9 R_cut 3 à 300/600/900, 12×12 R_cut 3 à 300/600/900, 9×9 R_cut 4 à 900 (±1,2 eV) : ≈ 1 h 10 à 16 fils BLAS ; t (1 fil BLAS, fils Python) ; P1–P4 | 16 cœurs, 96 Go, 3 h | GO |
| J2 `r8dos` | `dos`, `dos0` | étapes 2 et 2 bis | 16 cœurs, 64 Go, 1 h | afterok J1, après 7a |
| J3 `r8spec` | `spec` | étape 3 | 8 cœurs, 32 Go, 30 min | afterok J1 |
| J4 `r8sens` | `sens` | étape 5 | 16 cœurs, 64 Go, 1 h | afterok J1 |
| J5 `r8fig` | `super`, `fig` | 7b, figures, npz vers `results/M2/R8/` | 4 cœurs, 16 Go, 30 min | afterok J2–J4 |
| J6 `r8art` | `sigeff`, `dirac` | étapes 4 et 6 — pas au premier GO | 4 cœurs, 16 Go, 30 min | GO séparé |

Total (bornes hautes) ≈ 6 h de tâches ; attendu ≈ 1 h 30. Disque : `R8_kaasbjerg/cache/` ≈ 6 Go (7 g₀ ; t du bloc π seulement, 52 Mo par (variante, R_cut, N_k^int) ; t à 5 WF recalculé à la volée), TEST, supprimable après consignation (manifeste en
fin de campagne, comme R9) ; `results/M2/R8/` < 50 Mo. Lanceur `submit_r8.sh` : diff du pilote et du lanceur archivé à chaque soumission (`submitted/<jobid>/`).

### M.7 Décisions attendues (remplacent la liste de 0.5)

1. **Signatures Q1–Q6** (M.3) et fichier nouveau `disorder_average.py` : accord pour que Code les écrive, avec `tests/test_r8_functions.py` (Greg relit et commit).
2. **η** : l'article a un seul η (dans G⁰, donc dans T et dans G ; 50 meV Fig. 13, 25 meV Fig. 14). Proposé : η_t = 20 meV gelé dans t (cadre R8, config v2) et
   η_G = 15 / 25 / 50 meV dans G_k, la différence étant rapportée ; confirmer.
3. Courbe parfaite à η_G 15 meV : 300² comme l'article (première couronne à 53,5 meV ; 26,8 à 600²) + 1 200² en complément ?
4. Ṽ pour V₀ = 2Ṽ/A_cell (étape 6, article) : élément intrabande près de K comme l'article (R9 D.1, disque K, tel quel : 76,96 eV Å² π / 78,16 π* → V₀ = 29,23 /
   29,69 eV), ½ Re Tr à (K, K) (75,87 → 28,82 eV) ou norme de ligne (107,22 → 40,72 eV) ; en variante alignée, sur M2 exact ou sur V_loc − C_N·P_boîte.
   (Kaasbjerg : ∼70 eV Å² → V₀ ≈ +27 eV.)
5. Constantes de Ḡ₀ (éq. 48, étape 6) : nos bandes (ħv_F = 5,459 eV Å, A_cell = 5,266 Å², proposé) ou celles de l'article (v_F = 10⁶ m/s soit ħv_F = 6,582 eV Å,
   A_cell = 5,24–5,25 Å²).
6. Outil de 7a : `mutool trace` (proposé, aucune installation) ou PyMuPDF (installation dans le `.venv`).
7. `.gitignore` : exception pour `results/M2/R8/`.
8. `utils/lattice.build_k_path` (`nk=100` en dur) : corriger ou laisser.
9. Statut TEST au GO.

**STOP — mise à jour de la phase 0 terminée le 2026-09-28.** Rien n'est calculé ; rien n'est écrit dans `src/`, `tests/`, `scripts/` ni `results/` ; attente de
l'accord sur les signatures, de l'audit de l'image minimale, puis du GO.

## Réponses de Greg à M.7 (2026-09-28) — autorisés : 7a et écriture des fonctions ; pas de GO de calcul

Décisions (J1 et suivants attendent l'audit de l'image minimale puis le GO) :
1. Signatures Q1–Q5 acceptées, avec Q3 = inversion générale par lots seulement (pas de forme fermée 2 × 2) et Q6 = correction sur place de
   `utils/lattice.build_k_path` (même nom, points ∝ longueur cartésienne, coins une fois, argument `nk` respecté) au lieu d'un `kpath` dans
   `disorder_average.py`. **Code écrit tout** (Q1–Q5, `build_k_path`, `tests/test_r8_functions.py`) ; Greg relit et commit.
2. η : principal comme proposé (η_t = 20 meV gelé, η_G selon l'étape) + variante **`eta_unique`** (η_t = η_G : 50 meV pour la DOS, 25 meV pour A_k), 9×9 tel quel
   seulement, N_k^int 600 (deux g₀ de plus) ; elle entre en 2, 3 et 7b à côté de la principale (différence de protocole, pas d'ingrédient).
3. DOS parfaite : 300² (protocole de l'article) + 1 200² en complément.
4. V₀ = 2Ṽ/A_cell, Ṽ = élément intrabande près de K de R9 D.1 (π : 29,23 eV ; π* : 29,69 eV ; les deux rapportés) ; ½ Re Tr (28,82 eV) en second ;
   norme de ligne abandonnée ; variante alignée : valeur de R9 D si elle existe, sinon approximation (i), étiquetée.
5. Ḡ₀ (éq. 48) : nos constantes en principal (ħv_F = 5,459 eV Å, A_cell = 5,266 Å²), celles de l'article en sensibilité ; ajouter à l'étape 6 le pôle sans Λ,
   1/V₀ = Re Ḡ₀^W(ε) avec Ḡ₀^W = g₀ sur site p_z(A) de J1, et le Λ de l'éq. 48 qui reproduit Re Ḡ₀^W au pôle.
6. `mutool trace`.
7. Pas d'exception au `.gitignore` : statut TEST, **rien dans `results/`** ; sorties dans `R8_kaasbjerg/out/`, copie dans `article/R8_kaasbjerg/` (rapport, tables,
   figures, npz < 5 Mo). (Remplace « `results/M2/R8/` » partout plus haut.)
8. `build_k_path` corrigé sur place (voir 1).
9. Statut TEST.
- Rétablis à l'étape 1 (information, deux variantes) : ‖T̄^{(W)}_{AB}(K)‖, ‖T̄(K) − T̄(C₃K)‖, ‖T̄(K) − T̄(K′)‖.

## 7a — Kaasbjerg, Fig. 13 (haut, V_A) et Fig. 14 (colonne V_A) : extraction et lecture (2026-09-28 ; aucune de nos DOS n'est calculée)

`r8_driver.py extract` sur le nœud de connexion (quelques secondes) ; PDF `article/R8_kaasbjerg/ref/kaasbjerg_2020_prb101_045433.pdf` (md5
25e169d573441ffac17bd380cd77379c), `mutool trace` (MuPDF 1.22.0) des pages 14 et 15 ; sorties `out/7a/` (`7a_results.json`, `extract_log.txt`, courbes csv,
points blancs csv, profils de la colonne K npz) ; figure de contrôle `fig/7a_controle.{pdf,png}` (points extraits sur le rendu `pdftoppm` à 400 dpi :
superposition exacte des trois courbes et des 56 points à l'œil).

**Méthode.** Chemins en coordonnées de page (pt, y vers le bas). Cadre = fond blanc des axes ; graduations = traits d'encre de 0,216 pt posés sur les bords du
cadre ; valeurs des graduations lues sur le rendu et **vérifiées sur la forme des glyphes** du fichier vectoriel (chiffres dessinés en chemins, classés par
signature : 0, 1, 2, 5, point, moins) : les 20 étiquettes des trois panneaux concordent ; ordre de la légende de la Fig. 13 (tirets = pristine, bleu = « 0.1 », vert
= « 1.0 ») vérifié de même. Étalonnage affine par moindres carrés sur les cinq graduations de chaque axe :

| panneau | axe | échelle | résidu max des graduations |
|---|---|---|---|
| Fig. 13 haut | énergie | 0,016599 eV/pt | 6e-15 eV |
| Fig. 13 haut | DOS | 1,82584e-3 eV⁻¹/pt | 7,3e-7 eV⁻¹ |
| Fig. 14, V_A, 0,1 % | énergie | 0,025295 eV/pt | 2,0e-5 eV |
| Fig. 14, V_A, 1 % | énergie | 0,025295 eV/pt | 1,0e-5 eV |

**Fig. 13 (haut, V_A).** Trois polylignes (pristine 141 sommets, 0,1 % 148, 1 % 154) sur [−1,2 ; +1,2] eV ; sommets sur une grille de **10 meV** (plus petit pas
0,0100 eV ; pas jusqu'à 0,13 eV là où la courbe est droite : points alignés omis par le logiciel de tracé) ; épaisseur des traits 0,697 pt = 1,27e-3 eV⁻¹ (11,6 meV).
Interpolation linéaire au pas de 1 meV (`fig13_VA_interp_1meV.csv`). Δρ = ρ(c_i) − ρ_pristine **tel que tracé** (légende : pristine à η 15 meV sur 300², c_i à η 50 meV).

| grandeur (eV, eV⁻¹) | c_i = 0,1 % | c_i = 1 % |
|---|---|---|
| maximum de Δρ sur [−1, 0] : position | **−0,200** | **−0,210** |
| hauteur | 0,00126 | 0,0307 |
| largeur à mi-hauteur [bornes] | 0,163 [−0,293 ; −0,130] | 0,306 [−0,373 ; −0,067] |
| maximum local de ρ sur [−1, 0] | aucun | −0,27 (0,0602) |
| ρ(0) ; minimum de ρ (position) | 0,00546 ; 0,00546 (0,00) | 0,0128 ; 0,00987 (+0,05) |
| rms de Δρ sur [−1, +1] | 0,00184 | 0,0107 |

Pristine : minimum 0,00684 à 0,00 ; ρ(−1) = 0,1178, ρ(+1) = 0,1366 ; ajustement ρ = α|ε| + β : α = 0,1074 eV⁻² sur [−0,5 ; −0,1], 0,1147 sur [0,1 ; 0,5].
Pour l'alignement des unités avant 7b (arithmétique) : A_cell/(π(ħv_F)²) = 0,05625 eV⁻² par maille et par spin avec nos constantes (0,1125 en comptant le spin) ;
0,0385 (0,0770) avec v_F = 10⁶ m/s et A_cell = 5,24 Å² ; rapports des pentes mesurées à 0,05625 : 1,91 et 2,04. L'étiquette « DOS (eV⁻¹) » ne dit pas si le spin est
compté (décision avant 7b). Les positions sont sur la grille de 10 meV de l'article ; la hauteur à 0,1 % vaut l'épaisseur d'un trait.

**Fig. 14 (colonne V_A).**
- Axe k : une seule graduation (K, x = 114,018 pt), aucune valeur de k ; carte de couleur en polygones : 47,0 largeurs de colonne de 1,9004 pt (48 points k, K au
  24ᵉ depuis la gauche) × 240 lignes de 10,00 meV sur [−1,2 ; +1,2] eV. L'étendue en k (« Γ ← K → M ») n'est pas lisible sans hypothèse ; le chemin Γ–K–M
  complet n'est pas tracé. Tirets rouges (dispersion parfaite) : 2 branches, 67 points.
- **c_i = 1 %** : 56 points blancs (maxima de A_k), deux par colonne de −11 à +14 (un seul aux colonnes −12, 15, 16, 17), énergies sur la grille de 10 meV.
  **Colonne K : +0,010 et +0,110 eV → écart 0,100 eV** (texte : « ∼100 meV »). Branche inférieure : **saut de −0,45 à −0,16 eV** entre les colonnes −4 et −3
  (côté Γ) et **de −0,21 à −0,50 eV** entre les colonnes 4 et 5 (côté M) : aucun maximum dans ]−0,45 ; −0,16[ côté Γ ni dans ]−0,50 ; −0,21[ côté M ; branche
  supérieure sans saut > 0,15 eV. Pas des maxima de la branche inférieure entre colonnes voisines, hors sauts : médiane 0,100 eV côté Γ, 0,070 eV côté M. (Le
  texte de l'article parle d'un dédoublement de la bande de conduction ; les sauts relevés ici pour V_A sont sur la branche inférieure.)
- **c_i = 0,1 %** : aucun point blanc ; la colonne K est au gris le plus foncé (0,128, saturé) sur [−0,035 ; +0,065] eV (11 lignes) : **gap à 0,1 % non lisible**.
  À 1 % la saturation couvre [−0,035 ; +0,165] eV et contient les deux points blancs.

Rien d'autre n'est lu sur les figures. Positions de référence pour 7b : Δρ max −0,200 / −0,210 eV, hauteurs 0,00126 / 0,0307 eV⁻¹, largeurs 0,163 / 0,306 eV ;
gap à K (1 %) 0,100 eV ; sauts de la branche inférieure (1 %) ]−0,45 ; −0,16[ (Γ) et ]−0,50 ; −0,21[ (M).

## Fonctions écrites (2026-09-28, décision 1 : Code écrit tout ; non commité, Greg relit et commit)

| # | fichier | fonction | contenu |
|---|---|---|---|
| Q1 | `src/…/defects/many_body/disorder_average.py` (nouveau, 230 lignes) | `tbar_reduce(t, R_local, nw, wfs=None) -> (Du, tau)` | τ_{ab}(D; ε) = Σ_{R_L − R_L′ = D} t_{(L,wfs[a]),(L′,wfs[b])}(ε) ; table des différences de `local_tmatrix._diff_table` |
| Q2 | idem | `tbar_k(tau, Du, k, U=None, k_chunk=65536)` | T̄^{(W)}_k = Σ_D e^{−2πik·D} τ(D) ; U donné → U†T̄U (hors diagonale compris) |
| Q3 | idem | `green_k(Hk, Sigma, e, eta)` | G_k = [(e + iη)𝟙 − H_k − Σ_k]⁻¹, `np.linalg.inv` par lots (pas de forme fermée) |
| Q4 | idem | `dos_average(Hk, tau, Du, k, c_cell, egrid, eta, linear=False, e_chunk=16, k_chunk=32768)` | ρ = −(1/(πN_k)) Σ_k Im Tr G_k, Σ = c_cell T̄ (c_cell scalaire ou tableau : T̄ partagé entre concentrations) ; ρ₀ ; `linear` : −(1/(πN_k)) Σ_k Im Tr[g₀_k T̄_k g₀_k] ; boucles par paquets de k et d'énergies (mémoire) — `k_chunk` ajouté à la signature acceptée |
| Q5 | idem | `spectral_path(Hk, Tbar, c_cell, egrid, eta)`, `spectral_maxima(A, egrid, prominence=0.0, refine=True)` | A_k = −2 Im Tr G_k (éq. 28) ; maxima en ε à k fixé (`find_peaks`), sommet de la parabole par trois points |
| Q6 | `src/…/utils/lattice.py` | `build_k_path(high_sym_points, nk, B) -> (k, labels, idx, s)` | corrigé sur place : `nk` = nombre total de points (respecté), intervalles ∝ longueur cartésienne (plus grands restes, ≥ 1 par segment), coins une fois et exacts, s = abscisse cumulée ; **signature changée** (argument `B` ajouté, `s` renvoyé) — aucun appelant dans `src/`, `scripts/`, `tests/`, `notebooks/`, `memoire/`, `article/` |

Aucune autre ligne de `src/` modifiée ; `local_tmatrix`, `pole_criterion` (g₀, t, `cluster_ldos`) importés sans modification. Tests `tests/test_r8_functions.py`
(176 lignes, format de `tests/test_r9_functions.py`) : référence = résolvante directe d'un réseau p_z périodique de 6 × 6 mailles (`tb_models.graphene_pz_tb`,
bloc π), potentiel hermitien aléatoire sur 5 mailles, T = V[1 − G₀V]⁻¹ par inversion des matrices réelles 72 × 72 : (1) t local = bloc de T ; Q1 + Q2 = N⟨k|T|k⟩ en
13 k (dont K et hors grille) ; énergie seule ; base de bandes ; (2) `wfs=(3, 4)` sur la disposition à 5 WF = bloc π seul ; (3) Q3, Q4 (trois concentrations,
paquets de 2 énergies et 7 k) = inversion directe ; c = 0 ⇒ ρ = ρ₀ ; (4) terme `linear` = (1/π) Im Tr[t ∂_ε g₀] ; (5) Q5 = inversion directe sur un chemin Γ–K–M ;
sommet exact de pics paraboliques ; (6) `build_k_path` : Γ–K–M à 601 points → coins aux indices 0, 400, 600 (|ΓK| : |KM| = 2 : 1), |ΓK| = 4π/(3a), pas constant
par segment, `nk` = 7 → 0, 4, 6, erreur si `nk` < nombre de coins. Tolérance 1e-12 relative (1e-14 pour les identités exactes). `PYTHONPATH=src pytest
tests/test_r8_functions.py` : **6 passés** (1,5 s). md5 : `disorder_average.py` aa5c3fe2…, `test_r8_functions.py` f4ca6071…, `lattice.py` 0eb72093… ; pilote
`r8_driver.py` 08862e45… (sous-commande `extract` seule).

**STOP — 7a et fonctions faites le 2026-09-28.** Aucune de nos DOS n'est calculée ; rien dans `results/` ; J1 et suivants attendent l'audit de l'image minimale
puis le GO. Pour 7b, décision attendue : l'axe « DOS (eV⁻¹) » de la Fig. 13 compte-t-il le spin (pentes mesurées 0,107 / 0,115 eV⁻² contre 0,0563 par maille et
par spin avec nos constantes) ?
