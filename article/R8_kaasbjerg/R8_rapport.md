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
