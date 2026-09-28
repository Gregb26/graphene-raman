# Module électron-photon (série EM) — plan de travail M0–M4

Contexte pour le dépôt `graphene-raman`. Rédigé le 2026-09-28.

- **M0 à M4** : étapes du code que Greg écrit lui-même, avec Code en mode technicien.
- **EM1, EM2, …** : prompts de calcul exécutés par Code (cluster, QE, Wannier90).
- **P28** : modifications du texte du mémoire. Aucune n'est faite avant que tous les résultats EM soient validés.

**Ce plan est un guide, pas un cadre strict.** En implémentant, Greg ajuste les signatures et les structures quand ça aide. En cas d'écart, le code (`src/electron_defect_interaction/optics/velocity_operator.py` et `tests/test_velocity_operator.py`) fait foi, et EM.md est réaligné ensuite. Dernier alignement : 2026-09-28.

---

## 0. But

Calculer les éléments de matrice de vitesse interbandes v_mn^(μ)(k) par interpolation de Wannier, **avec la connexion de Berry complète** (Wang et al. 2006, Yates et al. 2007). Le couplage électron-photon du §2.5 du mémoire en découle :

g_mn^(λ)(k) = −e √(ħ / 2ε₀ωΩ_BvK) · e_λ · v_mn(k)

À ω fixé, le préfacteur est constant et absorbé dans la constante globale, comme au ch. 3 : **on calcule et on montre v, pas g**. La conductivité optique σ(ω) sert de validation physique (absorption universelle πα).

### Formule centrale

$$\hbar\,v^{(\mu)}_{mn}(\mathbf k)=\Big[V^{(\mathbf k)\dagger}\Big(\partial_\mu H(\mathbf k)+i\big[H(\mathbf k),A_\mu(\mathbf k)\big]\Big)V^{(\mathbf k)}\Big]_{mn}$$

avec

- ∂_μH(k) = Σ_R iR_μ e^{ik·R} H(R) / ndegen(R) ;
- A_μ(k) = Σ_R e^{ik·R} r_μ(R) / ndegen(R) ;
- r_μ,ij(R) = ⟨0i|r̂_μ|Rj⟩.

Forme équivalente, avec H̄ = V†∂_μH V et Ā = V†A_μV : ħv_mn = H̄_mn **+** i(ε_m − ε_n)Ā_mn (m = ligne ; vérifié numériquement le 2026-09-28 : écart 10⁻¹⁴ avec +, 60 eV·Å avec −). Wang et al. (2006) écrivent la même relation avec les indices dans l'autre ordre, v_nm = H̄_nm − i(ε_m − ε_n)Ā_nm ; recopiée avec mn, elle change de signe. Le terme de Berry est nul sur la diagonale et vaut ε_L en préfacteur sur l'anneau résonant. Le signe a été vérifié numériquement pour la convention e^{+ik·R} de Wannier90.

### Pourquoi on ne peut pas négliger le terme de Berry

C'est ce que fait l'éq. (2.5.7) du mémoire dans sa version actuelle. Estimation en liaisons fortes (t = 2.7 eV), sur l'anneau à 2.33 eV :

- |v_cv| est faux de −34 % à +27 % selon la direction ;
- le nœud de |e_x·v_cv|² tourne de 16° ;
- σ_xx/σ₀ = 1.23, σ_yy/σ₀ = 1.14 et σ_xy/σ₀ = 0.08 au lieu de 1.09, 1.09 et 0 : la symétrie C₃ est brisée ;
- déplacer un atome d'une maille fait passer σ_xx/σ₀ à 2.29.

Le résultat sans ce terme dépend donc du choix arbitraire de la maille d'attache des fonctions de Wannier.

### Potentiel non local

Avec des pseudo-potentiels NC, p̂ ≠ m_e v̂. Le couplage minimal passe par v̂ = (i/ħ)[Ĥ, r̂] (Ismail-Beigi, Chang et Louie 2001). La formule ci-dessus inclut automatiquement [V_NL, r̂], parce que H(R) est le hamiltonien KS complet.

---

## 1. Décisions verrouillées (2026-09-25)

1. **Wannierisation** : celle de la référence du ch. 4, maille unitaire, grille 27×27×1, **même `.chk`**, donc même jauge que le ch. 4.
2. **r(R) complet** : on ne se limite pas aux centres. « Centres seuls » sert de test de sensibilité.
3. **Hermitisation de r dans l'espace k** : A(k) → ½[A(k) + A(k)†], comme postw90 (`get_oper.F90`). C'est équivalent à r_ij(R) → ½[r_ij(R) + r_ji(−R)*], puisque tout R a son −R. La version en k est infaillible et donne le même objet que postw90.
4. **Pas de Wigner-Seitz** (wsvec.dat ignoré), ni pour H ni pour ∂H ni pour r. C'est la convention de la chaîne du ch. 4 pour H. L'effet est au plus de 6.8×10⁻⁵ Å sur r et de 1.4×10⁻⁵ eV sur H.
5. **Figure du §2.5** avec la courbe « sans Berry ». Celle-ci est définie comme l'éq. (2.5.7) telle qu'écrite, en jauge du réseau e^{ik·R}, avec les mailles d'attache choisies par Wannier90. La légende devra le dire.
6. **Texte** : on n'y touche qu'une fois tous les résultats validés et rentrés (P28).

---

## 2. Données (rapport EM1, 2026-09-25)

- **Référence**, en lecture seule : `graphene/qe/defects/unit_cell/27x27/`. Une copie suivie par git se trouve dans `wannier/27x27/wannier_tb.dat`.
- **Répertoire de travail EM1** : `graphene/qe/electron_photon/EM1_tb/`. Le `_tb.dat` régénéré par `restart = plot` est identique à l'octet près, hors ligne de date.
- **Binaire** : wannier90.x 3.1.0 du module `quantumespresso/7.5`. `~/bin/wannier90.x` (ifort) **plante** en `restart = plot`.
- **Wannierisation** :
  - num_wann = 5, num_bands = 20 ;
  - projections `C1:sp2;pz` / `C2:pz` ;
  - WF 1–3 = liaisons σ, centrées aux milieux des liaisons (étalement 0.611 Å²) ;
  - WF 4 = p_z de C1 (1.423660, 0, 0) Å ; WF 5 = p_z de C2 (2.847320, 0, 0) Å (étalement 0.967 Å²).
- **Réseau** (Å) : a₁ = (2.135490, −1.232926, 0), a₂ = (2.135490, 1.232926, 0), a₃ = (0, 0, 15.875316). K = (2/3, 1/3, 0) en coordonnées réduites.
- **Énergies** :
  - **E_D = −4.238895 eV** (dégénérescence à K : 1.8×10⁻⁷ eV) ;
  - dis_froz_max = −1.74 eV = E_D + 2.499 eV ;
  - dis_froz_min = dis_win_min = −25 eV.
- **R** : 741 vecteurs ; ndegen = 717 × 1 et 24 × 2 ; Σ 1/ndegen = 729.
- **Hermiticité** :
  - H hermitien à 10⁻¹⁵ eV ;
  - r **n'est pas hermitien** : max |r_ij(R) − r_ji(−R)*| = 2.5×10⁻³ Å (x) et 1.7×10⁻³ Å (y), maximum dans le bloc σ (R = (−1,0,0), i = 1, j = 3) ;
  - rapport de Frobenius anti-hermitien/hermitien : 1.2×10⁻³ ;
  - en k : max |r(k) − r(k)†| = 4.9×10⁻³ Å à K.
  - C'est l'erreur des différences finies de l'éq. (44) de WYSV06. Budget d'erreur sur ħv_cv à 2.33 eV : environ 0.01 eV·Å, soit ≲ 0.2 % de ħv_F.
- **Décroissance** : r hors diagonale décroît jusqu'à 4×10⁻⁵ Å au bord de la cellule de Wigner-Seitz, soit 2×10⁻⁴ du maximum.

### Format de `wannier_tb.dat` (Wannier90 3.1.0)

1. Ligne 1 : date. Lignes 2–4 : a₁, a₂, a₃ en Å. Ligne 5 : num_wann. Ligne 6 : nrpts.
2. ndegen : 15 entiers par ligne.
3. **Bloc H** : pour chaque R, une ligne vide, la ligne `R₁ R₂ R₃`, puis nW² lignes `j i Re Im`. **j varie le plus vite**, et la valeur est ⟨0j|H|Ri⟩ : **ligne = j, colonne = i**.
4. **Bloc r** : même structure, lignes `j i Re(x) Im(x) Re(y) Im(y) Re(z) Im(z)`, valeur ⟨0j|r|Ri⟩ en Å.
5. **Aucune valeur n'est divisée par ndegen dans le fichier.** La même règle s'applique aux deux blocs : X(k) = Σ_R e^{ik·R} X(R)/ndegen(R).
6. Le bloc R = 0 de r porte les centres sur sa diagonale.
7. Les composantes z de r viennent de b_z = 2π/c : elles sont moins précises, et **inutiles pour la lumière dans le plan**. On ignore v_z.

---

## 3. Conventions du module

- **Unités** : Å, eV, Å⁻¹ ; ħv en eV·Å ; σ en unités de σ₀ = e²/4ħ = 6.0853×10⁻⁵ S.
- **Vecteurs en 3D partout** (z = 0), comme dans le `_tb.dat`.
- **X_mn(R) = ⟨0m|X̂|Rn⟩ ; X(k) = Σ_R e^{+ik·R} X(R)/ndegen(R)**, en jauge du réseau (pas de τ dans les phases).
- **Appariement** : jauge du réseau pour H ⇔ r(R) complet, **centres compris**. En jauge atomique (phases e^{ik·(R+τ_j−τ_i)}), il faudrait d'abord retirer les centres de r(0), sinon ils sont comptés deux fois.
- **Une seule routine de Fourier** pour H, iR_μH et r : `fourier` (F4). Faire ou non de `Hwr_to_Hwk` (chaîne du ch. 4 : coordonnées réduites, `eigh` systématique, une vingtaine d'appelants) un emballage autour de `fourier` : décision en attente.
- **Conteneurs aux frontières, tableaux dans les calculs** : `WannierTB` (produit par le modèle et le lecteur) et `GridTB` (grille k) servent aux appelants ; les fonctions de calcul (`fourier`, `velocity`, Kubo) ne prennent que des tableaux, parce qu'elles servent aussi sur des blocs de la grille, sur k ± h e_μ, sur des anneaux et sur des chemins.
- **Occupations par l'énergie** (ε < μ), jamais par l'indice de bande : μ = 0 en M0, μ = E_D en M1+.
- **Dégénérescence à K** : n'évaluer aucun élément interbande exactement à K. On utilise des grilles décalées, et les anneaux ne passent pas par K.
- **Kubo**, partie absorptive, T = 0, spin inclus :

$$\frac{\sigma_{\alpha\beta}(\omega)}{\sigma_0}=\frac{8\pi}{\hbar\omega\,N_k A_{\text{cell}}}\sum_{\mathbf k}\ \sum_{\varepsilon_v<\mu<\varepsilon_c}\mathrm{Re}\big[(\hbar v^{\alpha}_{cv})^{*}\,\hbar v^{\beta}_{cv}\big]\,g_\eta(\hbar\omega-\varepsilon_c+\varepsilon_v)$$

  - g_η est une gaussienne normalisée de largeur η, en eV⁻¹ ; A_cell = |a₁ × a₂| en Å².
  - La somme porte sur toute la zone de Brillouin (les deux vallées).
  - Absorbance : A = σ/(ε₀c) → πα = 2.293 % quand σ = σ₀.

### Structure de données `WannierTB`

C'est ce que produit le modèle de M0 **et** le lecteur de M1.

| Champ | Forme | Unité |
|---|---|---|
| `lattice` | (3, 3), **colonnes** a₁, a₂, a₃ (`lattice[:, i]` = aᵢ, comme `B[:, i]` = bᵢ ; les listes de vecteurs k, R_cart sont en lignes) | Å |
| `R_int` | (nR, 3) entiers | — |
| `R_cart` | (nR, 3), `R_int @ lattice.T` | Å |
| `ndegen` | (nR,) | — |
| `H_R` | (nR, nW, nW) complexe | eV |
| `r_R` | (nR, 3, nW, nW) complexe | Å |
| `index` | dict, `tuple(R)` → ligne de R dans `R_int` | — |
| `minus` | (nR,) entiers, ligne de −R : `X_R[minus]` = X(−R) dans l'ordre de `R_int` | — |
| `t` | float, saut du modèle jouet | eV |
| `a_cc` | float, longueur de liaison du modèle jouet | Å |

- `R_cart`, `index` et `minus` se déduisent de `R_int` et `lattice` : ce sont des champs stockés, que chaque constructeur doit remplir (le lecteur F8 compris). Alternative envisagée : des `@property`, recalculées à la lecture.
- `t` et `a_cc` n'existent pas dans Wannier90 : le lecteur F8 leur donnera des valeurs par défaut.

### Structure `GridTB` (côté k)

Construite par `make_grid_tb(tb, N)` (F2 + F3), pour les appelants seulement.

| Champ | Forme | Unité |
|---|---|---|
| `B` | (3, 3), colonnes b₁, b₂, b₃ | Å⁻¹ |
| `k_red` | (N², 3), coordonnées réduites | — |
| `k_cart` | (N², 3), `k_red @ B.T` | Å⁻¹ |
| `K` | (3,), (2b₁ + b₂)/3 | Å⁻¹ |

---

## 4. M0 — Modèle de liaisons fortes synthétique

**But** : valider toute la chaîne après le lecteur (Fourier, vitesse, anneau, Kubo) sur un cas dont la réponse est connue indépendamment.

**F1 — `make_graphene_tb(t=2.7, a_cc=1.42, c=15.0, shift_B=(0,0,0))` → `WannierTB`**
- **Contenu :**
  - nW = 2 (A, B) ; a = √3·a_cc ; a₁ = a(1, 0, 0), a₂ = a(½, √3/2, 0), a₃ = (0, 0, c).
  - **Géométrie imposée** : les chiffres « sans Berry » en dépendent.
  - Avec S = {(0,0,0), (−1,0,0), (0,−1,0)} et L = `shift_B` :
    - H_AB(R) = −t pour R ∈ S − L ;
    - H_BA(R) = −t pour R ∈ L − S.
  - r(0) = diag(τ_A, τ_B + L), où L est converti en vecteur cartésien, avec τ_A = 0 et τ_B = (a₁+a₂)/3 ; ndegen = 1.
- **Déplacer B** : on garde fixes les vecteurs de saut δ = R + τ_B − τ_A. τ_B → τ_B + L impose **R → R − L** pour les ⟨0A|H|RB⟩ (et R → R + L pour les ⟨0B|H|RA⟩).
- **Vérifications :**
  - H_ij(R) = H_ji(−R)* ;
  - |R + τ_B − τ_A| = a_cc pour chaque saut, avec ou sans déplacement.

**F2 — `reciprocal(lattice)`** → B (3, 3), colonnes b_j, avec a_i·b_j = 2πδ_ij (vérifié par `test_reciprocal`, plus d'`assert` dans la fonction).

**F3 — `k_grid(B, N, shift=0.5)`** → `(k_red, k_cart, K)`
- k_red (N², 3) réduit et k_cart = k_red @ Bᵀ (N², 3) cartésien, avec k = ((i+s)/N) b₁ + ((j+s)/N) b₂ (ordre `'ij'` : j varie le plus vite) ; K = (2b₁ + b₂)/3.
- **`make_grid_tb(tb, N)`** emballe F2 et F3 dans un `GridTB`.
- **Vérification :** H(K) = 0.

**F4 — `fourier(X_R, R_cart, ndegen, k, deriv=None)`** → X_k (Nk, …)
- **Calcul :** X(k) = Σ_R e^{ik·R} X(R)/ndegen. Quand `deriv` est vrai, le facteur iR_μ est appliqué pour les trois μ à la fois : X_R doit être (nR, nW, nW) et la sortie est (Nk, 3, nW, nW). La fonction accepte des dimensions supplémentaires en fin de tableau (aplaties en (nR, M), un seul produit matriciel).
- **Usages :** H_k (Nk, nW, nW), dH_k (Nk, 3, nW, nW), A_k (Nk, 3, nW, nW).
- **`dagger(X)`** = X† et **`hermitize(A_k)`** = ½(A + A†), sur les deux derniers axes. `hermitize` est appliquée par l'appelant après `fourier` : sans effet en M0, mais testée.
- **Découpage :** traiter la grille par **blocs de k**, à l'extérieur de la fonction. La matrice de phases (Nk, nR) fait environ 38 Go pour 1800² points et 741 R. Prévoir des blocs de 1–2×10⁴ k en M1+ (environ 240 Mo).
- **Vérifications :**
  - H_k hermitien ;
  - H_AB(k) = −t Σ_{R∈S} e^{ik·R} analytique ;
  - dH_k égal à la différence finie centrée de H_k (h = 10⁻⁵ Å⁻¹), à 10⁻⁶ eV·Å près (absolu).

**F5 — `velocity(H_k, dH_k, A_k=None)`** → ε (Nk, nW), V (Nk, nW, nW), ħv (Nk, 3, nW, nW)
- **Calcul :** `eigh`, puis ħv = V†(dH + i[H, A])V. `A_k=None` donne la version sans Berry.
- **Mode « centres seuls » :** A_k d'un `centres_only(tb)` (F9, reporté à M1). En M0, il est identique au mode complet.
- **Vérifications :**
  - ħv hermitien à ~10⁻¹³ ;
  - ħv_nn = ∇ε_n par différences finies, loin de K ;
  - à q = 10⁻³ Å⁻¹ de K : |ħv^x_cv|² + |ħv^y_cv|² → (ħv_F)² et |ħ∇ε| → ħv_F, avec **ħv_F = 3ta_cc/2 = 5.751 eV·Å**.

**F6 — `ring(tb, K, hw, mu=0.0, ntheta=720, tol=1e-12, maxsteps=100)`** → k (nθ, 3), θ (nθ,), q (nθ,), q₀
- **Calcul :** bissection vectorisée le long des nθ rayons issus de K jusqu'à ε_c − ε_v = ħω (fonction `_gap`), avec c et v les bandes qui encadrent μ, choisies par l'énergie. nθ multiple de 6.
- **Réutilisation :** M3, panneau (a) de la figure, et liste de k d'EM2.
- **Vérifications à 2.33 eV, avec Berry :**
  - |ħv^x_cv|² s'annule à θ = 0 et π, c'est-à-dire q = k − K ∥ x ;
  - ⟨|ħv^x_cv|²⟩_θ / ((ħv_F)²/2) = **0.989** (déformation trigonale).
- **Vérifications sans Berry**, B en maille 0 :
  - plus de nœud à 0 ni à π (|ħv^x_cv|² = 2.05 eV²·Å²) ; deux nœuds à 193.6° et 343.4°, qui ne sont plus à 180° l'un de l'autre (le « vers 164° » d'origine est 343.4° − 180°) ;
  - |v_cv| sans/avec ∈ [0.664, 1.271] ;
  - ⟨|v_x|²⟩ sans/avec = 1.125.

**F7 — `kubo_accumulate(ε, ħv, ω, mu, eta)` → S (nω, 3, 3) ; `kubo_normalize(S, ω, N_k, A_cell)` → σ/σ₀**
- **Calcul :** S est la somme de la formule de Kubo sur un bloc de k, et la normalisation se fait à la fin (§3).

**Pilote — `sigma_on_grid(tb, N, ω, mu, eta, mode, chunk)`**
- Boucle sur les blocs : F4 → F5 → F7, puis normalisation.
- C'est **la** fonction de M4.

### Critères de sortie de M0

Grille 1800² décalée, η = 0.04 eV gaussienne, t = 2.7 eV, a_cc = 1.42 Å.

| ħω | Maille de B | Avec Berry (σ_xx, σ_yy, σ_xy)/σ₀ | Sans Berry |
|---|---|---|---|
| 1.0 eV | 0 | 1.0156, 1.0156, 0 | 1.0389, 1.0234, 0.0134 |
| 2.33 eV | 0 | 1.0938, 1.0938, 0 | 1.2266, 1.1381, 0.0767 |
| 1.0 eV | déplacée de a₁ (`shift_B = (1,0,0)`) | inchangé | 1.2248, 1.0234, 0.0403 |
| 2.33 eV | déplacée de a₁ | inchangé | 2.2896, 1.1381, 0.2301 |

- Avec la même grille et le même élargissement : les quatre chiffres. Sinon, un écart d'environ 10⁻³.
- σ_yy sans Berry est inchangé par le déplacement : L = a₁ est purement selon x. C'est un contrôle de plus.
- Ces valeurs ont été obtenues par deux implémentations indépendantes, en jauge du réseau et en jauge atomique.

### Tests pytest (`tests/test_velocity_operator.py`)

Lancer avec `.venv/bin/python -m pytest tests/test_velocity_operator.py -v`. Fixtures : `tb`, paramétrée indirectement sur `shift_B` = (0,0,0) et (1,0,0) (« deux jauges »), et `grid` = `make_grid_tb(tb, N = 100)`.

**Faits (F1–F6, 67 exécutions ; à partir de F6, les tests sont écrits par Code) :**

F1–F4 :
- `test_tb_hermitian` (deux jauges) : liste de R fermée sous R → −R, H_ij(R) = H_ji(−R)*, idem pour r, ndegen(−R) = ndegen(R).
- `test_tb_bond_length` (deux jauges) : trois sauts A → B, |δ| = a_cc.
- `test_tb_shift_same_bonds` : nR = 5 et 7, mêmes δ à l'ordre près.
- `test_reciprocal`.
- `test_fourier_analytic` : H_AB(k) = −t(1 + e^{−ik·a₁} + e^{−ik·a₂}).
- `test_dirac_point` (deux jauges) : H(K) = 0.
- `test_fourier_hermitian` (deux jauges) : H_k, A_k et dH_k hermitiens.
- `test_hermitize` : symétrise, idempotente, laisse A_k inchangé.
- `test_fourier_finite_difference` (deux jauges) : dH_k contre la différence finie centrée.

F5 (`velocity_from_tb` = fourier → hermitize → velocity) :
- `test_velocity_diagonalizes` (deux jauges) : V unitaire, V†HV = diag(ε), bandes croissantes.
- `test_velocity_hermitian` (deux jauges × avec/sans Berry) : ħv hermitien (4×10⁻¹⁵).
- `test_velocity_inputs_unchanged` : `velocity` ne modifie pas ses entrées.
- `test_velocity_diagonal_gradient` (deux jauges × avec/sans Berry) : diagonale = ∇ε par différences finies, gap > 1 eV.
- `test_fermi_velocity` (deux jauges) : ħv_F = 3ta_cc/2 à q = 10⁻³ Å⁻¹, interbande et vitesse de groupe.
- `test_velocity_gauge` : |ħv|² invariant sous shift_B avec Berry ; sans Berry, x change, y non.

F6 :
- `test_ring_on_shell` (deux jauges × 4 énergies) : |k − K| = q et ε_c − ε_v = ħω, recalculé sans `_gap`.
- `test_ring_dirac_limit` (deux jauges) : q = q₀ à 0.1 eV.
- `test_ring_symmetry` (deux jauges × 4 énergies) : C₃ et miroir θ → −θ.
- `test_ring_node` (deux jauges × 4 énergies) : |ħv^x_cv|² = 0 à θ = 0 et π avec Berry.
- `test_ring_isotropy` (deux jauges × 4 énergies) : ⟨|v^x_cv|²⟩ = ⟨|v^y_cv|²⟩, ⟨Re(v^x* v^y)⟩ = 0.
- `test_ring_average` (deux jauges) : rapport 1 à 0.1 eV (limite de Dirac) et 0.98909 à 2.33 eV.
- `test_ring_without_berry` : les chiffres « sans Berry » ci-dessus.

**À faire :**
- `test_kubo_reference` : les huit valeurs du tableau, à 10⁻⁴ près.
- `test_gauge_shift` : avec Berry inchangé ; sans Berry, σ_xx = 2.2896 à 2.33 eV.

### Notes d'implémentation (M0)

Les démonstrations et ordres de grandeur derrière le code et les tolérances des tests. Les docstrings y renvoient.

**Modèle et `shift_B`.** Le saut A(0) → B le long de la liaison δ aboutit à l'orbitale B de la maille R = τ_A + δ − τ_B, qui doit être un vecteur du réseau (vérifié). `shift_B` = L déclare que l'orbitale B de la maille 0 est l'atome situé en τ_B + L. Les liaisons ne changent pas, seules les étiquettes bougent : H_AB vit sur S − L et H_BA sur L − S, avec S = {(0,0,0), (−1,0,0), (0,−1,0)}. C'est le choix que fait Wannier90 quand il place un centre dans une maille voisine.

**Terme de Berry.** ħv̂ = i[Ĥ, r̂] écrit dans la base de Wannier donne ∂H + i[H, A]. ∂H seul traite chaque fonction de Wannier comme un point situé en R ; i[H, A] ajoute où elles se trouvent vraiment (centres, et éléments de position hors diagonale dans les données réelles). Dans la base des bandes, le terme vaut i(ε_m − ε_n)Ā_mn (§0) : il est nul sur la diagonale. Quand A ne contient que les centres τ_j (M0, ou mode « centres seuls »), ∂H + i[H, A] = U†(∂H_at)U avec U = diag(e^{−ik·τ_j}) et H_at le hamiltonien en jauge atomique (phases e^{ik·(R+τ_j−τ_i)}) : le terme de Berry transforme alors simplement la dérivée en celle de la jauge atomique. Vérifié numériquement à 10⁻⁹ près, soit la précision de la différence finie.

**Invariance de jauge (test_velocity_gauge).** Déplacer B de L donne H′(k) = U H(k) U† avec U = diag(1, e^{ik·L}). La dérivée de U ajoute iL_μ[P, H] à ∂H′ (P, projecteur sur B), et le centre déplacé τ_B + L ajoute exactement −iL_μ[P, H] à i[H′, A′]. Donc ∂H′ + i[H′, A′] = U(∂H + i[H, A])U†, V′ = UV, et ħv′ = ħv à la phase arbitraire d'`eigh` près : on compare les |ħv_mn|². Sans Berry, iL_μ[P, H] reste. Avec L = a₁, seule la composante x change.

**Anneau (F6).** Pour un cône parfait, l'anneau est un cercle de rayon q₀ = ħω/(2ħv_F). La déformation trigonale en fait un triangle arrondi : à 2.33 eV, q/q₀ = 1.10 vers les points M (θ = 0°, 120°, 240°) et 0.945 à 60°, 180° et 300°. Le long de chaque rayon, le gap croît de façon monotone de 0 (en K) jusqu'à 2t = 5.4 eV au point M (|KM| = 0.85 Å⁻¹). L'intervalle [0, 2q₀] contient donc exactement un croisement tant que ħω ≪ 2t (vérifié jusqu'à 2.54 eV sur 720 rayons). Au-delà de ~2t, l'anneau rejoint ses voisins : c'est la singularité de van Hove. Il faut environ 39 itérations pour descendre de 2q₀ ≈ 0.4 Å⁻¹ à tol = 10⁻¹² Å⁻¹, ce qui fixe le gap à ~2ħv_F·tol ≈ 10⁻¹¹ eV (mesuré : 5×10⁻¹²). On n'évalue jamais le gap en K même : les deux bandes y sont dégénérées à μ, et le décompte des bandes sous μ dépendrait de l'arrondi.

**Écarts mesurés et tolérances des tests (grille 100², décalée de ½) :**

| Grandeur | Écart mesuré | Tolérance |
|---|---|---|
| dH/dk contre différences finies | 10⁻⁹ à 10⁻⁸ eV·Å | atol 10⁻⁶ |
| ħv hermitien | 4×10⁻¹⁵ eV·Å | atol 10⁻¹² |
| diag ħv contre ∇ε, gap > 1 eV | 1.7×10⁻⁸ eV·Å | atol 10⁻⁷, rtol 0 |
| idem, deux points les plus proches de K et K′ (gap 0.098 eV) | 1.5×10⁻⁶ (erreur ∝ h²v_F/q²) | masqués |
| ħv_F à q = 10⁻³ Å⁻¹, point par point | 7.1×10⁻⁴ relatif (cos 3θ) | rtol 2×10⁻³ |
| ħv_F, moyenne sur 12 θ | 3.8×10⁻⁷ relatif | rtol 10⁻⁵ |
| \|ħv\|² sous shift_B, avec Berry | 3.7×10⁻¹³ eV²·Å² | atol 10⁻¹⁰ |
| anneau : ε_c − ε_v − ħω | 5×10⁻¹² eV | atol 10⁻¹⁰ |
| anneau à 0.1 eV : q/q₀ − 1 | ±3×10⁻³ | rtol 5×10⁻³ |

Pièges NumPy et pytest rencontrés en M0 : `memoire/EM/notes_numpy_pytest.md`.

---

## 5. M1 — Lecteur du `_tb.dat` et diagnostics des données réelles

**But** : remplacer la source des données sans toucher au reste de la chaîne.

**F8 — `read_w90_tb(path)`** → `WannierTB`
- Lecture selon le format du §2.
- **Attention à l'ordre des indices** : dans une ligne `j i …`, la valeur va dans `X_R[iR, j−1, i−1]` (ligne = j).
- On garde les valeurs brutes, avec la division par ndegen au même endroit que pour H dans `Hwr_to_Hwk`.
- **Vérifications :**
  - le bloc H coïncide avec `read_w90_HR` du `_hr.dat` (liste de R, ndegen, |ΔH| ≤ 5×10⁻⁷ eV, soit l'arrondi F12.6) ;
  - Σ 1/ndegen = 729 ;
  - diag r(0) = centres du .wout (4×10⁻⁷ Å) ;
  - **ε interpolées aux 729 points grossiers = `.eig`** pour les bandes de la fenêtre gelée (à ~10⁻⁶ eV près). C'est le test de la transformée avec ndegen ;
  - dégénérescence à K à E_D = −4.238895 eV.

**F9 — `centres_only(tb)`** → copie avec r_R réduit à la diagonale de R = 0.

**F10 — `hermiticity_report(tb, blocks)`** → normes anti-hermitiennes de r par bloc : σ (WF 1–3), p_z (WF 4–5), croisé σ/p_z.
- **À établir :** le défaut dans le bloc p_z seul. C'est le seul qui compte pour π→π*, et il est probablement bien plus petit que 2.5×10⁻³ Å.

**F11 — `symmetry_report(tb)`** → max |H_{σ,pz}(R)| et max |r^{x,y}_{σ,pz}(R)|.
- **Attendu :** du bruit numérique, par la symétrie miroir σ_h. Le 0.220 Å de EM1 entre σ et p_z doit être en z uniquement.

**Critère de sortie** : toutes les vérifications de F8, avec les normes de F10 et F11 consignées.

---

## 6. M2 — v(k) sur les données réelles

**But** : faire tourner F4 et F5, inchangés, sur le `_tb.dat` de référence, en trois modes (complet, centres seuls, sans Berry).

**F12 — `kpath(B, points, n)`** → chemin Γ–K–M–Γ. Sert à inspecter ε et ħv_nn.

**F13 — `fermi_velocity(tb, K, q=1e-3)`** → ħv_F moyenné sur les directions, par deux voies : |ħ∇ε| et (|ħv^x_cv|² + |ħv^y_cv|²)^{1/2}.

- **Vérifications :**
  - ħv hermitien à ~10⁻¹³ **après** l'hermitisation de A (sans hermitisation : défaut d'environ 10⁻² eV·Å, à noter) ;
  - diagonale = différences finies des ε interpolées ;
  - les deux voies de F13 concordent ;
  - ħv_F plausible pour la PBE (v_F ≈ 0.8–0.85×10⁶ m/s, soit ħv_F ≈ 5.3–5.6 eV·Å, à comparer à la littérature).
- **Sortie :** ħv_F (PBE) pour le texte.

**Critère de sortie** : les trois modes tournent, ħv(k) est hermitien et ħv_F est consigné.

---

## 7. M3 — Symétries, nœud et anneaux résonants

**But** : les contrôles physiques sur les données réelles, et les chiffres par anneau pour le texte.

**F14 — `pz_block(tb)`** → `WannierTB` réduit aux WF 4–5.
- **Test :** v_cv (π→π*) dans le plan est identique entre la matrice 5×5 et le bloc 2×2, à la taille près du bloc croisé de H (F11). La lumière dans le plan ne mélange pas σ et π.

**F6 (réutilisé) — anneaux** à 1.96, 2.33 et 2.54 eV (633, 532 et 488 nm) autour de K, avec c et v de part et d'autre de μ = E_D.

**F15 — `ring_stats(tb, rings, modes)`** → pour chaque anneau et chaque mode :
- ⟨|e_x·ħv_cv|²⟩_θ et ⟨|e_y·ħv_cv|²⟩_θ, en unités de (ħv_F)²/2 ;
- l'angle du nœud ;
- min et max de |v_cv| par rapport au mode complet.

- **Vérifications :**
  - mode complet : nœud à q ∥ e ;
  - mode complet : moyennes x et y égales (C₃) ;
  - l'écart centres seuls − complet est petit (valeur à consigner : c'est l'approximation de liaisons fortes) ;
  - le mode sans Berry déplace le nœud.

**F16 — `shift_home_cell(tb, j, L)`**, optionnel : test de jauge sur les données réelles. Pour la WF j déplacée de L (|R j′⟩ = |R+L, j⟩) :
- X′_{ij′}(R) = X_{ij}(R + L) pour i ≠ j ;
- X′_{j′i}(R) = X_{ji}(R − L) pour i ≠ j ;
- X′_{j′j′}(R) = X_{jj}(R), avec en plus **r′_{jj}(0) = r_{jj}(0) + L**.
- On diffuse 1/ndegen dans X **avant** le décalage : ndegen devient propre à chaque élément.
- **Attendu :** mode complet inchangé, mode sans Berry modifié.

**F17 — `ring_kpoints_crystal(rings, B)`** → liste de k en coordonnées réduites pour EM2 : 24 à 48 points sur l'anneau à 2.33 eV, plus quelques-uns à 1.96 et 2.54 eV. On écrit `em2_kpoints_crystal.txt`.

**Critère de sortie** : tableau de F15 (3 anneaux × 3 modes), test 5×5 = 2×2 réussi, liste de k d'EM2 écrite.

---

## 8. M4 — Conductivité σ(ω) et données de figure

**But** : la validation physique finale et les données de la figure du §2.5.

**Pilote (M0) — `sigma_on_grid`** sur les données réelles, μ = E_D, trois modes, ħω de 0.2 à 6 eV.
- **Convergence :**
  - pour résoudre l'anneau, il faut un pas de grille ≲ η/(2ħv_F), soit environ 3.6×10⁻³ Å⁻¹ pour η = 0.04 eV, donc N ≳ 800 ;
  - balayer N = 900, 1200, 1800 et η = 0.02, 0.04, 0.08 eV ;
  - retenir le couple stable à 10⁻³ près sur σ(ε_L).

**F18 — `frozen_window_limit(tb, grid, dis_froz_max, mu)`** → ħω_froz = min_k {ε_c − ε_v : ε_c > dis_froz_max}.
- Au-delà de cette énergie, les π* sortent de la fenêtre gelée : l'interpolation reste lisse mais n'est plus exacte aux points grossiers.
- Trace une ligne verticale sur la figure. Les trois énergies laser tombent largement en dessous, puisque leurs π* sont vers E_D + 1.0–1.3 eV.

**F19 — `map_around_K(tb, K, half_width, n)`** → carte de |e_x·ħv_cv|² et de Δε(k) autour de K, pour le panneau (a) avec les iso-contours Δε = ε_L.

- **Vérifications :**
  - mode complet : σ_xx = σ_yy et σ_xy = 0 ;
  - σ/σ₀ → 1 à basse énergie ;
  - montée vers la singularité de van Hove en M ;
  - le mode sans Berry brise l'isotropie ;
  - concordance avec postw90 `kubo` (EM2).
- **Sorties (npz)** : `em_sigma_{mode}_N{N}_eta{eta}.npz` (ω, σ_αβ, paramètres), `em_map_K.npz` (F19), `em_ring_stats.json` (F15), plus ħv_F et ħω_froz.

**Critère de sortie** : σ(ω) convergé et concordant avec postw90, chiffres consignés.

---

## 9. Répartition Greg / Code

| Élément | Auteur |
|---|---|
| F1, F4, F5, F7 (modèle, transformée, vitesse, Kubo) | **Greg** |
| F6, F13–F15 (anneau, v_F, statistiques) | Greg |
| F8 (lecteur `_tb.dat`), validé par son test contre `read_w90_HR` | Code (« écris-le ») ou Greg |
| Découpage et performance, échafaudage pytest, F16–F19, figures | Code |
| EM2, EM3 (cluster, QE, Wannier90) | Code |

Mode technicien : skill `technicien`, avec les mots-clés « explique » (par défaut), « indice », « montre », « écris-le » et « fin technicien ».

---

## 10. Après M4

### EM2 — références indépendantes (Code)

- **DFT directe** :
  - nscf sur les points de `em2_kpoints_crystal.txt` (même pseudo-potentiel, ecutwfc 100 Ry, nbnd 20, `nosym`), à partir de la densité de `defect_uc_dense_27` (scratch ou miroir `qe_tmp_backup/`) ;
  - puis `bands.x` avec `lp = .true.`, qui donne les éléments ⟨ψ_c|p̂|ψ_v⟩. Code documente les unités et vérifie dans la doc de 7.5 si le commutateur [V_NL, r] y est inclus ;
  - comparaison des |v_cv|² invariants de jauge, k par k. Écart attendu : quelques pour cent au plus, limité par les différences finies de r et le désenchevêtrement.
- **postw90 `kubo`** sur la copie d'EM1 :
  - réglages : `berry = true`, `berry_task = kubo`, `fermi_energy = E_D`, élargissement gaussien fixe de 0.04 eV, `berry_kmesh` égal au N retenu (noms des mots-clés à confirmer dans la doc 3.1.0) ;
  - postw90 donne σ en S/cm pour une cellule 3D : σ_2D = σ_3D · c avec c = 15.875316 Å, puis on divise par σ₀ ;
  - postw90 applique use_ws_distance, un écart négligeable (≤ 10⁻⁵).

### EM3 — figure et chiffres (Code)

- Figure à deux panneaux :
  - (a) carte de |e_x·v_cv|² autour de K avec les anneaux à 1.96, 2.33 et 2.54 eV ;
  - (b) σ(ω)/σ₀ : mode complet, postw90, sans Berry en pointillé, ligne ħω_froz.
- Tableau : ħv_F, σ(ε_L)/σ₀ aux trois énergies laser, ⟨|e·v_cv|²⟩ sur les anneaux, écart entre centres seuls et complet.

### P28 — texte du mémoire, une fois tout validé

- **§2.5.2** : remplacer (2.5.7) par la formule complète et reformuler la phrase « p̂ = m_e v̂ » (Ismail-Beigi 2001).
- **Nouvelle §2.5.3 « Résultats pour le graphène »**, en miroir de la §2.3.3, d'au plus 1.5 p. Elle mentionne :
  - la sélection π↔π* par σ_h ;
  - le budget d'erreur des différences finies (≲ 0.2 %) ;
  - le test centres seuls ;
  - la courbe sans Berry.
- **Passages à rouvrir.** Les numéros de page viennent du PDF du Projet, antérieur aux résumés. Faire un grep de « électron-photon », « hors du cadre », « à un près » et `g^{(\lambda)}` dans tout le mémoire, résumé, abstract et contributions compris.
  - p. 6, intro : « dont l'évaluation numérique est laissée hors du cadre… les deux chapitres suivants calculent tous les autres ingrédients ».
  - p. 37, §2.5.2 : « Nous ne calculons pas v_mn(k) explicitement… ».
  - p. 38, §2.6 : « Les deux derniers termes sont les ingrédients que ce mémoire calcule » → les trois.
  - p. 38, §2.7 : « explicitée mais laissée hors du cadre ».
  - p. 55, §3.2.3 : « son évaluation numérique est laissée hors du cadre ».
  - p. 127, §6.3 : « Elle demande d'abord d'évaluer le couplage électron-photon… » → le passage va dans la Synthèse, et la liste des corrections GW gagne le couplage électron-photon (|g^(λ)|² varie comme v_F²).
  - p. 129, paragraphe de clôture : « calculés, à un près ».
  - **À garder** : « (sauf la paire phonon-photon) » au §2.6, qui est physique.
- **Autour du texte** : symboles.tex (A_μ, r(R), σ₀, α), annexe D (chaîne électron-photon), Zotero : Wang 2006, Ismail-Beigi 2001, Grüneis 2003, Nair 2008, Mak 2008 (Yates 2007 y est déjà).

---

## 11. Statut

| Étape | Statut | Date |
|---|---|---|
| EM1 — inventaire et r(R) | fait (rapport `EM1_rapport.md`) | 2026-09-25 |
| M0 — modèle de liaisons fortes | en cours : F1–F6 faits et testés (67 exécutions), F7 (Kubo) et pilote à faire | 2026-09-28 |
| M1 — lecteur et diagnostics | à faire | |
| M2 — v(k) réel | à faire | |
| M3 — symétries et anneaux | à faire | |
| EM2 — DFT directe et postw90 | après M3 | |
| M4 — σ(ω) et données de figure | à faire | |
| EM3 — figure et chiffres | après M4 | |
| P28 — texte | après validation complète | |

---

## Références

- X. Wang, J. R. Yates, I. Souza, D. Vanderbilt, *Phys. Rev. B* **74**, 195118 (2006) — éq. (44), éléments de r ; connexion de Berry en interpolation de Wannier.
- J. R. Yates, X. Wang, D. Vanderbilt, I. Souza, *Phys. Rev. B* **75**, 195121 (2007) — vitesses et propriétés spectrales par interpolation de Wannier.
- N. Marzari, D. Vanderbilt, *Phys. Rev. B* **56**, 12847 (1997) — centres et différences finies.
- S. Ismail-Beigi, E. K. Chang, S. G. Louie, *Phys. Rev. Lett.* **87**, 087402 (2001) — couplage des potentiels non locaux au champ électromagnétique.
- A. Grüneis et al., *Phys. Rev. B* **67**, 165402 (2003) — absorption optique inhomogène autour de K (nœud).
- R. R. Nair et al., *Science* **320**, 1308 (2008) ; K. F. Mak et al., *Phys. Rev. Lett.* **101**, 196405 (2008) — absorption universelle πα.