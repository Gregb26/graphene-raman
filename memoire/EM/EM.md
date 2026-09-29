# Module électron-photon (série EM) — plan de travail M0–M4

Contexte pour le dépôt `graphene-raman`. Rédigé le 2026-09-28.

- **M0 à M4** : étapes du code que Greg écrit lui-même, avec Code en mode technicien.
- **EM1, EM2, …** : prompts de calcul exécutés par Code (cluster, QE, Wannier90).
- **P28** : modifications du texte du mémoire. Aucune n'est faite avant que tous les résultats EM soient validés.

**Ce plan est un guide, pas un cadre strict.** En implémentant, Greg ajuste les signatures et les structures quand ça aide. En cas d'écart, le code (`src/electron_defect_interaction/electron_photon/` et `tests/`) fait foi, et EM.md est réaligné ensuite. Dernier alignement : 2026-09-29.

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
  - en k : max |r(k) − r(k)†| = 4.9×10⁻³ Å à K, mais **1.5×10⁻² Å** sur 20 000 k tirés au hasard, et ce maximum est dans le bloc **p_z** (σ : 7×10⁻³ Å ; croisé : 10⁻¹⁰ Å). Mesure du 2026-09-28, à consolider par F10.
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

**F5 — `velocity(H_k, dH_k, A_k=None)`** → ε (Nk, nW), V (Nk, nW, nW), ħv (Nk, 3, nW, nW) ; **`compute_velocity(tb, k, mode='berry')`** → H_k, ε, V, ħv (la chaîne complète pour une liste de k)
- **Calcul :** `eigh`, puis ħv = V†(dH + i[H, A])V. `A_k=None` donne la version sans Berry.
- **Mode « centres seuls » :** A_k d'un `centres_only(tb)` (F9, fait en M1). En M0, il est identique au mode complet.
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

**F7 — `kubo_accumulate(eps, hv, hw, mu=0, eta=0.04)` → S (nω, 3, 3) ; `kubo_normalize(S, hw, nk, A_cell)` → σ/σ₀ ; `gaussian_eta(x, eta)`**
- **Calcul :** S est la somme de la formule de Kubo sur un bloc de k, et la normalisation se fait une seule fois à la fin (§3). Les transitions sont choisies par l'énergie et mises à plat en une liste, et la somme sur elles est un seul produit matriciel. `hw` peut être un scalaire.

**Pilote — `sigma_on_grid(tb, N, hw, mu=0, eta=0.04, mode='berry', chunk=10⁵)`**
- Boucle sur les blocs : `compute_velocity` (F4 → F5) → `kubo_accumulate`, puis normalisation avec N_k = N². `mode` vaut `'berry'` ou `'no_berry'`. Il n'y a pas de mode `'centres'` (décidé le 2026-09-28) : « centres seuls » est une propriété du modèle, pas de la chaîne, et on passe `centres_only(tb)` avec `'berry'`.
- Environ 3 s par appel à N = 1800 en M0. N = 300 est déjà à 3×10⁻⁴ près.
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

### Tests pytest (`tests/test_{tb_model,kgrid,velocity_operator,ring,kubo,diagnostics,wannier_io}.py`)

Lancer avec `.venv/bin/python -m pytest tests -v`. Fixtures partagées dans `tests/conftest.py` : `tb`, paramétrée indirectement sur `shift_B` = (0,0,0) et (1,0,0) (« deux jauges »), et `grid` = `make_grid_tb(tb, N = 100)`.

**Données réelles (M1 et après)** : tout ce que les tests savent du 27×27 (formes, ordre des WF, E_D, fenêtre gelée, valeurs figées) est dans `W90_REF` de `tests/conftest.py`, avec l'empreinte sha256 des quatre fichiers lus (fixture `w90_ref`). Si les données de référence changent, les tests sur données réelles (48 aujourd'hui) s'arrêtent tous sur « reference data changed », et les autres passent : on met alors à jour les empreintes et `W90_REF`, et rien d'autre.

**Faits (F1–F7 et pilote, 78 exécutions ; à partir de F6, les tests sont écrits par Code) :**

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

F5 (`compute_velocity` = fourier → hermitize → velocity, écrit par Greg ; utilisé aussi par `sigma_on_grid`) :
- `test_velocity_diagonalizes` (deux jauges) : V unitaire, V†HV = diag(ε), bandes croissantes.
- `test_velocity_hermitian` (deux jauges × avec/sans Berry) : ħv hermitien (4×10⁻¹⁵).
- `test_compute_velocity_modes` : la chaîne, modes `'berry'` et `'no_berry'`, `ValueError` sinon.
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

F7 et pilote :
- `test_gaussian_eta` : intégrale 1, second moment η² (η est l'écart-type), valeur au sommet.
- `test_kubo_accumulate_synthetic` : une transition fabriquée à la main contre la forme fermée g(ħω − gap)·Re(v_α* v_β) ; un k sans transition ; invariance par la phase de v_cv.
- `test_kubo_accumulate_no_transition` : un bloc sans transition donne S = 0.
- `test_sigma_dirac_limit` : σ/σ₀ = 1 à 0.3 eV (1.00145), isotrope, nul hors du plan.
- `test_kubo_reference` (4 cas : avec/sans Berry × deux jauges) : les 16 valeurs du tableau sur la grille 1800², à 10⁻⁴ près (écart mesuré 4.8×10⁻⁵ ; ~13 s en tout).
- `test_kubo_gauge_shift` (N = 300) : avec Berry, σ inchangé à 5×10⁻¹⁵ ; sans Berry, σ_xx change et σ_yy non.
- `test_sigma_on_grid_blocks_and_inputs` : indépendance vis-à-vis de `chunk`, `hw` scalaire, mode inconnu refusé (`ValueError`).

Les tests de F7 détectent trois bugs injectés dans une copie du code : normalisation par la taille du dernier bloc, gaussienne en largeur à mi-hauteur, paires (c, v) inversées.

Au total : 78 exécutions, ~16 s. Pour sauter le tableau à 1800² : `-k "not kubo_reference"`.

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

**Kubo (F7).** La formule de Kubo donne Re σ_αβ = (πe²g_s/(ωN_kA_cell)) Σ Re[v^α* v^β]_cv δ(ħω − ε_c + ε_v), avec v = (ħv)/ħ et g_s = 2. Divisée par σ₀ = e²/4ħ, elle laisse le préfacteur 8π/(ħωN_kA_cell). Un cône de Dirac parfait donne exactement 1 ; la déformation trigonale ajoute un terme en ω² (1.00145 à 0.3 eV, 1.0156 à 1 eV). Les conventions qui reproduisent le tableau :
- g_η est une gaussienne d'**écart-type** η, et non de largeur à mi-hauteur (2.355 η) ;
- N_k est le nombre de points de **toute** la grille, jamais celui d'un bloc ;
- A_cell = |a₁ × a₂| est l'aire **dans le plan** (5.2388 Å² en M0) ;
- la grille décalée de ½ couvre toute la zone de Brillouin, donc les deux vallées.

Pièges NumPy et pytest rencontrés en M0 : `memoire/EM/notes_numpy_pytest.md`.

---

## 5. M1 — Lecteur du `_tb.dat` et diagnostics des données réelles

**But** : remplacer la source des données sans toucher au reste de la chaîne.

**F8 — lecture du `_tb.dat`** → `WannierTB`
- **Répartition (décidée le 2026-09-28)** : la lecture du fichier reste dans `io/wannier_io.py` ; le `WannierTB` est construit dans `electron_photon/tb_model.py`.
- **Lecteur (fait, 2026-09-28)** : `read_w90_tb(path)` → `(HR, R, ndegen, rR, lattice)`, soit H(R) (nR, nW, nW) en eV, R entiers (nR, 3), ndegen (nR,), r(R) (nR, 3, nW, nW) en Å et `lattice` (3, 3) en colonnes.
  - Il remplace `read_w90_HR`, qui ne lisait que le bloc H. Tous les appelants de `src/` et `scripts/` sont passés à `read_w90_tb`.
  - `read_w90_HR` reste un emballage qui retourne `(HR, R, ndegen)`, pour les pilotes de campagne figés (`article/R4_quasi_lie`, `article/R9_controles`), leurs répertoires de travail sur rorqual et les notebooks.
  - Un seul analyseur (`_read_tb_blocks`) lit les sections H et r : même convention d'indices pour les deux. Dans une ligne `j i …`, la valeur va dans `X_R[iR, j−1, i−1]` (ligne = j).
  - Valeurs brutes : aucune division par ndegen, r non hermitisé, H vérifié hermitien (`check_hermicity_HR`).
  - Les 13 `_tb.dat` du dépôt ont leur bloc r ; sur les 13, H, R et ndegen sont exactement identiques (`np.array_equal`) à ceux de l'ancien `read_w90_HR`. Lecture du 27×27 : 0.07 s.
- **Construction (fait, Greg, 2026-09-28)** : `make_wannier_tb(path)` dans `tb_model.py` appelle `read_w90_tb` et remplit `R_cart`, `index` et `minus` comme `make_graphene_tb`. `t` = 2.7 eV et `a_cc` = 1.42 Å gardent les valeurs de M0.
  - Seul `ring` s'en sert, pour q₀ et l'intervalle [0, 2q₀]. Sur les données réelles à 2.33 eV : q/q₀ = 0.95–1.29, car le vrai ħv_F vaut **5.47 eV·Å** (v_F = 8.3×10⁵ m/s), contre 5.75 pour les valeurs par défaut.
  - Tests (`tests/test_tb_model.py`, fixtures `tb_w90` et `eig_w90` de `conftest.py`) :
    - champs dérivés cohérents ;
    - la chaîne cartésienne (`reciprocal` puis `compute_velocity`) reproduit le `.eig` à 1.3×10⁻⁵ eV, comme la transformée en coordonnées réduites ;
    - le `K` de `GridTB` est bien le point de Dirac du réseau réel, à E_D ;
    - A(K) brut non hermitien (2.3×10⁻³ Å en x, 4.9×10⁻³ en y), hermitien après `hermitize` ;
    - `ring` sur les données réelles, μ = E_D : sur la couche à 10⁻⁹ eV, entre les bandes p_z.
  - Tests de mutation : `R_cart` sans `.T` (5 échecs), `minus` = identité (1), t = 10 eV (1, l'anneau), r hermitisé en R (1).
- **Vérifications du lecteur** (`tests/test_wannier_io.py`, 9 tests, sur `wannier/27x27/` suivi par git) :
  - formes, réseau en colonnes (valeurs du §2) ;
  - Σ 1/ndegen = 729, 717 × 1 et 24 × 2, ndegen(−R) = ndegen(R) ;
  - H hermitien à 10⁻¹⁵ eV ; r **non** hermitien, 2.54×10⁻³ Å (x) et 1.65×10⁻³ Å (y), signe que le lecteur ne l'hermitise pas ;
  - diag r(0) = centres du `.wout` à 4×10⁻⁷ Å, partie imaginaire nulle : vérifie aussi l'entrelacement Re(x) Im(x) Re(y)… ;
  - **ordre ligne/colonne** : les valeurs propres n'y sont pas sensibles (H(R)ᵀ donne H(k)*, même spectre ; mesuré 4×10⁻¹⁴ eV). D'où un test géométrique : les trois plus grands |H_{45}(R)| = 2.909 eV sont en R = 0, −a₁, −a₂, à la distance a_cc = 1.42366 Å entre les centres ; avec les indices transposés, on aurait R = 0, +a₁, +a₂ à 3.77 Å ;
  - **ε interpolées aux 729 points grossiers = `.eig`** dans la fenêtre gelée (4 ou 5 bandes sous −1.74 eV) : écart mesuré **1.3×10⁻⁵ eV**, et non ~10⁻⁶. C'est le résidu de `use_ws_distance`, que Wannier90 applique et pas nous. Sans ndegen, l'écart monte à 7.7×10⁻⁵ eV seulement : à 27×27, ndegen = 2 ne touche que 24 R au bord, où H ~ 10⁻⁵ eV. Ce test ne distingue donc ndegen que d'un facteur ~6 ; c'est Σ 1/ndegen qui le teste vraiment ;
  - dégénérescence à K : 6×10⁻⁸ eV, à E_D à 1.3×10⁻⁶ eV près ;
  - `read_w90_HR` retourne exactement les trois premières sorties.
  - Tests de mutation du lecteur : indices permutés (1 échec, l'ordre ligne/colonne), Re/Im mal entrelacés (2), r hermitisé dans le lecteur (1), ndegen ignoré (3), réseau en lignes (2).
  - **Abandonné** : la comparaison avec un `_hr.dat`, dont il n'existe aucun exemplaire dans le dépôt. L'identité avec l'ancien lecteur et le test `.eig` la remplacent.

**F9 — `centres_only(tb)`** → copie avec r_R réduit à la diagonale de R = 0. **Fait (Greg, 2026-09-28)**, dans `tb_model.py`.
- `dataclasses.replace(tb, r_R=…)` : `tb` n'est pas modifié, et les autres champs sont partagés, pas copiés.
- A(k) = diag(τ_n), constant en k et hermitien sans `hermitize` : c'est l'approximation de liaisons fortes (Peierls). Elle jette les autres r_mn(R) : jusqu'à 3.4×10⁻² Å dans le bloc p_z, 0.1 Å dans le bloc σ, et 0.22 Å entre σ et p_z, mais en z seulement.
- On l'utilise sous la forme `compute_velocity(centres_only(tb), k, 'berry')`, sans mode dédié.
- Tests (`tests/test_tb_model.py`) :
  - identité en M0, dans les deux jauges ;
  - seule la diagonale de r(0) survit, l'entrée est intacte, les autres champs sont partagés ;
  - A(k) constant et hermitien ;
  - ε et ħv_nn inchangés ;
  - |ħv_mn|² égal à celui de la **jauge atomique**, calculé indépendamment (dérivée de H_at, sans terme de Berry), à 2×10⁻¹² eV²·Å² ; le mode complet en diffère jusqu'à ~30 eV²·Å².
  - Mutations : entrée modifiée sur place (3 échecs), `np.eye` oublié (2), diagonale placée dans la mauvaise ligne R (5).
- **Premier coup d'œil (2026-09-28)** sur les anneaux résonants, à consolider en M2 et M3 : ⟨|ħv_cv|²⟩ centres seuls / complet = 0.9976 à 0.5 eV, 0.9904 à 1 eV et **0.951 à 2.33 eV**. L'écart croît à peu près comme (ħω)² ; à l'énergie du laser, l'approximation de liaisons fortes coûte ~5 % sur |v_cv|².

**F10 — `hermiticity_report(tb, blocks, k)`** → défaut d'hermiticité de r par bloc, en R et en k. **Fait (Greg, 2026-09-28)**, dans le nouveau module `electron_photon/diagnostics.py` (rapports sur un modèle, sans le modifier) :
- `blocks` : dict nom → (lignes, colonnes), par exemple σ = WF 1–3, p_z = WF 4–5, croisé = (σ, p_z) ; `k` : points cartésiens (grille de `make_grid_tb`) ; retour : dict nom → `max_R` (3,), `frob_R`, `max_k` (3,).
- Défaut en R : D = r(R) − r(−R)† (= 2 × partie anti-hermitienne), avec `minus` et `dagger`. Défaut en k : A(k) − A(k)† sur A **brut** (sans `hermitize`).
- **Normalisation sans les centres** : la partie hermitienne contient les centres τ_n, qui dépendent de l'origine choisie ; un rapport de Frobenius qui les inclut est arbitraire. EM1 donne 1.21×10⁻³ avec les centres, **6.1×10⁻³** sans. On normalise par r − r_centres (`centres_only`).
- **Valeurs cibles (mesures du 2026-09-28, grille 60²)** :

  | Bloc | max_R x, y (Å) | frob_R (sans centres) | max_k x, y (Å) |
  |---|---|---|---|
  | σ | 2.54×10⁻³, 1.65×10⁻³ | 1.2×10⁻² | 6.9×10⁻³, 5.8×10⁻³ |
  | p_z | 1.22×10⁻³, 1.37×10⁻³ | 2.6×10⁻² | 1.28×10⁻², 1.42×10⁻² |
  | croisé | ~10⁻¹¹ | ~3×10⁻¹¹ | ~10⁻¹¹ |

  Sur la grille 100², max_k du bloc p_z monte à (1.34, 1.50)×10⁻² Å : un maximum en k est un échantillonnage, il croît avec la grille ; les valeurs en R sont exactes. Composante z : ~10⁻¹⁴ partout. En R, le bloc p_z est **moins** défectueux que σ en absolu, mais plus en relatif ; en k, c'est lui qui domine (somme cohérente sur R, maximum loin de K, en k_red ≈ (0.51, 0.89)).
- Sur le modèle M0 (centres seuls), S est nul : `frob_R` vaut `nan` (avec un `RuntimeWarning` de NumPy, 0/0 non traité explicitement).
- Tests (`tests/test_diagnostics.py`, 6) : M0 sans défaut et `frob_R` = nan (deux jauges) ; valeurs du tableau (σ, p_z, grille 60², rtol 10⁻³) ; bloc croisé au niveau du bruit ; rapport inchangé quand on déplace l'origine (centres + c). Mutations : S avec les centres (5 échecs), A hermitisé avant la mesure (2), indices appariés au lieu du sous-bloc (6), mauvais axes de réduction (4), transposée sans conjugaison (2).
- **Budget sur ħv_cv (pour M2)** : ce que `hermitize` jette vaut δħv_cv = i(ε_c − ε_v)(V†A_anti V)_cv. Sur l'anneau de 2.33 eV : max |δħv_cv| = **7.2×10⁻³ eV·Å, soit 0.13 % de ħv_F** (2.4×10⁻³ eV·Å à 1 eV). Le budget de §2 (≲ 0.2 %) tient ; l'estimation de 0.6 % notée plus tôt prenait le maximum sur toute la zone et oubliait le facteur ½ de la partie anti-hermitienne.

**F11 — `symmetry_report(tb, even, odd)`** → règles de sélection du miroir σ_h (z → −z), dans `diagnostics.py`. **Fait (Greg, 2026-09-28)**, avec l'aide `extract_block(X, rows, cols)` (sous-bloc des deux derniers axes), que `hermiticity_report` utilise aussi :
- `even` = WF paires sous σ_h (σ : [0, 1, 2]), `odd` = WF impaires (p_z : [3, 4]). H, x et y sont pairs, z est impair : un élément ⟨0m|O|Rn⟩ s'annule quand la parité totale (m, O, n) est impaire.
- Retour : `H_mixed` = max_R |H_{pair,impair}(R)| ; `r_mixed` (3,) = max_R |r^α_{pair,impair}(R)| (x, y interdits, z permis) ; `rz_same` = max_R |r^z| dans les blocs pair-pair et impair-impair (interdit).
- Rien à faire en k : la transformée est linéaire, un bloc nul en R l'est à tout k.
- **Valeurs cibles (27×27, 2026-09-28)** : `H_mixed` = 5.5×10⁻¹⁰ eV (pour des |H| jusqu'à 15 eV) ; `r_mixed` = (8.4×10⁻¹², 7.9×10⁻¹², **0.2204**) Å ; `rz_same` = 7.9×10⁻¹² Å sans les centres (5.9×10⁻¹¹ avec : c'est leur cote z ; le miroir est le plan du graphène, donc on mesure r^z sans les centres, comme en F10).
- Tests (`tests/test_diagnostics.py`, 5) : valeurs du 27×27 ; brisure détectée quand on ajoute à la main 10⁻³ à H_{σ,p_z}, à x_{σ,p_z} ou à z_{p_z,p_z} ; rapport inchangé quand le feuillet est déplacé en z₀ = 7.9 Å. C'est ce dernier test qui impose de retirer les centres : avec eux, `rz_same` vaudrait z₀. Mutations : centres gardés (1 échec), bloc pair-pair pour H (2), mauvais axes (2), composante x au lieu de z (2), blocs impair-impair oubliés (1).
- **Conséquence** : pour la lumière dans le plan, σ et π sont découplés à ~10⁻¹¹ près, dans H comme dans v_x, v_y. Ça justifie `pz_block` (F14) et le choix c, v = π*, π ; le 0.220 Å de EM1 est bien en z seulement.

**Critère de sortie** : toutes les vérifications de F8, avec les normes de F10 et F11 consignées.

---

## 6. M2 — v(k) sur les données réelles

**But** : faire tourner F4 et F5, inchangés, sur le `_tb.dat` de référence, en trois modes (complet, centres seuls, sans Berry).

**F12 — `kpath(B, points, n)`** → chemin Γ–K–M–Γ, dans `kgrid.py`. Sert à inspecter ε et ħv_nn. **Fait (Greg, 2026-09-28).**
- `points` : liste de (étiquette, k réduit) ; n points en tout, répartis au prorata des longueurs (pas uniforme en |k|). Retour : `k` (Nk, 3) cartésien, `x` (Nk,) distance cumulée en Å⁻¹ (abscisse des figures), `ticks` et `labels` des coins.
- Attendu (a = 2.4659 Å) : |ΓK| = 4π/3a = 1.6987 Å⁻¹, |KM| = 2π/3a = 0.8494, |MΓ| = 2π/(√3 a) = 1.4711.
- En K exactement, ħv_nn de π et π* n'est pas défini (bandes dégénérées) ; ε l'est.
- `utils/lattice.build_k_path` existe déjà (coordonnées réduites), mais écrase son argument `nk` par 100 : ne pas le réutiliser tel quel.
- Tests : `tests/test_kgrid.py` (longueurs des segments ; |Δk| = Δx à chaque pas, sans doublon, pas uniforme à 10 % près ; chaque sommet échantillonné une fois, à son tick) et, sur les données réelles, `test_velocity_path_gradient_real` (`tests/test_velocity_operator.py`) : ħv_nn projeté sur la direction du chemin = dε_n/dx à **4×10⁻⁴ eV·Å** (différences centrées, 3000 points), dans les deux modes. On exclut les points voisins d'un sommet et les croisements de bandes, où les bandes triées par énergie ont un coude (écart de 5 eV·Å sinon). Mutations : distance du segment suivant mal initialisée (le bug de la première version, 4 échecs), `endpoint=True` (4), dernier sommet oublié (4), même nombre de points par segment (3), `k_red @ B` sans transposée (4).

**F13 — `fermi_velocity(tb, K, mu, q=1e-3, ntheta=360, mode='berry')`** → ħv_F moyenné sur les directions, par deux voies : |ħ∇ε| (diagonale) et (|ħv^x_cv|² + |ħv^y_cv|²)^{1/2} (interbande). Dans `ring.py` (grandeurs autour de K). **Fait (Greg, 2026-09-28)** ; retour `{'intra_avg', 'inter_avg', 'pi', 'pi_star'}`, `mu` obligatoire (pas de constante des données dans `src`) :
- Cercle de rayon q autour de K ; v et c choisis par l'énergie (μ = E_D).
- Voie intrabande : |ħv_nn| dans le plan pour π et π*, moyennés sur θ, puis la moyenne des deux : l'asymétrie électron-trou est linéaire en q et s'annule dans la moyenne.
- Voie interbande : pour un cône de Dirac, H = ħv_F σ·q et ⟨c|σ|v⟩ = i ẑ×q̂ à une phase près, donc |ħv_cv| = ħv_F.
- **Valeurs cibles (27×27, 2026-09-28)**, q = 10⁻³ Å⁻¹ : π 5.4703, π* 5.4681, moyenne **5.4692 eV·Å** ; interbande 5.4692 (anisotropie ±0.13 %, ±0.27 % sans Berry). À q = 10⁻² : π 5.4803, π* 5.4583, moyenne 5.4693 (convergé) ; à q = 0.1 : distorsion trigonale de ±12 %.
- ħv_F = 5.469 eV·Å ↔ **v_F = 8.31×10⁵ m/s** (ħ = 6.582×10⁻¹⁶ eV·s).
- Piège rencontré : `.real` sur ħv_cv (élément hors diagonale, de phase arbitraire à chaque k) donnait 5.319 au lieu de 5.469 ; il faut |ħv_cv|². μ au-dessus de toutes les bandes lève `IndexError` (c = nW), pas l'assert : bruyant quand même.
- Tests (`tests/test_ring.py`, 6) : M0 (les deux voies = 3ta_cc/2, π = π*, deux jauges) ; valeurs réelles ; asymétrie π/π* ×10 de q = 10⁻³ à 10⁻², moyenne stable à 10⁻⁵ ; voie intrabande identique dans les trois modes, interbande à 10⁻⁴ ; μ qui coupe π* sur le cercle → assert. Mutations : `.real` sur l'interbande (3 échecs), intrabande = π seul (2), v et c inversés (1), assert retiré (1), rayon oublié (4).

- **Vérifications :**
  - ħv hermitien à ~10⁻¹³ **après** l'hermitisation de A (sans hermitisation : défaut d'environ 10⁻² eV·Å, à noter) ;
  - diagonale = différences finies des ε interpolées ;
  - les deux voies de F13 concordent ;
  - ħv_F plausible pour la PBE (v_F ≈ 0.8–0.85×10⁶ m/s, soit ħv_F ≈ 5.3–5.6 eV·Å, à comparer à la littérature).
- **Sortie :** ħv_F (PBE) pour le texte.

**Critère de sortie** : les trois modes tournent, ħv(k) est hermitien et ħv_F est consigné. **Rempli le 2026-09-28** :
- ħv hermitien à **1.4×10⁻¹³ eV·Å** dans les trois modes (500 k au hasard) ; avec A brut (sans `hermitize`), le défaut atteint **7×10⁻² eV·Å** sur toute la zone (et non ~10⁻²) : `test_velocity_hermitian_real` ;
- diagonale = dε/dx le long de Γ–K–M–Γ à 4×10⁻⁴ eV·Å (F12) ;
- les deux voies de F13 concordent à 10⁻⁵ ; **ħv_F (PBE) = 5.469 eV·Å, v_F = 8.31×10⁵ m/s**, dans la fourchette PBE attendue (5.3–5.6 eV·Å).

---

## 7. M3 — Symétries, nœud et anneaux résonants

**But** : les contrôles physiques sur les données réelles, et les chiffres par anneau pour le texte.

**F14 — `pz_block(tb, pz)`** → `WannierTB` réduit aux WF `pz` (liste d'indices en argument, pas de constante des données dans `src`) : `H_R` et `r_R` coupés au bloc pz × pz, centres compris ; liste des R et `ndegen` inchangés. Dans `tb_model.py`, comme `centres_only`. **Fait (Greg, 2026-09-29).**
- `extract_block` a été déplacé de `diagnostics.py` vers `tb_model.py` (module le plus bas ; l'inverse faisait un import circulaire) ; `diagnostics` l'importe de là. Il n'est pas exporté par le paquet : dans un notebook, l'importer depuis `tb_model`.
- **Valeurs (27×27, 2026-09-29)** : sur un cercle q = 0.1 autour de K, ε et |ħv^{x,y}_cv|² du 2×2 = ceux du 5×5 à **5×10⁻¹⁵ eV** et **1.5×10⁻¹³ eV²·Å²** ; `fermi_velocity` donne 5.46919 eV·Å et l'anneau de 2.33 eV ⟨|ħv_cv|²⟩ = **31.443 eV²·Å²**, identiques au 5×5. La lumière dans le plan ne mélange pas σ et π (miroir, F11) ; seul r^z entre σ et p_z est perdu.
- **Piège** : loin de K, π croise des bandes σ, et son indice dans le 5×5 trié change : comparer `eps5[:, 3:]` à `eps2` sur toute la zone donne 4.4 eV d'écart, alors que chaque bande du 2×2 coïncide avec *une* bande du 5×5 à 10⁻¹⁴. Comparer près de K, ou apparier par l'énergie.
- Tests (`tests/test_tb_model.py`, 7) : `extract_block` contre `np.ix_` (lignes ≠ colonnes) ; M0 inchangé par `pz_block(tb, [0, 1])` (deux jauges) ; formes, valeurs = bloc `np.ix_`, copies (pas de mémoire partagée), entrée intacte, champs R partagés ; cercle q = 0.1 (ε à 10⁻¹², |ħv_cv|² à 10⁻¹⁰) ; toute la zone, bandes appariées par l'énergie à 10⁻¹² (et l'appariement par indice échoue, > 1 eV) ; même `fermi_velocity` et même anneau que le 5×5. Mutations : r non coupé (3 échecs), indices appariés `H_R[..., pz, pz]` (6), entrée modifiée sur place (4), r centres seuls (3), colonnes = lignes dans `extract_block` (6, dont les tests de F10 et F11), bloc transposé (1).

**F6 (réutilisé) — anneaux** à 1.96, 2.33 et 2.54 eV (633, 532 et 488 nm) autour de K, avec c et v de part et d'autre de μ = E_D.

**F15 — `ring_stats(tb, K, mu, hws, ntheta=720)`** → dict ħω → variante (`'full'`, `'centres_only'`, `'no_berry'`) → `{'avg', 'node', 'ratio'}`. Dans `ring.py`. **Fait (Greg, 2026-09-29)** :
- `avg` : ⟨|e_x·ħv_cv|²⟩_θ et ⟨|e_y·ħv_cv|²⟩_θ, en unités de (ħv_F)²/2 (ħv_F calculé dedans par `fermi_velocity`) ;
- `node` : δ = angle du nœud − angle de e (e_x, e_y), en degrés, cherché dans la moitié d'anneau tournée vers +e (sinon `argmin` prend l'un ou l'autre des deux nœuds, de même profondeur) ;
- `ratio` : min et max sur θ de |ħv_cv| / |ħv_cv|_complet, normes dans le plan √(P_x + P_y), point par point sur le même anneau (un seul anneau par ħω pour les trois variantes).
- Pièges rencontrés : priorité de `%` et `*` dans le ramenage (`x % 2*np.pi` = `(x % 2)*np.pi`) ; indice d'`argmin` pris dans le tableau masqué ; norme calculée avec P² (P est déjà |·|²) ; `return` dans la boucle.
- Tests (`tests/test_ring.py`, 9 ; valeurs de référence dans `W90_REF.ring_stats`) : M0 (moyennes de `test_ring_average`, nœud de e_x à 0, de e_y à −4°, centres seuls = complet, deux jauges) ; invariance de jauge du mode complet et nombres sans Berry de `test_ring_without_berry` (δ_x = −16.5°, rapport 0.664–1.271) ; tableau ci-dessous (rtol 10⁻⁵) ; symétries (δ_y = 0, C₃, ⟨y⟩ sans Berry = centres à 10⁻¹², nœud de e_x de part et d'autre de q ∥ x) ; moyenne recalculée directement sur l'anneau ; nœud = minimum local trouvé indépendamment, profondeur < 10⁻⁴ ; ħω scalaire. Mutations : `tb` au lieu de `model` (2 échecs), P sans module au carré (8), priorité de `%` (7), `argmin` sur d (7), indice hors masque (6), norme avec P² (2), `return` dans la boucle (5), masque retiré (6), unité sans ½ (4), φ inversés (7).

- **Vérifications** (corrigées par les valeurs ci-dessous) :
  - mode complet : nœud **exactement** à q ∥ e seulement si e est le long d'une ligne miroir passant par K (ici e_y ; e_x dans M0, dont le réseau est tourné de 30°). Pour e_x, le gauchissement trigonal le décale de +6° à +8° (1.5° à 0.5 eV) ;
  - moyennes x et y égales (C₃) en modes complet et centres seuls, à 10⁻⁴ près ; **pas** sans Berry, qui n'est pas covariant sous C₃ ;
  - l'écart centres seuls − complet (c'est l'approximation de liaisons fortes) : −5 % sur ⟨|v_cv|²⟩ à 2.33 eV, presque uniforme en θ ;
  - le mode sans Berry déplace le nœud de e_x (+7° → −9.5°), pas celui de e_y (fixé par le miroir).
- **Valeurs (2026-09-29, 27×27, μ = E_D, 720 angles, un seul anneau par ħω pour les trois modes)**. Moyennes en unités de (ħv_F)²/2, ħv_F = 5.46919 eV·Å ; δ = angle du nœud − angle de e, cherché dans |δ| < 90° ; rapport = |ħv_cv|_mode / |ħv_cv|_complet point par point :

  | ħω (eV) | mode | ⟨x⟩ | ⟨y⟩ | δ_x | δ_y | rapport min–max | ⟨\|v\|²⟩ / complet (déduit de `avg`) |
  |---|---|---|---|---|---|---|---|
  | 1.96 | complet | 1.0375 | 1.0374 | +6.0° | 0 | 1 | 1 |
  | 1.96 | centres | 1.0006 | 1.0005 | +6.0° | 0 | 0.981–0.983 | 0.9644 |
  | 1.96 | sans Berry | 1.1307 | 1.0005 | −8.0° | 0 | 0.662–1.231 | 1.0271 |
  | 2.33 | complet | 1.0512 | 1.0511 | +7.0° | 0 | 1 | 1 |
  | 2.33 | centres | 0.9995 | 0.9994 | +7.0° | 0 | 0.973–0.976 | 0.9508 |
  | 2.33 | sans Berry | 1.1835 | 0.9994 | −9.5° | 0 | 0.578–1.283 | 1.0383 |
  | 2.54 | complet | 1.0594 | 1.0594 | +7.5° | 0 | 1 | 1 |
  | 2.54 | centres | 0.9984 | 0.9983 | +8.0° | 0 | 0.968–0.972 | 0.9424 |
  | 2.54 | sans Berry | 1.2170 | 0.9983 | −10.5° | 0 | 0.528–1.315 | 1.0455 |

  δ est limité par la grille (0.5° ; à 1440 angles : +7.3° et −9.5° à 2.33 eV). ⟨y⟩ sans Berry = ⟨y⟩ centres seuls : τ_B − τ_A est selon x dans ce réseau, donc le terme de Berry des centres n'a pas de composante y.

**F16 — `shift_home_cell(tb, j, L)`**, optionnel, **reporté (2026-09-29)** : l'invariance de jauge est déjà montrée sur M0 (`test_ring_stats_gauge`, `test_ring_without_berry`) ; sur les données réelles, il faudrait étendre la liste des R et un `ndegen` par élément. Test de jauge sur les données réelles. Pour la WF j déplacée de L (|R j′⟩ = |R+L, j⟩) :
- X′_{ij′}(R) = X_{ij}(R + L) pour i ≠ j ;
- X′_{j′i}(R) = X_{ji}(R − L) pour i ≠ j ;
- X′_{j′j′}(R) = X_{jj}(R), avec en plus **r′_{jj}(0) = r_{jj}(0) + L**.
- On diffuse 1/ndegen dans X **avant** le décalage : ndegen devient propre à chaque élément.
- **Attendu :** mode complet inchangé, mode sans Berry modifié.

**F17 — `ring_kpoints_crystal(tb, K, mu, npoints)`** → dict ħω → (k_red, θ) : points de `ring` en coordonnées cristallines, k_red,i = a_i·k/2π (`k @ tb.lattice / 2π`). Dans `ring.py`. **Fait (Greg, 2026-09-29).**
- Liste d'EM2 : `npoints = {2.33: 48, 1.96: 12, 2.54: 12}`, 72 k ; premier point à 2.33 eV (0.73959906, 0.40626573, 0) ; chaque anneau est centré sur K = (2/3, 1/3, 0) à 2×10⁻⁷ près (C₃ des données). Sur l'anneau de 2.33 eV, |ħv_cv| dans le plan va de 4.20 à 7.24 eV·Å : EM2 compare cette norme, sans nœud.
- Pas d'écriture de fichier dans `src` : le fichier `K_POINTS crystal` (poids compris) est écrit par le script de la campagne EM2, qui vérifie aussi que la cellule du nscf est le `unit_cell_cart` de `wannier.win`.
- Piège : `lattice.T` ou `B` à la place de `lattice` donnent K = (0.211, −0.455, 0) ou (4.33, 0, 0), sans erreur.
- Tests (`tests/test_ring.py`, 4) : M0 (retour en cartésien = points de `ring`, anneau centré sur K exactement, deux jauges) ; liste d'EM2 (clés, tailles, premier point, conversion contre un `solve` indépendant avec B, centre) ; chaque point sur la couche ħω à 10⁻¹⁰ eV. Mutations : `lattice.T` (4 échecs), sans 2π (4), `inv(lattice)` (4), `ntheta` oublié (3), même ħω partout (2), μ oublié (2 erreurs), θ inversé (3).

**Critère de sortie** : tableau de F15 (3 anneaux × 3 modes), test 5×5 = 2×2 réussi, liste de k d'EM2 écrite. **Rempli le 2026-09-29** (la liste est produite par F17 ; le fichier QE est écrit avec la campagne EM2).

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
| F8 : lecteur `read_w90_tb` dans `io/wannier_io.py` (fait, à la demande de Greg) ; construction du `WannierTB` | Code ; Greg |
| Découpage et performance, échafaudage pytest, F16, F18, F19, figures (F17 codée par Greg) | Code |
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
| M0 — modèle de liaisons fortes | **fait** : F1–F7 et pilote testés (78 exécutions), tableau de référence reproduit à 4.8×10⁻⁵ | 2026-09-28 |
| M1 — lecteur et diagnostics | **fait** : F8 (`read_w90_tb`, `make_wannier_tb`), F9 (`centres_only`), F10 (`hermiticity_report`), F11 (`symmetry_report`) ; 118 tests en tout | 2026-09-28 |
| M2 — v(k) réel | **fait** : F12 (`kpath`), F13 (`fermi_velocity`) ; ħv_F = 5.469 eV·Å ; 134 tests en tout | 2026-09-28 |
| M3 — symétries et anneaux | **fait** : F14 (`pz_block`), F15 (`ring_stats`), F17 (`ring_kpoints_crystal`) ; F16 reportée ; 160 tests en tout | 2026-09-29 |
| EM2 — DFT directe et postw90 | à préparer (liste de k prête) | |
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