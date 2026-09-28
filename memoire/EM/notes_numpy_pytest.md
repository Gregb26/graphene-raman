# Pièges NumPy et pytest rencontrés en M0

Aide-mémoire tiré de l'écriture de F1–F6 (2026-09-25 au 2026-09-28). Il reprend les explications qui se trouvaient dans les docstrings du module de M0 (aujourd'hui le paquet `electron_photon/`) et de ses tests.

## Les fonctions à deux arguments qui n'en prennent qu'un

Les fonctions NumPy élément par élément (`sqrt`, `arctan`, `exp`, `cos`…) ont presque toutes `out` en deuxième position. Leur passer deux tableaux ne lève pas d'erreur, mais le résultat est faux.

- `np.arctan(y, x)` écrit arctan(y) dans x. Il faut `np.arctan2(y, x)`.
- `np.sqrt(a**2, b**2)` renvoie |a|, écrit dans b². Il faut `np.sqrt(a**2 + b**2)` ou `np.hypot(a, b)`.
- `np.linalg.norm(x, y)` : le deuxième argument est `ord`, le type de norme. Pour une norme ligne par ligne, c'est `np.linalg.norm(X, axis=1)`.

## Formes et diffusion

- **Ajouter un axe de taille 1** : `H[:, None]`, de forme (Nk, 1, nW, nW), multiplie les trois composantes μ de A, de forme (Nk, 3, nW, nW). Même chose avec `q[:, None] * e` pour une distance par ligne.
- **`@` contre `*`** : `@` est le produit matriciel sur les deux derniers axes, `*` le produit élément par élément avec diffusion. `eps[:, None, :] @ np.eye(n)` ne fabrique pas de matrice diagonale ; `eps[:, :, None] * np.eye(n)` si.
- **`np.diag` ne connaît pas les piles** : sur un tableau (Nk, nW), il renvoie la diagonale de ce tableau rectangulaire. Pour une pile de matrices : `np.diagonal(X, axis1=-2, axis2=-1)`, dont l'axe diagonal passe en dernier.
- **`np.take_along_axis(eps, idx[:, None], axis=1)`** lit `eps[i, idx[i]]` pour chaque ligne. L'indice doit avoir autant de dimensions que `eps`, et le résultat, de forme (N, 1), se ramène à (N,) avec `[:, 0]`. Oublier ce `[:, 0]` fait diffuser en silence en (N, N) plus loin.
- **Remodelage** : `reshape` suit l'ordre C (dernier axe le plus rapide) à l'aplatissement comme à la reconstruction. Aplatir (Nk, 3, 3) en (3Nk, 3), puis remodeler le résultat en (Nk, 3, …), ne mélange aucun axe.
- **Un seul point k** se passe comme une liste d'un point : `K[None, :]`, de forme (1, 3), ou `k[[0]]`.
- **Tableau aux formes hétérogènes** : `np.array([tab_12, tab_12, 0])` lève « inhomogeneous shape (3,) + … ». Le premier nombre est le nombre d'éléments de la liste, ou de valeurs d'un tuple non déballé. Pour assembler des colonnes : `np.column_stack` avec `np.zeros_like(theta)` à la place du 0.

## Copies, tuples, booléens

- **`W = dH` ne copie pas** : `W += …` modifie aussi `dH`. Il faut `dH.copy()`. Pour tester qu'une fonction ne touche pas ses entrées : faire une copie avant l'appel, puis `np.array_equal`.
- **Un tuple de retour s'indexe** (`f(...)[-1]`) mais se déballe plus sûrement (`_, _, _, hv = f(...)`) : le déballage échoue bruyamment si le nombre de valeurs change.
- **`assert tableau > x`** lève « truth value of an array … is ambiguous ». Il faut `np.all(tableau > x)`, avec la comparaison **à l'intérieur** de `np.all` : `np.all(tab) > x` compare `True` (c'est-à-dire 1) à x.
- **Priorité des opérateurs** : `n % 3 and n % 2 == 0` se lit `(n % 3) and (n % 2 == 0)`. Or `720 % 3` vaut 0, donc faux. Il faut `n % 6 == 0`.
- **`(hi - lo).max`** sans parenthèses est la méthode elle-même, pas son résultat.
- **`a = 1, b = 0`** est une `SyntaxError` ; deux affectations sur une ligne se séparent par `;`.

## Indices et symétries

- **`np.roll(q, n)`** décale un tableau de façon cyclique. Sur 720 angles, une rotation de 120° correspond à `np.roll(q, 240)`.
- **Miroir θ → −θ** : l'indice de −θ_j est (−j) mod n, soit `(-np.arange(n)) % n`, ou `np.roll(q[::-1], 1)`. `q[::1]` ne renverse rien.
- **Minima locaux sur un cercle** : `(x < np.roll(x, 1)) & (x < np.roll(x, -1))`.
- **Paires ±R** : un dictionnaire `index = {tuple(R): iR}` et `minus = [index[tuple(-R)] for R in R_int]`. `X_R[minus]` liste alors X(−R) dans l'ordre de `R_int`.
- **Trier des vecteurs** : `d[np.argsort(clé)]` garde chaque ligne entière ; `np.sort(d, axis=0)` trie chaque colonne séparément et mélange les vecteurs.
- **`np.linalg.solve(lattice, x)`** donne les coordonnées réduites de x, c'est-à-dire les n tels que lattice @ n = x.

## Tolérances

- **`np.allclose(a, b)` teste |a − b| ≤ atol + rtol·|b|**, avec rtol = 1e-5 par défaut. Sur des valeurs de ~6, ça ajoute ~6×10⁻⁵ sans qu'on le voie. On passe `rtol=0` pour une tolérance absolue, et `atol=0` pour une tolérance relative.
- **`np.isclose`** compare deux scalaires et retourne un seul booléen.
- **Moyenne de grandeurs à phase arbitraire** : la phase de chaque vecteur propre d'`eigh` est arbitraire à chaque k. On moyenne |ħv_mn|², ou conj(x)·y pour le même élément, jamais les amplitudes complexes elles-mêmes.
- **Moyenne angulaire** : avec des θ équidistants et `endpoint=False`, la moyenne arithmétique est la moyenne (1/2π)∫dθ. Avec `endpoint=True`, θ = 2π serait un doublon de 0.
- **Graine aléatoire** : `rng = np.random.default_rng(0)`, puis `rng.random(shape)`.

## pytest

- **`tests/conftest.py`** est chargé automatiquement par pytest : les fixtures qui y sont définies servent à tous les fichiers de tests, sans import.
- **Un argument d'une fonction `test_…` est une fixture**, sauf s'il est paramétré. Une constante comme `N` ou `t` se met au niveau du module, pas en argument.
- **Fixture** : une fonction décorée `@pytest.fixture` ; pytest l'appelle et passe son résultat aux tests qui la nomment. Une fixture peut dépendre d'une autre, comme `grid(tb)`.
- **Paramétrisation indirecte** : `@pytest.mark.parametrize("tb", [(0,0,0), (1,0,0)], indirect=True)` envoie chaque valeur à la **fixture** `tb`, dans `request.param`. La fixture lit `getattr(request, "param", défaut)`, pour les tests non paramétrés.
- **Paramétrisation directe** : sans `indirect`, la valeur arrive telle quelle dans l'argument du test.
- **Deux décorateurs empilés** donnent le produit cartésien de leurs valeurs, avec des identifiants comme `[2.33-tb1]`.
- **Plusieurs noms à la fois** : `@pytest.mark.parametrize("hw, ratio", [(0.1, 1.0), (2.33, 0.98909)])` passe les couples ensemble.
- **Ne pas tester une fonction avec elle-même** : pour vérifier `ring`, on recalcule le gap par `velocity_from_tb`, pas par `_gap`.
- **Tester le test** : réintroduire le bug (par exemple retirer `.copy()` ou le terme de Berry) et vérifier que le test échoue.
- **Commandes** : `.venv/bin/python -m pytest tests -v` (ou un seul fichier, `tests/test_ring.py`), avec `-k mot` pour filtrer par nom.
