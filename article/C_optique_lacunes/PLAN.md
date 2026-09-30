# C — Conductivité optique du graphène avec lacunes (plan, 2026-09-29)

Statut : **plan**, aucun calcul. Perspective C d'`campagnes/EM/EM.md` §11 : **article, pas le mémoire** (Greg, 2026-09-29). Préalables, dans l'ordre :
le ménage du dépôt et la réappropriation de la chaîne matrice T par Greg, puis le mémoire, puis C0. Développement et tests sur le
portable, production sur rorqual. Bibliographie : `article/C_optique_lacunes/biblio/Optique du graphène avec lacunes.md` (notes de lecture dans `biblio/research_notes/`) (niche : aucun σ(ω) tiré d'une matrice T ab initio ;
aucune matrice T ab initio spin-résolue de la lacune). Mode prévu : Greg code (technicien), Code écrit tests, docs et calculs de production.

## But et formule

σ(ω) du graphène avec une concentration c de lacunes (par maille), bulle de Kubo-Greenwood sans vertex, Drude et interbande ensemble :

  σ_xx/σ₀ = (8 / π N_k A) ∫dε [f(ε) − f(ε + ħω)]/ħω · Σ_k Tr[ ħv^x ImG_k(ε) ħv^x ImG_k(ε + ħω) ]

- G_k(ε) = [ε + iη − H_k − Σ_k(ε)]⁻¹, Σ_k = c T̄_k (Kaasbjerg 2020, non autocohérente, exacte au premier ordre en c) ;
- **ImG = (G − G†)/2i, la partie anti-hermitienne de la matrice**, pas `np.imag` élément par élément (faux hors base propre) ;
- normalisation (dérivée, à vérifier numériquement en C1) : à c = 0 et η → 0, ImG → −π δ(ε − ε_n) et la formule redonne exactement
  8πS/(ħω N_k A) de `kubo.py` ;
- partie réelle seulement (absorption) ; Im σ demanderait aussi Re G (plus tard, comme B).

## Ce qui existe (inventaire du 2026-09-29)

- **Chaîne matrice T** : `campagnes/R/R8_kaasbjerg/r8_driver.py` enchaîne déjà V_loc (`prep_size` l.582), g₀ (`g0_make` l.726), t et τ(D)
  (`t_pi` l.777, `local_t_cache` + `tbar_reduce`), DOS (`dos_average` l.788), A_k (`tbar_k` + `spectral_path` l.1070) ;
  fonctions dans `defects/many_body/{local_tmatrix,disorder_average,pole_criterion}.py`. Production : 9×9, R_cut = 3 (29 mailles,
  dim 145), η = 0.02 eV, k_int = 300², c = 1 % (`config/production.json`).
- **Même wannierisation que le §2.5** : manifeste `wannier/27x27/wannier_manifest.json` (tb, u, u_dis, sha256 vérifiés) ; mêmes
  conventions e^{+ik·R}/ndegen des deux côtés ; `tbar_k` cohérent (D = R_L − R_L').
- **Données** : locales, seulement des produits (`results/M2/resonance_*.npz`, `campagnes/R/R8_kaasbjerg/out/` : DOS à 0.1 et 1 %, Σ(K)) ;
  sur rorqual : `M_dense_9x9` (3.4 Go), caches V_loc (`R6_production_corrigee/cache/Vloc_M2_9x9.npz`) et g₀. Aucun cache τ(D) local.
- **Banc d'essai liaisons fortes déjà là** : `tb_models.removed_site_vloc` (site p_z retiré sur le vrai H de Wannier) : la « lacune de
  modèle » passe par la même chaîne, ce qui répond au « qu'apporte l'ab initio ? » d'un rapporteur.
- **Vitesse** : `compute_velocity` donne ħv dans la base des bandes (V du H_k propre) ; G_k tourné dans cette base, la trace est
  invariante (G n'y est pas diagonal, sans importance).

## Étapes

- **C0 — Données et conventions** (Code) : rapatrier `Vloc_M2_9x9.npz` (petit) ou le recalculer ; recalculer g₀, t, τ(D) en local
  (R8 : prep 19 s, g₀ 1–8 min selon la grille) ; vérifier que la DOS obtenue redonne `out/dos/dos_9x9.npz`. Grilles : k décalée
  (électron-photon) contre MP non décalée (ch. 4) à harmoniser.
- **C1 — Bulle propre** (Greg) : fonction de bulle par blocs de k sur une grille d'énergie (pas ≲ largeur/4, 2.5 meV comme R8).
  Astuce : X_k(ε) = ħv^x ImG_k(ε) empilé en (n_ε, 25 n_k) ; Σ_k Tr[X(ε)X(ε')] = **un produit matriciel** (n_ε × n_ε), puis
  σ(ω) le long des diagonales ε' = ε + ħω avec les poids de Fermi. Tests : c = 0 contre `sigma_on_grid` (η → 0) ; Σ = −iΓ constant
  (cône : Drude de largeur 2Γ et bord à 2μ élargi, formes analytiques) ; indépendance vis-à-vis des blocs.
- **C2 — Σ de la lacune** (Greg + Code) : brancher Σ_k(ε) = c T̄_k(ε) ; contrôles : DOS de la même G = `dos_average` ; c → 0 ; règle de
  somme (poids spectral total de Re σ conservé à O(c) près).
- **C2b — SCTMA** (matrice T autocohérente, dans la foulée de C2) : à chaque énergie, point fixe Σ(ε) = c T̄[G_loc(ε; Σ)], où t est
  recalculé avec le propagateur local habillé G_loc(R, R'; ε) = (1/N) Σ_k e^{ik(R−R')} [ε + iη − H_k − Σ_k(ε)]⁻¹ au lieu de g₀
  (généraliser `local_green_batch` avec un Σ_k(ε)). Mélange linéaire, énergies indépendantes ; coût = n_iter × g₀ (production rorqual).
  Régularise la divergence en E_D jusqu'à Γ_c ~ √c (Ostrovsky 2006) : c'est la fenêtre IR des signatures des lacunes. Tests : c → 0
  redonne le non autocohérent ; site retiré sur le cône contre l'équation FSBA de Peres 2006. SCTMA n'est pas la CPA (défaut à
  plusieurs mailles) ; double comptage O(c²), du même ordre que les diagrammes croisés négligés. Concurrent : Kumari et al. 2026
  (SCTMA dérivée de la DFT, graphène intercalé Au, ARPES).
- **C3 — Physique** (Code) : c = 0.1 et 1 %, μ = E_D et E_D ± 0.3 eV, 300 K. Sorties : Δσ/σ₀ dans la fenêtre de Pauli par 10¹¹ cm⁻²
  (traduit en I_D/I_G, Cançado 2011), pic à ħω ≈ |μ − ε_rés|, largeur du Drude, bord à 2μ, pic de van Hove ; même chaîne avec
  `removed_site_vloc` (site retiré) : l'écart ab initio / modèle est le résultat de l'article.
- **C4 — Validité** : ne rien revendiquer sous Γ_c ~ 0.1–0.2 eV (Σ non autocohérente, vertex) ; limite ω → 0 contre un Boltzmann
  τ_tr de la même matrice T (mesure directe du vertex manquant) ; en option, le vertex en échelle exact dans l'espace de l'amas
  (noyau séparable : système linéaire de dimension (n_L n_w)², 58² = 3364 avec le bloc π) à ω → 0 et en quelques ω ;
  règle de somme f comme diagnostic (la bulle seule la viole à O(c)) ; convergence en R_cut (Re Σ ne l'est pas à R_cut = 3, NOTES_TGAMMA).
- **C5 — Lacune relaxée et spin (le plus neuf)** : la lacune relaxée est magnétique (liaison pendante σ après Jahn-Teller), donc elle
  se traite avec le spin : la maille vierge n'est pas magnétique, donc ψ_nk est commun aux deux spins et M^NL ne dépend
  pas du spin ; seul Ṽ^σ = V^σ_d − V_p change : `pp.x` avec `spin_component = 1, 2` sur la super-cellule nspin 2 (R1, relaxée), puis
  deux M^L → T^↑, T^↓ → Γ^↑, Γ^↓, séparation de la résonance (à comparer aux dizaines de meV de la STM, Zhang 2016), σ = σ^↑ + σ^↓.
  Colinéaire : pas de spin-flip ni de Kondo. Demande rorqual (pp.x + deux M^L) ; comparer nspin 1 et 2 relaxés à la même référence.

## Décisions (Greg, 2026-09-29)

1. Développement et validation sur le portable, production sur rorqual.
2. Lacune : 9×9 non relaxée, nspin 1 (production du ch. 4) pour C0–C4 ; la relaxée (R1) avec le spin, en C5.
3. η : probablement deux (η_t pour t, η_G pour G, comme R8), à confirmer en C1.
4. Autocohérence : SCTMA dans la foulée (C2b), après le non autocohérent qui en est la première itération.
5. Article, pas le mémoire.

## Pièges connus

- ImG matriciel (voir plus haut) ; grilles d'énergie de ε et ε + ħω communes ; t(ε) n'existe que sur E_D ± 3 eV (ħω ≲ 2.5 eV) ;
- Γ_T(E_D) converge lentement en k_int (+4.8 % de 300² à 600², NOTES_TGAMMA §6d) ; pôle du bloc σ à −0.81 eV (C16) ;
- base des bandes non définie en K (dégénérescence) : exclure K de la grille (grille décalée) ;
- DOS du ch. 4 par spin, σ/σ₀ du §2.5 spin compris ;
- scripts R10 du ch. 4 : `alignment_C` et `matrices_dir` importés mais absents de `config.py` commité (ImportError en local, 2026-09-29).
