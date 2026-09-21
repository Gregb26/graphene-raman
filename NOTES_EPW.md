# NOTES_EPW.md — volet électron-phonon (chapitre 5), état au 2026-09-21 (P2–P7, P9–P11 faits ; P18 chaîne degauss 0.02 Ry faite, §1g)

Calculs dans `graphene/qe/epw/` (voisin du dépôt, hors git). Ce fichier et `NOTES_EPW_REPERES.md` (ex-`graphene/qe/epw/CLAUDE.md`)
vivent dans le dépôt `ab-initio-defects` depuis le 2026-09-21 (liens symboliques aux anciens emplacements) ; scripts et résultats
versionnés dans le dépôt (`scripts/epw_*.py`, `scripts/submit_epw_*.sh`, `results/epw/`). Tout tourne en job SLURM
(compte rrg-cotemich-ac), jamais sur le nœud de connexion.

## 1. État de la chaîne 24k-24q

Grille grossière 24×24 k / 24×24 q, degauss 0.002 Ry, ecutwfc 100 Ry, assume_isolated '2D', QE 7.5 / EPW 6.0.

| étape | répertoire | job(s) | état |
|---|---|---|---|
| phonons DFPT 24×24 q (existants, degauss 0.002) | `24k-24q/phonons/` | antérieurs | terminés ; `_ph0/` 140 Go **en attente d'arbitrage** (voir §5) |
| `save/` (équivalent pp.py, `scripts/epw_pp_save.py`) | `phonons/save/` | 20637120 | terminé, 123 fichiers vérifiés (dvscf, dyn xml, phsave) ; `--check` validé sur 16k-8q (21/21 identiques) |
| q2r + matdyn (graphene.freq, graphene.ifc.xml) | `phonons/` | 20637181 | terminé |
| scf 24×24 + bands.x → `graphene.qe.bands` | `bands/` | 20637311-20637313 | terminé |
| nscf 24×24, nbnd 32, liste de 576 k (crystal) | `epw1/` | 20637121 | terminé |
| epw1 (wannierize, elph, epbwrite, epwwrite) | `epw1/` | 20637122 | **TIMEOUT 24 h** : OOM de la tâche 0 (5.2 Go, 4 Go/cœur) pendant l'écriture d'epmatwp, après 3 h 15 de calcul des g ; wannierisation et les 16 `.epb` (2 Go chacun) complets |
| epw1 relance (`epw1_restart.in` : epbread, wannierize=false via graphene.ukk ; 1 nœud, 16 tâches, 12 Go/cœur) | `epw1/` | 20830241 | terminé en 2 min (MaxRSS 5.66 Go) ; epmatwp 867 Mo, crystal/epwdata/dmedata/vmedata/wigner/selecq.fmt écrits |
| DFPT de référence en \|g\| : ph.x `electron_phonon='prt'`, k = Γ et k = K, q sur Γ–K–M–Γ (31 blocs, chemin en **cartésien 2π/a**) | `dfpt_g_G/`, `dfpt_g_K/` | 20637184…20637473 | terminés (31 blocs q chacun) |
| P1 EPW : band_plot (band.eig, phband.freq) ; prtgkk à k = Γ et k = K sur 474 q | `band_freq_interp/`, `epw_g_G/`, `epw_g_K/` | 20830301, 20830302, 20830303 | terminés (30–40 s chacun) ; liens symboliques vers `../epw1/` (fmt, ukk, epmatwp, graphene.save) |
| P1 post-traitement : `scripts/submit_epw_p1_post.sh 24k-24q 24k24q` (extraction prtgkk → `epw_g_g{ib}{jb}_{nu}.dat`, puis `epw_validate.py`) → `results/epw/validation_24k24q.npz` | dépôt | 20830496 | terminé ; résultats ci-dessous et dans `results/epw/logs/epw_p1_post_20830496.out` |

Dépendances : P1 post ← (20830301, 20830302, 20830303) ← 20830241. P2 ← verdict P1 (seuils : bandes ≤ 5 meV sur
±3 eV, phonons ≤ 1 meV, |g| ≤ 2 % par sous-espaces dégénérés). Rien de P2–P5 n'est lancé.

Résultats P1 (24k-24q, `results/epw/validation_24k24q.npz`, chiffres bruts) :
- Bandes EPW vs bands.x (décalage de référence −2.3849 eV retiré ; E_D bands.x −1.8540 eV, EPW −4.2389 eV) : σ1–σ3 max 0.80 / 3.80 / 4.11 meV sur tout le chemin (hors fenêtre ±3 eV) ; π : 1.43 meV chemin, 0.79 meV max et 0.40 meV rms sur |E−E_D| ≤ 3 eV ; π* : 138.5 meV chemin, 20.06 meV max et 2.86 meV rms sur ±3 eV (la fenêtre gelée s'arrête à E_D + 2.5 eV ; seuil demandé 5 meV).
- Phonons EPW vs matdyn (6 modes) : max 2.59 / 3.54 / 6.96 / 4.30 / 2.43 / 2.45 cm⁻¹ (0.32–0.86 meV), rms 0.58–1.36 cm⁻¹ ; seuil 1 meV.
- Décroissance (R ≥ 0.9 R_max = 34.2 Å) : H 3.5e-6, dynmat 9.0e-4, epmate 1.4e-4, epmatp 3.3e-4.
- |g| DFPT vs EPW (sommes gauge-invariantes par sous-espaces dégénérés, G > 20 meV) : k = K : 45 triplets, max |ΔG| 9.2 meV, rms 2.8 meV, rel. médiane 1.07 %, max 10.3 % (à s = 0.658, |ΔG| 3.1 meV) ; k = Γ : 120 triplets, rel. médiane 0.63 %, hors q = M : max 4.99 %, médiane 0.24 %, |ΔG| ≤ 2.3 meV ; **à q = M (s = 0.636) : 22 triplets aberrants (max 100 %, |ΔG| 222 meV)**.
- **Décision 2026-09-10 (utilisateur) : EPW est une méthode validée, pas de re-validation fine de |g|.** Les vérifications bandes/phonons/décroissance ci-dessus suffisent. Le test |g| garde ses deux limites connues (4 q appariés sur 30, chemin DFPT cartésien vs EPW crystal ; double comptage des groupes de modes à q = M) ; seule une figure de contrôle aux q déjà appariés, q = M exclu, a été produite (`figures/fig_epw_g_control`, optionnelle : 7 sommes à k = Γ, médiane 0,24 %, max 4,99 % ; 15 sommes à k = K, médiane 1,07 %, max 10,3 %).

Chiffres déjà disponibles (24k-24q) :
- Wannier EPW : spreads sp² 3 × 0.6106 Å², pz 2 × 0.9655 Å², total 3.7627 Å² (chapitre 4, 27×27 : 0.6109 / 0.9670, total 3.7668 Å²) ; centres identiques.
- E_F grille grossière EPW = −4.2351 eV ; E_D (chapitre 4, Wannier 27×27) = −4.2389 eV. Le « Fermi energy fine k-mesh = −0.449 eV » d'epw1 est un artefact de nkf = 1, à ignorer.
- Validation 16k-8q (ancienne, `results/epw/validation_16k8q.npz`) : bandes σ ≤ 4.4 meV, π 1.7 meV, π* 349 meV ; phonons ≤ 3.6 cm⁻¹ ; ratios de décroissance H 1.5e-5, dynmat 9.9e-4, epmate 2.7e-4, epmatp 1.3e-3. Son test |g| est invalide (chemin DFPT lu en cartésien alors qu'il avait été donné en cristal) : ne pas s'y référer.

## 1b. P2 elecselfen (2026-09-10)

Répertoires `24k-24q/epw2_selfen_<nkf>_T<T>_dg<degaussw>/` (liens vers ../epw1/{crystal,epwdata,dmedata,vmedata,wigner}.fmt, graphene.ukk, graphene.epmatwp, graphene.save ; **pas** de lien selecq.fmt, EPW l'écrase). Entrée : elecselfen, efermi_read, fermi_energy = −4.2389 (E_D), fsthick 3.5, temps, degaussw, ngaussw 0, mp_mesh_k, restart (pas 1000 q), etf_mem 1, nkf = nqf. 32 tâches (120²) / 64 tâches (240²), 4 Go/cœur.
Post-traitement : `scripts/epw_selfen_post.py --dir <rép> --tag <nkf>_dg<dg>` → `results/epw/selfen_<tag>_T<T>.npz` (module `src/electron_defect_interaction/electron_phonon/selfen.py`) ; `--table … --ref …` pour le tableau de convergence. Γ^ep = 2 Im Σ ; énergies d'epw.out relatives à E_F, remises en absolu une fois ; poids de la zone irréductible reconstruits (orbites du groupe 12 ops + renversement du temps, base réciproque, partagées entre les k listés : EPW en liste 2461 pour 120², 1261 orbites) ; états limités à |ε − E_F| ≤ fsthick ; Γ^ep(ε) = moyenne lorentzienne η = 0.02 eV (chapitre 4) ; médianes sur ε dans ±3 et ±1.2 eV.

| run | job | temps | médiane Γ^ep(ε) ±3 eV | ±1.2 eV | état |
|---|---|---|---|---|---|
| 120², 300 K, dg 0.02 | 20830878 | 6 min 41 | 55.872 meV | 21.583 meV | fait (premier elecselfen réussi → `_ph0` 24k-24q supprimé) |
| 120², 300 K, dg 0.01 | 20831014 | 6 min 46 | 55.689 (−0.33 %) | 21.405 (−0.82 %) | fait |
| 120², 300 K, dg 0.05 | 20831015 | 6 min 25 | 57.486 (+2.89 %) | 22.137 (+2.56 %) | fait |
| 240², 300 K, dg 0.02 (**production 300 K**) | 20831016 | 40 min 37 | 57.168 (+2.32 % vs 120²) | 21.882 (+1.39 %) | fait ; 9721 k irréductibles, 8826 états |
| 240², 10 K, dg 0.02 (**production 10 K**) | 20831356 | 36 min 50 | 55.204 (−3.44 % vs 300 K) | 21.010 (−3.99 %) | fait |

degaussw retenu : 0.02 eV ; grille retenue : 240² (120² → 240² : +2.32 % / +1.39 %, sous le critère 5 % ; pas de 360²). Valeur à E_D : 19.3 meV à 120², 31.5 meV à 240² — artefact de la moyenne lorentzienne dans une région de densité d'états nulle, quantifié en §1f point 4 ; pas un critère.
Γ^ed(c = 1 %) = c × Γ N_cells (`results/M/resonance_9x9.npz`, matrice T, R_cut 3, η 0.02) : médiane 23.436 meV (±3 eV), 10.844 meV (±1.2 eV) ; Γ^ed/Γ^ep(300 K) : médiane 0.36, min 0.088 à +1.87 eV (van Hove π*), max 2.01 à −0.755 eV ; égalité à −1.62 et −0.225 eV.
Figures (P5, `scripts/make_figures_epw.py --prod-tag 240_dg0.02`, PDF+PNG dans `figures/`) : `fig_epw_validation` (bandes + phonons), `fig_epw_gamma` ((a) convergence degaussw/grille, (b) 300 K et 10 K), `fig_epw_vs_ed` (Γ^ep 300 K et Γ^ed c = 1 %, même axe, log). Les trois figures ont été inspectées. À signaler sur fig_epw_gamma/fig_epw_vs_ed : bosse de Γ^ep à E_D (≈ 35 meV, largeur ≈ 0.1 eV) sur les deux grilles, zone à très peu d'états (moyenne lorentzienne mal échantillonnée) ; pics de van Hove à −2.55 et +1.88 eV. P3 (phonselfen) non demandé dans le plan révisé.

## 1c. P6 phonselfen (2026-09-11) — largeurs de phonons γ_qν, **prêt pour ch. 3/5**

Répertoires `24k-24q/epw6_phself_<q>_<nkf>_dg<σ>/` (q = K, G, path, zoom ; mêmes liens qu'en P2, pas de selecq.fmt). Entrée : phonselfen, elecselfen = .false., epwread, efermi_read + fermi_energy = −4.2389 (E_D), fsthick 3.5 (comme P2), temps = 10 300 (nstemp 2, un seul run), degaussw, ngaussw 0, **mp_mesh_k = .false.** (grille k complète : q quelconque n'a pas la symétrie du cristal), restart = .false., nkf1 = nkf2 = nkf, filqf. Fichiers q : `K_only.kpt`, `G_only.kpt`, `qpath_GKMG.kpt` (101 points, Γ/K/M exacts, 42/21/37 segments), `qzoom_GK.kpt` (23 points à |δq| ≤ 0.02–0.04 en fraction de segment autour de Γ et K). 64 tâches, 4 Go/cœur ; 1200² × 101 q : 6 min 33 (job 20912110), zoom 2 min (20912141).
Post-traitement : `scripts/epw_phself_post.py --dir <run> [<run2>] --tag <tag>` → `results/epw/phself_<tag>.npz` (module `electron_phonon/phself.py` : lecture de `linewidth.phself.<T>K` + coordonnées q d'epw.out, fusion et tri le long de Γ–K–M–Γ) ; `--table … --ref …` pour la convergence à K. **Convention (corrigée en P7, §1d) : gamma___ d'EPW est numériquement la LARGEUR TOTALE (FWHM) pour ce système ; HWHM = gamma___/2. npz : `gamma_epw` (brut), `gamma_fwhm` = gamma_epw, `gamma_hwhm` = gamma_epw/2.**

Convergence à K (γ(A1') = gamma___ d'EPW brut, meV, 10 K / 300 K ; réf. 1200², 0.02) :

| nkf | σ = 0.01 | σ = 0.02 | σ = 0.05 |
|---|---|---|---|
| 120² | 0 / 0 | 0 / 0 | 0.002 / 0.002 |
| 240² | 2.903 / 2.499 | 5.799 / 4.993 (+119 %) | 3.444 / 2.969 |
| 480² | 4.524 / 3.697 | 3.646 / 3.024 (+38 %) | 2.646 / 2.192 |
| 720² | 4.212 / 3.450 | 2.851 / 2.351 (+7.8 %) | 2.643 / 2.188 |
| 960² | 3.454 / 2.837 | 2.666 / 2.194 (+0.8 %) | 2.643 / 2.188 |
| 1200² | 2.969 / 2.438 (+12 %) | **2.645 / 2.177** | 2.644 / 2.188 (+0.5 %) |

Le couple (240², 0.02) demandé n'est PAS convergé (+119 %) ; à 120² γ = 0 : l'anneau résonant π→π* (rayon ħω/2v ≈ 0.009 Å⁻¹) est plus petit que le pas de grille (0.012 Å⁻¹ à 240²). Retenu : **1200², σ = 0.02 eV** (σ = 0.05 converge dès 480² à la même limite ; σ = 0.01 pas encore convergé à 1200²). Γ (E2g, HWHM 10 K) : 480² 1.135, 960² 1.324, 1200² 1.335 (σ 0.05 : 1.334) → même conclusion.

Production (1200², 0.02, `results/epw/phself_path_1200_dg0.02.npz`, 123 q = chemin + zoom) :

| mode | ω (meV / cm⁻¹) | T | gamma___ EPW = FWHM (meV / cm⁻¹) | HWHM = gamma___/2 (meV) |
|---|---|---|---|---|
| E2g (Γ, LO = TO exactement) | 162.69 / 1312.2 | 10 K | 1.3346 / 10.76 | 0.6673 |
| | | 300 K | 1.2235 / 9.87 | 0.6118 |
| A1' = mode le plus haut à K (branche 3 à degauss 0.002, branche 6 à 0.02 ; sélection par caractère dans les scripts depuis P18, jamais par index) | 120.53 / 972.1 | 10 K | 2.6448 / 21.33 | 1.3224 |
| | | 300 K | 2.1770 / 17.56 | 1.0885 |

(Tableau corrigé le 2026-09-11 après P7 ; l'ancienne lecture « gamma___ = HWHM » donnait des FWHM de 21.5 / 42.7 cm⁻¹, 2× la littérature.)

Chemin (300 K, valeurs gamma___ brutes = FWHM) : γ nul (< 3e-4 meV) partout sauf |q| ≲ ħω/v autour de Γ et K (anomalies de Kohn) ; autour de Γ, γ(E2g) monte de 1.22 meV à Γ à 2.85 meV au bord |q| = ħω/v ≈ 0.027 Å⁻¹ (s = 0.0068) puis tombe à 0 (canal interbande fermé) ; à q ≠ 0 les deux branches LO/TO diffèrent (3.18 vs 0.54 meV à |q| = 0.027 Å⁻¹, structure sin²θ/cos²θ du vertex E2g) ; à K, γ(A1') est maximal à K et s'annule pour |δq| ≳ 0.01 en s. Acoustiques à Γ : ω ≈ 0, γ = 0 (masqués dans la figure par ω < 5 meV). Aucune dépendance résiduelle au smearing au couple retenu (0.02 vs 0.05 : 0.5 %). γ diminue de 10 K à 300 K (−8 % à Γ, −18 % à K : facteurs de Fermi sur l'anneau).
Figure : `figures/fig_epw_phonselfen` ((a) γ_qν(300 K) sur Γ–K–M–Γ, six branches triées en fréquence, (b) E2g/A1' aux deux T), virgule décimale partout (appliquée aussi aux trois autres figures EPW, `fr()`/`virgule()` dans `make_figures_epw.py`).

**Anomalies à signaler (chiffres bruts) :**
1. Fréquences : ω(E2g, Γ) = 162.7 meV (1312 cm⁻¹) et ω(A1', K) = 120.5 meV (972 cm⁻¹) dans cette chaîne, contre ≈ 196 meV (1580 cm⁻¹) et ≈ 150–160 meV attendus : anomalies de Kohn adiabatiques exagérées par degauss 0.002 Ry (27 meV) sur 24×24 k (K sur la grille, dégénérescence π/π* à E_F). Décision arbitrée (degauss 0.002 conservé) : documenter, ne pas recalculer. Toute utilisation au ch. 3 doit citer ces ω.
2. Largeurs vs littérature : résolu en P7 (§1d) — gamma___ d'EPW est numériquement la FWHM ; FWHM 9.9–10.8 cm⁻¹ à Γ et 17.6–21.3 cm⁻¹ à K, conformes à la littérature (≈ 10–11 et ≈ 20–22 cm⁻¹). Le vertex ⟨D²⟩ est à 1–5 % de Piscanec 2004 (LDA).
3. Aucune bosse numérique ; le pic à |q| = ħω/v près de Γ est physique (fermeture du canal interbande).

## 1d. P7 — ⟨D²_Γ⟩, ⟨D²_K⟩ et localisation du facteur 2 (2026-09-11)

Script `scripts/epw_d2_extract.py` (aucun nouveau run EPW ; sortie `results/epw/d2_extract_24k24q.npz`). Conventions (documentées dans le script) : S_g = Σ_{i,j∈π,π*} Σ_ν |g_ij^ν(k = K, q)|² (somme de bloc, invariante de jauge ; ν ∈ LO+TO à q = Γ, A1' à q = K) ; D² = |g|² · 2Mω/ħ (calculé en SI : 5746.7 Å⁻² par eV de ħω et par eV² de g²) ; ⟨D²⟩ = moyenne par entrée = S_D/8 (Γ) ou S_D/4 (K). Identification avec le ⟨D²⟩_F de Piscanec appuyée par (i) l'identité TB ⟨D²_Γ⟩ = (9/4)(∂t/∂a)² = 45.6 pour ∂t/∂a = 4.5 eV/Å, qui est exactement S_D/8 dans la convention de déplacement d'EPW, (ii) la formule cône de Dirac ci-dessous qui redonne la FWHM ≈ 11 cm⁻¹ de la littérature à partir de 45.6 et v_F = 5.5 eV·Å. **Valeurs de référence citées de mémoire (Piscanec 2004 LDA 45.6 / 92.05 ; Lazzeri 2008 GW ≈ 62.8 / ≈ 193 eV²/Å²) : Zotero inaccessible depuis le cluster, à vérifier dans les PDF.**
v_F de notre chaîne : EPW (band.eig, |ε − E_D| ∈ [0.02, 0.3] eV) π 5.496, π* 5.431 → 5.464 eV·Å ; bands.x 5.31 eV·Å.

Comptabilité des facteurs (route 2, cône de Dirac, E_F = E_D, T → 0) : Im Π = π Σ_k w_k Σ_ij |g|² (f_i − f_j) δ(ε_j,k+q − ε_i,k − ħω) avec Σ_k w_k = 2 (spin, wkf = 2/N_k dans bzgrid.f90 l. 376, vérifié dans la source qe-7.5) ; seul π → π* contribue ; anneau 2v|k| = ħω, ∫d²k δ = πħω/(2v²) ; poids interbande moyen par mode = S^ν/4 (structure sin²θ/cos²θ, **confirmée** par γ_LO ≠ γ_TO à q ≠ 0 dans le zoom) ; vallées : 2 à q = Γ, 1 à q = K (k ∈ K' + K ≡ Γ). D'où γ_HWHM(E2g, par mode) = A_c ħω S_g(LO+TO)/(16 v²) et γ_HWHM(A1') = A_c ħω_K S_g(K)/(16 v²) ; ω s'élimine : γ = A_c ħ² S_D/(32 M v²). Intégrale vérifiée par force brute (grille 2400², gaussienne σ = 0.02 : 0.6685 vs 0.6684 meV analytique). N(E_F) imprimé par EPW est par spin (dosef/2 dans la source) et n'intervient pas.

| ⟨D²⟩ (eV²/Å², moyenne par entrée) | Γ | K | K/Γ |
|---|---|---|---|
| route 1a — g EPW (chaîne 24×24, ω 162.7 / 120.5 meV) | 43.6 (0.96) | 86.3 (0.94) | 1.98 |
| route 1b — g DFPT « prt » (grille 16×16, ω 194.2 / 159.8 meV) | 45.1 (0.99) | 91.6 (1.00) | 2.03 |
| route 2 — inversion de γ EPW (1200², σ 0.02, 10 K) | 87.0 (1.91) | 344.6 (3.74) | 3.96 |
| Piscanec 2004 (LDA, mémoire) | 45.6 | 92.05 | 2.02 |
| Lazzeri 2008 (GW, mémoire) | ≈ 62.8 | ≈ 193 | ≈ 3.1 |

Les deux jeux de g (EPW et DFPT direct, grilles k et ω différents) donnent le même vertex à 3–6 % : le point K sur la grille 24×24 (degauss 0.002 Ry) décale ω de 17–25 % mais pas D². Les routes 1 et 2 ne se recoupent PAS : γ_EPW = 1.996 × Im Π^{Dirac}(g EPW) à Γ et 3.992 × à K (entiers).

**Verdict (cas 2) : ⟨D²⟩ ≈ Piscanec (0.94–1.00), le vertex n'est pas renforcé ; le facteur 2 est un facteur de comptage entre la sortie gamma___ de phonselfen et Im Π (HWHM).** Lecture de la source qe-7.5 (selfen.f90 l. 1158–1191, bzgrid.f90 l. 376, printing.f90) : formule codée = π Σ wkf |g|²(f−f)δ avec Σwkf = 2 et g identique à celui de prtgkk (colonne 8 = epc_unsym) — aucun facteur explicite trouvé ; l'écart 2× est donc établi **empiriquement** (Γ : force brute + structure angulaire vérifiée). Conséquence : gamma___ = FWHM numériquement (10.8 cm⁻¹ à Γ, 21.3 à K à 10 K, = littérature). Le facteur 4 à K se décompose en 2 (commun) × 2 (vertex A1' entièrement interbande sur l'anneau, hypothèse de modèle cohérente avec γ_K/γ_Γ = 2 chez EPW et dans la littérature ; non testée). Correction appliquée : `phself.py` → `gamma_epw` (brut), `gamma_fwhm` = gamma_epw, `gamma_hwhm` = gamma_epw/2 ; npz, tableau §1c et `fig_epw_phonselfen` (axes « largeur totale ») re-générés. **Réserve levée le 2026-09-14 (§1f point 2)** : le prtgkk sur l'anneau confirme gamma___ = 2 Im Π aux deux points (rapport 1.997), et le vertex A1' est purement interbande sur l'anneau (S/2), ce qui explique le facteur 4 de la route 2 à K.

**État : prêt pour ch. 3/5** avec la convention gamma___ = FWHM ; ω(E2g, Γ) = 162.7 et ω(A1', K) = 120.5 meV restent les fréquences de la chaîne (artefact K-sur-grille du degauss 0.002 Ry, §1c), à citer telles quelles ou à remplacer par les DFPT « prt » (194.2 / 159.8 meV) si le chapitre le permet.

## 1e. P9 — fig_epw_decay : décroissance de H, D, g en représentation de Wannier (2026-09-11)

Données : `24k-24q/epw1/decay.{H,dynmat,epmate,epmatp}` (écrits par epw1 24k-24q, R en Å, valeurs en Ry telles qu'étiquetées par EPW ; déjà lus dans `results/epw/validation_24k24q.npz`, clés `decay_<q>_r/_v`) ; aucun recalcul. Figure : bloc 5 de `scripts/make_figures_epw.py` → `figures/fig_epw_decay.{pdf,png}` (2×2 comme fig_locality : (a) H(R_e), (b) D(R_p), (c) g(R_e), (d) g(R_p) ; y en eV = Ry × 13.605693 pour toutes les quantités, y compris D et g dont le « Ry » d'EPW désigne une matrice dynamique et un élément de matrice de Wannier ; axes x depuis 0 jusqu'au rayon WS ; pointillé = apothème WS).
Cellule de Wigner-Seitz de la supercellule 24×24 (a = 2.4659 Å) : apothème (demi-largeur) 29.59 Å, rayon aux sommets 34.17 Å = R_max des fichiers ; 601 vecteurs R, 57 distances distinctes.

| quantité | max (R = 0) | plancher (médiane R ≥ 0.9 R_max) | chute | 1 / 2 / 3 / 4 ordres atteints à |
|---|---|---|---|---|
| H(R_e) | 15.1 eV | 2.2e-5 eV (min 2.4e-6) | 5.8 ordres | 4.3 / 4.3 / 8.5 / 10.7 Å |
| D(R_p) | 17.4 eV | 1.6e-2 eV (min 1.1e-2) | 3.0 ordres | 2.5 / 4.9 / 17.3 / — Å |
| g(R_e) | 6.46 eV | 5.5e-4 eV (min 1.0e-4) | 4.1 ordres | 4.3 / 6.5 / 10.7 / 34.2 Å |
| g(R_p) | 6.46 eV | 2.0e-3 eV (min 1.2e-3) | 3.5 ordres | 4.3 / 6.5 / 11.3 / — Å |

Asymétrie de g : enveloppes identiques jusqu'à ≈ 5 Å (0.20 eV), puis g(R_p) décroît plus lentement et plafonne 3–4× plus haut que g(R_e) (à 10 Å : 1.3e-2 vs 8.8e-3 ; 15 Å : 3.5e-3 vs 1.2e-3 ; 25–30 Å : 2.2e-3 vs 7–8e-4 eV) — la portée côté phonon dépasse la portée côté électron ; le plancher de D (1.6e-2 eV, 3 ordres seulement) est du même type (partie non analytique / anomalies de Kohn de la matrice dynamique). Le dépôt LaTeX n'est pas sur le cluster : le PDF est livré dans `figures/` du dépôt de calcul (même chemin que les autres fig_epw_*), à copier.

## 1f. P11 — points fermés avant rédaction de §5.4 (2026-09-14 ; provenances fichier:ligne re-grep depuis les fichiers réels ; chemins relatifs à `graphene/qe/epw/24k-24q/`)

**1. Fonctionnelle et pseudopotentiel.** `abinit_processing/pseudo/C.upf` : ONCVPSP (D. R. Hamann, « scalar-relativistic version 3.3.0 08/16/2017 », l. 6), `pseudo_type="NC"` (l. 71), `relativistic="scalar"` (l. 72), **`functional="PBE"`** (l. 80) ; md5 34a24e64c0a39f27c6c36b90a16ac686. Chapitre 4 : même fichier (md5 identique dans `scratch/qe_tmp/defect_uc_dense_27/defect_uc_dense_27.save/C.upf`, `<pseudo_dir>` l. 24 et `<functional>PBE</functional>` l. 58 de son `data-file-schema.xml`). Même pseudo et même fonctionnelle (PBE) aux deux chapitres ; Piscanec 2004 est GGA-PBE + Troullier-Martins : la comparaison ⟨D²⟩ se fait à fonctionnelle égale, pseudos différents.

**2. Convention gamma___ — tranchée : gamma___ = 2 Im Π = largeur totale (FWHM).** Méthode : (i) lecture de la source qe-7.5 (P7, §1d) : le commentaire de `selfen.f90` dit « half width » et la formule codée est π Σ_k wkf |g|² (f−f) δ avec Σ wkf = 2 ; (ii) test direct `epw_g_ring/` (job 21063392, 25 s, prtgkk ; `ring.kpt` = 12 k sur deux anneaux autour de K, `GK.kpt` = q ∈ {Γ, K}) analysé par `scripts/epw_ring_check.py` → `results/epw/ring_check_24k24q.npz`. Sur l'anneau, |g₄₅|² interbande est lu directement (aucune hypothèse de structure angulaire) : E2g (q = Γ, modes 5+6) : somme de bloc S = 0.37275 eV² (constante à 1e-4), interbande ⟨|g₄₅|²⟩ = 0.09318 eV² = S/4.00 (fractions par angle 0.241–0.259) ; A1' (q = K, mode 3) : S = 0.49846, interbande 0.24922 = S/2.00 exactement, intrabande 2e-5 (vertex intervallée purement interbande sur l'anneau). Im Π du cône de Dirac avec ces |g|² (Σ wkf = 2, 2 vallées à Γ / 1 à K, ∫d²k δ(2vk−ħω) = πħω/(2v²), v = 5.464 eV·Å) : 0.6685 meV (Γ) et 1.3246 meV (K) ; gamma___ d'EPW = 1.3346 et 2.6448 meV → **rapport 1.997 aux deux points**. Conclusion unique : la valeur imprimée par phonselfen est 2 Im Π, i.e. la largeur totale, malgré le commentaire « half width » de la source ; la convention adoptée en P7 (`gamma_fwhm` = gamma___, `gamma_hwhm` = gamma___/2) est confirmée et la réserve de §1d est levée. Détail : les anneaux ont été placés à |ε − E_D| = 91 et 67 meV au lieu de ħω/2 = 81 et 60 meV (erreur de base réciproque dans le générateur, sans effet : S et la fraction interbande sont constants sur l'anneau).

**3. Fréquences : production vs DFPT direct.**

| mode | ω production 24×24 (degauss 0.002 Ry, K sur la grille k) | ω DFPT direct « prt » (réponse sur 16×16 k, K hors grille) | ⟨D²⟩ production (g EPW) | ⟨D²⟩ DFPT direct |
|---|---|---|---|---|
| E2g @ Γ | 1312.155 cm⁻¹ = 162.69 meV (`phonons/ph.out` l. 535–536 ; `phonons/graphene.freq` l. 3 ; EPW `epw6_phself_path_1200_dg0.02/epw.out` l. 184 sq.) | 1566.502 cm⁻¹ = 194.22 meV (`dfpt_g_K/ph.out` l. 445–446) | 43.6 eV²/Å² | 45.1 eV²/Å² |
| A1' @ K | 972.106 cm⁻¹ = 120.53 meV (`phonons/ph.out` l. 26863, q = (0, −2/3) ≡ K ; EPW l. 43 du chemin) | 1288.884 cm⁻¹ = 159.80 meV (`dfpt_g_K/ph.out` l. 66386, q = (0.5773503, −0.3333333) cartésien) | 86.3 eV²/Å² | 91.6 eV²/Å² |

Grilles : production `phonons/scf.in` l. 39–40 (24 24 1), `phonons/scf.out` l. 112 (61 k irréductibles), `phonons/ph.in` l. 10 (nq 24 24 1). DFPT direct : `dfpt_g_K/scf.in` l. 39–40 aussi 24 24 1 (`scf.out` l. 112 : 61 k, densité identique), mais le calcul de phonons « prt » utilise la liste `dfpt_g_K/nscf.in` l. 40–41 (`K_POINTS tpiba`, 31 points, pas (0.0360844, 0.0625) = |b|/16 : grille 16×16, K absent) — `nscf.out` l. 248 et `ph.out` l. 210 : 31 k (516 k+q) — pour la réponse DFPT. ⟨D²⟩ : `scripts/epw_d2_extract.py` (routes 1a/1b, moyenne par entrée), `results/epw/d2_extract_24k24q.npz`. Lecture : à Γ la production porte une anomalie de Kohn exagérée (−32 meV, K sur la grille avec 27 meV de smearing = terme de surface de Fermi fictif) ; à K, 16×16 sans point à K n'a pas d'anomalie du tout (1289 cm⁻¹) et la production l'exagère (972 cm⁻¹) ; le vertex D² est le même à 3–6 %. Recommandation pour §5.4 : citer les fréquences de production (ce sont celles des largeurs et de Σ) en documentant l'artefact, et donner les valeurs 16×16 comme borne « sans anomalie » ; ne pas présenter les 16×16 comme référence physique (l'anomalie réelle se situe entre les deux).

**4. Γ^ep(ε) près de E_D.** Chiffres (`results/epw/selfen_240_dg0.02_T300.npz` et `_T10`, `scripts/epw_selfen_post.py`) : la courbe moyennée présente à E_D une bosse de 18.9 meV (300 K ; 18.4 meV à 10 K) au-dessus des minima voisins (11.9 / 13.3 meV à ∓0.24 eV), largeur à mi-hauteur 0.185 eV, valeur au sommet 31.5 meV ; à 120² : 9.4 meV de haut, 0.60 eV de large. Or les états individuels près de E_D ont Γ^ep = 0.15 meV (ε = E_D, 300 K), 0.33 meV (±67 meV), 0.65 meV (±114 meV) et ≈ 0 à 10 K. Cause : la moyenne lorentzienne (η = 20 meV) là où la densité d'états s'annule — à 240² le premier état hors K est à ±67 meV avec un poids 6/N contre 2/N à K — est dominée par les queues des états lointains (|ε − E_D| ≳ 0.3 eV, Γ ≈ 10–30 meV), pas par les états locaux ; ce n'est pas l'échantillonnage de l'anneau ħω/2v (correction de la note §1c). Phrase citable : « Sur |ε − E_D| ≲ 0.1 eV, la moyenne lorentzienne (η = 20 meV) de Γ^ep(ε) présente une bosse de ≈ 19 meV de haut et 0.19 eV de large qui est un artefact de la moyenne dans une région de densité d'états nulle : les états individuels y ont Γ^ep < 1 meV (300 K), conformément au blocage de l'espace des phases. » Correction sans nouveau run : dans le post-traitement, tracer les états eux-mêmes ou une médiane par bande d'énergie au lieu de la moyenne lorentzienne près de E_D. Run court possible (non lancé, ≈ 5 min) : elecselfen avec `filkf` = ~100 k sur des anneaux autour de K (|ε − E_D| ≤ 0.3 eV) et nqf = 240², qui peuple la région et fait converger la moyenne ; il ne change pas les états.

**5. Paramètres de production (§5.4).**

| paramètre | valeur | fichier:ligne |
|---|---|---|
| maille (bohr), vide | a1 = (4.0354919061, −2.3298923383, 0), a2 = (4.0354919061, 2.3298923383, 0), c = 30 | `phonons/scf.in` l. 30–33 |
| atomes | C1 (1/3, 1/3, 0), C2 (2/3, 2/3, 0), deux espèces, même C.upf | `phonons/scf.in` l. 26–28, 35–37 |
| pseudo, fonctionnelle | ONCVPSP NC scalaire-relativiste, PBE | `abinit_processing/pseudo/C.upf` l. 6, 71–72, 80 |
| ecutwfc | 100 Ry | `phonons/scf.in` l. 13 ; `epw1/nscf.in` l. 14 |
| occupations | smearing Marzari-Vanderbilt, degauss 0.002 Ry | `phonons/scf.in` l. 14–16 |
| assume_isolated | '2D' | `phonons/scf.in` l. 18 |
| grille k SCF | 24×24×1 (61 irréductibles) | `phonons/scf.in` l. 39–40 ; `scf.out` l. 112 |
| conv_thr, mixing_beta | 1e-14, 0.3 | `phonons/scf.in` l. 22–23 |
| nbnd | 16 (scf), 32 (nscf EPW) | `phonons/scf.in` l. 17 ; `epw1/nscf.in` l. 18 |
| grille q DFPT | 24×24×1, ldisp, tr2_ph 1e-14, alpha_mix 0.3, asr | `phonons/ph.in` l. 8–12 |
| q2r / matdyn | zasr et asr 'crystal' | `phonons/q2r.in`, `phonons/matdyn.in` |
| nscf EPW | 576 k en liste crystal, calculation 'bands' | `epw1/nscf.in` l. 4, 40–41 |
| Wannier | nbndsub 5 ; proj C1:sp2;pz, C2:pz ; dis_win −25 / 15 eV ; dis_froz_max −1.74 eV ; num_iter 5000 ; conv_tol 1e-12, dis_num_iter 5000, dis_conv_tol 1e-12, dis_mix_ratio 0.5, guiding_centres | `epw1/epw1.in` l. 18–24, 29–34, 44 |
| epw1 | epbwrite puis relance epbread, wannierize .false., epwwrite | `epw1/epw1.in` l. 12–17 ; `epw1_restart.in` l. 12–17 |
| grilles grossières EPW | nk 24×24×1, nq 24×24×1 | `epw1/epw1.in` l. 55, 58 |
| elecselfen (production) | nkf = nqf = 240², mp_mesh_k, restart 1000, etf_mem 1 | `epw2_selfen_240_T300_dg0.02/epw.in` l. 9, 26–28, 30, 33 |
| phonselfen (production) | nkf 1200², mp_mesh_k .false., filqf `qpath_GKMG.kpt` (+ `qzoom_GK.kpt`), restart .false. | `epw6_phself_path_1200_dg0.02/epw.in` l. 27–28, 30, 33 |
| fsthick | 3.5 eV (les deux observables) | selfen l. 20 ; phself l. 21 |
| degaussw, ngaussw | 0.02 eV, 0 (gaussienne) | selfen l. 23–24 ; phself l. 24–25 |
| efermi_read, fermi_energy | .true., −4.2389 eV (= E_D) | selfen l. 18–19 ; phself l. 19–20 |
| températures | selfen : 300 K et 10 K (deux runs, `..._T300_...` l. 21–22 / `..._T10_...` l. 21–22) ; phself : temps = 10 300, nstemp 2 | phself l. 22–23 |
| ressources | epw1 16 tâches ; selfen 240² et phself 1200² 64 tâches, 4 Go/cœur | `epw1/submit.epw` l. 5, 7 ; `.../submit.epw` l. 5, 7, 16 |

Écarts fichiers ↔ notes : (a) §1c et §1d disaient « run DFPT prt sur 16×16 k » : précisément, sa densité SCF est 24×24 (61 k) et seule la réponse DFPT (nscf/ph) est sur la liste 16×16 (31 k) ; (b) §1c disait « demi-largeur, à reconfirmer » : levé (point 2) ; (c) la première ligne de commentaire de `epw6_phself_*/epw.in` dit encore « gamma_qnu = Im Pi (demi-largeur EPW) » — commentaire d'entrée sans effet, obsolète ; (d) §1c/§1d attribuaient la bosse à E_D à « l'échantillonnage de l'anneau » : corrigé (point 4). Aucun écart de valeur numérique entre les fichiers d'entrée et les tableaux §1–§1e.

## 2. Plan P0–P5 tel qu'arbitré

Cinq décisions (arbitrage du 2026-09-08) :
1. **degauss 0.002 Ry** (phonons existants conservés) ; l'écart avec le chapitre 4 (0.01 Ry) est documenté, pas recalculé.
2. **Grilles grossières 24k-24q** (pas 30k-30q).
3. **Température basse 10 K** (avec 300 K).
4. **Grilles fines 120² et 240²** ; 360² seulement si l'écart 120² → 240² dépasse 5 %.
5. **fsthick 3.5 eV**.

- P0 (fait) : save/, nscf 24×24 nbnd 32, epw1 avec dis_froz_max = E_D + 2.5 eV = −1.74 eV, dis_win −25/15 eV, projections C1:sp2;pz C2:pz, wdata comme 16k-8q.
- P1 (en cours de post-traitement) : bandes EPW vs bands.x, phonons EPW vs matdyn, ratios de décroissance, |g| DFPT vs EPW **par projection sur les sous-espaces dégénérés à k et k+q (jamais par index de bande)**, aux deux k (Γ, K).
- P2 : elecselfen, `efermi_read = .true.`, `fermi_energy = E_D`, fsthick 3.5, `temps = 300 10`, degaussw ∈ {0.01, 0.02, 0.05} eV, nkf = nqf ∈ {120², 240²} (360² conditionnel), mp_mesh_k, restart ; sortie Σ_nk(T, degaussw, grille).
- P3 : phonselfen (mêmes grilles).
- P4 : post-traitement : **Γ^ep = 2 Im Σ** (largeur totale, même convention que Γ^ed = −2 Im Σ du chapitre 4), npz tagués `units`, référence d'énergie E_D explicite, moyenne lorentzienne η = 0.02 eV comme au chapitre 4, aucun rescalage en sortie.
- P5 : figures (a)–(e) en français, dont Γ^ep contre Γ^ed à c = 1 %.

Après validation P1, proposer (pas exécuter) la liste des `_ph0` à archiver/supprimer (§5).

## 3. Conventions à respecter

- Figures : `figures/memoire.mplstyle`, texte en **français**, panneaux (a)(b)…, PDF + PNG, données lues uniquement depuis npz/save, scripts versionnés ; npz/CSV de production commis, jamais les logs ni les matrices.
- Γ^ep = 2 Im Σ (largeur totale) ; Γ^ed = −2 Im Σ ; les deux en eV (×1e3 → meV), jamais rescalés en sortie.
- Référence d'énergie : E_D (point de Dirac) explicite dans chaque npz (`E_D`), unités taguées (`units='eV'`) ; une seule conversion Ha → eV (HA2EV) au chargement.
- Chiffres bruts, sans conclusion sur la qualité Wannier/physique ; signaler les artefacts (bord de fenêtre, etc.).
- Chemin Γ–K–M–Γ dans ph.x 7.5 : **cartésien 2π/a** (pas de `q_in_cryst_coord` dans ph.x) : K = (0.5773503, −0.3333333), M = (0.2886751, −0.5). Dans EPW/W90 (crystal) : K = (2/3, 1/3), M = (1/2, 0). Le point (1/3, 1/3) n'est PAS K dans la maille QE à 60°.
- Format prtgkk : blocs `iq = … / ik = …`, colonnes ibnd jbnd imode enk enk+q omega |g_sym| |g| Re Im ; `scripts/epw_extract_gkk.py` reproduit octet pour octet les fichiers 16k-8q.
- epw1 multi-nœuds : la tâche 0 rassemble epmatwp ; prévoir ≥ 8 Go pour elle (OOM à 4 Go/cœur sur 24k-24q).
- Ne jamais regarder dans les répertoires d'autres utilisateurs ; l'utilisateur pousse les commits lui-même.

## 4. Écarts avec le chapitre 4 à documenter

| paramètre | chapitre 4 (M, t-matrice) | chapitre 5 (EPW) | vérification |
|---|---|---|---|
| degauss | 0.01 Ry | 0.002 Ry | scf comparés : E_D identique (−4.2389 eV), bandes 1–15 identiques < 0.01 meV, E_F −4.2199 (0.01) vs −4.2351 eV (0.002) |
| grille k de wannierisation | 27×27 (dense du 9×9 ; D = 24–32 selon la taille) | 24×24 | spreads 3.7668 vs 3.7627 Å² ; bandes EPW vs bands.x (P1) |
| nbnd | 20 (`diago_full_acc`) | 32 | fenêtre gelée identique −1.74 eV ; dis_win −25/15 identique |
| E_F / référence | E_D Wannier −4.2389 eV | efermi_read = E_D (P2) | |
| phonons | — | DFPT 24×24 q, degauss 0.002 | matdyn vs EPW (P1) |
| convention de largeur | Γ^ed = −2 Im Σ, η = 0.02 eV, grille 240², R_cut 3 | Γ^ep = 2 Im Σ, degaussw {0.01,0.02,0.05}, grilles 120²/240² | comparer à même η = 0.02 eV |

## 5. En attente de l'utilisateur

- `_ph0` supprimés le 2026-09-10 sur instruction : 30k-30q (334 Go ; **sans** save/ : ses dvscf sont perdus, les dyn*.xml restent, 438 Mo) et 24k-24q (140 Go ; save/ vérifié 123/123 identiques avant suppression). Restent : 16k-16q 30 Go, 12k-12q 9.8 Go, 16k-8q (non arbitrés).
- Commits à pousser (dépôt ab-initio-defects) : 14e07db, a9855af, cb94d4b, d2d0bac, d49f5bb.

## 6. Reprise (prochaine session)

0. P18 (2026-09-21) : chaîne 24k-24q_mv0.02 complète, §1g ; figures génériques = mv0.02, anciennes en `_mv0.002` ; A1' par caractère (P18b). Non arbitré : quelle chaîne le chapitre cite (0.002 documentée en §1c–1f, 0.02 en §1g) ; `_ph0` de 24k-24q_mv0.02 (140 + 56 + 56 Go) à archiver/supprimer après arbitrage.
1. Figures faites et inspectées (`scripts/make_figures_epw.py --prod-tag 240_dg0.02 --phself-tag path_1200_dg0.02`). Facteur 2 tranché empiriquement en P7 (§1d) ; vérifier les valeurs de Piscanec 2004 / Lazzeri 2008 dans Zotero et, si souhaité, lancer le prtgkk « anneau » de 30 s.
2. Commits à pousser (dépôt) : ba151b3, bc83f91 et le commit de retouche de figure ; vérifier `git log`.
3. Reste non arbitré : `_ph0` de 16k-16q, 12k-12q, 16k-8q ; P3 phonselfen si le plan le réintroduit.

## 1g. P18 — chaîne complète refaite à degauss 0.02 Ry (`24k-24q_mv0.02/`, 2026-09-20/21)

Chaîne `24k-24q_mv0.02/` calquée sur `24k-24q/` : **seul changement `degauss = 0.02` Ry** (Marzari-Vanderbilt) dans `phonons/scf.in` l. 17,
`bands/scf.in`, `epw1/nscf.in` l. 18, `dfpt_g_{G,K}/{scf,nscf}.in` ; `diff -r` vérifié (le reste des différences = commentaires de tête,
chemins des submit, logs vers fichiers, epw1 sur 1 nœud à 12 Go/cœur). `degaussw` EPW reste 0.02 eV. Rien d'écrasé dans `24k-24q/`.
Générateurs : `24k-24q_mv0.02/_p18_make_step2.sh` (répertoires EPW, liens vers ../epw1/, jamais selecq.fmt), `_p18_submit_step2b.sh`
(23 runs + post) ; post-traitements versionnés `scripts/submit_epw_p1_post.sh 24k-24q_mv0.02 24k24q_mv0.02` et
`scripts/submit_epw_p2_post_mv.sh`. Tags : `_mv0.02` sur tous les npz (`results/epw/*_mv0.02*.npz`).

**Jobs et coûts** (`sacct`) : scf 21482830 38 s (128 tâches) ; ph.x 21482831 **4 h 51 min** (128 tâches, `_ph0` **140 Go**) ;
pp_save 21482832 34 s (122 fichiers + phsave, `--check` 123/123 OK) ; q2r/matdyn 21482833 ; nscf 21482873 2 min 26 ; **epw1 21482874
4 h 58 min** (16 tâches, 1 nœud, 12 Go/cœur : pas d'OOM, epmatwp 867 Mo écrit du premier coup) ; bandes 21482875–77 ; band_plot /
prtgkk Γ / prtgkk K / anneau 21482890–93 (25–38 s) ; post P1 21483051 ; DFPT « prt » refait à 0.02 : `dfpt_g_G` 21483048–50 (ph.x 1 h 22,
`_ph0` 56 Go), `dfpt_g_K` 21482980–82 (ph.x 1 h 24, 56 Go) ; selfen 240² 300 K 21510232 36 min, 10 K 21510233 50 min, 120² 21510234
9 min ; phself chemin 21510235 7 min 44, zoom 21510236 6 min 40, table K 21510237–54 (18 × ≈ 10 min 30) ; post P2 21510255.
Incident : un `cd` relatif a soumis `dfpt_g_G` en double (21482972–77) ; chaîne annulée, répertoire nettoyé, resoumise (21483048–50),
`scf.out` vérifié (un en-tête, un JOB DONE, E_F −4.2094). Le post P1 initial (21482894) est mort de la dépendance annulée → 21483051.

**Correctif P18b — sélection des modes par caractère, jamais par index.** À 0.02 Ry, A1' à K redevient le mode le plus haut (branche 6,
1274.8 cm⁻¹) alors qu'il était la branche 3 (972 cm⁻¹) à 0.002. Les index codés en dur (`[…, 2]`, `(3,)`, « mode le plus proche de
120.5 meV ») ont été remplacés : `phself.mode_A1p` (A1' = plus grand γ à K), `phself.modes_E2g` (deux modes les plus hauts à Γ, dégénérescence
contrôlée à 0.5 meV) dans `epw_phself_post.py` (valeurs clés, table K), `epw_d2_extract.py` (route 1a : plus grande somme de bloc π/π* ;
route 2 : plus grand γ), `epw_ring_check.py` (plus grande somme de bloc), `make_figures_epw.py` (annotations, barres), `submit_epw_p2_post_mv.sh`.
`epw_extract_gkk.py` ne filtre pas ; `epw_validate.py` groupe les modes par fréquence. Non-régression sur 0.002 : mêmes chiffres qu'en §1c/§1d
(E2g 5+6, A1' mode 3, γ 1.3346/2.6448, ⟨D²⟩ 43.6/86.3, anneau 1.997). Le post P2 21510255 a tourné avec la copie Slurm du script prise
avant le correctif : sa ligne « anneau » a transmis γ_K = 0 ; l'anneau a été refait à la main (`epw_ring_check.py --gamma-G 1.3787
--gamma-K 2.7979 --v 5.4636 --tag 24k24q_mv0.02`), tout le reste du log était déjà avec les bons modes (« [modes 6] »).

**Tableau côte à côte** (chiffres bruts ; « écart » = mv0.02 relatif à mv0.002 ; ★ = écart > 5 % sur une quantité autre que ω et D²).

| quantité | mv0.002 (24k-24q) | mv0.02 (24k-24q_mv0.02) | écart | provenance |
|---|---|---|---|---|
| E_F scf 24×24 | −4.2351 eV | −4.2094 eV | +25.7 meV | `phonons/scf.out` l. 641 / l. 632 |
| E_D (EPW, gap min π/π* sur le chemin) | −4.238884 eV | −4.238898 eV | −0.01 meV | `validation_<tag>.npz` `bands_ED_epw` ; bands.x −1.8540 eV (offset −2.3849) les deux |
| Wannier : étalements sp² / p_z, Ω_tot, Ω_I (Å²) | 0.61056874 (×3) / 0.96547162, 0.96547259 ; 3.76265072 ; 3.020469979 | 0.61056861 (×3) / 0.96547110, 0.96547207 ; 3.76264925 ; 3.020468807 | −4e-5 % | `epw1/graphene.wout` « Final State » ; centres identiques (1.067745, ±0.616462), (2.135490, 0), (1.423660, 0), (2.847320, 0) Å |
| P1 bandes WF vs bands.x, max sur ±3 eV (rms) | π 0.79 (0.40), π* 20.06 (2.86) meV | π 0.81 (0.41), π* 20.06 (2.86) meV | — | log `epw_p1_post_<job>.out` ; bandes DFT elles-mêmes : ≤ 1 meV entre les deux chaînes, EPW ≤ 0.02 meV |
| P1 phonons EPW vs matdyn, max\|Δω\| modes 1–6 (cm⁻¹) | 2.59 / 3.54 / 6.96 / 4.30 / 2.43 / 2.45 | 4.68 / 0.73 / 0.52 / 0.38 / 0.21 / 0.30 | ★ modes 2–6 ÷5 à ÷13, mode 1 +81 % | `validation_<tag>.npz` `ph_err` |
| planchers de décroissance (médiane R ≥ 0.9 R_max, eV) H / C / g(R_e) / g(R_p) | 2.18e-5 / 1.57e-2 / 5.52e-4 / 1.96e-3 (chutes 5.8 / 3.0 / 4.1 / 3.5 ordres) | 2.18e-5 / **5.97e-4** / 5.52e-4 / **5.52e-4** (5.8 / 4.5 / 4.1 / 4.1 ordres) | ★ C −96 %, g(R_p) −72 % | `make_figures_epw.py` bloc 5 (log) ; ratios `decay_*_ratio` : dynmat 9.04e-4 → 3.70e-5, epmatp 3.34e-4 → 1.43e-4 |
| P1 \|g\| k = Γ (sommes gauge-invariantes, G > 20 meV) | 30 q : médiane 0.63 %, max 100 % (q = M) ; hors M : 7 sommes, méd. 0.24 %, max 4.99 % | 30 q : médiane 2.26 %, max 100 % (q = M) ; hors M : 9 sommes, méd. 0.34 %, **max 53.3 %** | ★ | `validation_<tag>.npz` `g_G` ; les deux sommes > 5 % sont à s = 0.614 (q voisin de M), états initiaux σ à −11.87 et −7.30 eV (absolus), groupes de modes 172.5 et 58.5 meV : 52.1 → 26.5 et 63.5 → 29.7 meV ; même famille que le double comptage à q = M |
| P1 \|g\| k = K | 45 triplets : méd. 1.07 %, max 10.3 % ; contrôle 15 sommes | 54 triplets : méd. 0.42 %, max 10.18 % ; contrôle 18 sommes | — | idem `g_K` |
| **contrôle \|g\| retenu (fig_epw_g_control, 2026-09-21)** : toutes les sommes G > 20 meV, **q = M (s = 0.634) et son voisin s = 0.614 exclus** (même regroupement ambigu des modes), sans filtre fsthick | k = Γ : 7 sommes, méd. 0.24 %, max 4.99 % (s = 0.044), \|ΔG\| max 2.32 meV ; k = K : 15 sommes, méd. 1.07 %, max 10.31 % (s = 0.658), \|ΔG\| max 8.45 meV | k = Γ : **7 sommes, méd. 0.23 %, max 4.97 %** (s = 0.044), \|ΔG\| max 2.31 meV ; k = K : **18 sommes, méd. 0.42 %, max 10.18 %** (s = 0.658), \|ΔG\| max 22.70 meV | — | `make_figures_epw.py --control` (bloc « contrôle \|g\| », log) ; la figure n'a plus ni sur-titre ni légende (2026-09-21, cohérence avec le texte) : l'exclusion est documentée ici et dans le log du script ; q retenus : Γ s ∈ {0, 0.044, 0.658}, K s ∈ {0, 0.044, 0.132, 0.658} |
| **\|g\| restreint à la fenêtre fsthick** (complément 2026-09-21 : états initial ET final à \|ε − E_D\| ≤ 3.5 eV, G > 20 meV), k = Γ | 4 sommes, toutes à q = M : méd. 9.18 %, max 45.6 % ; hors M : **aucune** | 5 sommes : méd. 18.0 %, max 53.3 % ; hors M : **1 somme** (s = 0.614, σ à −3.06 eV → π* à +1.66 eV, groupe de modes 58.5 meV : DFPT 63.5 vs EPW 29.7 meV = 53.3 %) | — | `epw_validate.py --fsthick 3.5` (clés `g_G_win_*`, `g_fsthick`) ; sans seuil G : 12 / 12 sommes dans la fenêtre (hors M 2 / 4, G < 1 meV sauf celle-ci) |
| idem, k = K | 1 somme : 1.07 % (s = 0.132, π/π*(K) → π* à +1.56 eV, ω 186 meV, 428.4 vs 423.8 meV) | 1 somme : 1.44 % (166.2 vs 163.8 meV, ω 199 meV) | — | clés `g_K_win_*` ; sans seuil G : 3 / 4 sommes. **Restriction abandonnée pour la figure le 2026-09-21 (une somme par k) : les clés `_win_` restent dans les npz sans être utilisées ; le contrôle retenu est la ligne ci-dessus** |
| ω(E2g, Γ) matdyn / EPW | 1312.16 cm⁻¹ = 162.69 meV | **1550.51 cm⁻¹ = 192.24 meV** | +18.2 % | `phonons/graphene.freq` l. 3 ; phself npz `omega` |
| ω(A1', K) matdyn | 972.11 cm⁻¹ = 120.53 meV (branche 3) | **1274.81 cm⁻¹ = 158.06 meV (branche 6)** | +31.1 % | `graphene.freq` (q = (0.577350, −0.333334) l. 252–253) ; six modes à K : [532.66, 532.66, 972.11, 996.66, 1208.54, 1208.54] → [532.62, 532.62, 996.41, 1216.38, 1216.38, 1274.81] |
| ω(M) LO / TO matdyn | 1343.37 / 1393.81 | 1343.18 / 1393.63 | −0.01 % | `graphene.freq` l. 502–503 |
| DFPT « prt » 16×16 (K hors grille) : ω(E2g,Γ), ω(A1',K), ω(M) LO/TO | 1566.50 ; 1288.88 (mode 6) ; 1343.29 / 1393.76 | 1567.49 ; 1289.16 (mode 6) ; 1343.20 / 1393.68 | +0.06 % / +0.02 % | `dfpt_path_freq_<tag>.npz` (`scripts/epw_dfpt_path_freq.py`, blocs « Diagonalizing » de `dfpt_g_K/ph.out`) |
| ⟨D²⟩ Γ / K, route 1a (g EPW) | 43.6 / 86.3 eV²/Å² (S_g 0.37272 / —) | 45.0 / 91.4 (S_g 0.32586 / 0.40240 eV²) | +3.2 % / +5.9 % | `d2_extract_<tag>.npz` ; D² ∝ S_g·ω : S_g baisse de 12–19 %, ω monte de 18–31 % |
| ⟨D²⟩ route 1b (DFPT prt) | 45.1 / 91.6 | 45.1 / 91.7 | — | idem ; EPW/DFPT : S_D 0.997 / 0.997 |
| v_F EPW (π / π*) | 5.496 / 5.431 → 5.464 eV·Å | 5.496 / 5.431 → 5.464 | — | `d2_extract` log |
| Γ^ep(ε) médiane ±3 / ±1.2 eV (courbe lorentzienne η 0.02), 120² 300 K, dg 0.02 (dg 0.01 ; dg 0.05) | 55.872 / 21.583 meV (55.689 / 21.405 ; 57.486 / 22.137) | 55.390 / 21.288 (55.234 / 20.900 ; 56.592 / 21.782) | −0.9 / −1.4 % (−0.8 / −2.4 ; −1.6 / −1.6 %) | `selfen_120_dg{0.01,0.02,0.05}_mv0.02_T300.npz` (dg 0.01 / 0.05 refaits en complément le 2026-09-21, jobs 21511876 / 21511878, 8 et 7 min) ; table de convergence mv0.02 relative à 120² dg 0.02 (`epw_selfen_post.py --table`) : dg 0.01 −0.28 / −1.82 %, dg 0.05 +2.17 / +2.32 %, 240² 300 K +2.83 / +1.71 %, 240² 10 K −0.72 / −2.43 % (mv0.002 : −0.33 / −0.82, +2.89 / +2.56, +2.32 / +1.39, −1.2 / −2.6 %) ; même conclusion : dg 0.02, 240² |
| idem 240² 300 K (**production**) | 57.168 / 21.882 | **56.956 / 21.653** | −0.37 % / −1.05 % | `selfen_240_dg0.02_mv0.02_T300.npz` ; 120² → 240² : +2.32 / +1.39 % → +2.83 / +1.71 % |
| idem 240² 10 K | 55.204 / 21.010 | 54.991 / 20.772 | −0.39 % / −1.13 % | `selfen_240_dg0.02_mv0.02_T10.npz` |
| médiane par état, ±3 / ±1.2 eV, 240² 300 K | — | 101.120 / 25.016 meV | — | log post P2 |
| van Hove de Γ^ep (240², 300 K) | 151.8 meV à −2.550 eV ; 274.0 à +1.880 | 151.4 à −2.550 ; 272.3 à +1.880 | −0.3 / −0.6 % | `selfen_240_dg0.02*_T300.npz` `Gamma_e` |
| bosse à E_D (240², 300 K) : valeur, minima ∓0.24 eV, hauteur, largeur mi-hauteur | 31.49 ; 11.92 / 13.26 ; 18.9 meV ; 0.185 eV | 31.14 ; 10.94 / 12.30 ; 19.5 meV ; 0.190 eV | hauteur +3 % | idem (10 K : 30.54 / 18.4 → 30.25 / 19.1 ; 120² : 19.33 / 9.4 / 0.595 → 18.87 / 10.1 / 0.505) |
| γ(E2g, Γ) FWHM = gamma___, 10 K / 300 K | 1.3346 / 1.2235 meV (10.76 / 9.87 cm⁻¹) | **1.3787 / 1.3123** (11.12 / 10.58) | ★ +3.3 % / **+7.3 %** | `phself_path_1200_dg0.02_mv0.02.npz` `key_values` ; HWHM = /2 |
| γ(A1', K) FWHM, 10 K / 300 K | 2.6448 / 2.1770 (21.33 / 17.56) | **2.7979 / 2.5444** (22.57 / 20.52) | ★ **+5.8 % / +16.9 %** | idem ; contrôle demandé ≈ 1.33 / 2.64 : bons modes (ω s'élimine à v_F et S_g fixés ; S_g change de −12/−19 %, la variation de γ vient des facteurs de Fermi et du vertex) |
| convergence γ(A1', K) gamma___ 10 K / 300 K (meV), σ = 0.01 / 0.02 / 0.05 | §1c (réf. 1200², 0.02 : 2.6448 / 2.1770 ; 960² +0.8 % ; 0.05 +0.5 %) | 120² : 0 / 0 ; 0 / 0 ; 0.093 / 0.092 — 240² : 0.048 / 0.041 ; 1.772 / 1.526 (−36.7 %) ; 2.714 / 2.378 — 480² : 0.188 / 0.175 ; 1.940 / 1.782 (−30.7 %) ; 2.797 / 2.535 — 720² : 4.512 / 4.114 ; 2.988 / 2.721 (+6.8 %) ; 2.798 / 2.536 — 960² : 2.464 / 2.229 ; 2.787 / 2.533 (−0.40 %) ; 2.798 / 2.536 — **1200² : 2.665 / 2.432 (−4.7 %) ; 2.798 / 2.544 ; 2.798 / 2.536 (+0.01 %)** | même conclusion : 1200², σ 0.02 retenu (960² à −0.4 %, σ 0.05 à +0.01 %) ; à 240² l'écart est −37 % au lieu de +119 % | log post P2 (table `--ref 1200_dg0.02`) |
| anneau (P11) : S(Γ), interbande ; S(K), interbande ; Im Π ; gamma___/Im Π | 0.37275, S/4.00 ; 0.49846, S/2.00 ; 0.6685 / 1.3246 meV ; **1.997 / 1.997** | 0.32589, S/4.00 ; 0.40234, S/2.00 ; 0.6907 / 1.4023 meV ; **1.996 / 1.995** | — | `ring_check_24k24q_mv0.02.npz` (refait à la main, voir P18b) ; énergies d'anneau ±91 et ±67 meV |
| Γ^ed/Γ^ep (c = 1 %, 300 K, ±3 eV) : médiane ; min ; max ; croisements | 0.360 ; 0.088 à +1.875 eV ; 2.006 à −0.755 ; −1.617, −0.224 eV | 0.361 ; 0.088 à +1.875 ; 2.011 à −0.755 ; −1.619, −0.215 | — | `ed_vs_ep_<tag>.npz` (`scripts/epw_ed_vs_ep.py`, Γ^ed interpolé sur la grille EPW : médiane Γ^ed 23.497 meV, §1b citait 23.436 sur sa propre grille) |

Lecture (constats) : ω(E2g, Γ) et ω(A1', K) remontent à 1550.5 et 1274.8 cm⁻¹, à 1.1 % des DFPT « prt » hors grille (1567.5 / 1289.2) ; il
subsiste un adoucissement à K de 14 cm⁻¹ et à Γ de 17 cm⁻¹. La branche A1' redevient la plus haute à K. Le vertex ⟨D²⟩ passe de 0.96/0.94 à
0.99/0.99 de Piscanec (identique au DFPT direct à 0.3 %). Les largeurs γ restent à la limite Dirac-cône × 2 (anneau 1.996/1.995) ; leurs
valeurs montent de 3–17 % (facteurs de Fermi : la baisse 10 K → 300 K passe de −8 %/−18 % à −5 %/−9 %). Γ^ep(ε) et le rapport Γ^ed/Γ^ep
sont inchangés à ≈ 1 %. Les planchers de C(R_p) et g(R_p) baissent d'un ordre de grandeur (partie non analytique de la matrice dynamique
réduite avec les anomalies). Signalés > 5 % (★) : γ à 300 K (+7.3 % Γ, +16.9 % K) et γ(A1') à 10 K (+5.8 %), max de la table \|g\| à k = Γ
hors M (53 %, une seule somme σ → π* à s = 0.614, voisin de M ; en excluant aussi s = 0.614 le contrôle retenu donne Γ max 4.97 %, K max 10.18 %, comme à 0.002), erreurs phonons EPW/matdyn (mode 1 +81 %, modes 2–6 divisées par 5–13), planchers de
décroissance C et g(R_p).

**Figures (étape 3, point décimal — `fr()` sans virgule, `virgule()` neutralisée, 2026-09-21)** : anciennes renommées `fig_epw_*_mv0.002.{pdf,png}`
(git mv, non supprimées) ; nouvelles sous les noms génériques `fig_epw_{validation,gamma,vs_ed,g_control,phonselfen,decay}` (chaîne mv0.02,
`make_figures_epw.py --prod-tag 240_dg0.02_mv0.02 --phself-tag path_1200_dg0.02_mv0.02 --val-tag 24k24q_mv0.02 --sel-suffix _mv0.02 --control
--kohn-val-tags 24k24q,24k24q_mv0.02 --dfpt-tag 24k24q`) ; nouvelle `fig_epw_kohn_degauss` (branches ν = 5, 6 sur Γ–K–M–Γ : matdyn 0.002,
matdyn 0.02, DFPT direct 16×16 aux 31 q du chemin). Les huit figures ont été inspectées. Remarque : sur fig_epw_kohn_degauss, à 0.002 la
branche A1' à K est la 3 (972 cm⁻¹) : elle n'apparaît pas dans ν = 5, 6 pour cette chaîne (les deux branches tracées y sont le doublet E' à
1208.5) ; à 0.02 la branche 6 porte l'anomalie A1' (158.1 meV) et la 5 le E' (150.8). fig_epw_gamma (a) a les quatre courbes (120² σ 0.01 / 0.02 / 0.05, 240² σ 0.02) depuis le complément
du 2026-09-21. Tous les npz `_mv0.02` et `dfpt_path_freq_*`, `ed_vs_ep_*` sont commis ; logs non.

**Étiquettes des figures (2026-09-21, cohérence avec le texte du mémoire, `scripts/make_figures_epw.py`)** : fig_epw_validation —
légendes « DFT / EPW » et « DFPT / EPW » sans parenthèses, axe $\hbar\omega_{\nu\mathbf{q}}$ ; fig_epw_g_control — titres
$\mathbf{k} = \Gamma$ et $\mathbf{k} = K$ seulement, plus de légende ni de sur-titre ; fig_epw_decay — « $H$ : Hamiltonien »,
indices $ij$ (au lieu de $mn$) dans $H_{ij}$ et $g_{ij}^{(\alpha\mu)}$ ; fig_epw_kohn_degauss — titre « Anomalies de Kohn en
fonction de l'élargissement de Marzari-Vanderbilt $\sigma$ », légende « $N\times N$ $\mathbf{k}$, $\sigma$ = X Ry » (24×24 à 0.002
et 0.02 Ry pour matdyn, 16×16 à 0.002 Ry pour le DFPT direct ; options `--kohn-sigmas`, `--dfpt-sigma` remplacent `--kohn-labels`),
axe $\hbar\omega_{\nu\mathbf{q}}$ ; fig_epw_phonselfen — axe $\gamma_{\nu\mathbf{q}}$. Données et chiffres inchangés.
