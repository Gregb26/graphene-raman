# La matrice T ab initio comble un vide optique étroit

Aucun travail publié trouvé ne calcule σ(ω) du graphène avec des lacunes diluées à partir d'une matrice T ab initio. La perspective C occupe donc une niche réelle mais étroite, et chacune de ses briques a un précédent direct. Les modèles de liaisons fortes ont déjà prédit les signatures optiques des diffuseurs résonants : un fond infrarouge (IR) sous 2μ jusqu'à ≈ 0,4 σ₀ à 0,5 % d'impuretés, un pic à ħω ≈ μ, un pic à ħω ≈ t, et un Drude effacé sous μ_c ≈ 0,15 eV à 0,4 % de lacunes. Dans ces modèles, la lacune est toujours un site retiré et le courant est celui de Peierls. Côté ab initio, Kaasbjerg (PRB 101, 045433) construit la matrice T DFT de la lacune de graphène, mais en dc seulement et avec un potentiel indépendant du spin. Lu, Zhou et Bernardi ont déjà la construction M^L + M^NL et son interpolation de Wannier, au niveau de Born. Chan et al. (Nano Lett. 2026) ont la seule optique tirée d'une matrice T ab initio, pour des excitons de MoS₂. Nous n'avons trouvé aucune matrice T ab initio spin-résolue pour la lacune magnétique, ce qui confirme l'intuition de Greg sous réserve d'une recherche de citations incomplète. C'est l'angle le plus neuf, mais aussi le plus fragile : le moment π dépend des choix numériques, et la séparation d'échange agit surtout sous ~0,1–0,2 eV, là où Σ = cT̄ non autocohérente n'est plus contrôlée. Côté expérience, les lacunes n'expliquent pas le résidu de 0,2–0,3 σ₀ sous 2E_F du graphène de bonne qualité, et aucune mesure IR ne suit une densité de lacunes calibrée par I_D/I_G. Verdict : un article de niveau PRB est réaliste si la chaîne est validée contre le Kubo exact en liaisons fortes, si la fenêtre de validité est bornée, et si les résultats sont donnés par unité de densité et traduits en I_D/I_G.

*Légende de fiabilité.* Les étiquettes suivantes indiquent jusqu'où chaque affirmation a été vérifiée :

| Étiquette | Signification |
|---|---|
| [texte] | Contenu lu dans le texte intégral (arXiv, ar5iv ou PMC). |
| [résumé] | Contenu lu dans le résumé seulement. |
| [méta] | Seule la notice bibliographique a été contrôlée (Crossref, arXiv, OpenAlex). |
| [inféré] | Raisonnement des notes ou de ce rapport, sans source directe. |
| [?] | Non vérifié. |

La recherche a été faite le 2026-09-29. Les listes « cité par » de Google Scholar n'ont pas pu être interrogées.

## Les modèles ont déjà cartographié les signatures optiques des lacunes

Depuis 2006, deux familles de travaux calculent σ(ω) du graphène avec lacunes ou diffuseurs résonants.

La première famille insère une self-énergie locale, indépendante de k et **autocohérente**, dans une bulle de Kubo. Il y a deux variantes :
- **FSBA pour les lacunes.** Peres, Guinea et Castro Neto utilisent l'approximation de Born complète autocohérente. Leur bulle ne contient que des produits Im G_AA, et aucun vertex n'est discuté ([Peres 2006](https://arxiv.org/abs/cond-mat/0512091)).
- **CPA pour les diffuseurs unitaires.** Peres, Stauber et Castro Neto la combinent à la RPA pour les impuretés chargées et aux phonons ([Peres EPL 2008](https://arxiv.org/abs/0803.2816) ; [Stauber 2008](https://arxiv.org/abs/0809.2578)).

C'est l'analogue formel le plus proche du schéma « bulle + Σ = cT̄ » de Greg, à une différence près : leur Σ est autocohérente et la sienne ne l'est pas.

La seconde famille résout Kubo **exactement en espace réel** sur des réseaux π de liaisons fortes. Yuan et al. utilisent la propagation temporelle sur 8192² sites à 300 K ([Yuan 2010](https://arxiv.org/abs/1007.3930) ; [Yuan 2011](https://arxiv.org/abs/1109.3485)). Cysne et al. obtiennent σ(ω) complexe par développement de Chebyshev, à 0,4 % de lacunes compensées ([Cysne 2016](https://arxiv.org/abs/1608.04368)). Ces calculs contiennent implicitement le vertex, les diagrammes croisés et la localisation. Ce sont donc les références naturelles pour tester une bulle non autocohérente.

Les signatures prédites sont précises [texte] :
- **Fond dans la fenêtre bloquée.** Les diffuseurs résonants (lacunes, H) laissent un fond dans 0 < ħω < 2μ qui **atteint ≈ 0,4 σ₀ à 0,5 %** et croît quand le dopage baisse. Le désordre non résonant le ramène presque à zéro.
- **Pic à ħω ≈ μ**, dû aux transitions entre la bande d'impureté et le niveau de Fermi.
- **Nouveau pic à ħω ≈ t ≈ 2,7 eV** à 1–10 % d'impuretés. Il vient des transitions entre les états de milieu de bande et les singularités de van Hove, et croît avec la concentration. Le pic à 2t s'étale en même temps ([Yuan 2011](https://arxiv.org/abs/1109.3485)).
- **Lien avec l'expérience.** Yuan et al. reproduisent le fond de Li et al. avec ≈ 0,25 % d'impuretés résonantes. Ils notent toutefois que des amas de potentiel gaussiens (flaques) donnent un fond semblable. **Le fond IR seul n'identifie donc pas le type de défaut.**
- **Drude à 0,4 % de lacunes.** Cysne et al. trouvent ħ/τ ≈ 0,07 eV à μ = 0,5 eV. Le seuil de Fermi est effacé dès μ = 0,2 eV, et **le Drude est complètement lavé sous μ_c ≈ 0,15 eV** ([Cysne 2016](https://arxiv.org/abs/1608.04368)).
- **À μ = 0,** la FSBA prévoit un maximum de basse fréquence qui suit presque une loi en ω/√n_i ([Peres 2006](https://arxiv.org/abs/cond-mat/0512091)).

Ces modèles n'ont pas ce que la perspective C apporte :
- **La lacune est un site p_z retiré.** Le seul intrant DFT est l'ajustement des adatomes de Wehling et al. (V ≈ 2t, ε_d ≈ −t/16), repris par Yuan 2011 ([Wehling 2010](https://arxiv.org/abs/1003.0609)).
- **Le courant est celui de Peierls**, J = −(ie/ħ) Σ t_ij (r_j − r_i) c†_i c_j. La position y est diagonale, sans dipôle intra-maille ni connexion de Berry ([Yuan 2010](https://arxiv.org/abs/1007.3930)).
- **Le vertex n'est jamais quantifié.** Aucun article trouvé ne chiffre la correction de vertex à σ(ω) pour des diffuseurs résonants, ni l'écart entre bulle et Kubo exact.

Deux travaux de modèles méritent d'être cités. Shubnyi et al. traitent par matrice T de Fano–Anderson l'élargissement du Drude par diffusion résonante ([Shubnyi 2019](https://arxiv.org/abs/1903.10363)). Skrypnyk et Loktev montrent que la conductivité minimale se fixe à **l'énergie de la résonance, pas au point de Dirac** ([Skrypnyk–Loktev 2010](https://arxiv.org/abs/1004.4606)). La position ab initio de la résonance a donc des conséquences physiques directes.

## Chaque brique ab initio existe, mais jamais l'assemblage optique

Le travail le plus proche est la théorie atomistique de matrice T de Kaasbjerg ([Kaasbjerg 2020](https://arxiv.org/abs/1911.00530)) [texte, via ar5iv]. Elle traite précisément la lacune de carbone du graphène, ainsi que la substitution N et des défauts de dichalcogénures de métaux de transition (TMD). Sa méthode :
- GPAW, LDA, PAW, base LCAO DZP ;
- super-cellule 11×11 avec 10 Å de vide ;
- potentiel de défaut égal à la différence entre les potentiels cristallins défectueux et vierge, corrections atomiques PAW incluses (l'analogue de M^NL) ;
- **Σ_k(ε) = N_i T_kk(ε), non autocohérente**, exacte au premier ordre en concentration.

L'article montre l'effondrement de l'approximation de Born en 2D : dans les TMD, Born surestime les taux jusqu'à plusieurs ordres de grandeur. Mais il s'arrête à la DOS, à la largeur de raie et à la conductivité dc de Boltzmann. **Il n'y a ni optique, ni vertex, et le potentiel V_i(r) ⊗ ŝ₀ est indépendant du spin.**

L'article évite explicitement l'étape Wannier grâce à sa base fixe. La route « Wannier + amas de 29 mailles » de Greg est donc une variante méthodologique du même objet, pas une nouveauté en soi [inféré]. Semantic Scholar recense 12 citations de Kaasbjerg 2020, et aucune ne calcule σ(ω) du graphène ([Semantic Scholar](https://api.semanticscholar.org/graph/v1/paper/DOI:10.1103/PhysRevB.101.045433/citations)).

Les éléments de matrice ne sont pas neufs non plus. La construction de M^L via ΔṼ_L(q+G) avec les fonctions d'onde de la maille primitive, plus un terme KB non local, à partir de V_KS(d) − V_KS(p) dans Quantum ESPRESSO, est celle de Lu, Zhou et Bernardi ([Lu 2019](https://arxiv.org/abs/1901.03449)) [texte]. L'interpolation de Wannier de M est celle de Lu, Park, Zhou et Bernardi ([Lu 2020](https://arxiv.org/abs/1910.14516)). Ces deux travaux restent au niveau de la règle d'or, pour Si et Cu, en dc. Pour l'article de 2020, le niveau Born est [inféré] du résumé.

La seule matrice T ab initio qui atteint un spectre optique est **excitonique**. Chan, Haber, Naik, Qiu et da Jornada projettent les éléments électron-défaut sur les enveloppes BSE, pour la lacune de soufre de MoS₂. Ils obtiennent l'absorption et la photoluminescence ([Chan 2026](https://arxiv.org/abs/2505.15523)) [texte]. Il n'y a là ni bulle de porteurs libres ni canal de Drude.

Plusieurs autres familles touchent au sujet sans le couvrir :
- **Spectres DFT de super-cellules défectueuses.** Lacunes mono- à tétra- en 5×5 et 6×6 ([Singh 2013](https://arxiv.org/abs/1311.2213)), hydrogène ordonné ([Putz 2014](https://arxiv.org/abs/1309.1016)), Kubo–Greenwood LCAO avec dépliement ([Lee 2020](https://arxiv.org/abs/2006.15894)). Tous travaillent à des concentrations ordonnées de l'ordre du pour cent, sans durée de vie due au désordre.
- **Liaisons fortes paramétrées par DFT + Kubo en espace réel.** Cette communauté reste en dc ([Lherbier 2012](https://arxiv.org/abs/1204.4574) ; [Leconte 2011](https://arxiv.org/abs/1111.3566) ; [Saloriutta 2012](https://arxiv.org/abs/1211.7170)).
- **Absorption par porteurs libres assistée par impuretés, ab initio.** Elle existe au niveau de Born pour des impuretés coulombiennes dans les oxydes ([Peelaers 2015](https://doi.org/10.1103/PhysRevB.92.235201)).
- **Matrice T autocohérente dérivée de la DFT.** Un préprint de juillet 2026 enchaîne DFT, liaisons fortes et SCTMA pour le graphène intercalé Au, avec comparaison ARPES mais sans optique ([Kumari 2026](https://arxiv.org/abs/2607.28296)). La moyenne sur le désordre ab initio du graphène est donc un sujet actif.

| Brique | Ce qui existe | Ce qui manque | Fiabilité |
|---|---|---|---|
| Signatures σ(ω) des lacunes (modèles) | Yuan 2011, Cysne 2016 (Kubo exact, liaisons fortes) ; Peres 2006/2008 (FSBA/CPA + bulle) | Lacune ab initio, courant au-delà de Peierls, vertex chiffré | [texte] |
| Matrice T DFT de la lacune de graphène | Kaasbjerg 2020 (GPAW-LCAO, Σ = N_i T_kk, Boltzmann dc) | Fréquence finie, interbande, spin | [texte] |
| M^L + M^NL (QE) et interpolation Wannier | Lu–Zhou–Bernardi 2019 ; Lu–Park–Zhou–Bernardi 2020 (Born, Si/Cu) | Au-delà de Born, 2D, graphène | [texte] / Born 2020 [inféré] |
| Matrice T ab initio → spectre optique | Chan et al. 2026 (excitons, lacune S de MoS₂) | Porteurs libres, Drude + interbande, semi-métal | [texte] |
| σ(ω) DFT de graphène défectueux | Singh 2013, Putz 2014, Lee 2020 (super-cellules ordonnées) | Régime dilué, moyenne sur le désordre | [texte] / [résumé] |
| Moyenne sur le désordre dérivée de la DFT | Kumari et al. 2026 (SCTMA, Au intercalé, ARPES) | Optique, lacune | préprint [texte] |
| Diffusion spin-résolue dans le graphène | Kochan 2014, Soriano 2011 (hydrogène, modèles ajustés) ; Fedorov 2013 (adatomes non magnétiques) | Lacune, ab initio sans ajustement, optique | [texte] ; méthode Fedorov [?] |
| Vitesse de Wannier avec connexion de Berry dans un σ(ω) désordonné | Rien trouvé ; les modèles utilisent Peierls | Tout | [inféré] ; la littérature des cristaux propres n'a pas été couverte |
| σ(ω) IR mesurée à densité de lacunes calibrée | Rien trouvé ; Chen 2009 (dc), Lee 2016 (hydrogénation), Whelan 2024 (THz, joints de grains) | IR sous 2E_F en fonction de I_D/I_G | absence constatée, pas prouvée |

## La lacune magnétique n'a pas de matrice T ab initio spin-résolue

Tous les calculs DFT polarisés en spin trouvent la monolacune magnétique. Environ 1 μB vient de la liaison pendante σ, qui subsiste après la distorsion Jahn–Teller (géométrie 5-9). L'état π quasi localisé ajoute une part variable, de sorte que **les totaux publiés vont de 1,04 à 2,0 μB**. Les fonctionnelles hybrides donnent 2 μB entiers ([Valencia–Caldas 2017](https://arxiv.org/abs/1611.08246)).

La part π est très sensible aux choix de calcul :
- **Concentration.** Le moment vaut 1,12 à 1,53 μB selon la concentration ([Yazyev–Helm 2007](https://arxiv.org/abs/cond-mat/0610638)). Palacios et Ynduráin trouvent que **les moments π étendus s'annulent à toute concentration expérimentalement pertinente** ([Palacios–Ynduráin 2012](https://arxiv.org/abs/1203.6485)).
- **Smearing.** Le moment reste sensible à la largeur de smearing même avec des milliers de points k ([Wang–Pantelides 2012](https://doi.org/10.1103/PhysRevB.86.165438)) [résumé].
- **Géométrie.** Une structure non plane, non magnétique, ne coûte que ~50 meV ([Padmanabhan–Nanda 2016](https://arxiv.org/abs/1605.03921)).
- **État de spin.** Un calcul multiréférence prédit un triplet ([Casartelli 2013](https://arxiv.org/abs/1303.1924)).

Côté expérience :
- La magnétométrie de graphène irradié mesure un **spin ½ de 0,1 à 0,4 μB par lacune**, sans ordre jusqu'aux températures de l'hélium liquide ([Nair 2012](https://arxiv.org/abs/1111.3775)). Le moment a une double origine, σ et π, en parts à peu près égales ([Nair 2013](https://arxiv.org/abs/1301.7611)).
- La STM voit l'état Vπ **scindé de quelques dizaines de meV**, avec un g effectif ≈ 40 ([Zhang 2016](https://arxiv.org/abs/1604.06542)).
- Un effet Kondo à T_K = 30–90 K a été revendiqué ([Chen 2011](https://arxiv.org/abs/1004.3373)), puis contesté au profit d'une interaction électron-électron ([Jobst–Weber 2012](https://doi.org/10.1038/nphys2297)). La STM montre ensuite un écrantage Kondo contrôlé par la grille et la courbure ([Jiang 2018](https://doi.org/10.1038/s41467-018-04812-6)).

La littérature sur la diffusion dépendante du spin traite presque exclusivement **l'hydrogène**, avec des modèles ajustés à la DFT :
- Kochan, Gmitra et Fabian écrivent que les moments « sur des diffuseurs résonants, comme les lacunes et les adatomes » sont des points chauds de spin. Ils ne calculent pourtant que l'hydrogène. Leur échange est ajusté (J_h = −0,82 eV), puis remplacé par un J générique de −0,4 eV ([Kochan 2014](https://arxiv.org/abs/1306.0230)) [texte].
- Soriano et al. font du Kubo en champ moyen de Hubbard, là encore pour l'hydrogène ([Soriano 2011](https://arxiv.org/abs/1105.1005)).
- Fedorov et al. sont entièrement ab initio, mais pour des adatomes non magnétiques avec couplage spin-orbite ([Fedorov 2013](https://doi.org/10.1103/PhysRevLett.110.156602)).
- Les pipelines de Kaasbjerg (ŝ₀), de Bernardi et de Chan ne sont pas magnétiques.

**Nous n'avons trouvé aucune matrice T ab initio spin-résolue de la lacune, ni aucune résistivité à deux courants ou σ(ω) spin-résolue.** L'intuition de Greg est confirmée, avec la réserve d'une recherche par mots-clés.

Le coût marginal serait faible. Dans pp.x, `spin_component = 1` ou `2` avec `plot_num = 1` donne directement V^↑ et V^↓ ([QE INPUT_PP](https://www.quantum-espresso.org/Doc/INPUT_PP.html)). Les super-cellules 9×9 relaxées en nspin = 1 et 2 existent déjà (campagne R1 du dépôt).

Les limites sont toutefois nettes [inféré] :
- **Canaux.** Un calcul colinéaire donne deux canaux qui conservent le spin. Il n'y a ni spin-flip (les points chauds de Kochan exigent un opérateur s·S), ni Kondo, ni réponse à la question spin 1 contre spin ½.
- **Nature de l'écart ↑/↓.** Avec une référence vierge commune, Ṽ^↑ − Ṽ^↓ se réduit au potentiel d'échange-corrélation de magnétisation.
- **Domaine d'énergie.** Surtout, la séparation d'échange de dizaines de meV modifie Γ(ε) près de la résonance. Aux énergies optiques de 1–2 eV, l'écart ↑/↓ est probablement faible ; c'est une hypothèse à tester.

## Les lacunes seules n'expliquent pas le résidu expérimental sous 2E_F

Toutes les mesures IR à grille sur monocouche trouvent un blocage de Pauli incomplet :
- **Li et al.** (exfolié, 45 K) estiment le résidu **jusqu'à 0,3 σ₀, presque indépendant de la grille**. Le bord à 2E_F est large de **~1400 cm⁻¹ (≈ 174 meV)**, contre ~500 cm⁻¹ attendus de la température ([Li 2008](https://arxiv.org/abs/0807.3780)) [texte, version arXiv].
- **Horng et al.** trouvent un poids de Drude **inférieur de 20 à 45 %** à la valeur théorique. La règle de somme relie ce déficit exactement au résidu interbande. Les auteurs jugent que ni la diffusion unitaire pure ni les impuretés chargées ne suffisent ([Horng 2011](https://arxiv.org/abs/1007.4623)).
- **Yan et al.** (CVD sur quartz) voient l'absorption passer de 1,9 % à 0,4 % et un poids de Drude inférieur de 28 à 35 %. Le déficit varie d'un échantillon à l'autre, et certains échantillons n'en ont pas. L'inhomogénéité de dopage vue en Raman élargit le bord ([Yan 2011](https://arxiv.org/abs/1111.3714)).
- **Mécanisme concurrent.** Sur substrat polaire, l'absorption assistée par phonons peut atteindre **20–25 % de σ₀** dans la fenêtre bloquée ([Scharf 2013](https://doi.org/10.1103/PhysRevB.87.035414)) [résumé].

Point décisif : **Li et al. ont testé la théorie à diffuseurs unitaires** (« edge defects, cracks, vacancies »). Ils l'ont jugée trop faible et affectée de la mauvaise dépendance en grille, car elle prédit un résidu qui décroît avec la tension ([Li 2008](https://arxiv.org/abs/0807.3780)). Cette théorie supposait n_i = 4×10⁻⁵ diffuseurs unitaires par carbone ([Peres EPL 2008](https://arxiv.org/abs/0803.2816) ; [Stauber 2008](https://arxiv.org/abs/0809.2578)).

La calibration Raman rend l'argument quantitatif. Dans le régime L_D ≥ 10 nm, **n_D (cm⁻²) = (1,8 ± 0,5)×10²² λ_L⁻⁴ (I_D/I_G)**. Le rapport I_D/I_G est non monotone et passe par un maximum vers L_D ≈ 3 nm ([Cançado 2011](https://arxiv.org/abs/1105.0175)). Deux corrections s'ajoutent :
- I_D/I_D′ ≈ 7 signale des défauts de type lacune, contre ≈ 13 pour les sp³ ([Eckmann 2012](https://arxiv.org/abs/1207.2058)).
- L'intensité D baisse avec le dopage ([Bruna 2014](https://doi.org/10.1021/nn502676g)).

Le tableau ci-dessous est [inféré] : il découle de la formule de Cançado à 514 nm, avec 3,82×10¹⁵ C/cm².

| Fraction de lacunes x | n_D (cm⁻²) | L_D (nm) | I_D/I_G | Régime |
|---|---|---|---|---|
| 7×10⁻⁷ | 2,7×10⁹ | ~110 | ~0,01 | Pic D de ~1 % de G, typique des bons échantillons ([Ni 2010](https://doi.org/10.1021/nl101399r)) |
| 4×10⁻⁵ (théorie comparée par Li) | 1,5×10¹¹ | ~14 | ~0,6 | Pic D bien visible |
| 10⁻⁴ | 3,8×10¹¹ | ~9 | ~1,5 | Limite de validité de la formule |
| 2,5×10⁻³ (ajustement de Yuan sur Li) | 9,5×10¹² | ~1,8 | — | Au-delà du maximum : défauts massifs |
| 1/162 (super-cellule 9×9) | 2,4×10¹³ | ~1,2 | — | Stade 2, formule invalide |

Il faudrait donc **environ 10³ fois plus de lacunes** que ce qu'autorise un pic D de 1 % pour produire le fond que Yuan ajuste sur Li [inféré]. Un σ(ω) de lacunes diluées n'« expliquera » pas le résidu du graphène exfolié de qualité. Il peut en revanche **borner** la contribution des lacunes, et prédire le signal du graphène irradié.

Côté ω → 0, le meilleur ancrage est Chen et al. Sur du graphène irradié Ne⁺/He⁺ à 500 eV, ils trouvent **1/μ linéaire en dose**. Cette pente correspond à un rayon de lacune de 2,3–2,9 Å dans la formule des états de milieu de bande, et elle est 4 fois plus forte que pour des impuretés chargées ([Chen 2009](https://arxiv.org/abs/0903.2602)).

Les autres mesures sur défauts contrôlés sont rares :
- **Hydrogénation.** Lee et al. voient l'absorption interbande de Dirac et le pic du point M **diminuer** avec la couverture en H, avec un IR lointain non-Drude ([Lee 2016](https://doi.org/10.1016/j.carbon.2016.03.008)) [résumé]. La diminution de poids spectral compte donc autant que l'élargissement.
- **Pic UV.** Le pic est un exciton de Fano à 4,62 eV, 0,6 eV sous l'énergie GW ([Mak 2011](https://arxiv.org/abs/1012.2922)). Une DFT sans BSE le place mal ; il faut comparer des variations relatives.
- **Aucune expérience ne mesure σ(ω) IR dans la fenêtre de Pauli en fonction d'une densité de lacunes calibrée par Raman.** Le seul couplage Raman–conductivité optique sur un même échantillon concerne le THz et les joints de grains ([Whelan 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10850153/)).

## Quatre objections de rapporteur, et leurs parades

| Objection | Pourquoi elle porte | Parade proposée |
|---|---|---|
| **Vertex absent** | Pour un diffuseur court-portée, le vertex en échelle est nul dans la SCTMA ([Ostrovsky 2006](https://arxiv.org/abs/cond-mat/0609617)) [texte]. Un argument C₃ l'annule aussi pour une matrice T d'onde s sur site [inféré]. Mais une matrice T sur un amas de 29 mailles a des canaux anisotropes, donc un vertex a priori non nul, qui croît vers 1–3 eV [inféré]. La bulle donne 1/τ et non 1/τ_tr. | Comparer la limite ω → 0 de la bulle au Boltzmann τ_tr calculé avec la même matrice T (comme Kaasbjerg) : c'est une mesure directe du vertex en dc. |
| **Non-autocohérence près de E_D** | Im Σ ∝ c/(\|E\| ln²(D/\|E\|)) diverge ([Peres 2006](https://arxiv.org/abs/cond-mat/0512091)). L'autocohérence régularise sous Γ_c ∝ √c, à des logarithmes près ([Ostrovsky 2006](https://arxiv.org/abs/cond-mat/0609617)). Estimation Γ_c ~ Δ√(c/ln(1/c)) ≈ 0,1 eV à 0,1 % et ≈ 0,2 eV à 0,4 % [inféré, à vérifier sur l'éq. (56) d'Ostrovsky], cohérente avec μ_c ≈ 0,15 eV de Cysne. Plus bas, modes zéro et classes BDI/Gade échappent à toute matrice T à une impureté ([Häfner 2014](https://arxiv.org/abs/1404.6138) ; [Ostrovsky 2010](https://arxiv.org/abs/1006.3299)). | Ne rien revendiquer quantitativement pour \|μ\|, ħω ≲ Γ_c. Faire tourner la lacune de liaisons fortes dans la même chaîne et la comparer au Kubo exact de Yuan/Cysne (ou KITE, [João 2020](https://arxiv.org/abs/1910.05194)). Cela mesure le vertex et l'autocohérence d'un coup. |
| **Magnétisme de la lacune** | Le moment π, qui porte la résonance, est la quantité DFT la plus fragile (smearing, fonctionnelle, concentration, planéité). Le choix nspin déplace la résonance, et la conductivité minimale suit la résonance ([Skrypnyk–Loktev 2010](https://arxiv.org/abs/1004.4606)). L'expérience mesure 0,1–0,4 μB et un possible écrantage Kondo. | Présenter nspin = 1 et nspin = 2 comme deux bornes. Tester la sensibilité au smearing et à N. Comparer la séparation ↑/↓ de la résonance aux dizaines de meV de la STM. |
| **Concentration des super-cellules contre échantillons réels** | 1/162 en 9×9, contre 10⁻⁶–10⁻⁴ dans les échantillons. Comme c entre linéairement dans Σ, la vraie question est la convergence de Ṽ_ed : images, relaxation, et artefact ΔE_F de 0,2–0,9 eV pour N ≠ 3m (fait du dépôt). S'y ajoute la disparition du moment π à faible concentration ([Palacios–Ynduráin 2012](https://arxiv.org/abs/1203.6485)). | Converger la matrice T sur la famille N = 6, 9, 12. Publier Δσ/σ₀ par 10¹¹ cm⁻² de lacunes, traduit en I_D/I_G. |

Deux objections secondaires suivront. La première : **« qu'apporte l'ab initio par rapport à Yuan 2011 ? »**. La réponse doit être une différence chiffrée entre la lacune ab initio et le site retiré, calculée avec la même bulle : position de la résonance, canaux σ, asymétrie électron-trou. La vitesse avec connexion de Berry est un plus de précision, pas un argument de nouveauté à elle seule [inféré]. La seconde est la **pertinence expérimentale**, qu'on traite en présentant le calcul comme une prédiction pour le graphène irradié et comme une borne sur le résidu de Li.

## Verdict : une niche PRB réelle, à condition de la borner

**La niche existe et elle est publiable, vraisemblablement dans Phys. Rev. B ou Phys. Rev. Materials [inféré], mais pas comme « nouvelle méthode ».** Chaque brique a un précédent direct, et un rapporteur citera Kaasbjerg, Lu–Bernardi, Yuan et Chan dès la première page. Ce qui est défendable, « à notre connaissance », est la combinaison : le premier σ(ω) (Drude + interbande) du graphène avec lacunes diluées tiré d'une matrice T ab initio non ajustée et d'une vitesse de Wannier complète.

Trois conditions en font un résultat, et la première n'est pas négociable :
1. **Banc d'essai en liaisons fortes.** Il répond à la fois au vertex et à l'autocohérence.
2. **Fenêtre de validité explicite**, |μ|, ħω ≫ Γ_c.
3. **Résultats par unité de densité**, traduits en I_D/I_G.

Le contenu physique qui porterait l'article est l'écart ab initio / site retiré sur le fond IR, le pic à ħω ≈ μ et le pic à ħω ≈ t. Si cet écart est faible, l'article se réduit à une validation.

**L'angle de la matrice T spin-polarisée est le plus original et le moins cher.** On n'y trouve aucun précédent, Kaasbjerg utilise ŝ₀, et les potentiels V^↑ et V^↓ sont déjà accessibles. Mais il est physiquement plus solide pour **les largeurs Γ^σ(ε) et la résistivité à deux courants** que pour le σ(ω) IR, parce que son effet se concentre là où Σ = cT̄ est la moins fiable. La recommandation se décompose ainsi :
- Dans l'article optique, présenter nspin = 1 et nspin = 2 comme bornes.
- Réserver la matrice T spin-résolue (Γ^↑, Γ^↓, séparation de la résonance comparée à la STM) à une section dédiée ou à une courte lettre compagnon.
- Avant toute revendication de priorité, balayer les listes « cité par » de Kaasbjerg 2020, Kochan 2014 et Yazyev–Helm 2007, et chercher du côté KKR (groupe Mertig), l'endroit le plus probable d'un précédent caché.
- Garder en tête le risque de concurrence [inféré] : Kaasbjerg cite les potentiels dépendants du spin comme extension naturelle, le groupe da Jornada maîtrise la matrice T optique, et le groupe Bena fait de la SCTMA dérivée de la DFT sur le graphène.

## Vérifié, inféré, incertain

**Vérifié sur le texte :**
- l'absence d'optique et de spin chez Kaasbjerg 2020 ;
- les signatures de Yuan 2011 et Cysne 2016 ;
- le test des diffuseurs unitaires par Li 2008, et leur résidu de 0,3 σ₀ (version arXiv, la version publiée est modifiée) ;
- la formule de Cançado ;
- les chiffres de Nair 2012 et de Chen 2009 ;
- le fait que Kochan 2014 traite l'hydrogène et non la lacune.

**Inféré :**
- l'estimation de Γ_c (dérivation propre, à confronter à Ostrovsky) ;
- l'argument de symétrie sur le vertex ;
- le tableau x → I_D/I_G et l'écart de ~10³ ;
- le niveau Born de Lu 2020 ;
- l'hypothèse d'un faible écart ↑/↓ aux énergies optiques ;
- la cible de revue.

**Incertain ou non vérifié :**
- la complétude de la recherche de citations (Semantic Scholar seulement, Google Scholar inaccessible) ;
- une éventuelle étude KKR spin-polarisée de la lacune de graphène ;
- la méthode de Fedorov 2013 ;
- l'annulation des vertex en CPA mono-site attribuée à Velický (1969) ;
- la largeur de bord de Yan 2011 (« 2001 cm⁻¹ » dans le texte extrait, probablement ~2000 cm⁻¹) ;
- les constantes de Lucchese 2010 ;
- le contenu optique éventuel de la revue de Peres (RMP 2010) et de Skrypnyk–Loktev 2006/2007.

Une incohérence de citation est aussi à noter : la réf. 11 de Li 2008 imprime PRB 78, 085418 avec le préprint arXiv:0803.2816, qui est en réalité l'article EPL 84, 38002.

## Références principales et statut de vérification

| Référence | Identifiant | Statut |
|---|---|---|
| Peres, Guinea, Castro Neto, PRB 73, 125411 (2006) | [10.1103/PhysRevB.73.125411](https://doi.org/10.1103/PhysRevB.73.125411) ; arXiv:cond-mat/0512091 | [texte] |
| Peres, Stauber, Castro Neto, EPL 84, 38002 (2008) | [10.1209/0295-5075/84/38002](https://doi.org/10.1209/0295-5075/84/38002) ; arXiv:0803.2816 | [texte] |
| Stauber, Peres, Castro Neto, PRB 78, 085418 (2008) | [10.1103/PhysRevB.78.085418](https://doi.org/10.1103/PhysRevB.78.085418) ; arXiv:0809.2578 | [texte] |
| Yuan, De Raedt, Katsnelson, PRB 82, 115448 (2010) | [10.1103/PhysRevB.82.115448](https://doi.org/10.1103/PhysRevB.82.115448) ; arXiv:1007.3930 | [texte] |
| Yuan, Roldán, De Raedt, Katsnelson, PRB 84, 195418 (2011) | [10.1103/PhysRevB.84.195418](https://doi.org/10.1103/PhysRevB.84.195418) ; arXiv:1109.3485 | [texte] |
| Cysne et al., PRB 94, 235405 (2016) | [10.1103/PhysRevB.94.235405](https://doi.org/10.1103/PhysRevB.94.235405) ; arXiv:1608.04368 | [texte] |
| Ostrovsky, Gornyi, Mirlin, PRB 74, 235443 (2006) | [10.1103/PhysRevB.74.235443](https://doi.org/10.1103/PhysRevB.74.235443) ; arXiv:cond-mat/0609617 | [texte] (préfacteur de Γ_η mal extrait) |
| Skrypnyk, Loktev, PRB 82, 085436 (2010) | [10.1103/PhysRevB.82.085436](https://doi.org/10.1103/PhysRevB.82.085436) ; arXiv:1004.4606 | [texte] |
| Shubnyi et al., PRB 99, 235421 (2019) | [10.1103/PhysRevB.99.235421](https://doi.org/10.1103/PhysRevB.99.235421) ; arXiv:1903.10363 | [résumé] |
| Wehling et al., PRL 105, 056802 (2010) | [10.1103/PhysRevLett.105.056802](https://doi.org/10.1103/PhysRevLett.105.056802) ; arXiv:1003.0609 | [texte] |
| Kaasbjerg, PRB 101, 045433 (2020) | [10.1103/PhysRevB.101.045433](https://doi.org/10.1103/PhysRevB.101.045433) ; arXiv:1911.00530 | [texte] (ar5iv) |
| Lu, Zhou, Bernardi, PR Materials 3, 033804 (2019) | [10.1103/PhysRevMaterials.3.033804](https://doi.org/10.1103/PhysRevMaterials.3.033804) ; arXiv:1901.03449 | [texte] |
| Lu, Park, Zhou, Bernardi, npj Comput. Mater. 6, 17 (2020) | [10.1038/s41524-020-0284-y](https://doi.org/10.1038/s41524-020-0284-y) ; arXiv:1910.14516 | [résumé] |
| Chan, Haber, Naik, Qiu, da Jornada, Nano Lett. 26, 961 (2026) | [10.1021/acs.nanolett.5c04479](https://doi.org/10.1021/acs.nanolett.5c04479) ; arXiv:2505.15523 | [texte] (PMC) |
| Kumari et al., préprint (2026) | [arXiv:2607.28296](https://arxiv.org/abs/2607.28296) | [résumé], préprint |
| Putz, Gmitra, Fabian, PRB 89, 035437 (2014) | [10.1103/PhysRevB.89.035437](https://doi.org/10.1103/PhysRevB.89.035437) ; arXiv:1309.1016 | [résumé] |
| Singh et al., APL 102, 023101 (2013) | [10.1063/1.4781382](https://doi.org/10.1063/1.4781382) ; arXiv:1311.2213 | [résumé] |
| Lherbier et al., PRB 86, 075402 (2012) | [10.1103/PhysRevB.86.075402](https://doi.org/10.1103/PhysRevB.86.075402) ; arXiv:1204.4574 | [résumé] |
| Yazyev, Helm, PRB 75, 125408 (2007) | [10.1103/PhysRevB.75.125408](https://doi.org/10.1103/PhysRevB.75.125408) ; arXiv:cond-mat/0610638 | [résumé] |
| Palacios, Ynduráin, PRB 85, 245443 (2012) | [10.1103/PhysRevB.85.245443](https://doi.org/10.1103/PhysRevB.85.245443) ; arXiv:1203.6485 | [résumé] |
| Wang, Pantelides, PRB 86, 165438 (2012) | [10.1103/PhysRevB.86.165438](https://doi.org/10.1103/PhysRevB.86.165438) | [résumé] |
| Valencia, Caldas, PRB 96, 125431 (2017) | [10.1103/PhysRevB.96.125431](https://doi.org/10.1103/PhysRevB.96.125431) ; arXiv:1611.08246 | [résumé] |
| Kochan, Gmitra, Fabian, PRL 112, 116602 (2014) | [10.1103/PhysRevLett.112.116602](https://doi.org/10.1103/PhysRevLett.112.116602) ; arXiv:1306.0230 | [texte] |
| Soriano et al., PRL 107, 016602 (2011) | [10.1103/PhysRevLett.107.016602](https://doi.org/10.1103/PhysRevLett.107.016602) ; arXiv:1105.1005 | [résumé] |
| Fedorov et al., PRL 110, 156602 (2013) | [10.1103/PhysRevLett.110.156602](https://doi.org/10.1103/PhysRevLett.110.156602) | [méta] ; méthode [?] |
| Nair et al., Nat. Phys. 8, 199 (2012) | [10.1038/nphys2183](https://doi.org/10.1038/nphys2183) ; arXiv:1111.3775 | [texte] |
| Nair et al., Nat. Commun. 4, 2010 (2013) | [10.1038/ncomms3010](https://doi.org/10.1038/ncomms3010) ; arXiv:1301.7611 | [résumé] (une note le laisse non vérifié) |
| Zhang et al., PRL 117, 166801 (2016) | [10.1103/PhysRevLett.117.166801](https://doi.org/10.1103/PhysRevLett.117.166801) ; arXiv:1604.06542 | [résumé] |
| Li et al., Nat. Phys. 4, 532 (2008) | [10.1038/nphys989](https://doi.org/10.1038/nphys989) ; arXiv:0807.3780 | [texte] (arXiv, la version publiée diffère) |
| Horng et al., PRB 83, 165113 (2011) | [10.1103/PhysRevB.83.165113](https://doi.org/10.1103/PhysRevB.83.165113) ; arXiv:1007.4623 | [texte] |
| Yan et al., ACS Nano 5, 9854 (2011) | [10.1021/nn203506n](https://doi.org/10.1021/nn203506n) ; arXiv:1111.3714 | [texte] |
| Chen, Cullen, Jang, Fuhrer, Williams, PRL 102, 236805 (2009) | [10.1103/PhysRevLett.102.236805](https://doi.org/10.1103/PhysRevLett.102.236805) ; arXiv:0903.2602 | [texte] |
| Cançado et al., Nano Lett. 11, 3190 (2011) | [10.1021/nl201432g](https://doi.org/10.1021/nl201432g) ; arXiv:1105.0175 | [texte] |
| Eckmann et al., Nano Lett. 12, 3925 (2012) | [10.1021/nl300901a](https://doi.org/10.1021/nl300901a) ; arXiv:1207.2058 | [résumé] |
| Lee et al., Carbon 103, 109 (2016) | [10.1016/j.carbon.2016.03.008](https://doi.org/10.1016/j.carbon.2016.03.008) | [résumé] |
| Mak, Shan, Heinz, PRL 106, 046401 (2011) | [10.1103/PhysRevLett.106.046401](https://doi.org/10.1103/PhysRevLett.106.046401) ; arXiv:1012.2922 | [texte] |

## Conclusion

Le paradoxe central de la perspective C est que la physique qui rend l'optique des lacunes intéressante se situe toute sous ~0,1–0,2 eV. C'est le cas de la résonance au point de Dirac, du magnétisme, du fond IR qui grandit quand le dopage baisse, et du pic à ħω ≈ μ. C'est justement là que Σ = cT̄ non autocohérente et la bulle sans vertex sont le moins contrôlées. À l'inverse, le régime où la méthode est solide (μ ≳ 0,2–0,3 eV, visible) est celui où les signatures des lacunes sont faibles. L'amélioration la plus rentable n'est donc probablement pas le vertex mais une self-énergie autocohérente de type FSBA/SCTMA. Le résultat le plus défendable, à court terme, est une **borne quantitative ancrée sur I_D/I_G** plutôt qu'une explication du résidu expérimental [inféré].

La matrice T spin-résolue suit la même logique. Elle est neuve et presque gratuite avec les calculs R1 existants, mais sa valeur scientifique tient d'abord aux largeurs Γ^↑, Γ^↓ et à la séparation de la résonance, comparables directement à la STM, avant l'optique. Trois groupes actifs (Kaasbjerg, da Jornada, Bena) sont à une extension de ce terrain. La fenêtre de priorité est réelle, mais elle n'est probablement pas longue [inféré].
