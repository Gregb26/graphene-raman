# R4 — Diagnostic de l'état quasi-lié de la lacune (préparation article, hors mémoire)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 8 figures de `fig/` (hors mémoire) et `submit_r4.sh`. `r4_driver.py` est **gardé uniquement parce que `../R5_base_vs_M/r5_driver.py` l'importe**, lui-même importé par `../R7_tailles_3m/r7_driver.py` pour deux figures du mémoire ; il partira avec eux quand ces figures seront refaites. Gardés : rapport, `sections/`, json et npz de `prep/ d1/ d3/ d4/ d5/ d6/ r/`. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R4_quasi_lie/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

> **Mise au style du 2026-09-30 (étape F0)** : `r4_driver.py`, gardé seulement parce que R7 l'importait, est retiré. Les fonctions de tracé ont été extraites, à contenu identique, dans `scripts/fig/make_figures_controles.py`, qui lit les json/npz de cette campagne et écrit `figures/electron_defect/fig_*` ; rien d'autre ne change. Dernier commit qui contient ces fichiers : `d9a3588`.

## Fiche R4 — État quasi-lié de la lacune 9×9 : vérité DFT contre la chaîne M → Wannier → T
- **Question** : où la DFT place-t-elle les états de la lacune (π et σ) de la 9×9, et quelle étape de la chaîne M → Wannier → matrice T s'en écarte ? Le facteur α de D3 sert au diagnostic ; ce n'est jamais un paramètre.
- **Pourquoi** : les rapports ne donnent pas d'autre raison que ce but (non consigné). La phase 0 retrouve les valeurs de production à confronter à la DFT : min |det|/max = 2,09e-4 à −2,530 eV, pic de −Im T̄ à −1,29 eV, p_z–p_z sur site = 6,617 eV.
- **Méthode** : post-traitement seulement (pp.x sur R1 en J0, aucun pw.x). Fenêtre ε − E_D ∈ [−3, +1] eV, parité miroir ⟨σ_h⟩, poids w₂ dans un disque de 2 Å (seuil 0,0812). D6 : banc liaisons fortes (jouet t = 2,7 eV, puis vrai H(R)). D1 : vérité DFT (alignement Lu 2019). D4 : escalier QE → (a1) Bloch 16 bandes, H = diag ε + M/81 → (a2)/(a3) M dense 16/20 bandes → (b) 5 WF → (c-all)/(c-3) H(R) + M_W repliés. D2 : éléments de M. D3 : critère de pôle avec α·M_loc, par bloc de parité. D5 : ΔV^L pour N = 5…12. R : lecture de code pour les géométries relaxées.
  Jobs : J0 21796633/4, J1 21796852, J2 21796853, J3 21796854, J4 21797745 → J6 21797985, J5a 21797930, J5b 21797931.
- **Résultat** :
  - D1 (QE, E_D = −4,23847 eV) : doublet σ localisé à +0,101 eV (w₂ 0,713). États π localisés à −1,759 ×2 (0,114), −0,737 (0,246) et +0,269 eV (0,102). Sur R1 relaxée, le σ est à −0,002 eV (nspin1) et à −0,850 eV (nspin2 ↑).
  - D4 : la grosse marche est QE → (a1). Le π passe de −0,737 à −1,350 eV (Δ −0,61) et le σ de +0,101 à −2,812 eV (Δ −2,91). De (a1) à (c-3), ce sont ensuite les mêmes états localisés à ≤ 0,3 eV près. Portes de repliement : 4,9e-14 eV et 2,4e-12 eV (avec k exacts).
  - D3 : porte de régression retrouvée (2,0908e-4 à −2,530 eV). Ce minimum est un doublet σ : vecteur propre à 100 % σ, dont 0,976 sur les sp² de l'atome retiré. Le bloc π n'approche jamais zéro : min |λ| = 0,775 (300²) / 0,789 (600²).
  - D6 : le jouet retrouve les racines attendues à ≤ 0,0039 eV (1200²). U_c = 5,94 eV (jouet) contre 5,67 eV (vrai H(R)).
  - D2 : p_z–p_z sur site = 6,6166 = 0,3113 (M^L) + 6,3053 (M^NL) eV. D5 : ⟨ΔV^L⟩_3D vaut 28,6 meV en 5×5, 9,1 en 9×9 et 5,1 meV en 12×12.
- **Interprétation / décision** : le rapport donne les chiffres bruts, sans interprétation physique. Constat retenu : l'écart apparaît dès la marche (a1), donc avant Wannier et la matrice T. D'où la campagne R5 : la cause est-elle la base tronquée ou M ? L'erratum R5 (2026-09-25) corrige deux points :
  - (i) les w₂ de (a2)–(c-3) étaient faux (index d'ondes planes construit avec k81 au lieu des k du dense). Après correction, la paire σ est localisée dans toutes les variantes (w₂ 0,254–0,490).
  - (ii) R5 A.2 : dans les M de production v1 (`results/M/`), M^L est en norme super-cellule, donc N_cells = 81 fois trop petit face à M^NL. Tous les chiffres R4 qui passent par M (D2, D3, D4, bloc π de D6.2) sont dans cette convention fausse. La correction a été confiée à R6.
- **Où sont les données** : `sections/tables.md` et `sections/D1.md`…`D6.md`, `R.md` ; `prep/prep_results.json`, `prep/M_NL_coarse_new_9x9.json` ; `d1/d1_results.json`, `d1/d1_states.npz`, `d1/projwfc_inputs/*/projwfc.in` ; `d3/d3_results.json`, `d3/d3_curves.npz` ; `d4/d4_results.json`, `d4/d4_spectra.npz` ; `d5/d5_results.json`, `d5/d5_lines.npz` ; `d6/d6_results.json`, `d6/d6_curves.npz` ; `r/r_results.json`. Hors dépôt : `cache/` (3.4G) et `prep/M_coarse_*.npy`.

---

- But : savoir où la DFT place les états de la lacune (π et σ) de la 9×9 et quelle étape de la chaîne
  M → Wannier → matrice T s'en écarte. Le facteur α de D3 est un diagnostic, jamais un paramètre.
- Prompt d'origine : R4 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) D6 → D1 → D4 → D2 → D3 → D5 → R → rapport → STOP.
- Statut : **TEST** (post-traitement seulement ; aucun calcul pw.x ; données de production en lecture seule).
- Dates : phase 0 le 2026-09-25 (matin) ; GO reçu le 2026-09-25 ; jobs J0–J6 exécutés et rapport terminé le 2026-09-25 (STOP).
- Répertoire de travail (celui-ci, hors dépôt, règle 5 de CLAUDE.md) : `graphene/qe/defects/R4_quasi_lie/`, à côté de
  `super_cell/9x9/` (données du ch. 4) et de `super_cell_relaxed/9x9/` (R1) qu'il prolonge. Copie versionnée prévue :
  `graphene-raman/campagnes/R/R4_quasi_lie/` (rapport `R4_rapport.md`, pilote, tables, figures ; jamais les npz > 5 Mo ni les `slurm-*`).
- Contenu : `R4_rapport.md` (rapport complet : phase 0, méthodes, D6, D1, D4, D2, D3, D5, R), pilote unique `r4_driver.py`
  (sous-commandes prep, d6, d1, d5, r, d4, g0, d3, tables, figs) + `submit_r4.sh`, résultats `prep/ d6/ d1/ d4/ d3/ d5/ r/` (json, npz),
  `sections/` (rédaction + `tables.md` généré), `fig/`, caches `cache/` (3.4G, non versionnés), `slurm-*`, `JOBID`, `failed_runs/`.
- Jobs : J0 21796633/4 (pp.x, dans R1), J1 21796852, J2 21796853, J3 21796854, J4 21797745 → 21797985 (J6, avec d1), J5a 21797930, J5b 21797931.
- Aucune suppression ; tout nettoyage passe par un manifeste et un GO séparé. Git en lecture seulement.
