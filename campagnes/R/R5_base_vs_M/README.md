# R5 — Base ou M ? Diagnostic de la marche (a1) de R4 ; vérité DFT selon la taille (préparation article, hors mémoire)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 5 figures de `fig/` (hors mémoire), `submit_r5.sh`, `submit_r5_m128.sh` et `r5_deltav_pw.py` (promu dans `src/…/defects/deltav_pw.py`). `r5_driver.py`, `r5_sc_projection.py`, `r5_basis_diagnostics.py`, `r5_alignment_ext.py` sont **gardés uniquement parce que `../R7_tailles_3m/r7_driver.py` les importe** pour deux figures du mémoire ; ils partiront avec lui quand ces figures seront refaites. Gardés : rapport, tables et json/npz de `a/ b/ c/`, inputs du nscf 128 bandes. Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R5_base_vs_M/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

> **Mise au style du 2026-09-30 (étape F0)** : `r5_driver.py`, `r5_sc_projection.py`, `r5_basis_diagnostics.py`, `r5_alignment_ext.py`, gardés seulement parce que R7 les importait, sont retirés (les copies promues dans `src/` restent). Les fonctions de tracé ont été extraites, à contenu identique, dans `scripts/fig/make_figures_controles.py`, qui lit les json/npz de cette campagne et écrit `figures/electron_defect/fig_*` ; rien d'autre ne change. Dernier commit qui contient ces fichiers : `d9a3588`.

## Fiche R5 — Base ou M ? Origine de la marche (a1) de R4 et vérité DFT selon la taille
- **Question** : l'écart entre H_p + M en ondes planes (16 bandes × 81 k : π à −1,350 eV, pas de doublet σ) et QE (−0,737 et +0,101 eV) vient-il de la troncature de la base repliée ou de la matrice M elle-même ? Et où la DFT place-t-elle l'état quasi-lié pour N = 5…12 ?
- **Pourquoi** : R4 (D4) a situé l'essentiel de l'écart à la marche QE → (a1), c'est-à-dire avant Wannier et la matrice T.
- **Méthode** : post-traitement, plus un seul calcul QE (le nscf 128 bandes, lancé sur GO séparé).
  - A (données existantes) : recouvrements c_nk = ⟨nk|Ψ_QE⟩, poids hors base δ, quotient de Rayleigh et résidu (A.1) ; ΔV^L et ΔV^NL appliqués directement à des états purs, comparés à M/81 (A.2) ; entrelacement (A.3) ; w₂ corrigés (A.4) ; moyenne de M^L soustraite (A.5).
  - C : protocole D1 de R4 pour N = 5…12, ajustements ε_∞ + a/N^p par famille (N = 3m / N ≠ 3m).
  - B : nscf 128 bandes, M à 128 bandes, puis (a1) pour n ∈ {16…128}.
  Jobs : J1 A 21811069, J2 C 21811223, J3 nscf 21808943, J4 M 128 bandes 21810702, J5 B 21810855.
- **Résultat** :
  - A.2 : ⟨nk|ΔV^L|n′k′⟩ calculé directement vaut 81,0000 × M^L/81 sur toutes les paires ; pour M^NL le rapport est 1,0000 (écart 3,4e-15 eV). Dans `M_ed_9x9` (juin) comme dans `M_dense_9x9` (septembre), M^L est en norme super-cellule et M^NL en norme maille, alors que les sidecars déclarent `unit_cell` pour les deux.
  - A.1 : la base n'est pas en cause. δ ≤ 0,0009 pour les 28 états π (base 16) et 0,036 / 0,024 pour la paire σ. Pour l'état π 320, R − ε_QE = −306 meV en convention de production, contre +29 meV avec M^L × 81 (poids 0,999 sur l'état à −0,727).
  - A.4 : avec M^L × 81 (variante diagnostique), π à −0,727 contre −0,737 (+10 meV), +0,263 contre +0,269, −1,775 contre −1,759 eV. La paire σ reste loin : +1,48 (16 bandes) / +1,07 eV (20 bandes) contre +0,101. Après correction du bug d'index, w₂(a2) = w₂(a1) à 2,7e-4.
  - B, à n = 128 : en convention de production, π à −1,431 eV, σ à −2,940 eV, et aucun σ à moins de 0,3 eV de +0,101. Avec M^L × 81 : π à −0,760 (−23 meV), σ à +0,116 (+15 meV), et un σ à moins de 0,3 eV dès n = 48.
  - C : le doublet σ a w₂ = 0,702–0,713 à toutes les tailles. Famille 3m : π à −1,065 (6), −0,737 (9) et −0,551 eV (12) ; l'ajustement en 1/N donne ε_∞ = −0,046 eV (rms 0,007). Famille N ≠ 3m : π entre −0,328 et −0,495 eV, dans le gap à Γ.
- **Interprétation / décision** : le rapport ne tranche pas explicitement (« rapporté sans conclure »). Mais la B.3 montre que δ tend vers 0 avec n alors que l'écart reste dans H_in : il vient donc de M (le facteur 81 entre M^L et M^NL), pas de la base. Décision : lancer la campagne R6 de correction, avec M2 = N_cells·M^L + M^NL dans `results/M2/`, `results/M/` gelé, et une porte A.2 à ≤ 1e-6 eV. R5 n'a touché à aucun fichier de production ; il a ajouté l'erratum à `R4_rapport.md`. Sort du `.save` et des M à 128 bandes : non consigné.
- **Où sont les données** : `a/A_tables.md`, `a/a_results.json`, `a/a_overlaps.npz`, `a/a_qe_all.npz`, `a/a_spectra.npz` ; `b/B_tables.md`, `b/b1_results.json`, `b/b1_parity.npz`, `b/b2_results.json`, `b/b3_results.json`, `b/b3_overlaps.npz`, `b/M_{L,NL,ed}_9x9_nb128.json` (sidecars), `b/nscf_nb128.in`, `b/nscf_nb128.diff`, `b/nscf_nb128.out`, `b/submit_nscf_nb128.sh` ; `c/C_tables.md`, `c/c_results.json` ; `r5_log.txt`. Hors dépôt : `b/M_*_9x9_nb128.npy` (3 × 1,72 Go) et le `.save` 128 bandes `qe_tmp/R5_uc9x9_nb128/` (1,59 Go).

---

- But : savoir si l'écart entre H_p + M en ondes planes (16 bandes × 81 k : état π à −1,350 eV, pas de doublet σ) et QE
  (−0,737 et +0,101 eV) vient de la troncature de la base repliée (bandes) ou de la matrice M elle-même, et où la DFT
  place l'état quasi-lié quand la taille de super-cellule varie (N = 5…12).
- Prompt d'origine : R5 (Greg, 2026-09-25). Ordre : phase 0 → STOP → (GO) A → C → (GO nscf) B → rapport → STOP.
- Statut : **TEST** (post-traitement ; un seul calcul QE autorisé, le nscf 128 bandes de B, sur GO séparé, dans un outdir neuf ;
  données de production en lecture seule ; routines de production intouchées ; git en lecture seulement).
- Dates : phase 0 le 2026-09-25 ; GO reçu le 2026-09-25 ; jobs J1–J5 exécutés et rapport terminé le 2026-09-25 (STOP).
- Répertoire de travail (celui-ci, hors dépôt, règle 5 de CLAUDE.md) : `graphene/qe/defects/R5_base_vs_M/`, à côté de
  `R4_quasi_lie/` qu'il prolonge. Copie versionnée prévue : `graphene-raman/campagnes/R/R5_base_vs_M/` (rapport, pilote, tables,
  figures ; jamais les `.npy` de M, les `slurm-*`, `JOBID`, ni les fichiers > 5 Mo).
- Contenu : `R5_rapport.md` (phase 0, méthodes, A, C, B, fichiers, manifestes), modules de campagne `r5_sc_projection.py`,
  `r5_deltav_pw.py`, `r5_basis_diagnostics.py`, `r5_alignment_ext.py` (Greg a refusé l'écriture dans `src/`), pilote `r5_driver.py`,
  `submit_r5.sh`, `submit_r5_m128.sh`, `b/nscf_nb128.in` + `.diff` + `submit_nscf_nb128.sh` + `nscf_nb128.out`, résultats `a/ c/ b/`
  (json, npz, tables ; `b/M_*_nb128.npy` 3 × 1,72 Go non versionnés), `fig/`, `r5_log.txt`, `JOBID`, `slurm-*`.
- Jobs : J3 nscf 21808943 ; J4 M 128 bandes 21810702 ; J1 A 21811069 (tables par `atables`) ; J2 C 21811223 ; J5 B 21810855.
- Constat principal (A.2) : M^L de production en norme super-cellule, M^NL en norme maille (facteur 81) ; production non touchée.
- Aucune suppression ; tout nettoyage passe par un manifeste et un GO séparé.
