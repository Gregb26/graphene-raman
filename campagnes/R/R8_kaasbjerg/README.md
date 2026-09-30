# R8 — Reproduction des figures 13 et 14 de Kaasbjerg, PRB 101, 045433 (2020)

> **Ménage du 2026-09-30 (décision de Greg : on garde les résultats, rapports, tables et ce qui régénère les figures du mémoire ; les pilotes de diagnostic et les figures hors mémoire sont retirés).** Retirés de cette copie : les 3 figures hors mémoire (`7a_controle`, `dos_c_rho0_300`, `superposition_rho0_300` ; leurs données restent dans `out/`) et `ref/kaasbjerg_2020_prb101_045433.pdf` (article sous droits d'auteur). `r8_driver.py` et `submit_r8.sh` sont gardés : ils produisent les cinq figures du mémoire. Gardés : rapport, `out/` (json, npz < 5 Mo), `fig/` (5 figures). Tout ce qui est retiré reste dans l'historique git ; dernier commit qui le contient : `32efd9e` (ex. `git show 32efd9e:article/R8_kaasbjerg/<chemin> > <fichier>`). La fiche ci-dessous résume la campagne ; le rapport reste la référence.

## Fiche R8 — Reproduction des Fig. 13 et 14 de Kaasbjerg (PRB 101, 045433)
- **Question** : avec M2 de production et la matrice t locale (bloc π), au protocole de Kaasbjerg (Σ_k = c_i T̄_kk, c_i par maille), retrouve-t-on sa DOS moyennée sur le désordre (Fig. 13, c_i = 0,1 % et 1 %) et sa fonction spectrale Γ–K–M (Fig. 14, c_i = 1 %) pour des lacunes sur un seul sous-réseau (V_A) ?
- **Pourquoi** : non consigné explicitement (le rapport donne des « chiffres bruts, sans interprétation »). Le README présente R8 comme la reproduction d'une référence publiée avec la chaîne de production. Les figures servent au mémoire (étapes 1–5 et 7) et à l'article (étape 6).
- **Méthode** : extraction vectorielle des figures du PDF (7a, grille de 10 meV). T̄_k est tiré de t (R_cut 3, N_k^int 900, η_t 20 meV), puis G_k = [(ε + iη_G) − H_k − c_i T̄_k]⁻¹. On calcule la DOS en 9×9 (lacune A ; 300² et 600², η_G 50 meV) et A_k = −2 Im Tr G_k sur 601 points (η_G 25 meV), en trois variantes : tel quel, aligné plateau, eta_unique. Suivent les sensibilités (étape 5 ; N_k^int, R_cut, 12×12, η_G), Σ^eff à K (étape 4) et le modèle de Dirac (étape 6), ces deux derniers sur le nœud de connexion.
  - Jobs : 22024568 `r8g0`, 22026235 `r8gate`, 22026236 `r8dos`, 22026237 `r8spec`, 22026238 `r8sens`, 22026239 `r8fig` (tous COMPLETED, 2026-09-29).
- **Résultat** :
  - Portes P1–P4 : PASS (P1 ≤ 8,0e-16 ; P2 5,3e-15 ; P3 4,6e-14 ; P4 12/12 et 2/2 identiques).
  - Maximum de ρ − ρ₀ (9×9, 600², tel quel) : −0,1625 eV (0,1 %) et −0,1650 eV (1 %) ; hauteurs 0,00283 / 0,02676 états/eV/maille/spin ; largeurs 0,191 / 0,201 eV. Aligné plateau : −0,1475 / −0,1500 eV. Kaasbjerg : −0,200 / −0,210 eV, 0,00126 / 0,0307 eV⁻¹, 0,163 / 0,306 eV.
  - Superposition (7b) : écart de position (nous − Kaasbjerg) de +0,0375 à +0,0600 eV. Sur [−1, +1] eV, rms ×1 de 0,0333 à 0,0351 contre ×2 de 0,0040 à 0,0081 (ρ₀ 300²). Avec ρ₀ 1 200², le rms ×2 tel quel vaut 0,00296 (0,1 %) et 0,00643 (1 %).
  - A_K à 1 % (tel quel) : maxima à +0,0076 et +0,1252 eV, écart 0,1176 eV (aligné 0,1182), contre +0,010 / +0,110 eV (écart 0,100 eV) chez Kaasbjerg. Sauts de la branche inférieure : −0,4357 → −0,1325 (côté Γ) et −0,1330 → −0,4300 eV (côté M), contre ]−0,45 ; −0,16[ et ]−0,50 ; −0,21[. À 0,1 % : un seul maximum (+0,0091 eV).
  - Sensibilités (1 %, position du maximum) : de −0,1725 à −0,1500 eV tel quel, de −0,1525 à −0,1325 eV aligné. Étape 4 : solutions à K en +0,0074 ; +0,0821 ; +0,1215 eV (tel quel). Étape 6 : pôle sans Λ à −0,1870 eV (Ṽ 76,96 eV Å²), Λ reproduisant 4,8 / 21,5 eV.
- **Interprétation / décision** : les écarts à Kaasbjerg sont rapportés mais non investigués (amendement du 28 sept.). Le rms ×2 est 5 à 8 fois plus petit que le rms ×1. Greg a pourtant décidé le 2026-09-29 de garder les DOS en états/eV/maille/par spin ; le ×2 reste à titre d'information. Les figures de DOS sont refaites avec ρ₀ 1 200² (η_G 15 meV). La portée du mémoire (2026-09-30) couvre les étapes 1–5 et 7 (`dos_c`, `spectral_GKM`, `sigma_K`, `sensibilites`, `superposition`) ; l'étape 6 reste pour l'article. Statut TEST, rien dans `results/`.
- **Où sont les données** (relatif à `campagnes/R/R8_kaasbjerg/`) :
  - Rapport : `R8_rapport.md` (tables).
  - Portes : `out/gate/gate_results.json`, `res_<taille>_<variante>_nk<N>.npz`.
  - DOS : `out/dos/dos_results.json`, `dos_9x9.npz`.
  - Spectre : `out/spec/spec_results.json` (`spectral_GKM_9x9.npz`, 14 Mo, non copié).
  - Sensibilités : `out/sens/sens_results.json`, `sens_dos.npz`.
  - Superposition : `out/fig/fig_results.json` (+ `fig_results_rho0_300.json`).
  - Extraction 7a : `out/7a/7a_results.json`, `fig14_K_profiles.npz`, csv.
  - Étapes 4 et 6 : `out/sigeff/sigeff_results.json`, `sigma_K_9x9.npz` ; `out/dirac/dirac_results.json`, `dirac_9x9.npz`.
  - Préparation : `out/prep/prep_9x9.json`, `prep_12x12.json`.
  - Figures du mémoire : `fig/` (`r8_driver.py fig`, `sigeff` pour `sigma_K`).

---

- But : DOS (c_i = 0,1 % et 1 %, Fig. 13) et fonction spectrale Γ–K–M (c_i = 1 %, Fig. 14) du graphène avec lacunes sur un seul
  sous-réseau, avec M2 de production et la matrice t locale (bloc π), protocole de Kaasbjerg : Σ_k = c_i T̄_kk, c_i par MAILLE.
- Prompt : R8 (Greg, 2026-09-27) + ajout du 27 sept. + amendement du 28 sept. après R9 (N_k^int 900, variantes tel quel / aligné plateau, 2 bis, 7a avant les DOS).
- Statut : **TEST** (post-traitement seul ; résultats consignés dans `R8_rapport.md`, `out/` et `fig/` ; rien dans `results/`).
- Date : 2026-09-27 (phase 0) ; mise à jour du 2026-09-28 (amendement après R9).
- Lecture seule : `results/M2/`, `wannier/27x27`, `wannier/24x24`, json de R9 (jamais `R9_controles/cache/`) ; g₀ recalculé. Aucun calcul QE. Git en lecture.
- Copie versionnée : `campagnes/R/R8_kaasbjerg/` (rapport, pilote, tables, figures ; jamais npz > 5 Mo ni slurm).
- État : 1er GO (étapes 1, 2, 2 bis, 3, 5, 7b) et GO article (étapes 4, 6) exécutés le 2026-09-29 ; résultats dans `R8_rapport.md`, `out/`, `fig/` ; STOP.
- Portée du mémoire (2026-09-30) : étapes 1–5 et 7 ; figures `dos_c`, `spectral_GKM`, `sigma_K`, `sensibilites`, `superposition` ; l'étape 6 (Dirac) reste pour l'article.
