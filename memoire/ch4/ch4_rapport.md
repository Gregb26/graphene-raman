# Chiffres du ch. 4 — rapport de session (2026-09-30)

Prompt de Greg (2026-09-30) : table v1 → final, régions d'alignement, NOTES_TGAMMA au final, nombres de `défauts.tex`. Pas de campagne R,
aucun job. Ordre : phase 0 → STOP → GO (parties 1, 2, 3, 4, décisions D1–D5 par défaut, précisions 1.b, 1.g, 3) → rapport → STOP.
Dépôt `graphene-raman`, HEAD `66602ed` (arbre propre au départ). Emplacement unique `memoire/ch4/`, un seul pilote `ch4_chiffres.py`.
Chiffres bruts, sans interprétation. Rien dans `src/`, `scripts/`, `config/`, `results/`, `article/`, `figures/` ; seul fichier modifié hors
`memoire/ch4/` : `NOTES_TGAMMA.md`. Rien n'est commité par Code.

## Phase 0 (rien n'écrit, dans le chat)

- 0.1 : R6 `table_v1_v2.md` 126 lignes, R10 `table_v2_plateau.md` 192 lignes, aucun doublon ; règle d'appariement = grandeur normalisée + réécritures
  d'étiquettes (liste dans `table_v1_final.md`, section « Règles ») ; 78 appariées à l'identique, 117 après réécritures.
- 0.2 : plan des modifications de `NOTES_TGAMMA.md` (l.10–16 results_dir, l.65 alignement, §2, C9–C16, ligne N_k^int, §4 en-tête, §6e, §8, E_res).
- 0.3 : lignes de `CLAUDE.md` fausses après R10 (l.243, l.255, l.315, table de classement) → `CLAUDE_md.diff` (proposé, non appliqué).
- 0.4 : porte PASS (13 C_N redonnés au bit).
- 0.5 : sources localisées pour a–l et pour les grandeurs retirées ; constats : « (c-3) avec Kumagai–Oba » n'existe pas (R10 A.3 aligné = site unique) ;
  « max|M| hors k = k′ » = R9 R.2, pas R10 C.3 ; R8 a trois variantes ; `défauts.tex` est dans le dépôt (`memoire/défauts.tex`, commit 66602ed).
- Décisions (défaut appliqué au GO) : D1 npz finaux pour les lignes v1 seulement (4 lignes : 48, 51, 77, 80 de R6 ; Friedel par bloc depuis
  `resonance_criteria_9x9.npz`) ; D2 tab:rcut_M : deux valeurs, différence sur la première ; D3 les 8 lignes C14 ±25 en « retiré » ; D4 a, c depuis
  les `scf.in` ; D5 `défauts.tex` lu dans le dépôt, lecture seule.

## Porte de la partie 2 (`check`)

Règle « distance vraie ≥ 0,75 r_max, moyenne des `shift10` » sur `article/R10_plateau/a/profiles_<S>.npz` → les 13 C_N de `config/production.json`
avec un écart de **0,0 eV** (seuil 1e-9) ; n et rms identiques à `a1_results.json` (rms = écart-type de population). md5 `a1_results.json`
3d98bbd244eba8e1c1c645dcf675432b = `alignment.source_md5`. md5 des sources d'en-tête : `table_v1_v2.md` f8972c4de041dbd79e70b2899aa4f941 ;
`table_v2_plateau.md` 807bb352fcb989f01d94481955fe4c10 ; `results/M2_plateau/MD5SUMS_2026-09-30.txt` 9c29daedb9c417ba795bb39c0ed71bee
(`md5sum -c` : 33/33 OK). Réseau : a = |a₁|/N = 2,465852 Å pour les 13 `scf.in` (CELL_PARAMETERS bohr) = `wannier.wout` a_1 ;
c = 30 bohr = 15,875316 Å = `wannier.wout` a_3 = `c6_results.json` c_A. r_max = a_sc/√3 (9×9 : 12,8129 Å) ; N·a/2 = 11,0963 Å = (√3/2) r_max.

## Partie 1 — `table_v1_final.md`

Comptes : appariées 117 (dont 9 lignes critère/Friedel de R6 appariées aux 3 lignes combinées de R10), v1 seulement 9, final seulement 81, total
**207 lignes**. Statuts mécaniques : remplacé 85, nouveau 87, retiré (décision) 24 (E_res × 14, C14 67 meV / ±25 × 9, det −2,530 × 2 — les lignes
R6 l.89–90 —, convention intensive × 2 ; le statut mécanique est conservé après « ; »), inchangé 6, sans équivalent final 5 (Born/T min et max,
Γ_T à c = 0,1 %, ħ/Γ : pas de scalaire dans les npz finaux ; non-régression nbnd 16 → 20 : non rejouable).
Écarts entre les deux tables sources (rapportés dans le fichier, non corrigés) : 31, dont 20 différences de colonne « fichier source » (R6
« level1_summary.csv / level2_summary.csv » vs R10 « level1_summary.csv » ; R10 ajoute « c/c2_results.json » aux lignes tab:rcut_M), 3 différences
de forme de la colonne v2 (R6 l.87 « 3132.872 » vs R10 « 3132.87 » ; l.88 « 1.638 » vs « +1.638 » ; l.90/93/96 λ_min complexe vs |λ|), et les
cellules Friedel de R10 réduites à la première valeur (R6 donne bande / Lloyd / ±3 eV ; le final complet est lu dans le npz, D1).
Chiffres nouveaux a–l : chaque section nomme la source, l'alignement et l'E_D. Précisions du GO appliquées : b — ligne « (c-3) aligné » = modèle
de R10 A.3 avec C_9 à site unique (−24,65 meV), étiquetée ainsi, aucune valeur Kumagai–Oba ; QE donné deux fois (R6 2.1 site unique, R10 A.2
Kumagai–Oba) ; E_D et décalage rigide (`rigid_shift_fit`) par marche ; g — C.3 de R10 (« convention intensive », écart relatif de max|M| des M
denses complets par famille : 3m 2,257e-2, non-3m 2,057e-2) et R9 R.2 (`offdiag_intensive`, trois variantes) rapportés sous leur vrai nom ;
h — R8 `eta_unique` signalée dans la table, non utilisée.

## Partie 2 — `alignement_regions.md`, `.npz`, `.pdf/.png`

13 tailles, régions (A) ≥ 0,75 r_max, (B) ≥ N·a/2 (sites à < 1e-6 Å du cercle inclus et comptés : 3 pour N pair, 0 pour N impair), (C) ≥ min(N·a/2, c/2)
(= (B) pour 5×5 et 6×6 ; rayon c/2 = 7,9377 Å pour N ≥ 7), (D) site unique publié et site vraiment le plus loin. Chiffres bruts (meV) :
9×9 (A) 53 sites, −25,1437, rms 9,0509 ; (B) 20, −31,6899, rms 10,9424 ; (C) 83, −24,9918, rms 9,5205 ; (D) −24,6514 / −45,8805 ;
(B) − (A) = −6,5462 = 0,72 rms(A). Sur les 13 tailles, |(B) − (A)|/rms(A) va de 0,095 (10×10) à 0,866 (27×27). Figure au style du mémoire
(`figures/memoire.mplstyle`, palette `scripts/_palette.py`), moyenne ± rms contre N, familles 3m (cercles) / non-3m (carrés) ; non installée dans `figures/`.

## Partie 3 — `NOTES_TGAMMA.md` (racine) modifié en place ; diff complet ci-dessous et dans `NOTES_TGAMMA.diff`

- En-tête : paragraphe « Mise à jour R10 (2026-09-30) — état final » (results_dir, matrices_dir, nomenclature, règle « valeur en tête = final,
  [v1 : …] », renvoi à `memoire/ch4/`) ; phrase du paragraphe R6 mise au passé.
- §1 : ligne « Super-cellule de référence » (matrices brutes, alignement en base de Wannier) ; ligne « Alignement du potentiel » **remplacée**
  (définition de Kumagai–Oba, région, sphères, rms = incertitude, symbole ΔV_PA^(N), les 13 C_N de la config, application M_W(R,R) − C_N·𝕀₅ sur les N²
  mailles de la boîte, approximation (i) exacte à 0,06 meV à 9×9 (R9 A.2, `R9_rapport.md` l.402), images de Wigner-Seitz, et la phrase
  « la production n'a jamais soustrait la moyenne de ΔV (`subtract_mean=False`, `scripts/compute_M.py:118`) ») ; ligne « Coût » (R10 C).
- §2 : titre et les 10 lignes en valeurs finales, sources `results/M2_plateau/…` ; E_res marqué « retiré du ch. 4 (R10 B.1) », valeurs gardées ;
  grandeurs non publiées par R10 (Born/T min/max, écart ponctuel de Friedel, Γ_T à c = 0,1 % scalaires) dites telles quelles.
- §3 : C9 (E_res retiré), C10, C11, C12, C13, C14 (**nouvelle définition** : ±9,05 meV autour du V_loc aligné), C15, C16, ligne N_k^int.
  Pourcentages marqués « arithmétique », vérifiés par script : C11 0,37 / 0,18 / 0,43 / 1,37 % ; C12 3,60 / 1,83 / 8,11 % ; N_k^int −0,09 / −0,03 / +4,82 %.
- §4 : en-tête (renvois ; point 4 clos). §5 inchangé. §6 : nouvelle sous-section 6e (final, `nkint_check_9x9.csv`, deux tables, écarts arithmétiques).
- §7 inchangé. Nouveau §8 « R10 » : décision, code, portes, valeur non alignée → finale pour Re Σ (R_cut 0…4) et Kaasbjerg Fig. 3 (D.1 K′ en tête,
  K, D.2 ½ Tr), Fig. 13–14 (R8, deux variantes, valeurs lues), constats ouverts (C.6, nscf 27×27 7 avertissements, modèle replié aligné à site unique,
  points de R6 encore ouverts).
- `notes` : diff de 180 lignes ; contrôle de 12 valeurs finales (csv/npz → texte) : 12/12 retrouvées.
- Le mot « Lu » n'apparaissait pas dans le fichier ; aucune valeur non alignée n'y figure hors §8 (et le §6d historique, gardé).

## Partie 4 — `defauts_nombres.md` (`memoire/défauts.tex`, 749 lignes, md5 dans l'en-tête, lecture seule)

553 nombres relevés (hors commentaires, \label, \ref, \eqref, \cite, \includegraphics, années, indices/exposants) ; 399 entiers courts (< 100, sans
décimale) non comparés ; **87 appariés, 67 non appariés** ; 59 lignes de la table principale (sur 207) jamais citées. Le mémoire cite les valeurs v1
(ex. tab:L_NL 0.421 / 3.100 / 6.5 → lignes v1 ; médianes 2524, 2509, 2597 → v1 arrondies ; « quasi-annulation à ∼ −2,5 eV », « partie imaginaire
minimale à −1,3 eV » → lignes retirées/remplacées) ; les non appariés sont surtout des paramètres (300², 100 Ry, 30 bohr, fenêtres Wannier) et les
médianes de niveau 1 à R_cut 0–2 des tailles autres que 9×9 (absentes des tables sources). Un nombre du texte peut s'apparier à plusieurs lignes par
coïncidence (0.5, 0.1, 3.10) : correspondances mécaniques, pas des identifications.

## Erreurs et incohérences d'archives trouvées (rapportées, non corrigées)

1. `table_v2_plateau.md` : colonne « fichier source » des lignes de niveau 1 sans `level2_summary.csv` (R6 l'avait) ; v2 de R6 l.87 (3132.872) et
   l.88 (« 1.638 ») recopiés avec une décimale de moins / un signe « + » ; cellules critère/Friedel combinées (λ en module, Friedel réduit à une valeur).
2. `NOTES_TGAMMA.md` l.6 : « Chemins relatifs à `ab-initio-defects/` » (ancien nom, déjà signalé l.16) — laissé.
3. Prompt : « (c-3) avec l'alignement de Kumagai–Oba = modèle aligné N = 9 de R10 A.3 » (b) alors que ce modèle utilise C_9 à site unique (e) ;
   « max|M| hors k = k′ (R10 C.3) » (g) alors que C.3 est la convention intensive sur M complets et R9 R.2 le test hors diagonale — tranché par les
   précisions du GO.
4. `défauts.tex` dans le dépôt (commit 66602ed) alors que le prompt le voulait hors dépôt (D5 : lu là, lecture seule).
5. R8 `spec_results.json` : `K_gap_eV = None` à c_i = 0,1 % (un seul maximum) ; la table l'écrit « — (un seul maximum) ».

## Fichiers et état du dépôt

- Créés : `memoire/ch4/{ch4_chiffres.py, table_v1_final.md, alignement_regions.md, alignement_regions.npz, alignement_regions.pdf,
  alignement_regions.png, NOTES_TGAMMA.diff, defauts_nombres.md, CLAUDE_md.diff, README.md, ch4_rapport.md}` (non suivis).
- Modifié : `NOTES_TGAMMA.md` (racine). Proposé, non appliqué : `CLAUDE_md.diff`.
- Non touchés : `src/`, `scripts/`, `config/`, `results/`, `article/`, `figures/`, `memoire/défauts.tex`. Aucune suppression.

**STOP — Greg relit et commit.**

## Annexe A — diff complet de `NOTES_TGAMMA.md` (`git diff`)

```diff
diff --git a/NOTES_TGAMMA.md b/NOTES_TGAMMA.md
index 77b0550..c678547 100644
--- a/NOTES_TGAMMA.md
+++ b/NOTES_TGAMMA.md
@@ -10,11 +10,22 @@ code fait réellement. **Les jugements sont laissés au lecteur ; ce document co
 **Mise à jour R6 (2026-09-27).** La campagne R6 a corrigé la normalisation de la partie locale de M : M^L de production était en
 norme super-cellule, trop petit d'un facteur N_cells = N² (constat R5-A.2). Les M corrigés (v2 : M2 = N_cells·M^L + M^NL, sidecar
 `M_normalization = v2`) et tous les produits sont dans `results/M2/` (`config/production.json` : `results_dir`) ; `results/M/` est gelé
-(README). Les §1, §2, §3 et §6 donnent maintenant les valeurs v2, les valeurs v1 entre crochets « [v1 : …] » ; les §4 et §5 sont l'état
+(README). Après R6, les §1, §2, §3 et §6 donnaient les valeurs v2 ; les §4 et §5 sont l'état
 du 2026-09-18 (v1), inchangés ; le §7 résume R6. Rapport : `graphene/qe/defects/R6_production_corrigee/R6_rapport.md` (copie
 `article/R6_production_corrigee/`) ; table complète ancien → nouveau : `article/R6_production_corrigee/etape3/table_v1_v2.md`. Le dépôt
 s'appelle `graphene-raman` depuis le 2026-09-24 (ex `ab-initio-defects`) ; les provenances `fichier:ligne` datent du 2026-09-18.
 
+**Mise à jour R10 (2026-09-30) — état final.** La production du chapitre 4 est `results/M2_plateau/` (`config/production.json` :
+`results_dir`) : matrices M2 inchangées (`results/M2/`, clé `matrices_dir` ; `results/M/` et `results/M2/` gelés, `results_dir_frozen`) et
+**alignement de Kumagai–Oba** du potentiel de défaut (bloc `alignment`, symbole du mémoire ΔV_PA^(N) = C_N ; définition au §1, ligne
+« Alignement du potentiel », et au §8). Nomenclature : « non aligné » (= « tel quel » des archives, `results/M2`), « alignement à site unique »
+(= « Lu »), « alignement de Kumagai–Oba » (= « plateau (i) »), « final » = v2 + Kumagai–Oba. Dans ce document, **la valeur en tête est la
+valeur finale** ; entre crochets, seulement « [v1 : …] » (ce que dit le mémoire actuel) ; la valeur non alignée n'apparaît qu'au §8
+(alignement), pour Re Σ et pour la grandeur de Kaasbjerg (Fig. 3) près de K. Table complète v1 → final, chiffres nouveaux et nomenclature :
+`memoire/ch4/table_v1_final.md` (pilote `memoire/ch4/ch4_chiffres.py`) ; régions d'alignement : `memoire/ch4/alignement_regions.md`.
+Rapport R10 : `graphene/qe/defects/R10_plateau/R10_rapport.md` (copie `article/R10_plateau/`) ; table v2 → final :
+`article/R10_plateau/c/table_v2_plateau.md`.
+
 ## 0. Ce que le script fait, dans l'ordre (chaîne de production, taille de référence 9×9)
 
 Scripts : `scripts/compute_spectral_wannier.py` (carte niveau 1 : R_cut × grille × η), `scripts/resonance_metrics.py`
@@ -46,7 +57,7 @@ sur une grille d'énergie de pas η/ne_per_eta (`resonance_metrics.py:55`) →
 
 | paramètre | valeur | fichier:ligne / preuve |
 |---|---|---|
-| Super-cellule de référence, M | 9×9 (N_cells = 81), lacune sur le sous-réseau A, s_red = (13/27, 13/27) ; M dense zero-padded p = 3, D = 27, forme (20, 729, 20, 729), Hartree, norme `unit_cell` pour M^L et M^NL (`M_normalization = v2`, R6 : M2 = N_cells·M^L + M^NL) | `config/production.json:7` ; `results/M2/M_dense_9x9.json` (sidecar) [v1 : `results/M/M_dense_9x9.json`, M^L en norme super-cellule] |
+| Super-cellule de référence, M | 9×9 (N_cells = 81), lacune sur le sous-réseau A, s_red = (13/27, 13/27) ; M dense zero-padded p = 3, D = 27, forme (20, 729, 20, 729), Hartree, norme `unit_cell` pour M^L et M^NL (`M_normalization = v2`, R6 : M2 = N_cells·M^L + M^NL) ; matrices brutes, l'alignement de Kumagai–Oba est appliqué en base de Wannier (ligne « Alignement du potentiel », §8) | `config/production.json` (`matrices_dir`, `alignment`) ; `results/M2/M_dense_9x9.json` (sidecar) [v1 : `results/M/M_dense_9x9.json`, M^L en norme super-cellule] |
 | Cellule primitive dense (nscf) | grille 27×27×1 = 729 k, `nosym` + `noinv` (grille complète), nbnd = 20, `diago_full_acc`, ecutwfc 50 (unités Ha du XML, soit 100 Ry), ecutrho 200, smearing m-v degauss 5,0e-3 (Ha du XML), conv_thr 5e-15, `assume_isolated = 2D` | `scratch/qe_tmp/defect_uc_dense_27/defect_uc_dense_27.save/data-file-schema.xml` (`<nk>`, `<nbnd>`, `<nosym>`, `<noinv>`, `<ecutwfc>`, `<smearing>`, `<assume_isolated>`) ; `config/production.json:16` (nbnd_dense 20) |
 | Wannierisation 27×27 | num_wann = 5 (C1 : sp², p_z ; C2 : p_z), 20 bandes, fenêtre externe [−25, 15] eV, fenêtre gelée [−25, −1,74] eV (= E_D + 2,5 eV), Ω_I = 3,0246, Ω_D = 0,0136, Ω_OD = 0,7285, Ω_tot = 3,7668 Å² ; étalements 0,611 (×3, sp²) et 0,967 Å² (×2, p_z) ; run_id `baa17b88b1b69e51` | `wannier/27x27/wannier.wout` l. 125 (grille), 863 (num_wann), 1012–1013 (fenêtres), 2553–2555 (Ω), « Final State » ; `config/production.json` clé `wannier` ; `wannier/27x27/wannier_manifest.json` (sha256 de tb / u / u_dis) |
 | Maille | a₁ = (2,135490, −1,232926, 0), a₂ = (2,135490, 1,232926, 0) Å → a = 2,4659 Å, \|b\| = 4π/(√3 a) = 2,9423 Å⁻¹ ; c = 15,875 Å | `wannier/27x27/wannier.wout` l. 89–92 |
@@ -62,26 +73,26 @@ sur une grille d'énergie de pas η/ne_per_eta (`resonance_metrics.py:55`) →
 | Grille de ρ₀ | 900×900 (`--rho0-grid`), ρ₀ par cellule et par spin ; `rho0_240` (sur la grille de sortie) aussi stocké | `resonance_metrics.py:18, 80–82` |
 | Concentration | c = 1 % par cellule pour ρ_dis = ρ₀ + c δρ ; c = 0,1 % pour la comparaison Kaasbjerg | `config/production.json` clé `defect_concentration_for_dos` ; `resonance_metrics.py:83` ; `resonance_criteria.py:19` |
 | Born | t_B = V + V g₀ V (deux premiers termes) → Γ_Born = −2 Im ⟨V g₀ V⟩ (le terme du 1er ordre est réel) | `resonance_metrics.py:60` |
-| Alignement du potentiel | décalage G ≈ 0 = moyenne diagonale de M^L = **5 427,3 meV** (v2, norme unit_cell) [v1 : 67,0 meV] ; Γ_T recalculé avec M − ⟨M^L⟩·1 ; R6 : variante V_loc + C·1, C = ±25 meV (C14) | `resonance_metrics.py:25–26, 42, 61` ; log `results/M2/logs/post_res9x9_21862372.out` « [align] » |
+| Alignement du potentiel (ΔV_PA^(N) = C_N) | **Alignement de Kumagai–Oba** (Kumagai et Oba, PRB 89, 195205, 2014 ; région élargie, `memoire/ch4/alignement_regions.md`) : C_N = moyenne des potentiels de site ΔV (sphères de 1,0 Å, `alignment.atom_sphere_shifts`) sur les atomes à distance vraie ≥ 0,75 r_max de la lacune (image minimale) ; rms de ces potentiels = incertitude ; **9×9 : C_9 = −25,1437 meV, rms 9,05 meV, 53 atomes** (13 tailles : 5×5 −57,5626 ; 6×6 −50,5091 ; 7×7 −26,9143 ; 8×8 −14,4103 ; 9×9 −25,1437 ; 10×10 −9,3417 ; 11×11 −9,5507 ; 12×12 −18,6896 ; 15×15 −16,2185 ; 18×18 −15,0621 ; 21×21 −14,2343 ; 24×24 −13,4815 ; 27×27 −13,0058 meV). Application en base de Wannier : M_W(R, R) − C_N·𝕀₅ sur les N² mailles de la boîte (approximation (i), exacte à 0,06 meV à 9×9, R9 A.2 ; `local_tmatrix.defect_mwr`) ; étiquettes hors grille par images de Wigner-Seitz (`ws_images`, `ws_phase`). **La production n'a jamais soustrait la moyenne de ΔV** (`subtract_mean=False`, `scripts/compute_M.py:118`). C14 = ±rms (±9,05 meV) autour du V_loc aligné (§3) [v1 : décalage G ≈ 0 = moyenne diagonale de M^L, 67,0 meV, Γ_T recalculé avec M − ⟨M^L⟩·1] | `config/production.json` bloc `alignment` (`C_N_eV`, `source` = `article/R10_plateau/a/a1_results.json`, md5 3d98bbd2…) ; `defects/alignment.py` (`atom_sphere_shifts`, `true_min_image_dist`) ; `local_tmatrix.py` (`defect_mwr`) ; `wannier_interpolation.py` (`ws_images`, `ws_phase`, `Mwr_to_Mwk(ws=)`) ; ligne « [align] » des journaux `results/M2_plateau/logs/` ; R10 A.1 |
 | δρ | δρ(ε) = (1/π) Im Tr[t(ε) dg₀/dε] par défaut et par spin ; contrôle Lloyd δρ = −(1/π) d/dε Im ln det[1 − g₀V] | `resonance_metrics.py:62` ; `resonance_criteria.py:57–66` |
 | T̄ en K | paire π (bandes 3, 4 de 5), trace/2 de ⟨n K\| t \|n′ K⟩ | `resonance_metrics.py:86–90` |
 | Test d'or | η = 0,10 eV, k_int = k_out = grille grossière, R_local = boîte MP-duale complète (exact) ; seuil rel < 1e-8 et Γ ≥ −1e-8 | `scripts/test_local_tmatrix_real.py:20, 53, 59, 61` |
-| Coût (16 cœurs) | v2 (R6) : carte niveau 1 9×9 (45 combinaisons, R_cut 0–4) 53 min ; `resonance_metrics` 7 min 49 ; `resonance_criteria` (trois blocs, règle de somme par bloc) 6 h 27 ; test d'or 5×5 dense 1 h 52 [v1 : carte (36 combinaisons) 1 h 15 ; 7 min 39 ; 1 h 34 ; test ne_per_eta 3 min 35] | `sacct` jobs 21857276, 21862372, 21852238 [v1 : 20294202, 20421370, 20414406, 20684497] |
+| Coût (16 cœurs) | final (R10 C) : cartes niveau 1 `specwd_<S>` 1 h 17 / 2 h 01 / 2 h 47 / 2 h 38 / **1 h 01** / 1 h 11 (5, 6, 7, 8, 9, 12 ; R_cut 0–4) ; `post_res9` (métriques + critères 9×9) 39 min 45 ; `post_res6` 6 h 25 ; `post_res12` 2 h 20 ; `nkint` 11 min 33 ; `resigma` 2 min 47 + 8 min 14 ; C14 3 min 36 ; test d'or 5×5 dense non rejoué (D12 : recopié de R6, 1 h 52) [v1 : carte (36 combinaisons) 1 h 15 ; 7 min 39 ; 1 h 34 ; test ne_per_eta 3 min 35] | R10_rapport.md, table « Rapport — C » (jobs 22058839–44, 22058850–52, 22058845–47, 22058854) [v1 : 20294202, 20421370, 20414406, 20684497] |
 
-## 2. Valeurs de production (9×9, R_cut 3, 240², η 0,02, N_k^int 300²) — v2 (R6, `results/M2/`)
+## 2. Valeurs de production (9×9, R_cut 3, 240², η 0,02, N_k^int 300²) — final (R10, `results/M2_plateau/`)
 
 | quantité | valeur | source |
 |---|---|---|
-| médiane de \|Γ\|·N_cells sur les 41 266 états (carte niveau 1) | **3 132,60 meV** (Γ par défaut, ×81 = intensif) [v1 : 2 473,55] | `results/M2/level1_summary.csv` l. « 9x9,3,240,0.02 » ; `specwd_9x9_prod.npz` (job 21857276, 2026-09-26) |
-| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 3 740,17 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 177 204,74 meV ; Born/T médian 46,489 (min 1,106, max 211,72) [v1 : 2 343,59 ; 7 405,83 ; 3,302 (0,677, 16,45)] | log `results/M2/logs/post_res9x9_21862372.out` |
-| Re Σ médian, R_cut 3 | −26,65 meV ; \|Re Σ\|/Γ médian 0,125 [v1 : 722,3 ; 0,227] | `results/M2/m_rcut_resigma.csv` (reconstruit par R6 avec les définitions v1, `R6_production_corrigee/r6_m_rcut_resigma.py`) |
-| position de résonance E_res − E_D (argmax des états à ±1,5 eV) | −0,175 eV (R_cut 3, 240², 0,02 ; de même η 0,01 et R_cut 0, 2, 4) ; −0,202 (η 0,05) ; −0,227 (R_cut 1) : valeurs discrètes [v1 : −1,238 ; −1,183 (η 0,01) ; −1,292 (R_cut 2 et 4)] | `level1_summary.csv` ; `m_rcut_resigma.csv` |
-| pics des courbes (ε − E_D) | Γ_T : −0,180 eV ; Γ_Born : +1,695 eV ; Γ_T/ρ₀ : −0,175 eV ; δρ : −0,815 eV ; ρ_dis : +1,710 eV ; \|T̄\| : −0,172 ; −Im T̄ : −0,177 ; min \|Re T̄\| : −2,447 ; Re T̄(E_D) = 11,122 eV, Im T̄(E_D) = −1,796 eV ; zéros de Re T̄ : −2,447, −0,222, +1,638 [v1 : −1,24 ; +1,695 ; −0,015 ; −2,53 ; −2,53 ; −0,905 ; −1,29 ; −2,145 ; 2,521 / −0,091 ; aucun zéro] | `results/M2/resonance_9x9.npz` (`peak_*`, `Tbar_zero_crossings`) |
-| critère det / valeur propre | matrice complète (dim 145) : min de \|det[1 − Vg₀]\|/max = 1,32e-4 à −0,812 eV ; min_i \|λ_i\| = 0,0019 au même point (λ = +0,0000 + 0,0019 i) ; autres minima locaux −0,630 (7,5e-4), puis −0,255 … −0,132 eV (4,7e-4 – 5,7e-4). Bloc σ (dim 87) : même zéro (3,48e-4 ; λ = 0,0019 i). Bloc π (dim 58) : aucun zéro, min \|λ\| = 0,397 à −0,170 eV (\|det\|/max 2,33e-2) [v1 : 2,09e-4 à −2,530 eV ; 0,0108 (λ = +0,0003 + 0,0108 i) ; secondaires −2,19 … −1,97 eV] | log `post_res9x9_21862372.out` ; `results/M2/resonance_criteria_9x9.npz` (clés `_pi`, `_sigma`) |
-| règle de somme de Friedel | ∫δρ sur [−25,28, 8,66] eV (6 790 points, pas 5 meV) = **−1,0007** états (Tr[t g₀′]) et −1,0007 (Lloyd) ; bloc π −0,9980, bloc σ −0,0028 ; dans ±3 eV : +0,669 ; écart ponctuel max entre les deux formules 0,630 états/eV ; cumul aux bords : −3,2545 (−3 eV), −2,5862 (+3 eV) [v1 : −0,0569 ; +1,782 ; 0,614 ; −2,190 / −0,409] | idem |
-| Γ_T à c = 0,1 % sur ±1 eV | min 2,25 meV (+1,00 eV), max 37,84 meV (−0,18 eV), 7,79 meV à E_D ; ħ/Γ à ∓0,3 eV : 24 / 230 fs [v1 : 0,63 (+0,24) ; 5,55 (−1,00) ; 1,27 ; 419 / 1 025 fs] | idem |
-| localité de M_W (9×9 dense) | ‖M_W(0, 0)‖ = 41,879 eV ; p_z–p_z sur le site de la lacune 31,521 eV, p_z de l'autre sous-réseau 0,615 eV ; ‖M_W(R, 0)‖ par distance : \|R\| = a : 1,963 / 1,501 / 1,388 ; √3 a : 0,639 / 0,274 / 0,137 / 0,122 ; 2a : 0,287 / 0,278 eV [v1 : 9,365 ; 6,617 ; 0,044 ; a : 0,658 / 0,460 / 0,437 ; √3 a : 0,165 / 0,079 / 0,035 / 0,029 ; 2a : 0,082 / 0,081 — le texte du 2026-09-18 rangeait 0,460 dans la 2ᵉ couronne] | `results/M2/mwr_locality.npz` (`9x9_dense_dist`, `9x9_dense_w`) |
-| recentrage | R_d = [4, 4, 0] sur la boîte 27×27 (cohérent avec s_red = 13/27 = 4·3 + 1, p = 3) ; inchangé | log `results/M2/logs/resigma_21857279.out` « [recenter] » |
+| médiane de \|Γ\|·N_cells sur les 41 266 états (carte niveau 1) | **3 189,01 meV** (Γ par défaut, ×81 = intensif) [v1 : 2 473,55] | `results/M2_plateau/level1_summary.csv` l. « 9x9,3,240,0.02 » (3189.0136) ; `specwd_9x9_prod.npz` (job 22058843, 2026-09-29) |
+| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 3 771,46 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 1,7944e+05 meV ; Born/T médian 45,846 (min et max non publiés par R10) [v1 : 2 343,59 ; 7 405,83 ; 3,302 (0,677, 16,45)] | `article/R10_plateau/c/table_v2_plateau.md` l.101–103 (`resonance_9x9.npz`) |
+| Re Σ médian, R_cut 3 | **+605,08 meV** ; \|Re Σ\|/Γ médian 0,2173 [v1 : 722,3 ; 0,227] (valeur non alignée : §8) | `results/M2_plateau/m_rcut_resigma.csv` (med_ReSigma_meV 605.0832…, med_absReSigma_over_Gamma 0.21732…) |
+| position de résonance E_res − E_D (argmax des états à ±1,5 eV) — **retiré du ch. 4 (R10 B.1 : argmax discret sur une couronne de la grille de sortie ; gardé ici pour la traçabilité)** | −0,175 eV (R_cut 2, 3, 4 à η 0,01/0,02 ; R_cut 0 à η 0,01/0,02 ; −0,181 à R_cut 0 et 4, η 0,05) ; −0,227 (R_cut 1) ; contre la grille de sortie (N_k^int 900) : −0,134 / −0,175 / −0,175 / −0,172 eV à 120² / 240² / 480² / 960² (non aligné : −0,227 / −0,202 / −0,200 / −0,191) [v1 : −1,238 ; −1,183 (η 0,01) ; −1,292 (R_cut 2 et 4)] | `results/M2_plateau/level1_summary.csv` (argmax_E_minus_ED_eV) ; `article/R10_plateau/b/b_results.json` B1 |
+| pics des courbes (ε − E_D) | Γ_T : −0,180 eV ; Γ_Born : +1,695 eV ; Γ_T/ρ₀ : −0,170 eV ; δρ : −0,787 eV ; ρ_dis : +1,710 eV ; \|T̄\| : −0,170 ; −Im T̄ : −0,172 ; min \|Re T̄\| : −2,140 ; Re T̄(E_D) = 12,624 eV, Im T̄(E_D) = −2,119 eV ; zéros de Re T̄ : −2,1425, −2,130, −2,095, −2,070, −2,0475, −2,005, −2,000, −0,210 [v1 : −1,24 ; +1,695 ; −0,015 ; −2,53 ; −2,53 ; −0,905 ; −1,29 ; −2,145 ; 2,521 / −0,091 ; aucun zéro] | `results/M2_plateau/resonance_9x9.npz` (`peak_*`, `ReTbar_at_ED`, `ImTbar_at_ED`, `Tbar_zero_crossings`) |
+| critère det / valeur propre | matrice complète (dim 145) : min de \|det[1 − Vg₀]\|/max = 1,261e-4 à −0,785 eV ; min_i \|λ_i\| = 0,0019 à −0,787 eV. Bloc σ (dim 87) : 4,379e-4 à −0,7875 eV, \|λ\| = 0,00186. Bloc π (dim 58) : aucun zéro, min \|det\|/max 1,706e-2 à −0,170 eV, min \|λ\| = 0,3440 à −0,127 eV [v1 : 2,09e-4 à −2,530 eV ; 0,0108 (λ = +0,0003 + 0,0108 i) ; secondaires −2,19 … −1,97 eV] | `article/R10_plateau/c/table_v2_plateau.md` l.115–117 ; `results/M2_plateau/resonance_criteria_9x9.npz` (`sigma_flag_at`, `sigma_flag_minlam`, `sigma_flag_det_rel`, `dim_pi`, `dim_sigma`) |
+| règle de somme de Friedel | ∫δρ sur toute la bande = **−1,0005** états (Tr[t g₀′]) et −1,0005 (Lloyd) ; bloc π −0,9981, bloc σ −0,0024 ; dans ±3 eV : +0,702 (π −0,256, σ +0,958) ; écart ponctuel max et cumuls aux bords non publiés par R10 [v1 : −0,0569 ; +1,782 ; 0,614 ; −2,190 / −0,409] | `results/M2_plateau/resonance_criteria_9x9.npz` (`sumrule`, `sumrule_lloyd`, `sumrule_window`, `*_pi`, `*_sigma`) |
+| Γ_T à c = 0,1 % sur ±1 eV | courbe finale dans le npz (`x_c`, `Gamma_c`, c_compare 0,001) ; min / max / valeur à E_D et ħ/Γ non publiés par R10 [v1 : 0,63 (+0,24) ; 5,55 (−1,00) ; 1,27 ; 419 / 1 025 fs] | `results/M2_plateau/resonance_criteria_9x9.npz` |
+| localité de M_W (9×9 dense) | ‖M_W(0, 0)‖ = 41,927 eV ; p_z–p_z sur le site de la lacune 31,546 eV, p_z de l'autre sous-réseau 0,640 eV ; ‖M_W(R, 0)‖ hors site inchangés par l'alignement (seule la diagonale sur site change : \|R\| = a : 1,963 / 1,388 / 1,501 ; √3 a : 0,137 / 0,639 / 0,122 …) ; abscisse maximale 15,588 a (images de Wigner-Seitz, contre 22,517 a avec les étiquettes brutes) [v1 : 9,365 ; 6,617 ; 0,044 ; a : 0,658 / 0,460 / 0,437 ; √3 a : 0,165 / 0,079 / 0,035 / 0,029 ; 2a : 0,082 / 0,081 — le texte du 2026-09-18 rangeait 0,460 dans la 2ᵉ couronne] | `results/M2_plateau/mwr_locality.npz` (`9x9_dense_onsite_norm`, `9x9_dense_onsite_pzvac`, `9x9_dense_onsite_pzB`, `9x9_dense_dist`, `9x9_dense_w`) ; table R10 l.72–74 |
+| recentrage | R_d = [4, 4, 0] sur la boîte 27×27 (cohérent avec s_red = 13/27 = 4·3 + 1, p = 3) ; inchangé ; C_N soustrait sur les 81 mailles de la boîte (étiquettes brutes de `Mwk_to_Mwr`) puis recentrage (`defect_mwr`) | lignes « [align] » et « [recenter] » des journaux `results/M2_plateau/logs/` ; R10 C.1 |
 
 ## 3. Contrôles réellement implémentés (à confronter aux dix du mémoire)
 
@@ -98,21 +109,21 @@ Type : **P** = porte bloquante (le code refuse, `raise`), **T** = test PASS/FAIL
 | C6 | **Test d'or** : matrice t locale = matrice T dense (`compute_T`) sur la même grille grossière, même sous-espace 5 WF | T, bloquant avant production | `scripts/test_local_tmatrix_real.py` ; `scripts/test_local_tmatrix.py` (synthétique) ; `scripts/test_local_rcut.py` (support tronqué ⇒ écart) | rel < 1e-8 et Γ ≥ −1e-8 | v2 : PASS 5×5 dense 1,80e-13 (M2, job 21852238, 2026-09-26, `results/M2/logs/r6golden_21852238.out`) ; 6×6 et 12×12 non rejoués avec M2 [v1 : 5×5 dense 2,3e-14 (2026-09-05, après le correctif d'unités), 6×6 6,06e-9, 12×12 3,68e-9 (2026-09-07) ; logs `golden_dense_20294198/20444392/20435998.out`] |
 | C7 | Positivité Γ_nk ≥ 0 | P | `local_tmatrix.py:269–271` ; `rcut_resigma.py:55` | min Γ < −1e-8 ⇒ `AssertionError` | v2 : min Γ_loc = +0,2646 eV (golden 5×5 dense M2) ; jamais déclenché [v1 : +0,26 eV] |
 | C8 | g₀ batché = g₀ de référence (restructuration exacte) et évaluation « point d'énergie le plus proche » vs état par état | T | `scripts/test_local_green_batch.py:15` | rel < 1e-12 | ≤ 1,21e-14 pour R_cut 0–3 (rejoué 2026-09-18, grille 90², 7 énergies ; ne dépend pas de M) ; nearest-grid vs exact : 5,4e-3 avec ne_per_eta 8 (cas synthétique) |
-| C9 | Pas de la grille d'énergie ne_per_eta | C | `rcut_resigma.py --npe` (9e0d1a9, 2026-09-09) | — | v1 seulement (non rejoué avec M2) : ne_per_eta 2 → 32 : médiane Γ 2 476,19 → 2 474,33 meV (0,08 %), E_res inchangé ; log `npe_test_20684497.out` |
-| C10 | R_cut (support de V_loc) | C | carte niveau 1 R_cut 0–3 (0–4 depuis R6) ; `rcut_resigma.py` R_cut 0–4 | plateau ≤ 5 % (énoncé `compute_spectral_wannier.py:118`) | v2, 9×9 : 3 009,7 / 3 223,7 / 3 170,2 / 3 132,6 / 3 147,5 meV (R_cut 0…4) ; écart médian par état à R_cut 4 : 13,1 / 4,1 / 4,2 / **1,15 %** ; Re Σ médian 17,7 / 682,5 / 479,2 / −26,6 / −410,3 meV (`results/M2/m_rcut_resigma.csv`) [v1 : 2 596,9 / 2 488,9 / 2 466,8 / 2 473,5 / 2 468,0 ; 10,0 / 4,6 / 1,5 / 0,63 % ; Re Σ 515 / 694 / 723 / 722 / 720] |
-| C11 | Grille de sortie × η (plateau conjoint) | C | carte niveau 1 (grilles 60/120/240, η 0,05/0,02/0,01) | ≤ 5 % quand η/2 et grille ×2 | v2, 9×9, R_cut 3 : 120² → 240² : 0,49 % (η 0,02), 0,70 % (η 0,01) ; η 0,02 → 0,01 à 240² : 0,10 % ; 0,05 → 0,02 : 1,42 % ; `results/M2/level1_summary.csv` [v1 : 0,07 % ; 0,3 % ; 0,001 % ; 0,74 %] |
-| C12 | Taille de super-cellule N (niveau 2, familles N mod 3) | C | `scripts/level2_families.py`, `results/M2/level2_summary.csv`, `level2_families.csv` | N ≥ 7 (config l. 6) | v2, 5/6/7/8/9/12 : 3 288,9 / 3 156,0 / 3 052,0 / 3 020,4 / 3 132,6 / 3 162,8 meV ; 7–9 : étendue 3,66 % ; famille 3m (6, 9, 12) : 0,96 % ; non-3m (5, 7, 8) : 8,60 % [v1 : 2 524,3 / 2 509,2 / 2 487,5 / 2 478,5 / 2 473,5 / 2 458,7 ; 0,56 % ; 2,0 % ; 1,84 %] |
-| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | v2 : Born/T médian 46,5 sur ±3 eV (min 1,11, max 211,7) [v1 : 3,30] |
-| C14 | Indépendance à l'alignement du potentiel (composante G ≈ 0) | K | `resonance_metrics.py:25–26, 42, 61` ; R6 : `--shift-L-meV` (V_loc + C·1) | — | v2 : (i) M − ⟨M^L⟩_diag·1, décalage 5 427,3 meV ⇒ écart relatif max 5,2e-2 (médian 1,1e-2) sur Γ_T ; (ii) C14 redéfini par R6 (C = ±25 meV ajouté uniformément à ΔV^L) : médiane des états 3 132,87 → 3 187,96 (+25 meV) / 3 158,54 (−25 meV), E_res −0,175 → −0,175 / −0,227 eV, max rel \|ΔΓ_T\| 0,192 / 0,163 (`results/M2/resonance_9x9_shiftL.npz`) [v1 : décalage 67,0 meV ⇒ 6,4e-4 (médian 1,8e-4)] |
-| C15 | Règle de somme de Friedel, deux formules (Tr[t g₀′] vs Lloyd) | K | `resonance_criteria.py:57–75` ; R6 : par bloc (`--blocks full,pi,sigma`) | — | v2 : −1,0007 / −1,0007 états sur toute la bande (π −0,998, σ −0,003) ; +0,669 dans ±3 eV [v1 : −0,0569 / −0,0569 ; +1,78] |
-| C16 | Critère de résonance (det, valeur propre minimale) et position du pic | K | `resonance_criteria.py:41–55` ; `compute_spectral_wannier.py:112` ; R6 : par bloc, `--flag-eV` | — | v2 : minimum global à −0,812 eV (λ = 0,0019 i), porté par le bloc σ ; bloc π sans zéro (min \|λ\| 0,397 à −0,170 eV) ; E_res(argmax Γ) = −0,175 eV [v1 : minimum unique à −2,530 eV ; E_res −1,24 eV] |
+| C9 | Pas de la grille d'énergie ne_per_eta | C | `rcut_resigma.py --npe` (9e0d1a9, 2026-09-09) | — | v1 seulement (non rejoué avec M2 ni avec l'alignement) : ne_per_eta 2 → 32 : médiane Γ 2 476,19 → 2 474,33 meV (0,08 %), E_res inchangé (E_res retiré, R10 B.1) ; log `npe_test_20684497.out` |
+| C10 | R_cut (support de V_loc) | C | carte niveau 1 R_cut 0–3 (0–4 depuis R6) ; `rcut_resigma.py` R_cut 0–4 | plateau ≤ 5 % (énoncé `compute_spectral_wannier.py:118`) | final, 9×9 : 3 013,96 / 3 273,25 / 3 222,08 / **3 189,01** / 3 200,00 meV (R_cut 0…4) ; écart médian par état à R_cut 4 : 9,79 / 4,52 / 1,98 / **0,66 %** ; Re Σ médian 29,5 / 761,8 / 737,2 / 605,1 / 696,7 meV (`results/M2_plateau/m_rcut_resigma.csv`, `rel_med_dGamma`, `med_ReSigma_meV`) [v1 : 2 596,9 / 2 488,9 / 2 466,8 / 2 473,5 / 2 468,0 ; 10,0 / 4,6 / 1,5 / 0,63 % ; Re Σ 515 / 694 / 723 / 722 / 720] |
+| C11 | Grille de sortie × η (plateau conjoint) | C | carte niveau 1 (grilles 60/120/240, η 0,05/0,02/0,01) | ≤ 5 % quand η/2 et grille ×2 | final, 9×9, R_cut 3 (`results/M2_plateau/level1_summary.csv`) : 120² → 240² : 3 200,71 → 3 189,01 (η 0,02 ; 0,37 %, arithmétique), 3 208,42 → 3 202,75 (η 0,01 ; 0,18 %) ; η 0,02 → 0,01 à 240² : 3 189,01 → 3 202,75 (0,43 %) ; 0,05 → 0,02 : 3 232,61 → 3 189,01 (1,37 %) [v1 : 0,07 % ; 0,3 % ; 0,001 % ; 0,74 %] |
+| C12 | Taille de super-cellule N (niveau 2, familles N mod 3) | C | `scripts/level2_families.py`, `results/M2_plateau/level2_summary.csv`, `level2_families.csv` | N ≥ 7 (config l. 6) | final, 5/6/7/8/9/12 : 2 997,14 / 3 149,15 / 2 942,85 / 2 974,18 / 3 189,01 / 3 264,51 meV ; (max − min)/moyenne, arithmétique (= R9 clôture R.4) : famille 3m (6, 9, 12) 3,60 %, non-3m (5, 7, 8) 1,83 %, 7–9 : 8,11 % [v1 : 2 524,3 / 2 509,2 / 2 487,5 / 2 478,5 / 2 473,5 / 2 458,7 ; 0,56 % ; 2,0 % ; 1,84 %] |
+| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | final : Born/T médian 45,846 sur ±3 eV (min et max non publiés) [v1 : 3,30] |
+| C14 | Sensibilité à l'alignement du potentiel : **C = ±rms du plateau de Kumagai–Oba (±9,05 meV à 9×9) ajouté uniformément au V_loc aligné** (R10 D10) | K | `resonance_metrics.py --shift-L-meV 9.05,-9.05` (V_loc + C·1 sur la boîte) | — | final : médiane des états 3 188,35 → 3 237,29 (+9,05 meV) / 3 157,07 (−9,05 meV) ; pic de la courbe Γ_T −0,180 eV inchangé ; écart relatif de la courbe Γ_T : max 6,36e-2 / 6,10e-2, médian 1,37e-2 / 1,24e-2 (`results/M2_plateau/resonance_9x9_shiftL.npz` ; R10 C.1). Retirés : la variante « M − ⟨M^L⟩·1 » (décalage 5 427,3 meV en v2) et C = ±25 meV autour du V_loc non aligné (R6 3.5) [v1 : décalage 67,0 meV ⇒ 6,4e-4 (médian 1,8e-4)] |
+| C15 | Règle de somme de Friedel, deux formules (Tr[t g₀′] vs Lloyd) | K | `resonance_criteria.py:57–75` ; R6 : par bloc (`--blocks full,pi,sigma`) | — | final : −1,0005 / −1,0005 états sur toute la bande (π −0,9981, σ −0,0024) ; +0,702 dans ±3 eV [v1 : −0,0569 / −0,0569 ; +1,78] |
+| C16 | Critère de résonance (det, valeur propre minimale) et position du pic | K | `resonance_criteria.py:41–55` ; R6 : par bloc, `--flag-eV` | — | final : minimum global à −0,785 eV (det) / −0,787 eV (\|λ\| = 0,0019), porté par le bloc σ ; bloc π sans zéro (min \|λ\| 0,3440 à −0,127 eV) ; pic de la courbe Γ_T −0,180 eV ; E_res(argmax Γ) retiré (R10 B.1 ; valeur −0,175 eV) [v1 : minimum unique à −2,530 eV ; E_res −1,24 eV] |
 | C17 | K sur la grille de sortie et dégénérescence π/π* en K | K | `resonance_metrics.py:86` ; `config/production.json` clé `K_red` | — | k_out[38480] = (2/3, 1/3), E(π) = E(π*) = E_D |
 | C18 | **Porte A.2 (R6)** : ΔV^L (grille) et ΔV^NL (projecteur KB de l'atome retiré) appliqués directement aux états de Bloch purs repliés, comparés à M2/N_cells | P | `scripts/gate_M_normalization.py` (`defects/deltav_pw.py`, `wavefunctions/sc_projection.py`) | max \|écart\| ≤ 1e-6 eV, sinon refus | 14 M2 OK (8 grossiers ≤ 3,3e-14 eV, 6 denses ≤ 1,2e-8 eV ; job 21820491, 2026-09-25) ; M v1 refusé (M^L = N_cells × direct, rapport 81,000000) |
-| — | **Grille interne N_k^int** | C (depuis P13) | `rcut_resigma.py --nk-int` ; `submit_nkint_check.sh` ; `nkint_check_post.py` | 5 % | v2, 300 → 600 : médiane −0,10 %, Γ_T(E_D) +4,81 % [v1 : −0,14 % ; +5,36 %] ; §6 |
+| — | **Grille interne N_k^int** | C (depuis P13) | `rcut_resigma.py --nk-int` ; `submit_nkint_check.sh` ; `nkint_check_post.py` | 5 % | final, 300 → 600 : médiane −0,09 %, Re Σ médian −0,03 %, Γ_T(E_D) +4,82 % (arithmétique sur `results/M2_plateau/nkint_check_9x9.csv`) [v1 : −0,14 % ; +5,36 %] ; §6e |
 
 ## 4. Points ouverts constatés le 2026-09-18
 
-(État v1 du 2026-09-18, inchangé ; valeurs v2 aux §2, §3 et §6, points R6 au §7.)
+(État v1 du 2026-09-18, inchangé ; valeurs finales aux §2, §3 et §6e, points R6 au §7, R10 au §8. Le point 4 (E_res) est clos par le retrait d'E_res du ch. 4, R10 B.1.)
 
 1. **N_k^int n'a jamais été varié.** Toute la production (cartes niveau 1, niveau 2, résonance, critères, Re Σ, test
 ne_per_eta) utilise 300² = 90 000 points ; c'est aussi le défaut codé en dur de `scattering_rate` /
@@ -235,6 +246,30 @@ Constats (seuil 5 %) : médiane de la fenêtre 300 → 600 : −0,10 % ; Γ_T(E_
 +6,23 % d'une valeur de −25 meV (écart absolu 1,56 meV) ; zone ±0,3 eV : +2,11 % (Γ), −0,47 % (Re Σ) ; E_res : −0,238 / −0,175 / −0,202 /
 −0,202 eV.
 
+### 6e. Même balayage, final (R10 C, job 22058845 `nkint`, 2026-09-29, 11 min 33)
+
+`scripts/submit_nkint_check.sh 9x9`, mêmes options, V_loc aligné (Kumagai–Oba, C_9 = −25,1437 meV) ; sources `results/M2_plateau/resigma_9x9_rc3_nk{150,300,450,600}.npz`,
+post-traitement `nkint_check_post.py` (C1f, job 22058855) → `results/M2_plateau/nkint_check_9x9.csv`. Écarts relatifs à nk_int = 600 (arithmétique sur le csv).
+
+| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV, retiré) | Γ_T(E_D) (meV, 4 états à K) | écart |
+|---|---|---|---|---|---|---|---|
+| 150 | 3 210,02 | +0,57 % | 607,78 | +0,42 % | −0,238 | 5 662,16 | +38,43 % |
+| 300 (production) | 3 189,01 | −0,09 % | 605,08 | −0,03 % | −0,175 | 4 287,11 | +4,82 % |
+| 450 | 3 196,11 | +0,13 % | 604,19 | −0,18 % | −0,175 | 4 120,38 | +0,74 % |
+| 600 | 3 191,93 | 0 | 605,25 | 0 | −0,175 | 4 090,11 | 0 |
+
+États à |ε − E_D| ≤ 0,3 eV (280 états ; colonnes `z_medG`, `z_medRe`, `z_E_res` du csv) :
+
+| nk_int | médiane Γ N_cells (meV) | écart | Re Σ médian (meV) | écart | E_res − E_D (eV, retiré) |
+|---|---|---|---|---|---|
+| 150 | 4 288,84 | +41,12 % | 6 056,58 | +2,03 % | −0,238 |
+| 300 (production) | 3 045,28 | +0,20 % | 5 858,01 | −1,31 % | −0,175 |
+| 450 | 3 079,75 | +1,34 % | 5 942,08 | +0,10 % | −0,175 |
+| 600 | 3 039,06 | 0 | 5 935,94 | 0 | −0,175 |
+
+Constats (seuil 5 %) : médiane de la fenêtre 300 → 600 : −0,09 % ; Γ_T(E_D) 300 → 600 : +4,82 % (v2 : +4,81 % ; v1 : +5,36 %) ; Re Σ médian 300 → 600 : −0,03 % ;
+zone ±0,3 eV : +0,20 % (Γ), −1,31 % (Re Σ). Le run nk 300 reproduit la production finale (3 189,01 meV, `level1_summary.csv`, `m_rcut_resigma.csv`).
+
 ## 7. R6 (2026-09-25 → 2026-09-27) : passage à M2
 
 - **Cause** : dans tous les M de production, M^L venait d'états de Bloch normalisés sur la super-cellule et M^NL sur la maille (R5-A.2) :
@@ -257,3 +292,39 @@ Constats (seuil 5 %) : médiane de la fenêtre 300 → 600 : −0,10 % ; Γ_T(E_
      négatif en v2 (Re Σ médian à R_cut 4 : −410,3 meV).
   6. Points du §4 à relire avec les valeurs v2 : §4.3 (médiane des états 3 132,6 meV, de la courbe 3 740,2 meV), §4.4 (E_res −0,175 /
      −0,202 / −0,227 eV ; critère det −0,812 eV), §4.6 (Friedel −1,0007 ; +0,669 dans ±3 eV).
+
+
+## 8. R10 (2026-09-29 → 2026-09-30) : alignement de Kumagai–Oba, production finale `results/M2_plateau/`
+
+- **Décision** (Greg, R10 étape G) : le potentiel de défaut est aligné par l'alignement de Kumagai–Oba (§1, ligne « Alignement du potentiel » ;
+  région, sphères de 1,0 Å, rms = incertitude ; 13 tailles, `config/production.json` bloc `alignment`) ; l'alignement à site unique est abandonné
+  (les six tailles 5…12 étaient « sans plateau » au critère de R9 A.1, max|écart| 14–234 meV ; R10 D14 : aucun critère d'arrêt). E_res est retiré
+  du ch. 4 (R10 B.1 : argmax discret, saute entre couronnes de la grille de sortie, −0,227/−0,202/−0,200/−0,191 eV non aligné et
+  −0,134/−0,175/−0,175/−0,172 eV final à 120²…960², `article/R10_plateau/b/b_results.json`).
+- **Code** (commits c61c118, ea91d61, 23ee3bb) : `wannier_interpolation.ws_images`, `ws_phase`, `Mwr_to_Mwk(_pairs)(ws=)` (étiquettes hors grille par
+  images de Wigner-Seitz) ; `local_tmatrix.defect_mwr(Mbk, U, U_dis, k, MP, n_box, C_N)` (M_W(R, R) − C_N·𝕀₅ sur les N² mailles de la boîte, étiquettes
+  brutes, puis recentrage) ; `config.matrices_dir` / `alignment_C` ; `results_dir = results/M2_plateau`, `matrices_dir = results/M2`, `results_dir_frozen`.
+- **Portes** : C.0 (C_N = 0 redonne R9 à 0,0 ; V_loc aligné = R9 au bit ; `Mwr_to_Mwk(ws)` sur la grille MP = sans ws à 4e-16) ; D6 d'`analyze_M.py`
+  2,9e-15 / 1,3e-15 / 1,4e-10 ; P-b2, P-c2 de l'audit de l'image minimale (seuil porté à 1e-10 eV par Greg, 16/16 profils OK) ; test d'or 5×5 dense non
+  rejoué (D12, recopié de R6, 1,80e-13). `results/M2_plateau/MD5SUMS_2026-09-30.txt` (33 fichiers).
+- **Valeur non alignée → finale** (les deux seules grandeurs pour lesquelles ce document garde la valeur non alignée) : Re Σ médian, R_cut 3 :
+  −26,65 → **+605,08 meV** (R_cut 0…4 : 17,67 / 682,54 / 479,23 / −26,65 / −410,29 → 29,47 / 761,79 / 737,20 / 605,08 / 696,74 ;
+  `m_rcut_resigma.csv` des deux répertoires ; colonnes Δ et C_9 × Δn_L dans `memoire/ch4/table_v1_final.md`, section j). Grandeur de Kaasbjerg (Fig. 3,
+  R10 C.2, `article/R10_plateau/c/c2_results.json`), eV Å², moyenne sur le disque de 0,1471 Å⁻¹, images de Wigner-Seitz : **K′ (insensible à
+  l'alignement)** valence 81,177 → 81,168, conduction 82,305 → 82,297 ; **K** valence 76,956 → **81,605**, conduction 78,168 → **82,821** ; D.2 (bloc
+  2×2 π/π* à (K, K)) : ½ Tr 75,87 → **86,59** (partie locale 59,53 → 70,25 ; non locale 16,34, inchangée). Médiane des états 3 132,60 → 3 189,01 meV
+  (table complète : `article/R10_plateau/c/table_v2_plateau.md`).
+- **Fig. 13–14 de Kaasbjerg (R8, `article/R8_kaasbjerg/out/{dos,spec,7a}/*.json`)** : position du maximum de ρ − ρ₀ (600², η_G 50 meV) à
+  c_i = 0,1 % / 1 % : non aligné −0,1625 / −0,1650 eV, final −0,1475 / −0,1500 eV ; gap à K (A_k, 1 %, η_G 25 meV) : non aligné 0,1176 eV, final
+  0,1182 eV ; valeurs lues sur l'article : maximum à −0,20 / −0,21 eV, points blancs à K +0,010 / +0,110 (gap 0,100 eV). (R8 a une troisième variante,
+  η_t = η_G, signalée dans `memoire/ch4/table_v1_final.md`, non utilisée.)
+- **Constats ouverts (non jugés)** :
+  1. C.6 : ⟨ΔV⟩_3D des 18×18, 24×24, 27×27 anormal (−8,48 / −62,74 / −99,91 meV contre +1…+29 meV pour 5…15 et 21) ; l'anomalie est dans le vide
+     (|z| > 5 Å : −30,4 / −167,5 / −252,9 meV en moyenne), la part du feuillet (|z| ≤ 3 Å) est régulière (+2,905 / +1,605 / +1,261 meV) ; SCF tous
+     convergés ; cause inconnue (`article/R10_plateau/c/c6_results.json`).
+  2. nscf dense 27×27 de production (`defect_uc_dense_27`, source de M_dense_9x9 et de la wannierisation) : 7 avertissements « c_bands: 1 eigenvalues
+     not converged » (R5 phase 0, `article/R5_base_vs_M/R5_rapport.md` l.51 et l.202).
+  3. Le modèle replié « aligné » de R10 A.3 utilise C_9 à site unique (−24,65 meV), pas C_9 de Kumagai–Oba (−25,14 meV) ; aucun modèle replié n'a été
+     recalculé avec l'alignement final.
+  4. Les points 1–2 et 5 de « Constats ouverts » de R6 (§7) restent ouverts ; le point 1 (convention intensive six tailles, 7,7e-2) est remplacé par
+     deux lignes par famille (3m 2,3e-2, non-3m 2,1e-2, OK ; R10 C.3).
```

## Annexe B — diff proposé de `CLAUDE.md` (non appliqué)

```diff
--- CLAUDE.md	2026-09-30 08:40:44.000000000 -0400
+++ CLAUDE.md (proposé)	2026-09-30 09:50:15.643129401 -0400
@@ -240,11 +240,15 @@
 - Taille : `figure.figsize` du style (6.5 × 3.6 po) pour une figure pleine largeur ; deux panneaux
   côte à côte = largeur 6.5 po, hauteur ajustée. Sauvegarde en PDF (vectoriel, pour LaTeX) et PNG
   (prévisualisation) dans `figures/`.
-- Données : lues uniquement dans le `results_dir` de `config/production.json` (`results/M2/*.npz` depuis R6 ; `results/M/` gelé, v1)
+- Données : lues uniquement dans le `results_dir` de `config/production.json` (`results/M2_plateau/` depuis R10, 2026-09-30 ;
+  matrices M2 dans `results/M2/` = `matrices_dir` ; `results/M/` (v1) et `results/M2/` (v2 non aligné) gelés, `results_dir_frozen`)
   ou dans les `.save`, jamais dans les logs ;
   les scripts de figures sont versionnés, les npz/CSV de production aussi (pas les logs).
 - Unités : M est stocké en Hartree et converti en eV UNE fois via
   `matrix_io.load_M_checked(..., units=matrix_io.EV)` ; Γ en meV, énergies relatives à $E_D$.
+- Alignement du potentiel de défaut : bloc `alignment` de `config/production.json` (C_N = ΔV_PA^(N), alignement de Kumagai–Oba,
+  13 tailles), appliqué en base de Wannier par `local_tmatrix.defect_mwr` ; nomenclature (non aligné / site unique / Kumagai–Oba /
+  final) et table v1 → final des chiffres du ch. 4 : `memoire/ch4/` (`table_v1_final.md`, `NOTES_TGAMMA.md` §1 et §8).
 
 ## Échantillonnage Γ des super-cellules : N mod 3 (fait physique à retenir)
 
@@ -252,7 +256,7 @@
 uniquement si N est un multiple de 3 : seuls 6×6, 9×9, 12×12 incluent les états de Dirac dans le SCF.
 5×5, 7×7, 8×8 (et 10×10, 11×11) partagent un artefact d'échantillonnage : ΔE_F = E_F(d) − E_F(p) de
 0.2 à 0.9 eV à la création de la lacune, absent pour N = 3m (|ΔE_F| < 3 meV). Les deux familles se
-comparent séparément ; la famille de convergence honnête est N = 6, 9, 12. Voir results/M2/sampling_table.csv (identique à results/M/, ne dépend pas de M).
+comparent séparément ; la famille de convergence honnête est N = 6, 9, 12. Voir results/M2_plateau/sampling_table.csv (identique à results/M2/ et results/M/, ne dépend pas de M).
 
 ## Sous-réseau de la lacune (A pour 5/7/8/9, B pour 6/12) — équivalence par symétrie
 
@@ -312,7 +316,7 @@
 6. Le scratch reste la copie de travail (les `outdir` des `.in` et les liens
    `data/` y pointent) ; on ne réécrit jamais les `outdir`/`prefix` d'un run terminé.
 
-Classement au 2026-09-23 (EM1 ajouté le 2026-09-25, M4_sigma, EM2 et EM3 le 2026-09-29) :
+Classement au 2026-09-23 (EM1 ajouté le 2026-09-25, M4_sigma, EM2 et EM3 le 2026-09-29, ch4 le 2026-09-30) :
 
 | Campagne | Répertoire | Statut | Miroir |
 |---|---|---|---|
@@ -326,3 +330,4 @@
 | EM2 références indépendantes (DFT directe aux anneaux + postw90 `kubo`) | `graphene/qe/electron_photon/EM2/` (copie `memoire/EM/EM2/`) | PRODUCTION (faite le 2026-09-29, rapport `EM2_rapport.md`) | `qe_tmp_backup/em2_bands_27/` (84 k, md5 OK 2026-09-29) |
 | EM3 figure et chiffres du §2.5 (anneau de 2.33 eV, postw90 `transl_inv`, tableau) | `memoire/EM/EM3/` (seul emplacement : calcul local < 1 s) | PRODUCTION | aucun `.save` ; tout se relance avec `make_em3_data.py` |
 | B σ(ω) complexe 27×27 (perspective B d'EM.md, hors mémoire) | `memoire/EM/B_sigma_complex/` (seul emplacement : calcul local de 3 min) | PRODUCTION | aucun `.save` ; tout se relance avec `b_prod.py` |
+| Chiffres du ch. 4 : table v1 → final, régions d'alignement, NOTES_TGAMMA au final, nombres de défauts.tex (2026-09-30) | `memoire/ch4/` (seul emplacement : calcul local de quelques secondes, pas de répertoire de travail hors dépôt) | PRODUCTION | aucun `.save` ; tout se relance avec `ch4_chiffres.py all` |
```

## Partie 5 — compléments (2026-09-30, après le commit 90d2022 de Greg)

Même cadre : chiffres copiés ou lus sur les fichiers de `results/`, colonnes arithmétiques seulement, aucune autre modification. Fichiers modifiés :
`memoire/ch4/ch4_chiffres.py` (fonctions `compl_5_1` … `compl_5_4`, `complements`, appariement du tex refait), `table_v1_final.md` (section
« Compléments » ajoutée en fin de fichier ; la table principale et les sections a–l sont inchangées, seul l'en-tête porte le nouveau HEAD),
`defauts_nombres.md` (refait), `NOTES_TGAMMA.md` (trois lignes), `README.md`, ce rapport ; nouveau : `NOTES_TGAMMA_partie5.diff`.
`alignement_regions.*` inchangés (contenu de HEAD). `défauts.tex` non modifié (md5 70140bc016e0c2d6f73cfa1235815129, identique à la partie 4).

- **5.1 tab:rcut_M** : `sv_mismatch_pi_blocks` (π–π*) et `diag_max_dM_over_max` (diag) des trois `m_rcut_convergence.csv`, 9×9 et 12×12, R_cut 0…6, en %
  (fraction × 100 par décalage de la virgule). Exemple 9×9 R_cut 3 : π–π* 2,5718 (v1) / 8,8431 (non aligné) / 2,5737 (final) ; diag 2,0915 / 7,2103 / 1,7074.
  **Lignes en double : aucune** dans les trois csv (14 lignes, 14 couples (taille, R_cut) distincts chacun ; le csv v1 est écrit en deux blocs, R_cut 0–3
  puis 4–6). Ce qui est égal dans le csv final, ce sont des colonnes : `max_dM_over_maxM` = `diag_max_dM_over_max` pour 12×12 R_cut 0–3 et 9×9 R_cut 0 ;
  `diag_max_dM_over_max` = `diag_abs_mismatch` sur toutes les lignes des trois csv.
- **5.2 niveau 1** (`défauts.tex` l.664–670) : 6 tailles × R_cut 0…4 à 240², η 0,02, depuis les trois `level1_summary.csv` (30 lignes). Constat :
  `results/M/level1_summary.csv` n'a aucune ligne R_cut 4 ; le point 9×9, R_cut 4 du mémoire (2468) est dans `results/M/m_rcut_resigma.csv`
  (med_Gamma_meV 2468.003090813311), rapporté à côté.
- **5.3 tab:échantillonnage** : ΔE_F et `Ved_radial_1.42A_meV` de `sampling_table.csv` (md5 égaux dans les trois répertoires), 8 tailles ; colonne arithmétique
  aligné = non aligné − C_N (C_N de `config/production.json`) : 6 +58,7091 ; 9 +68,4437 ; 12 +69,7896 ; 5 −306,8374 ; 7 −373,1857 ; 8 −255,4897 ;
  10 −288,3583 ; 11 −244,3493 meV.
- **5.4 scalaires lus sur les courbes** (définitions de `r6_compare_v1_v2.py`, rstats et cstats) : **porte PASS**, les 8 valeurs publiées (l.73, 74, 98, 99 ;
  v1 et non aligné) sont redonnées à la dernière décimale. Final : Born/T min 1,050, max 202,851 (±3 eV) ; Γ_T à c = 0,1 % sur ±1 eV : min 2,48 meV (+1,00 eV),
  max 40,58 meV (−0,18 eV), 8,47 meV à E_D ; ħ/Γ à ∓0,3 eV : 26 / 212 fs. Dans la table principale, les lignes R6 l.73, 74, 98, 99 restent « sans équivalent
  final » (non modifiées) ; les valeurs sont dans la section « Compléments ».
- **5.5 `defauts_nombres.md`** : 606 nombres (la notation scientifique et les seuils « 10^{-12} » comptent maintenant pour un nombre entier ; les tirets
  « 2--32 » ne donnent plus de signe) ; **appariés 234, non appariés 63, entiers courts hors tabular non comparés 309** ; lignes de la table principale jamais
  citées 22/207 (partie 4 : 553 nombres, 87 / 67 / 399, 59/207). L'index d'appariement comprend la table principale et les compléments 5.1–5.4. Non appariés
  restants (section dédiée du fichier) : paramètres (300², 100 Ry, 2.466 Å, centres et étalements de Wannier, a_CC 1.42, A_cell 5.27, E_D −4.239),
  23.9 et 34.2 eV Å² (carte de M̃_π, absents des tables sources), nombres de cellules 49 et 113, seuils des tests (10^{-12}, 10^{-5}, −10^{-8}), test d'or
  grossier 5×10^{-14} / 5.1×10^{-14} / 6.1×10^{-9} / 3.7×10^{-9}, erreurs n_e (8×10^{-4}, 5×10^{-3}), grille interne 1.4×10^{-3}, ordres de grandeur du
  texte (10^9, 10^{23}, 10^{14}, 2×10^{13}), 4.3 %, −1.05 et −0.70 eV. Les entiers des `tabular` (N, N mod 3, R_cut) s'apparient souvent par coïncidence.
- **NOTES_TGAMMA.md** : seules les grandeurs de 5.4 y figuraient déjà comme « non publiées » : §2 ligne « médiane de la courbe » (Born/T min et max), §2 ligne
  « Γ_T à c = 0,1 % » (min, max, E_D, ħ/Γ), §3 C13. tab:rcut_M et tab:échantillonnage n'y figurent pas ; les médianes de niveau 1 de 5.2 y sont déjà en valeurs
  finales (C10, C12). Contrôle `notes` : 15/15 valeurs finales retrouvées. Diff ci-dessous (annexe C) et dans `NOTES_TGAMMA_partie5.diff`.

**STOP — Greg relit et commit.**

## Annexe C — diff de `NOTES_TGAMMA.md`, partie 5 (`git diff` contre 90d2022)

```diff
diff --git a/NOTES_TGAMMA.md b/NOTES_TGAMMA.md
index c678547..e973270 100644
--- a/NOTES_TGAMMA.md
+++ b/NOTES_TGAMMA.md
@@ -84,13 +84,13 @@ sur une grille d'énergie de pas η/ne_per_eta (`resonance_metrics.py:55`) →
 | quantité | valeur | source |
 |---|---|---|
 | médiane de \|Γ\|·N_cells sur les 41 266 états (carte niveau 1) | **3 189,01 meV** (Γ par défaut, ×81 = intensif) [v1 : 2 473,55] | `results/M2_plateau/level1_summary.csv` l. « 9x9,3,240,0.02 » (3189.0136) ; `specwd_9x9_prod.npz` (job 22058843, 2026-09-29) |
-| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 3 771,46 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 1,7944e+05 meV ; Born/T médian 45,846 (min et max non publiés par R10) [v1 : 2 343,59 ; 7 405,83 ; 3,302 (0,677, 16,45)] | `article/R10_plateau/c/table_v2_plateau.md` l.101–103 (`resonance_9x9.npz`) |
+| médiane de la **courbe** Γ_T(ε) (moyenne lorentzienne) sur ±3 eV | 3 771,46 meV — **autre médiane** que la précédente (courbe vs états) ; Γ_Born : 1,7944e+05 meV ; Born/T médian 45,846 (min 1,050, max 202,851 sur ±3 eV, lus sur les courbes `Gamma_Born`/`Gamma_T` du npz) [v1 : 2 343,59 ; 7 405,83 ; 3,302 (0,677, 16,45)] | `article/R10_plateau/c/table_v2_plateau.md` l.101–103 (`resonance_9x9.npz`) ; min/max : `memoire/ch4/table_v1_final.md`, compléments 5.4 (porte : relecture v1 et non alignée = valeurs publiées) |
 | Re Σ médian, R_cut 3 | **+605,08 meV** ; \|Re Σ\|/Γ médian 0,2173 [v1 : 722,3 ; 0,227] (valeur non alignée : §8) | `results/M2_plateau/m_rcut_resigma.csv` (med_ReSigma_meV 605.0832…, med_absReSigma_over_Gamma 0.21732…) |
 | position de résonance E_res − E_D (argmax des états à ±1,5 eV) — **retiré du ch. 4 (R10 B.1 : argmax discret sur une couronne de la grille de sortie ; gardé ici pour la traçabilité)** | −0,175 eV (R_cut 2, 3, 4 à η 0,01/0,02 ; R_cut 0 à η 0,01/0,02 ; −0,181 à R_cut 0 et 4, η 0,05) ; −0,227 (R_cut 1) ; contre la grille de sortie (N_k^int 900) : −0,134 / −0,175 / −0,175 / −0,172 eV à 120² / 240² / 480² / 960² (non aligné : −0,227 / −0,202 / −0,200 / −0,191) [v1 : −1,238 ; −1,183 (η 0,01) ; −1,292 (R_cut 2 et 4)] | `results/M2_plateau/level1_summary.csv` (argmax_E_minus_ED_eV) ; `article/R10_plateau/b/b_results.json` B1 |
 | pics des courbes (ε − E_D) | Γ_T : −0,180 eV ; Γ_Born : +1,695 eV ; Γ_T/ρ₀ : −0,170 eV ; δρ : −0,787 eV ; ρ_dis : +1,710 eV ; \|T̄\| : −0,170 ; −Im T̄ : −0,172 ; min \|Re T̄\| : −2,140 ; Re T̄(E_D) = 12,624 eV, Im T̄(E_D) = −2,119 eV ; zéros de Re T̄ : −2,1425, −2,130, −2,095, −2,070, −2,0475, −2,005, −2,000, −0,210 [v1 : −1,24 ; +1,695 ; −0,015 ; −2,53 ; −2,53 ; −0,905 ; −1,29 ; −2,145 ; 2,521 / −0,091 ; aucun zéro] | `results/M2_plateau/resonance_9x9.npz` (`peak_*`, `ReTbar_at_ED`, `ImTbar_at_ED`, `Tbar_zero_crossings`) |
 | critère det / valeur propre | matrice complète (dim 145) : min de \|det[1 − Vg₀]\|/max = 1,261e-4 à −0,785 eV ; min_i \|λ_i\| = 0,0019 à −0,787 eV. Bloc σ (dim 87) : 4,379e-4 à −0,7875 eV, \|λ\| = 0,00186. Bloc π (dim 58) : aucun zéro, min \|det\|/max 1,706e-2 à −0,170 eV, min \|λ\| = 0,3440 à −0,127 eV [v1 : 2,09e-4 à −2,530 eV ; 0,0108 (λ = +0,0003 + 0,0108 i) ; secondaires −2,19 … −1,97 eV] | `article/R10_plateau/c/table_v2_plateau.md` l.115–117 ; `results/M2_plateau/resonance_criteria_9x9.npz` (`sigma_flag_at`, `sigma_flag_minlam`, `sigma_flag_det_rel`, `dim_pi`, `dim_sigma`) |
 | règle de somme de Friedel | ∫δρ sur toute la bande = **−1,0005** états (Tr[t g₀′]) et −1,0005 (Lloyd) ; bloc π −0,9981, bloc σ −0,0024 ; dans ±3 eV : +0,702 (π −0,256, σ +0,958) ; écart ponctuel max et cumuls aux bords non publiés par R10 [v1 : −0,0569 ; +1,782 ; 0,614 ; −2,190 / −0,409] | `results/M2_plateau/resonance_criteria_9x9.npz` (`sumrule`, `sumrule_lloyd`, `sumrule_window`, `*_pi`, `*_sigma`) |
-| Γ_T à c = 0,1 % sur ±1 eV | courbe finale dans le npz (`x_c`, `Gamma_c`, c_compare 0,001) ; min / max / valeur à E_D et ħ/Γ non publiés par R10 [v1 : 0,63 (+0,24) ; 5,55 (−1,00) ; 1,27 ; 419 / 1 025 fs] | `results/M2_plateau/resonance_criteria_9x9.npz` |
+| Γ_T à c = 0,1 % sur ±1 eV | min 2,48 meV (+1,00 eV), max 40,58 meV (−0,18 eV), 8,47 meV à E_D ; ħ/Γ à ∓0,3 eV : 26 / 212 fs (lus sur la courbe `x_c`, `Gamma_c` du npz, c_compare 0,001, mêmes définitions que R6) [v1 : 0,63 (+0,24) ; 5,55 (−1,00) ; 1,27 ; 419 / 1 025 fs] | `results/M2_plateau/resonance_criteria_9x9.npz` ; `memoire/ch4/table_v1_final.md`, compléments 5.4 (porte : relecture v1 et non alignée = valeurs publiées) |
 | localité de M_W (9×9 dense) | ‖M_W(0, 0)‖ = 41,927 eV ; p_z–p_z sur le site de la lacune 31,546 eV, p_z de l'autre sous-réseau 0,640 eV ; ‖M_W(R, 0)‖ hors site inchangés par l'alignement (seule la diagonale sur site change : \|R\| = a : 1,963 / 1,388 / 1,501 ; √3 a : 0,137 / 0,639 / 0,122 …) ; abscisse maximale 15,588 a (images de Wigner-Seitz, contre 22,517 a avec les étiquettes brutes) [v1 : 9,365 ; 6,617 ; 0,044 ; a : 0,658 / 0,460 / 0,437 ; √3 a : 0,165 / 0,079 / 0,035 / 0,029 ; 2a : 0,082 / 0,081 — le texte du 2026-09-18 rangeait 0,460 dans la 2ᵉ couronne] | `results/M2_plateau/mwr_locality.npz` (`9x9_dense_onsite_norm`, `9x9_dense_onsite_pzvac`, `9x9_dense_onsite_pzB`, `9x9_dense_dist`, `9x9_dense_w`) ; table R10 l.72–74 |
 | recentrage | R_d = [4, 4, 0] sur la boîte 27×27 (cohérent avec s_red = 13/27 = 4·3 + 1, p = 3) ; inchangé ; C_N soustrait sur les 81 mailles de la boîte (étiquettes brutes de `Mwk_to_Mwr`) puis recentrage (`defect_mwr`) | lignes « [align] » et « [recenter] » des journaux `results/M2_plateau/logs/` ; R10 C.1 |
 
@@ -113,7 +113,7 @@ Type : **P** = porte bloquante (le code refuse, `raise`), **T** = test PASS/FAIL
 | C10 | R_cut (support de V_loc) | C | carte niveau 1 R_cut 0–3 (0–4 depuis R6) ; `rcut_resigma.py` R_cut 0–4 | plateau ≤ 5 % (énoncé `compute_spectral_wannier.py:118`) | final, 9×9 : 3 013,96 / 3 273,25 / 3 222,08 / **3 189,01** / 3 200,00 meV (R_cut 0…4) ; écart médian par état à R_cut 4 : 9,79 / 4,52 / 1,98 / **0,66 %** ; Re Σ médian 29,5 / 761,8 / 737,2 / 605,1 / 696,7 meV (`results/M2_plateau/m_rcut_resigma.csv`, `rel_med_dGamma`, `med_ReSigma_meV`) [v1 : 2 596,9 / 2 488,9 / 2 466,8 / 2 473,5 / 2 468,0 ; 10,0 / 4,6 / 1,5 / 0,63 % ; Re Σ 515 / 694 / 723 / 722 / 720] |
 | C11 | Grille de sortie × η (plateau conjoint) | C | carte niveau 1 (grilles 60/120/240, η 0,05/0,02/0,01) | ≤ 5 % quand η/2 et grille ×2 | final, 9×9, R_cut 3 (`results/M2_plateau/level1_summary.csv`) : 120² → 240² : 3 200,71 → 3 189,01 (η 0,02 ; 0,37 %, arithmétique), 3 208,42 → 3 202,75 (η 0,01 ; 0,18 %) ; η 0,02 → 0,01 à 240² : 3 189,01 → 3 202,75 (0,43 %) ; 0,05 → 0,02 : 3 232,61 → 3 189,01 (1,37 %) [v1 : 0,07 % ; 0,3 % ; 0,001 % ; 0,74 %] |
 | C12 | Taille de super-cellule N (niveau 2, familles N mod 3) | C | `scripts/level2_families.py`, `results/M2_plateau/level2_summary.csv`, `level2_families.csv` | N ≥ 7 (config l. 6) | final, 5/6/7/8/9/12 : 2 997,14 / 3 149,15 / 2 942,85 / 2 974,18 / 3 189,01 / 3 264,51 meV ; (max − min)/moyenne, arithmétique (= R9 clôture R.4) : famille 3m (6, 9, 12) 3,60 %, non-3m (5, 7, 8) 1,83 %, 7–9 : 8,11 % [v1 : 2 524,3 / 2 509,2 / 2 487,5 / 2 478,5 / 2 473,5 / 2 458,7 ; 0,56 % ; 2,0 % ; 1,84 %] |
-| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | final : Born/T médian 45,846 sur ±3 eV (min et max non publiés) [v1 : 3,30] |
+| C13 | Born vs matrice T | K | `resonance_metrics.py:60` | — | final : Born/T médian 45,846 sur ±3 eV (min 1,050, max 202,851 ; compléments 5.4) [v1 : 3,30] |
 | C14 | Sensibilité à l'alignement du potentiel : **C = ±rms du plateau de Kumagai–Oba (±9,05 meV à 9×9) ajouté uniformément au V_loc aligné** (R10 D10) | K | `resonance_metrics.py --shift-L-meV 9.05,-9.05` (V_loc + C·1 sur la boîte) | — | final : médiane des états 3 188,35 → 3 237,29 (+9,05 meV) / 3 157,07 (−9,05 meV) ; pic de la courbe Γ_T −0,180 eV inchangé ; écart relatif de la courbe Γ_T : max 6,36e-2 / 6,10e-2, médian 1,37e-2 / 1,24e-2 (`results/M2_plateau/resonance_9x9_shiftL.npz` ; R10 C.1). Retirés : la variante « M − ⟨M^L⟩·1 » (décalage 5 427,3 meV en v2) et C = ±25 meV autour du V_loc non aligné (R6 3.5) [v1 : décalage 67,0 meV ⇒ 6,4e-4 (médian 1,8e-4)] |
 | C15 | Règle de somme de Friedel, deux formules (Tr[t g₀′] vs Lloyd) | K | `resonance_criteria.py:57–75` ; R6 : par bloc (`--blocks full,pi,sigma`) | — | final : −1,0005 / −1,0005 états sur toute la bande (π −0,9981, σ −0,0024) ; +0,702 dans ±3 eV [v1 : −0,0569 / −0,0569 ; +1,78] |
 | C16 | Critère de résonance (det, valeur propre minimale) et position du pic | K | `resonance_criteria.py:41–55` ; R6 : par bloc, `--flag-eV` | — | final : minimum global à −0,785 eV (det) / −0,787 eV (\|λ\| = 0,0019), porté par le bloc σ ; bloc π sans zéro (min \|λ\| 0,3440 à −0,127 eV) ; pic de la courbe Γ_T −0,180 eV ; E_res(argmax Γ) retiré (R10 B.1 ; valeur −0,175 eV) [v1 : minimum unique à −2,530 eV ; E_res −1,24 eV] |
```
