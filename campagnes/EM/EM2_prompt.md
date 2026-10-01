# EM2 — références indépendantes du couplage électron-photon (prompt pour la session Code sur rorqual)

Série EM (couplage électron-photon, §2.5 du mémoire), après EM1. Rédigé le 2026-09-29.

## Contexte

Le module `electron_photon` du dépôt calcule ħv = V†(∂H + i[H, A])V par interpolation de Wannier sur la wannierisation 27×27
d'EM1. Avant de l'écrire dans le mémoire, on le confronte à deux références indépendantes :
- **A. DFT directe** : les éléments de matrice de p̂ de QE (`bands.x`, `lp = .true.`) aux points k des anneaux résonants ;
- **B. postw90 `kubo`** : σ(ω) de Wannier90 lui-même, comparé à notre σ(ω) de production (M4).

Lis d'abord `CLAUDE.md` puis `campagnes/EM/EM.md` (§2 données, §8 M4, §10 EM2), et le rapport d'EM1
(`campagnes/EM/EM1_tb/EM1_rapport.md`, chemins `$P` et `$S` en tête).

## Règles

- **Aucun `sbatch`/`srun` sans ma confirmation explicite** de la commande et des ressources. Des calculs tournent déjà :
  ne touche à aucun job existant, ni à aucun `outdir` d'un run terminé (règle 6).
- La référence `$P/graphene/qe/defects/unit_cell/27x27/` et le `.save` de `defect_uc_dense_27` sont en **LECTURE SEULE**.
- Ni commit ni push. Ne modifie ni `src/`, ni `tests/`, ni `campagnes/EM/EM.md` : consigne tout dans le rapport, la session locale réalignera EM.md.
- Pas de répertoires de collègues.
- Campagne selon CLAUDE.md : répertoire de travail `$P/graphene/qe/electron_photon/EM2/` (à côté d'EM1_tb), README de dix lignes,
  statut **PRODUCTION**, copie versionnée `campagnes/EM/EM2/` par `cp -p` (inputs, scripts, rapport, analyses ; pas de `.save`, `slurm-*`,
  ni de fichier > 5 Mo). Le `.save` du nscf (72 k, petit) est miroité dans `qe_tmp_backup/` avec md5 en fin de campagne.
- Binaires : QE 7.5 et Wannier90 3.1.0 du module (`module restore qe`), ceux de la référence ; pas `~/bin/wannier90.x` (SIGSEGV, EM1).

## Étape 0 — Inventaire (lecture seule)

1. Le dépôt sur rorqual doit contenir `ring_kpoints_crystal` (F17) et `campagnes/EM/M4_sigma/` (production M4).
   Sinon, **arrête-toi et demande-moi de faire le `git pull`** (ne le fais pas toi-même).
2. Relève dans la référence : `nscf.in` (ecutwfc, nbnd, pseudo, cellule, prefix, outdir), `wannier.win`, et la présence de
   `wannier.mmn`, `.chk`, `.eig` (postw90 a besoin du `.mmn` pour r(R) : vérifie).
3. **Cellule** : les vecteurs de `CELL_PARAMETERS` du nscf doivent être ceux de `unit_cell_cart` de `wannier.win`
   (4.0354919061, ∓2.3298923383, 0 bohr ; même ordre). Sinon, les coordonnées cristallines de la liste ne veulent pas dire la même
   chose : arrête-toi et rapporte.

## A. DFT directe aux anneaux

1. **Liste de k** : `ring_kpoints_crystal(tb, K, E_D, {2.33: 48, 1.96: 12, 2.54: 12})` sur `wannier/27x27/wannier_tb.dat`,
   `K = make_grid_tb(tb, 3).K`, `E_D = -4.238895` eV. 72 points ; premier point à 2.33 eV : (0.73959906, 0.40626573, 0).
   Écris `em2_kpoints_crystal.txt` (`K_POINTS crystal`, poids 1, au moins 10 décimales) et garde à côté ħω et θ de chaque point.
2. **Run QE** : `calculation = 'bands'` (le potentiel est celui de la référence, donc la même référence d'énergie que le `.eig`
   de Wannier90), mêmes paramètres que le `nscf.in` de la référence sauf la liste de k, `nbnd = 20`, `nosym`, `noinv`.
   `outdir` NEUF sur le scratch : copie `charge-density*` et `data-file-schema.xml` du `.save` de référence, n'écris jamais dedans.
3. **`bands.x`** avec `lp = .true.` et `filp`. Lis `PP/src/write_p_avg.f90` de QE 7.5 et consigne : format du fichier, unités de p,
   moyenne sur les états dégénérés, et surtout **si le commutateur [V_NL, r] est inclus**. Sans lui, v = p/m diffère de la vraie
   vitesse de quelques pour cent : c'est à chiffrer, pas à corriger.
   Si lp refuse un calcul avec smearing, ou si valence/conduction y sont ambiguës, relance le run
   bands avec occupations = 'fixed' (nelec inchangé) et dis-le.
4. **Comparaison, k par k** : π, π* choisis par l'énergie autour de E_D (π = bande 4 de QE près de K, à vérifier).
   - ε_π, ε_π* de QE contre `compute_velocity(tb, k, 'berry')` aux mêmes k : écart attendu de l'ordre du meV ;
   - **|ħv_cv| dans le plan** (√(|ħv^x_cv|² + |ħv^y_cv|²), invariant de jauge et sans nœud ; 4.20 à 7.24 eV·Å sur l'anneau de
     2.33 eV) : QE contre les trois variantes du dépôt (complet, `centres_only(tb)`, `no_berry`). Attendu : complet à quelques pour
     cent au plus ; l'écart des deux autres est ce qu'on veut montrer.
   - Contrôle d'unités : près de K, |ħv_cv| doit tendre vers ħv_F = 5.469 eV·Å.
      - Contrôle d'unités, indépendant de Wannier : ajoute à la liste 12 k sur un cercle de rayon
     q = 0.005 Å⁻¹ autour de K. Sur ce cercle, ħv_F^DFT = ⟨Δε/(2q)⟩_θ (pentes des valeurs propres QE),
     et |ħv_cv| → ħv_F (limite de Dirac). Convertis p avec les unités lues dans write_p_avg.f90 et
     compare à ħv_F^DFT : accord à 10⁻³ près = p = m_e v, commutateur non local inclus ; sinon, l'écart
     chiffre [V_NL, r], et on le rapporte sans corriger. Compare aussi ħv_F^DFT à 5.469 eV·Å (Wannier).

## B. postw90 `kubo`

1. Sur une copie d'EM1 (`.win`, `.chk`, `.eig`, `.mmn` en lien symbolique si gros) dans `EM2/postw90/`.
2. Réglages (noms à confirmer dans la doc 3.1.0) : `berry = true`, `berry_task = kubo`, `fermi_energy = -4.238895`,
   `kubo_freq_min = 0.2`, `kubo_freq_max = 6.0`, `kubo_freq_step = 0.01`, élargissement gaussien fixe (`kubo_adpt_smr = false`).
   **Largeur** : notre η = 0.04 eV est l'**écart-type** de la gaussienne. Wannier90 utilise, je crois, exp(−x²)/√π avec
   x = ΔE/largeur, donc largeur = √2 η = 0.05657 eV : **vérifie dans le code** (`utility_w0gauss`) avant de fixer la valeur.
      kubo_eigval_max au-dessus de toutes les bandes : vérifie son défaut dans la doc 3.1.0 ; s'il
   dépend de dis_froz_max, il couperait des transitions que notre calcul garde.
3. berry_kmesh = 1201 1201 1 (pas un multiple de 3 : K hors grille, comme notre grille décalée ;
   la convergence en N de M4 rend l'écart avec 1200 négligeable). Essai à 301 × 301 pour le temps
   et la mémoire, puis propose-moi les ressources du vrai calcul.
   
4. Conversion : postw90 donne σ en S/cm pour la cellule 3D ; σ_2D = σ_3D × 100 × c, c = 15.875316×10⁻¹⁰ m, puis σ/σ₀ avec
   σ₀ = e²/(4ħ) = 6.0853×10⁻⁵ S. Contrôle : σ/σ₀ → 1 à 0.2 eV.
5. Comparaison avec `campagnes/EM/M4_sigma/em_sigma_full_N1200_eta0.04.npz` (même grille de ħω) : σ_xx, σ_yy, σ_xy sur toute la
   courbe. Repères M4 : σ_xx/σ₀ = 1.2637, 1.4056, 1.5115 à 1.96, 2.33, 2.54 eV ; 1.0026 à 0.2 eV ; pic de van Hove à 4.05 eV
   (6.9). postw90 applique `use_ws_distance`, nous non (écart attendu ≤ 10⁻⁵ ; notre grille est décalée, la sienne centrée en Γ).

## C. Rapport

`EM2_rapport.md` dans le répertoire de travail (puis `cp -p` vers `campagnes/EM/EM2/`) : chemins, commandes, jobs, md5, unités
et conventions vérifiées (p de QE, largeur de postw90), tableaux des écarts (A : par anneau, min/moy/max ; B : aux trois lasers
et sur toute la courbe), figures de contrôle éventuelles, et ce qui reste à décider. Termine par un résumé de dix lignes
que je rapporterai à la session locale.
