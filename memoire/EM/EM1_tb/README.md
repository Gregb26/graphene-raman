# EM1 — éléments de position r(R) de la wannierisation de référence (ch. 4)

- But : produire `wannier_tb.dat` (H(R) et ⟨0i|r|Rj⟩) à partir du `.chk` de la wannierisation
  de référence du chapitre 4 (maille unitaire, grille 27×27×1, `graphene/qe/defects/unit_cell/27x27/`,
  LECTURE SEULE), par `restart = plot`, et vérifier jauge, centres, hermiticité, décroissance, format.
- Prompt d'origine : EM1 (série EM, couplage électron-photon, §2.5 du mémoire ; distincte des P et des R).
- Statut : PRODUCTION (2026-09-25). Aucun `.save` produit ; le `.save` du nscf d'origine est déjà
  miroité (`qe_tmp_backup/defect_uc_dense_27/`, md5 du 2026-09-17).
- Binaire : Wannier90 3.1.0 du module `quantumespresso/7.5` (`module restore qe`), celui de la référence.
- Contenu : `wannier.win` (copie + `restart = plot`), `wannier.chk`/`wannier.eig` (copies), `ref/` (liens vers
  la référence), sorties `wannier_tb.dat`, `wannier_hr.dat`, `wannier_wsvec.dat`, `wannier_u*.mat`,
  `wannier.bvec`, `wannier.wout`, `run.log`, `em1_check.py` + sorties `.txt`, `MD5SUMS_EM1_2026-09-25.txt`,
  rapport `EM1_rapport.md`.
- Copie versionnée : `memoire/EM1_tb/` dans le dépôt (créée le 2026-09-25, filtrée : sans `_tb.dat` (déjà suivi sous
  `wannier/27x27/`), `_hr.dat`, `.chk`, `.eig`, `.wout`, `.bvec`, `_u*.mat`).
- Suite : EM2 (nscf / bands.x, interpolation v_mn^(μ)(k) par Greg). Rien du mémoire touché.
