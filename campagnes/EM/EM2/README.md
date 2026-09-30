# EM2 — références indépendantes du couplage électron-photon (série EM, §2.5 du mémoire)

- But : confronter ħv = V†(∂H + i[H, A])V (module `electron_photon`, wannierisation 27×27 d'EM1) à (A) la DFT directe,
  ⟨c|p̂|v⟩ de QE 7.5 (`bands.x`, `lp = .true.`) sur 72 k des anneaux laser + 12 k d'un cercle de Dirac, et (B) σ(ω) de postw90 3.1.0 `kubo` contre M4.
- Prompt d'origine : EM2 (`campagnes/EM/EM2_prompt.md`, 2026-09-29) ; exécuté par Code sur rorqual, sans GO par étape (levé par Greg).
- Statut : PRODUCTION (2026-09-29). Référence `defects/unit_cell/27x27/` et `.save` `defect_uc_dense_27` lus seulement.
- A : `em2_kpoints.py` → `em2_kpoints_crystal.txt` ; `bands.in` (diff : `bands.in.diff`), `bands_pp.in`, `submit_bands.sh` (job 22034896),
  outdir scratch `qe_tmp/em2_bands_27` (miroir `qe_tmp_backup/em2_bands_27`) ; `em2_A_compare.py` → `em2_A_summary.txt`, `em2_A_table.txt`, `em2_A.npz` ; figure `em2_A_anneaux.{pdf,png}`.
- B : `postw90/` (`.chk`, `.eig` copiés d'EM1, `.mmn` en lien vers la référence ; `mkwin.sh`, `submit_postw90.sh` ; runs `k301*`, `k1201`,
  `k1201_ti`, `k1201_ti_nows`) ; `em2_B_ours.py` (+ `submit_B_ours.sh`, `em2_B_ours.log`, `.npz`) ; `em2_B_compare.py` → `em2_B_compare_*.txt`,
  `em2_postw90_sigma.npz` (k1201, lu par `scripts/fig/make_figures_em.py`) ; `em2_B_decompose.py` → `em2_B_decomposition.txt` ; figure `em2_B_sigma.{pdf,png}`.
- Rapport : `EM2_rapport.md` ; md5 : `MD5SUMS_EM2_2026-09-29.txt`. Copie versionnée : `campagnes/EM/EM2/` par `./em2_sync_copy.sh` (`cp -p`, liste explicite).
