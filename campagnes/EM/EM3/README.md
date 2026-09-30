# EM3 — figure et chiffres du couplage électron-photon (série EM, §2.5 du mémoire)

- But : données de la figure du §2.5 (|ħv_cv|(θ) sur l'anneau de 2.33 eV avec la DFT directe ; σ_xx(ω)/σ₀ avec postw90 `transl_inv`) et tableau des chiffres.
- Prompt d'origine : EM3 (Greg, 2026-09-29), session locale ; décisions 7–10 d'EM.md (§1).
- Statut : PRODUCTION (2026-09-29), calcul local (< 1 s), aucun `.save` ; seul emplacement (pas de répertoire de travail hors dépôt).
- `make_em3_data.py` → `em2_postw90_sigma_ti.npz` (run `k1201_ti` d'EM2, lecteur `read_postw90` d'`EM2/em2_B_compare.py`, g_s = 2),
  `em3_ring_2p33.npz` (`ring`, 720 angles, trois variantes), `em_table.md`, `em_table.csv`, `make_em3_data.log` (quatre contrôles bloquants).
- Figure : `scripts/fig/make_figures_em.py` → `figures/fig_em_coupling.{pdf,png}` (lit aussi `EM2/em2_A.npz` et `M4_sigma/`).
- Rapport : `EM3_rapport.md` (contrôles, figure, légende proposée pour P28, écarts au prompt).
- Relancer : `.venv/bin/python campagnes/EM/EM3/make_em3_data.py && .venv/bin/python scripts/fig/make_figures_em.py`.
- Suite : P28 (texte du mémoire), sur décision de Greg.
