# M4_sigma — σ(ω), cartes autour de K et chiffres des anneaux (série EM, §2.5 du mémoire)

- But : données de la figure et du tableau du §2.5 sur la wannierisation 27×27 (`wannier/27x27/`) : σ(ω)/σ₀ des trois variantes
  de la vitesse (complète, centres seuls, sans Berry), cartes Δε et |ħv_cv|² autour de K, `ring_stats`, ħv_F et ħω_froz.
- Prompt d'origine : M4 de `memoire/EM/EM.md` (§8), pilote puis production lancés en local par Code à la demande de Greg (2026-09-29).
- Statut : PRODUCTION (2026-09-29), calcul local, aucun `.save` ; tout est ici et se relance en ~4 min.
- `pilote/` : balayage N = 900, 1200, 1800 × η = 0.02, 0.04, 0.08 eV (`sweep.py`, `sweep.log`, npz), qui fixe N = 1200, η = 0.04 eV ;
  tableau de convergence complet dans `pilote/convergence.txt` (`convergence.py`).
- `m4_prod.py` → `em_sigma_{full,centres_only,no_berry}_N1200_eta0.04.npz`, `em_map_K.npz`, `em_ring_stats.json`,
  `em_scalars.json` (ħv_F, ħω_froz, σ aux lasers, provenance : commit et sha256 des modules), `m4_prod.log`.
- `make_table.py` → `em_table.csv` (format long : σ, ⟨|ħv_cv|²⟩, nœuds, rapports, par laser et variante) et `em_table.tex`
  (`tabular` booktabs pour le mémoire, généré : ne pas éditer à la main).
- Relancer : `.venv/bin/python memoire/EM/M4_sigma/m4_prod.py` (fonctions du paquet seulement).
- Suite : concordance avec postw90 `kubo` (EM2, rorqual), figure (EM3).
