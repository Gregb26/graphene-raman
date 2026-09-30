# B_sigma_complex — σ(ω) complexe des données 27×27 (perspective B d'EM.md §11)

- But : Re et Im de σ_xx(ω)/σ₀ (vitesse complète, noyau `gaussian_complex`, η = 0.04 eV, N = 1200, ħω = 0.20–6.00 eV) pour quatre cas :
  non dopé à T = 0 et 300 K, dopé n et p (μ = E_D ± 0.3 eV) à 300 K. Hors mémoire tant que Greg ne le décide pas ; prépare la perspective C.
- Prompt d'origine : perspective B (Greg, 2026-09-29), fonctions de Greg dans `kubo.py` (argument `kernel`), calcul lancé en local par Code.
- Statut : PRODUCTION (2026-09-29), calcul local, aucun `.save` ; seul emplacement.
- `b_prod.py` → `b_sigma_complex_N1200_eta0.04.npz` (hw, cases, mu, kT, sigma (4, nw, 3, 3) complexe, provenance), `b_prod.log`.
  La boucle partagée (une diagonalisation pour les quatre cas) est vérifiée contre `sigma_on_grid` à N = 120 (1e-13).
- `b_figure.py` → `b_sigma_complex.{pdf,png}` : Re/Im non dopé (signature de Kramers-Kronig autour de van Hove), et σ(μ) − σ(E_D)
  en dopage n et p contre le cône de Dirac à 300 K (quadrature, même noyau).
- Relancer : `.venv/bin/python memoire/EM/B_sigma_complex/b_prod.py && .venv/bin/python memoire/EM/B_sigma_complex/b_figure.py`.
- À savoir : non dopé à T = 0, Im σ converge seulement en 1/N (singularité w/Δ au point de Dirac, EM.md §11 B) ; les cas à 300 K sont convergés.
