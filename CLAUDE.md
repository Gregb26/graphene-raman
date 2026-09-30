# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project & scientific goal

Compute the **electron–defect scattering matrix** M_mn(k',k) for a point defect in a host
crystal (graphene is the test system), from first-principles DFT outputs. The matrix is the
sum of a local and a non-local part:

- **M^L** (local): `M^L_mn(k',k) = Σ_{G,G''} Ṽ_ed^L(q+G) · C*_mk'(G+G'') · C_nk(G'')`, where
  `Ṽ_ed = V_defective_supercell − V_pristine_supercell` is the defect potential.
- **M^NL** (non-local): KB-separable, prefactor **(4π)²/Ω_uc** (Ω_uc = **unit-cell** volume, NOT
  supercell — this is the equation's convention and what the code implements). Built from
  per-(atom, l, i, m_l) factors `B = C* · F_l(|k+G|) · Y_lm · e^{∓iK·τ}` contracted with the KB
  energies E_li.

The matrix feeds downstream physics (e.g. carrier lifetimes, which need a sum over bands), and
is **Wannier-interpolated** from a coarse k-grid onto a fine grid / k-path.

The project was **ported from ABINIT to Quantum ESPRESSO**; ABINIT support has been removed.
Work is conducted in French; code and docstrings are in English.

## Repository layout (2026-09-30)

- `src/graphene_raman/` — the package (renamed from `electron_defect_interaction` on 2026-09-30; see "Module structure").
  `scripts/{m,t,epw,fig,validation,slurm}/` — production, analysis, figure and validation scripts, SLURM launchers
  (index `scripts/README.md`). `tests/` — pytest suite. `config/production.json` — frozen production parameters.
- `results/M2_plateau/` — **live** production products (`results_dir`, R10: v2 + Kumagai–Oba alignment); `results/M2/` — the M
  matrices (`matrices_dir`; the `.npy` exist on rorqual only) and the frozen unaligned v2 products; `results/M/` — frozen v1;
  `results/epw/` — chapter 5; `results/electron/`, `results/phonon/` — chapter 2 data (QE convergence, bands, DOS, matdyn). `results/wannier/<D>x<D>/` — tracked Wannier90 outputs per grid. `figures/` — thesis figures.
- `campagnes/<série>/<campagne>/` — versioned copies of every computation campaign, whatever its destination (thesis or article):
  series `R/` (R1…R10, chain M/T and DFT of the vacancy), `EM/` (electron–photon), `M/` (`ch4/`, numbers of chapter 4); index in
  `campagnes/README.md` (status live/archive, destination). `article/C_optique_lacunes/` — plan and bibliography of the future article.
- `notes/` — `NOTES_TGAMMA.md` (t/Γ chain, ch. 4), `NOTES_EPW.md`, `NOTES_EPW_REPERES.md` (ch. 5). `admin/` — closed housekeeping
  logs (`CLEANUP.md`, `REECRITURE_HISTORIQUE_2026-09-29.md`: table of rewritten commit hashes).
- `INVENTAIRE_2026-09-16.md` (root) — full map of the repo written on 2026-09-30 (role and status of every script, the t-matrix chain
  function by function in §3, what is broken or obsolete in §11, the reorganisation plan in §13). Read it before moving anything.
- The thesis (LaTeX) is a **separate repo**, `../msc-graphene-raman-defects` (`memoire/chapitres/*.tex`, `PROVENANCE.md` = register
  figure/table/number → source in this repo). Override its location with `$MSC_THESIS`.

## Tech stack

- **DFT**: Quantum ESPRESSO (pw.x). Local potential via `pp.x` (`plot_num=1`, filplot — there is
  **no HDF5 output from pp.x**). Pseudopotentials: UPF v2.
- **Wannier functions**: Wannier90 (run via `pw2wannier90.x` + `wannier90.x` — outside this repo).
- **Python** package `graphene_raman` (src layout, editable install). Dependencies are declared in
  `pyproject.toml` (numpy, scipy, h5py, matplotlib, tqdm; extras `mpi` = mpi4py, `test` = pytest, `campaigns` = ase, spglib);
  there is no `requirements.txt`. `mpi4py` is imported at module load by `defects.local_R`, `local_G`, `non_local`.
- **Cluster**: `rorqual` (Alliance/Calcul Québec, SLURM). Large supercells run there with the MPI
  drivers. `mpi4py` is built against **MPICH** → launch with the **MPICH** `mpirun`/`srun`, never
  the OpenMPI launcher (OpenMPI `mpirun` gives a `PMI_Init` failure).

## Commands

```bash
# Editable install (already done in .venv; redo it if the repo folder is renamed or moved: the .pth stores the absolute path)
.venv/bin/python -m pip install -e . --no-deps --no-build-isolation

# ONE test launcher: pytest (config in pyproject.toml, testpaths = tests). 208 tests: EM series, matrix_io, R8/R9/R10 functions
# (synthetic) + the chain M / chain T validation scripts of scripts/validation/, wrapped by tests/test_scripts_{M,tmatrix}.py
# (each script keeps its PASS/FAIL and exit code; the test asserts on it). Markers: slow (> ~10 s: the two real-data kubo
# tests, local_green_batch, ks_reconstruction, zero-pad non-regression), needs_data (local data/graphene, skipped when absent),
# cluster (M matrices / scratch .save, skipped elsewhere).
.venv/bin/python -m pytest tests                       # .venv/bin/pytest has a dead shebang: always use `python -m pytest`
.venv/bin/python -m pytest tests -m "not slow"         # ~1 min
# The validation scripts still run standalone (from the repo root; data paths in scripts/validation/_paths.py, root override EDI_DATA):
.venv/bin/python scripts/validation/test_ks_reconstruction.py     # core M = M^L + M^NL pipeline (local 5x5 data)
.venv/bin/python scripts/validation/test_wannier.py               # Wannier interpolation pipeline
.venv/bin/python scripts/validation/test_zero_pad_dense.py        # zero-pad densification of M^L (exact)
.venv/bin/python scripts/validation/test_local_tmatrix.py         # synthetic golden test of the local t-matrix
.venv/bin/python scripts/validation/test_local_rcut.py            # R_cut, extract_V_loc, mwr_locality (synthetic)
.venv/bin/python scripts/validation/test_local_green_batch.py     # local_green_batch / scattering_rate_fast (results/wannier/27x27)
.venv/bin/python scripts/validation/validate_wannier_bands.py     # Wannier vs DFT bands (coarse grid) + figures
.venv/bin/python scripts/validation/compare_bands_qe.py           # Wannier vs DFT along a k-path (needs bands.dat)
.venv/bin/python scripts/validation/compare_bands_w90_qe.py       # Wannier90 .dat vs QE bands.dat (argparse)
#    Cluster only (scratch paths, M matrices): test_local_tmatrix_real.py (real golden test), test_pad_vs_full_supercell.py.

# Repo root: never hard-code the folder path (renamed ab-initio-defects -> graphene-raman on 2026-09-24).
#   Python: derive it from __file__ (config.ROOT); shell: PROJ=${GRAPHENE_RAMAN:-$(git -C "${SLURM_SUBMIT_DIR:-$PWD}" rev-parse --show-toplevel)}
#   (sbatch copies the script to the spool, so $0/BASH_SOURCE are useless there). Scripts outside the repo (graphene/qe/...)
#   use ${GRAPHENE_RAMAN:-$PROJECTS/graphene-raman}; both variables are exported in ~/.bashrc on rorqual (not locally).
#   Wannier manifests store paths relative to the manifest directory.

# Production of M (rorqual, SLURM; MPICH srun). Each stage is a fresh process; outputs go to matrices_dir (results/M2).
sbatch scripts/slurm/submit_M.sh                             # coarse M: compute_M.py --stage ml (MPI) -> nl -> combine
sbatch scripts/slurm/submit_M_dense.sh 9x9                   # dense M on the (p*N)^2 grid: compute_M_dense_stages.py ml -> nl -> combine
# Chain T (rorqual): compute_spectral_wannier.py, rcut_resigma.py, resonance_metrics.py, resonance_criteria.py, ... through
# scripts/submit_*.sh and campagnes/R/R10_plateau/submit_r10.sh. Figures (local): scripts/make_figures{,_memoire,_epw,_em}.py.
```

There is no build/lint step beyond the editable install. Scripts are run from the repo root (several use relative
paths such as `results/wannier/27x27` or `results/…` through `config`). The former `scripts/run.py` and
`scripts/compute_M_cluster.py` no longer exist (replaced by `compute_M.py` and `compute_M_dense_stages.py`).

## Critical physics conventions (get these wrong and results are silently off)

- **Units**: Hartree atomic units throughout (Ha, Bohr). The QE XML (`data-file-schema.xml`) is
  written in Hartree (energies, ecutwfc) — no Ry→Ha conversion for those.
- **pp.x potential is in Rydberg**: `qe_io.get_pot` multiplies by 0.5 (Ry→Ha) by default. It
  returns the array as `(nr3, nr2, nr1)`, so callers do `.transpose(2,1,0)` to get `[ix,iy,iz]`.
- **UPF energies (PP_DIJ, PP_LOCAL) are in Rydberg** → `read_upf` converts to Hartree.
- **Plane-wave normalization**: `Σ_G |C_nk(G)|² = 1`.
- **Miller indices / G-vectors**: signed-integer FFT convention; map to FFT grid with
  `utils.fft_utils.map_G_to_fft_grid` (uses `np.rint`, NOT truncation).
- **M index convention** (all M arrays): `M[bra_band, k', ket_band, k]`.
- **Unit-cell wavefunctions must be on the full k-grid** (use `nosym` + `noinv` in pw.x) so the
  coarse k-grid is a complete Γ-centered Monkhorst-Pack grid.
- **Wannier Fourier-transform conventions**:
  - Wannier90 writes `H(R)` in **eV**; QE eigenvalues are in **Hartree** → multiply eps by
    `HA_TO_EV = 27.211386245988` before comparing.
  - Interpolation **must** divide by the Wigner-Seitz degeneracies:
    `H(k) = Σ_R e^{2πi k·R} H(R) / ndegen(R)`. Omitting `ndegen` degrades the coarse-grid
    Wannier-vs-DFT agreement by ~10⁴×.
  - The M and H Fourier transforms use the **plain** MP-dual R box (a simple N1×N2×N3
    parallelepiped). `use_ws_distance` (W90's per-element R-shifts by Wannier centres) is
    **deliberately not applied**, so H and M stay mutually consistent; this leaves a benign ~3 meV
    residual vs `V†εV` on the coarse grid. Since R10, `wannier_interpolation.ws_images` / `ws_phase`
    (Wigner-Seitz images of the R labels, weight 1/n_tie) exist for **off-grid** k points; they are used by
    `m_rcut_convergence.py` and `r10_driver.py` only — the local t-matrix never needs them.
  - graphene high-symmetry path convention: **K = (2/3, 1/3, 0)**, **M = (1/2, 0, 0)**, Γ = (0,0,0).

## Module structure (`src/graphene_raman/`)

- **config.py** — single loader of `config/production.json`: `load_production`, `results_dir` (products), `matrices_dir`
  (M files), `epw_dir` (`results/epw`, chapter 5 products), `electron_dir`, `phonon_dir` (`results/{electron,phonon}`, chapter 2 data), `alignment_C(cfg, size)` (C_N in eV), `dense_paths(cfg, size)` (dense
  `.save`, M file, Wannier dir; the scratch path is rorqual's), `ROOT`, `HA2EV` (the only definition of the Ha→eV factor). **Every
  access to results goes through these functions**, never through a literal `results/...` path.
- **io/matrix_io.py** — `save_M` / `load_M_checked` / `read_manifest`: every M file has a JSON sidecar (Bloch norm, units,
  `M_normalization`); `load_M_checked(..., require_bloch_norm=UNIT_CELL, units=EV, require_normalization=M_NORM_V2)` is the only
  Ha→eV conversion and refuses anything that is not v2.
- **io/wannier_provenance.py** — sha256 gauge gate: `write_wannier_manifest`, `load_wannier_checked(manifest)`.
- **io/qe_gamma_io.py** — Γ-only supercell wavefunctions (gate A.2, R10). The projwfc.x reader (`projwfc_io`, R4 only) was removed on 2026-09-30.
- **io/qe_io.py** — the QE I/O backend. All compute functions take an `io=` module; pass `qe_io`.
  Functions take a `prefix.save/` dir: `get_C_nk`, `get_G_red`, `get_k_red`, `get_A_volume`,
  `get_B_volume`, `get_ecut`, `get_x_red`, `get_eigenvalues`, `get_ngfft`, `get_pot`.
- **io/pseudo_io.py** — `read_upf` (QE UPF, the default `pseudo_reader`), `fq_from_fr` (Hankel
  transform of the radial projectors). The ABINIT `.psp8` reader was removed on 2026-09-30.
- **io/wannier_io.py** — Wannier90 readers: `read_w90_mat` (U / U_dis `.mat`, asserts
  unitarity/isometry), `read_w90_tb` (`_tb.dat`, returns `(HR, R, ndegen, rR, lattice)`: H(R) in eV,
  position operator r(R) (nR, 3, nw, nw) in Angstrom, lattice in columns), `read_w90_HR` (former name, kept as a
  wrapper returning `(HR, R, ndegen)` for frozen campaign drivers), `read_w90_hr` (`_hr.dat`, returns
  `(HR, R, ndegen)`), `check_hermicity_HR`. Tests: `tests/test_wannier_io.py` on `results/wannier/27x27/`.
- **defects/local_R.py** — M^L in **real space**: `compute_ML_R` (serial, dense BLAS, fast for
  moderate supercells but O(D³)), `prep_realspace_inputs` + `compute_ML_R_mpi` (grid-distributed,
  Allreduce, coarse grid), `compute_ML_R_mpi_shared` (**production kernel of the dense M**: zero-padding, node-shared u_nk),
  `fourier_resample`. Builds folded unit-cell Bloch parts via `wavefunctions/fold_wfk_to_sc`. Production never subtracts the
  mean of ΔV (`subtract_mean=False`; the default of `compute_ML_R` is `True`).
- **defects/local_G.py** — M^L in **reciprocal space** (`compute_ML_G*`, `zero_pad_potential`). **Historical**: replaced in
  production by `compute_ML_R_mpi_shared`; kept as the reference of `test_zero_pad_dense.py` / `test_pad_vs_full_supercell.py`.
- **defects/alignment.py** — potential alignment: `vacancy_site`, `far_atom_alignment`, `true_min_image_dist`,
  `atom_sphere_shifts` (source of the C_N; the values are frozen in `config/production.json`).
- **defects/deltav_pw.py**, **wavefunctions/sc_projection.py** — gate A.2 (ΔV applied to pure Bloch states vs M/N_cells),
  used by `scripts/m/gate_M_normalization.py`.
- **defects/non_local.py** — M^NL: `compute_M_NL` (serial), `compute_M_NL_mpi` (distributes the bra
  k'-index, Allreduce); helpers `build_K_vectors`, `compute_phase`, `compute_angular_part`. The
  pseudopotential path argument is `pseudo_path`; the reader is `pseudo_reader` (default `read_upf`).
- **wannier/wannier_hamiltonian.py** — `Hwr_to_Hwk(Hwr, Rw, k, ndegen=None)` → (Hwk, eigenvalues,
  eigenvectors); the eigenvectors are the per-k unitary to the smooth Bloch gauge.
- **wannier/wannier_interpolation.py** — the interpolation chain: `Mbk_to_Mwk` (V†MV, V=U_dis·U) →
  `Mwk_to_Mwr` (double FT) → `Mwr_to_Mwk` (inverse FT to fine grid) → `Mwk_to_Mbk` (back-rotate via
  H(k)). Top-level `wannier_interpolate(M, k_coarse, k_fine, wannier_tb, u_path, u_dis_path)` →
  `(M_bk_fine, E_fine)` (used by `test_wannier.py` only). Helpers `_match_kpoint_order`, `_infer_mp_grid`; `Mwr_to_Mwk_pairs`
  (R9), `ws_images`, `ws_phase` (R10).
- **wannier/supercell_fold.py** — supercell Hamiltonian folded at Γ, LDOS from eigenpairs (R4–R9 campaigns).
- **defects/many_body/local_tmatrix.py** — **core of chain T** (M → M_W(R,R') → V_loc → g₀ → t → Γ): `defect_mwr` (gauge rotation,
  double FT, subtraction of C_N, recentering), `recenter_mwr`, `mwr_locality`, `extract_V_loc`, `mp_grid`, `local_green_batch`,
  `local_t`, `scattering_rate` (exact), `scattering_rate_fast` (production), `cluster_ldos`. R_cut is a norm on the **reduced R
  labels**, not a Cartesian distance.
- **defects/many_body/disorder_average.py** — disorder-averaged DOS and A_k (Kaasbjerg PRB 101, 045433): `tbar_reduce`, `tbar_k`,
  `dos_average`, `spectral_path`, `sigma_eff`, `dirac_*`; only user: `campagnes/R/R8_kaasbjerg/r8_driver.py`.
- **defects/many_body/pole_criterion.py** — resonance criterion (det, λ_min), `local_t_cache`. **tb_models.py** — synthetic
  5-WF graphene bench (tests, R4).
- **electron_phonon/{phself,selfen}.py** — post-processing of EPW outputs (Γ^ep = 2 Im Σ), chapter 5.
- **wavefunctions/wfk.py** — `compute_psi_nk` (real-space ψ from C_nk on the FFT grid).
- **wavefunctions/fold_wfk_to_sc.py** — `compute_psi_nk_fold_sc` (unfold unit-cell ψ onto the
  supercell grid).
- **utils/** — `lattice` (`red_to_cart`, `build_k_path`), `planewaves` (`mask_invalid_G`),
  `fft_utils` (`map_G_to_fft_grid`, `fft_grid_from_G_red`).
- **plotting/palette.py** + **plotting/memoire.mplstyle** — thesis style and palette: `use_style()` (style found from `__file__`,
  never from the cwd), `save(fig, stem, outdir)` (PDF + PNG with fixed metadata, byte-reproducible), colours `NAVY`, `ORANGE`, …,
  `CMAP_SEQ`, `CMAP_DIV`. Every figure script and campaign driver imports it (no `sys.path` to `scripts/`).
- **defects/many_body/single_defect.py** — historical dense Bloch T-matrix (`compute_T`, `compute_G0`, `compute_G`, M in
  supercell norm). Not used in production, but it is the **reference of the golden tests** (`test_local_tmatrix*.py`,
  `test_local_rcut.py`): do not delete.
- **electron_photon/** — electron-photon coupling (EM series, in progress), conventions in its
  `__init__.py`: `tb_model` (`WannierTB`, graphene toy model, `make_wannier_tb(path)` from a real `_tb.dat`, `centres_only(tb)` for the tight-binding
  approximation of r, `pz_block(tb, pz)` for the p_z-only model, `extract_block` (not exported)), `kgrid` (reciprocal lattice, k grids),
  `velocity_operator` (single Fourier routine `fourier`, velocity with Berry connection, whole chain
  `compute_velocity(tb, k, mode)`), `ring`
  (resonant k points around K, `fermi_velocity`, `ring_stats` for the three velocity variants on the laser rings, `ring_kpoints_crystal` for the EM2 k list, `map_around_K` for the maps of the figure), `kubo` (σ(ω)/σ₀, driver `sigma_on_grid`, doping and temperature through `mu` and `kT`, complex σ through `kernel` = `gaussian_complex` or `lorentzian_complex`, Dirac references `kubo_doped_finite_T_analytical` and `kubo_doped_complex_analytical`), `diagnostics` (reports on a model:
  `hermiticity_report`, `symmetry_report`, `frozen_window_limit`). Tests (pytest) in
  `tests/test_{tb_model,kgrid,velocity_operator,ring,kubo,diagnostics,wannier_io}.py` (the other pytest files,
  `test_matrix_io_units.py` and `test_r{8,9,10}_functions.py`, cover chain M/T functions), shared fixtures in `tests/conftest.py` (real data:
  `tb_w90`, `eig_w90`, `w90_ref` on the tracked `results/wannier/27x27/`; every data-specific value and the sha256 of the
  files are in `W90_REF`/`W90_SHA256` of `conftest.py`).
  Plan and conventions: `campagnes/EM/EM.md`.

Index convention everywhere: `M[bra_band, k', ket_band, k]`, shape `(nband, nk, nband, nk)`.

## Data / supercells (local mirror of the cluster runs)

Under `data/graphene/` (this whole dir is git-ignored; real files locally, symlinks to the scratch on rorqual via
`scripts/m/link_data.sh`):

- `unit_cell/qe/defect_unit_cell_{5x5,11x11,12x12}.save/` — unit cells on the full MP grid (`wfc*.hdf5`, `C.upf`,
  `data-file-schema.xml`). Local naming is `defect_unit_cell_<N>.save` (`scripts/validation/_paths.py`); on the cluster it is
  `defect_<N>.save` (`link_data.sh`, production scripts).
- `supercell/qe/defect_5x5_p.save/` (pristine, `Vks_5x5_p`) and `defect_5x5_d.save/` (defective, `Vks_5x5_d`) — the only
  supercell pair available locally (plus the `Vks_11x11_{p,d}` potentials); it feeds the local validation scripts.
- **No M matrix is usable locally**: the 43 M2 files (`results/M2/M_*.npy` + sidecars, listed in `results/M2/MD5SUMS_2026-09-25.txt`)
  are on rorqual only; `results/M/M_ed_*.npy` (v1, git-ignored) may be present locally but are frozen. Everything that needs M
  (chain T, real golden test) runs on the cluster; locally one works from the npz/csv products of `results/M2_plateau/`.
- Production sizes: supercells 5×5 … 12×12 (all pristine + defective, Γ only), dense unit-cell grids 24–32 (`config["dense"]`);
  alignment constants for 13 sizes up to 27×27.

The same `.save` directories live on `rorqual`; large supercells are computed there. Exact cluster
paths are not inspectable from this checkout — confirm names before launching.

## Validation tests & state

- **scripts/validation/test_ks_reconstruction.py** — Test A reconstructs H_mn(k)=T+⟨ψ|V_loc|ψ⟩+V^NL and checks
  it equals diag(eps) from the XML (validates qe_io, kinetic term, get_pot, UPF projectors); Test B
  is the null-defect check (defect = pristine ⇒ M=0). **PASS** (reconstruction ~1e-8 Ha).
- **scripts/validation/test_wannier.py** — 5 tests: parser properties (U unitary, U_dis isometry+projector),
  Wannier-gauge FT round-trip (exact), full pipeline gauge-invariant spectrum, fine-grid smoke
  test, and a round-trip on a real M (skipped locally: no M file, no `_hr.dat`). Runs on the local 11×11 data.
- **scripts/validation/validate_wannier_bands.py** — Wannier-interpolated bands vs QE DFT eigenvalues on the
  coarse 5×5 grid (with/without ndegen contrast) + Γ–K–M–Γ Dirac-cone figure. **PASS**.
- **scripts/validation/compare_bands_qe.py** — Wannier vs DFT along a continuous k-path from a `bands.x`
  `bands.dat`. Median agreement ~21 meV. Aligns each band structure on its own Dirac point.

- **scripts/validation/test_local_tmatrix.py**, **test_local_rcut.py** (synthetic), **test_local_green_batch.py** (H of `results/wannier/27x27`),
  **test_local_tmatrix_real.py** (real golden test: local t = dense `compute_T` × N_cells; cluster only) — chain T.
- **tests/** (pytest, 208 tests) — EM series, `matrix_io`, the R8/R9/R10 functions of chain T on synthetic data, and the
  standalone scripts above through `tests/test_scripts_M.py` (zero padding, KS reconstruction, Wannier pipeline, pad vs full
  supercell) and `tests/test_scripts_tmatrix.py` (golden tests, `local_green_batch`, real golden test). Known local hiccup:
  `test_zero_pad_non_regression` (dense M^L on the 5x5 pair, > 10 min) once died with an MPICH/libfabric error
  (`OFI poll failed`, network interface) on the laptop — an environment issue, not a physics failure.

The standalone scripts take their data paths from `scripts/validation/_paths.py` or module constants; adjust them to the available
`.save` names.

## Known pitfalls & documented bugs

- **MPI launcher**: `mpi4py` is MPICH-built. Use the MPICH `mpirun`/`srun`; the OpenMPI launcher
  fails with `PMI_Init`.
- **eV vs Hartree**: Wannier90 H(R) is in eV; QE eigenvalues in Hartree (×27.211386245988).
- **ndegen**: must divide `H(R)/ndegen(R)` in the interpolation (see conventions above).
- **`_infer_mp_grid`**: `get_k_red` carries float noise (a k of 0 can come back as ~1−1e-15);
  fold with `np.mod(np.round(k,6),1.0)` then treat ~1 as 0 before counting the grid (already done).
- **`map_G_to_fft_grid`**: use `np.rint` (not int truncation), which dropped indices.
- **bands.x vs nscf reference**: a separate `bands.x` run and the nscf run that seeded Wannier90
  differ by a **rigid ~2.4 eV energy-reference offset** on the low bands. Align on the graphene
  Dirac point at K (mean of sorted bands [3:5]), NOT on the XML `fermi_energy`.
- **band-comparison metric**: compare each Wannier band to the **nearest** DFT band; the
  "lowest-nw sorted" comparison mis-pairs at band crossings (σ* dipping below π*) → spurious eV errors.
- **pp.x**: `plot_num=1` is the total potential (with XC); `plot_num=11` is without XC. Use filplot
  (no HDF5).
- **M^NL prefactor** uses **Ω_uc** (unit cell), not Ω_sc.

## What NOT to do

- **Do not `sbatch`/`srun` on the cluster without explicit confirmation.** Cluster jobs consume
  shared allocation; always confirm the command and resources first.
- **Do not browse or read colleagues' directories** on the cluster (other users' `scratch`/`home`).
  Stay within this project's paths.
- **Do not commit anything: Greg makes every commit himself** (2026-09-30). Use git read-only (`status`, `log`, `diff`,
  `show`, `grep`); edit files in the working tree and list what changed. A tracked file is moved with a plain `mv`.
- **Never add `.npy` or `.save` data to the repo.** `data/`, `jobs/`, `*.npy`, `*.save/` are git-ignored, and `results/` is
  ignored except for white-listed npz/csv products (see `.gitignore`): computed matrices, wavefunctions, potentials and `.save`
  directories stay out.
- **Do not reintroduce ABINIT** (netCDF4, `.psp8`, `WFK.nc`, `abinit_io`) into the main pipeline;
  the port to QE is complete.

## Typical workflow for a new supercell

1. **DFT (QE)**: relax/scf the defective and pristine supercells; run the unit cell scf/nscf on the
   full Γ-centered MP grid (`nosym`, `noinv`) matching the supercell folding.
2. **Potentials**: `pp.x` with `plot_num=1` (filplot) for both the pristine and defective
   supercells → `Vks_*` files.
3. **(For interpolation) Wannier90**: `pw2wannier90.x` + `wannier90.x` on the unit cell to produce
   `wannier_u.mat`, `wannier_u_dis.mat`, `wannier_tb.dat` (and `_hr.dat`).
4. **Matrix**: `scripts/m/compute_M.py` (coarse grid; stages `ml` (MPI), `nl`, `combine`; launcher `submit_M.sh`) or
   `scripts/m/compute_M_dense_stages.py` (dense (pN)² grid by zero-padding; launcher `submit_M_dense.sh`), pointing at the
   `.save` dirs, the `Vks_*` potentials and `C.upf`. Output is `M[bra_band, k', ket_band, k]` in Hartree with its JSON
   sidecar, v2 normalization (M = N_cells·M^L(supercell norm) + M^NL; gate: `scripts/m/gate_M_normalization.py`), in
   `matrices_dir`.
5. **Validate**: `test_ks_reconstruction.py` (sanity on the unit cell / null defect). For bands,
   `validate_wannier_bands.py` and, if a `bands.x` `bands.dat` is available, `compare_bands_qe.py`.
6. **Chain T** (the production use of M): `local_tmatrix.defect_mwr` → `extract_V_loc` → `local_green_batch` → `local_t` →
   `scattering_rate_fast`, driven by `compute_spectral_wannier.py`, `rcut_resigma.py`, `resonance_metrics.py`,
   `resonance_criteria.py`; parameters and controls in `notes/NOTES_TGAMMA.md`. (`wannier_interpolate` onto a fine grid of
   band-basis M exists but is only exercised by `test_wannier.py`.)

## Figures du mémoire (conventions obligatoires)

- Style : `memoire.mplstyle` dans le paquet, chargé par `graphene_raman.plotting.palette.use_style()` dans tout script de figure
  (jamais par un chemin relatif au cwd). Palette fixe dans `graphene_raman.plotting.palette`
  (2026-09-21) : **principale = bleu marine #000080** (`\definecolor{darkblue}{rgb}{0,0,0.5}`, couleur des
  hyperliens du mémoire), portée par la grandeur de production de chaque figure (9×9, matrice T, EPW,
  240² à 300 K, chaîne degauss 0.02) ; orange #eb6834 pour le contraste (Born, DFPT direct, Γ^ed au ch. 5),
  vert, jaune, rose, bleu ciel ensuite ; références en gris. Cartes : `CMAP_SEQ` (blanc → marine),
  `CMAP_DIV` (marine ↔ blanc ↔ orange foncé). Le cycler du style reprend `CYCLE` dans le même ordre.
  Une teinte par catégorie, jamais recyclée ; une séquence = une seule teinte.
- Tout le texte des figures est en FRANÇAIS (titres d'axes, légendes, annotations, titres de
  panneaux). `text.usetex` est actif : symboles en LaTeX, unités entre parenthèses. Exemples :
  « Énergie $\varepsilon - E_D$ (eV) », « Taux d'amortissement $\Gamma$ (meV) »,
  « Rayon de coupure $R_\text{cut}$ (mailles) », « Densité d'états (états/eV/cellule) »,
  « Grille $k$ », « Matrice $T$ », « approximation de Born ». Pas de mélange anglais/français
  dans une même figure. Les panneaux sont étiquetés (a), (b), …
- Taille : `figure.figsize` du style (6.5 × 4.0 po) pour une figure pleine largeur ; deux panneaux
  côte à côte = largeur 6.5 po, hauteur ajustée. Sauvegarde en PDF (vectoriel, pour LaTeX) et PNG
  (prévisualisation) dans `figures/<chapitre>/` (depuis le 2026-09-30 : `electron/` (ch. 2), `electron_defect/` (ch. 4),
  `electron_phonon/` (ch. 5, et `fig_phfreq_phdos` du ch. 2), `electron_photon/` (§2.5) ; défauts `--outdir` des sept
  `scripts/fig/make_figures*.py`), par `palette.save` : PDF et PNG **reproductibles au bit** (métadonnées fixées), ce qui est le
  test de recette de la copie vers le dépôt du mémoire (régénérer, comparer les md5). Les 20 figures incluses : 16 par
  `make_figures{,_memoire,_epw}.py` + `fig_em_coupling` (`make_figures_em.py`) + `fig_kb_pseudo_C` (`make_figures_electron.py`,
  redessinée le 2026-09-30) + `fig_electron_convergence`, `fig_ebands_edos` (même script, ex `qe_pp`, données `results/electron/*.dat`)
  + `fig_phfreq_phdos` (`make_figures_phonons.py`, données `results/phonon/`, même contenu que `fig_epw_phonons` de `make_figures_epw.py`). Les 9 figures de
  contrôle du ch. 4 retenues (R7 : `fig_size_3m`, `fig_localized_3m` ; R9 : `fig_resonance_vs_nkint{,_plateau}`, `fig_rcut_aligned`,
  `fig_folded_vs_R7` ; R10 : `fig_offset_profiles_13`, `fig_levels_vs_invN`, `fig_kaasbjerg_plateau_ws`) sont produites par
  `make_figures_controles.py` depuis les json/npz de `campagnes/R/` (fonctions extraites des pilotes, retirés) ; les 5 figures R8 par
  `campagnes/R/R8_kaasbjerg/r8_driver.py fig|sigeff` (données `out/`).
- Données : lues uniquement dans le `results_dir` de `config/production.json` (`results/M2_plateau/` depuis R10, 2026-09-30 ;
  matrices M2 dans `results/M2/` = `matrices_dir` ; `results/M/` (v1) et `results/M2/` (v2 non aligné) gelés, `results_dir_frozen`)
  ou dans les `.save`, jamais dans les logs ;
  les scripts de figures sont versionnés, les npz/CSV de production aussi (pas les logs).
- Unités : M est stocké en Hartree et converti en eV UNE fois via
  `matrix_io.load_M_checked(..., units=matrix_io.EV)` ; Γ en meV, énergies relatives à $E_D$.
- Alignement du potentiel de défaut : bloc `alignment` de `config/production.json` (C_N = ΔV_PA^(N), alignement de Kumagai–Oba,
  13 tailles), appliqué en base de Wannier par `local_tmatrix.defect_mwr` ; nomenclature (non aligné / site unique / Kumagai–Oba /
  final) et table v1 → final des chiffres du ch. 4 : `campagnes/M/ch4/` (`table_v1_final.md`) et `notes/NOTES_TGAMMA.md` §1 et §8.

## Échantillonnage Γ des super-cellules : N mod 3 (fait physique à retenir)

Les super-cellules N×N sont calculées avec le seul point Γ. Le point K de la maille se replie sur Γ
uniquement si N est un multiple de 3 : seuls 6×6, 9×9, 12×12 incluent les états de Dirac dans le SCF.
5×5, 7×7, 8×8 (et 10×10, 11×11) partagent un artefact d'échantillonnage : ΔE_F = E_F(d) − E_F(p) de
0.2 à 0.9 eV à la création de la lacune, absent pour N = 3m (|ΔE_F| < 3 meV). Les deux familles se
comparent séparément ; la famille de convergence honnête est N = 6, 9, 12. Voir results/M2_plateau/sampling_table.csv (identique à results/M2/ et results/M/, ne dépend pas de M).

## Sous-réseau de la lacune (A pour 5/7/8/9, B pour 6/12) — équivalence par symétrie

Les super-cellules 6×6 et 12×12 existantes ont la lacune sur le sous-réseau B (coordonnées de maille
(2/3, 2/3)) ; 5×5, 7×7, 8×8, 9×9 sur A ((1/3, 1/3)). Pour une lacune isolée non relaxée, A et B sont
reliés par l'inversion (ou le miroir) du réseau en nid d'abeille : mêmes Γ, mêmes M à une permutation
près des fonctions de Wannier pz(A) ↔ pz(B) et à une rotation près de la ZB. Aucun run « A » n'est
refait ; les scripts rapportent l'on-site pz–pz du sous-réseau de la lacune, et chaque manifest de M
porte `vacancy_sublattice` (scripts/m/tag_vacancy_sublattice.py). Tous les runs de super-cellule (5–12)
utilisent `assume_isolated='2D'` (vérifié dans les scf.out : « running with the 2D cutoff »).

## Chapitre 5 — couplage électron-phonon (EPW)

Les calculs EPW vivent hors dépôt dans `graphene/qe/epw/` (voisin du dépôt) ; l'état courant des jobs est dans
`notes/NOTES_EPW.md` et les repères durables dans `notes/NOTES_EPW_REPERES.md` (dans `notes/` depuis le 2026-09-30 ; sur rorqual,
les liens symboliques `graphene/qe/epw/NOTES_EPW.md` et `graphene/qe/epw/CLAUDE.md` doivent pointer vers ces deux fichiers) : grilles 24k-24q à
degauss 0.002 et 0.02 Ry, décisions arbitrées, conventions Γ^ep = 2 Im Σ, chemin q cartésien dans ph.x, piège OOM d'epw1,
sélection A1'/E2g par caractère jamais par index. Le volet t/Γ (chapitre 4) a son pendant dans `notes/NOTES_TGAMMA.md`.
Scripts versionnés ici : `scripts/epw_*.py` (pp_save, extract_gkk, validate, selfen_post, phself_post, d2_extract, ring_check,
dfpt_path_freq, phdos_extract, ed_vs_ep), `scripts/submit_epw_p{1_post,2_post_mv}.sh`, `scripts/fig/make_figures_epw.py` (à lancer
avec les options de production de `submit_post.sh`, sinon l'ancienne chaîne 0.002 est tracée) ; post-traitement dans
`electron_phonon/` ; résultats dans `results/epw/` (npz commis, logs non). **Production du ch. 5 = chaîne `24k-24q_mv0.02`
(degauss 0.02 Ry ; tranché par Greg le 2026-09-30)** ; la chaîne 0.002 (`NOTES_EPW.md` §1c–1f, npz sans suffixe) est l'ancienne.

## Couplage électron-photon (série EM, mémoire §2.5)

**Pour tout ce qui touche au couplage électron-photon, lire d'abord `campagnes/EM/EM.md`** : but, formule
ħv = V†(∂H + i[H, A])V, décisions verrouillées, format du `_tb.dat`, conventions du module, plan M0–M4
(fonctions F1–F19, vérifications, critères de sortie avec les valeurs de référence σ/σ₀) et statut. EM.md est
un guide, pas un cadre : en cas d'écart le code fait foi, et on réaligne EM.md (ne pas « corriger » le code vers le plan). Le code
est `src/graphene_raman/electron_photon/` (un module par responsabilité : les fonctions de M1–M4 vont
dans le module de leur rôle, pas dans un module par étape ; la lecture de fichiers reste dans `io/`) ; les campagnes de calcul EM vont sous `campagnes/EM/<campagne>/`.
Les étapes M0–M4 sont codées par Greg lui-même en mode technicien (skill `technicien`) : n'écrire ni ne
modifier son code sans « écris-le » ou demande explicite. Données de production : `campagnes/EM/M4_sigma/` (calcul local,
`m4_prod.py`, pilote de convergence dans `pilote/`) ; figure : `scripts/fig/make_figures_em.py` → `figures/fig_em_coupling.{pdf,png}` ;
prompt d'EM2 : `campagnes/EM/EM2_prompt.md` ; données de la figure (anneau de 2.33 eV, postw90 `transl_inv`) et tableau des chiffres : `campagnes/EM/EM3/`.

## Campagnes de calcul (règle du 2026-09-17, `admin/CLEANUP.md` ; précisée le 2026-09-23)

1. **Une campagne = un répertoire + un README (ou rapport) de dix lignes** : but, prompt
   d'origine (P/R), statut `TEST` ou `PRODUCTION`, date. Le répertoire vit à côté de ce
   qu'il prolonge (ex. `graphene/qe/defects/super_cell_relaxed/9x9/` pour R1).
2. **Le statut est décidé au lancement.** Un `TEST` dont le résultat est consigné
   (rapport, table du mémoire, npz dans `results/`) a ses `.save`, `.wfcN` et `outdir`
   supprimables sans arbitrage.
3. **`PRODUCTION`** : `.save` miroité vers `graphene/qe/qe_tmp_backup/<même chemin que
   qe_tmp/>` avec md5 en fin de campagne (`rsync -a -r --no-o --no-g --open-noatime`,
   puis `md5sum -c`), `.wfcN` en vrac supprimés, **rien d'unique sur le scratch** (purge
   Alliance : 60 jours sans accès ni modification, sans préavis fiable).
5. **Deux emplacements par campagne, jamais trois.** Le **répertoire de travail** est
   hors dépôt, à côté de ce qu'il prolonge (`graphene/qe/defects/super_cell_relaxed/9x9/`
   pour R1) : il contient TOUT — inputs (`*.in`, `submit.*`), scripts, rapports, sorties
   (`*.out`, `slurm-*`, `JOBID`, projwfc/pdos), et les `.in` y portent les `outdir` scratch.
   La **copie versionnée** est dans le dépôt sous `campagnes/<série>/<campagne>/` (séries `R/`, `EM/`, `M/` ; depuis le
   2026-09-30, ex-`article/R*` et ex-`memoire/{EM,ch4}` ; l'index `campagnes/README.md` porte le statut et la destination) : inputs, `submit.*`, scripts, rapports, analyses `.txt`,
   README, et les `.out` de pw.x s'ils font moins de ~5 Mo ; jamais `.save`, `slurm-*`,
   `JOBID`, sorties projwfc/pdos ni fichiers > 5 Mo (les lister dans le rapport). On édite
   dans le répertoire de travail puis on resynchronise la copie (`cp -p`), pas l'inverse.
6. Le scratch reste la copie de travail (les `outdir` des `.in` et les liens
   `data/` y pointent) ; on ne réécrit jamais les `outdir`/`prefix` d'un run terminé.

Classement au 2026-09-23 (EM1 ajouté le 2026-09-25, M4_sigma, EM2 et EM3 le 2026-09-29, ch4 et R2–R10 le 2026-09-30 ; tout sous
`campagnes/` depuis le 2026-09-30). Règle du ménage du 2026-09-30 : les pilotes de diagnostic des campagnes closes sont retirés
(fiche en tête du README, `git show 32efd9e:article/<campagne>/…` pour les retrouver) ; on garde rapports, tables, données, inputs, et
ce qui régénère les figures du mémoire :

| Campagne | Répertoire | Statut | Miroir |
|---|---|---|---|
| R1 relaxation lacune 9×9, Γ, nspin 1/2 | `graphene/qe/defects/super_cell_relaxed/9x9/nspin{1,2}` | PRODUCTION | `qe_tmp_backup/vacancy_relaxed/nspin{1,2}` (md5 OK 2026-09-22, revérifié 2026-09-23) |
| R1b contrôles 3×3×1 | `graphene/qe/defects/super_cell_relaxed/9x9/k3x3` | TEST consigné (R1b_rapport.md) | aucun ; `.save` scratch supprimables (manifeste 7) |
| Chaîne M (SCF supercellules, NSCF denses, mailles unitaires) | `graphene/qe/defects/{super_cell,unit_cell}` | PRODUCTION | `qe_tmp_backup/` (md5 OK 2026-09-17) |
| EPW 24k-24q et 24k-24q_mv0.02 | `graphene/qe/epw/` | PRODUCTION (outdir dans le projet) | — |
| EPW grilles test avril 2026 | `graphene/qe/epw/{36k-30q,30k-24q,24k-12q,16k-16q,16k-12q,12k-12q}` | TEST consigné (NOTES_EPW) | — (étage 1 du ménage) |
| EM1 éléments de position r(R) de la wannierisation 27×27 (mémoire §2.5, `restart = plot`) | `graphene/qe/electron_photon/EM1_tb/` (copie `campagnes/EM/EM1_tb/`) | PRODUCTION | aucun `.save` produit ; nscf déjà miroité `qe_tmp_backup/defect_uc_dense_27` (md5 2026-09-17) |
| M4_sigma σ(ω), cartes et chiffres des anneaux (mémoire §2.5, calcul local) | `campagnes/EM/M4_sigma/` (seul emplacement : calcul local de ~4 min, sans répertoire de travail hors dépôt) | PRODUCTION | aucun `.save` ; tout se relance avec `m4_prod.py` |
| EM2 références indépendantes (DFT directe aux anneaux + postw90 `kubo`) | `graphene/qe/electron_photon/EM2/` (copie `campagnes/EM/EM2/`) | PRODUCTION (faite le 2026-09-29, rapport `EM2_rapport.md`) | `qe_tmp_backup/em2_bands_27/` (84 k, md5 OK 2026-09-29) |
| EM3 figure et chiffres du §2.5 (anneau de 2.33 eV, postw90 `transl_inv`, tableau) | `campagnes/EM/EM3/` (seul emplacement : calcul local < 1 s) | PRODUCTION | aucun `.save` ; tout se relance avec `make_em3_data.py` |
| B σ(ω) complexe 27×27 (perspective B d'EM.md, hors mémoire) | `campagnes/EM/B_sigma_complex/` (seul emplacement : calcul local de 3 min) | PRODUCTION | aucun `.save` ; tout se relance avec `b_prod.py` |
| Chiffres du ch. 4 : table v1 → final, régions d'alignement, NOTES_TGAMMA au final, nombres de défauts.tex (2026-09-30) | `campagnes/M/ch4/` (seul emplacement : calcul local de quelques secondes, pas de répertoire de travail hors dépôt) | PRODUCTION | aucun `.save` ; tout se relance avec `ch4_chiffres.py all` |
| R2 relaxations de la lacune 5×5…12×12 | `graphene/qe/defects/super_cell_relaxed/series/` (copie `campagnes/R/R2_size_series/`) | PRODUCTION (2026-09-24) | `qe_tmp_backup/vacancy_relaxed/series/` (`MD5SUMS_series_2026-09-24.txt`) |
| R4 états π/σ de la lacune 9×9, DFT vs chaîne M → Wannier → T | `graphene/qe/defects/R4_quasi_lie/` (copie `campagnes/R/R4_quasi_lie/`) | TEST consigné (R4_rapport.md ; M en v1, erratum R5) | aucun `.save` |
| R5 écart H_p + M vs QE : base ou M ? (constat A.2, facteur N_cells) | `graphene/qe/defects/R5_base_vs_M/` (copie `campagnes/R/R5_base_vs_M/`) | TEST consigné (R5_rapport.md) | — |
| R6 M2 = N_cells·M^L + M^NL, porte A.2, régénération du ch. 4 | `graphene/qe/defects/R6_production_corrigee/` (copie `campagnes/R/R6_production_corrigee/`) | PRODUCTION (2026-09-25→27) ; produits `results/M2/` gelés depuis R10 | matrices M2 sur `/project` (`results/M2/MD5SUMS_2026-09-25.txt`), reconstructibles par `assemble_M2.py` |
| R7 (+R7c) scf 15×15…27×27, famille 3m ; relaxations nspin 1 | `graphene/qe/defects/R7_tailles_3m/` (copie `campagnes/R/R7_tailles_3m/`) | PRODUCTION (2026-09-25/26) | `qe_tmp_backup/` (wfc 15/18, md5) |
| R8 DOS et A_k moyennés sur le désordre (Kaasbjerg Fig. 13/14) | `campagnes/R/R8_kaasbjerg/` (`r8_driver.py`, `out/`, `fig/`) | TEST ; **cinq figures vont au ch. 4** (`dos_c`, `spectral_GKM`, `sigma_K`, `sensibilites`, `superposition`) | aucun `.save` |
| R9 contrôles : C_N, E_res vs N_k^int, chaîne repliée, Kaasbjerg | `campagnes/R/R9_controles/` | TEST consigné, clos le 2026-09-29 (R9_rapport.md) | — |
| R10 base unique du ch. 4 : alignement de Kumagai–Oba (13 tailles), rejeu de la production | `graphene/qe/defects/R10_plateau/` (copie `campagnes/R/R10_plateau/`) | TEST dont les sorties C sont la **production du ch. 4** (`results/M2_plateau/`), clos le 2026-09-30 | aucun `.save` ; `results/M2_plateau/MD5SUMS_2026-09-30.txt` |
| C σ(ω) du graphène avec lacunes (perspective, article) | `article/C_optique_lacunes/` (`PLAN.md`, `biblio/`) | plan (2026-09-29), après le mémoire | — |
