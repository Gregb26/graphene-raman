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

## Tech stack

- **DFT**: Quantum ESPRESSO (pw.x). Local potential via `pp.x` (`plot_num=1`, filplot — there is
  **no HDF5 output from pp.x**). Pseudopotentials: UPF v2.
- **Wannier functions**: Wannier90 (run via `pw2wannier90.x` + `wannier90.x` — outside this repo).
- **Python** package `electron_defect_interaction` (src layout, editable install). Core deps:
  numpy, scipy, h5py, mpi4py, matplotlib, tqdm.
- **Cluster**: `rorqual` (Alliance/Calcul Québec, SLURM). Large supercells run there with the MPI
  drivers. `mpi4py` is built against **MPICH** → launch with the **MPICH** `mpirun`/`srun`, never
  the OpenMPI launcher (OpenMPI `mpirun` gives a `PMI_Init` failure).

## Commands

```bash
# Editable install (already done in .venv)
pip install -e .

# Run the venv python
.venv/bin/python <script>

# Validation tests are standalone scripts: they print PASS/FAIL and exit 0/1 (no pytest).
# Hard-coded data paths live in scripts/_paths.py (override the data root with EDI_DATA).
.venv/bin/python scripts/test_ks_reconstruction.py     # core M = M^L + M^NL pipeline
.venv/bin/python scripts/test_wannier.py               # Wannier interpolation pipeline
.venv/bin/python scripts/test_zero_pad_dense.py        # zero-pad densification of M^L (exact)
.venv/bin/python scripts/validate_wannier_bands.py     # Wannier vs DFT bands (coarse grid) + figures
.venv/bin/python scripts/compare_bands_qe.py           # Wannier vs DFT along a k-path (needs bands.dat)
.venv/bin/python scripts/compare_bands_w90_qe.py       # Wannier90 .dat vs QE bands.dat (argparse)

# Serial driver (local)
.venv/bin/python scripts/run.py                        # serial M^L(real) + M^NL

# MPI / cluster driver (SLURM or local mpirun): --method real|reciprocal, --bands "all"|csv.
# Use the MPICH mpirun/srun (mpi4py is MPICH-built). Works locally too: mpirun -n N python ...
srun -n 256 .venv/bin/python scripts/compute_M_cluster.py \
    --uc UC.save --sc-p SCp.save --sc-d SCd.save \
    --pot-p SCp.save/Vks_p --pot-d SCd.save/Vks_d --upf UC.save/C.upf --out M_ed.npy
```

There is no build/lint step beyond the editable install. To "run a single test", run the
corresponding `scripts/*.py` directly; each is self-contained with hard-coded data paths in
`main()` that you adjust to the available `.save` directories.

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
    residual vs `V†εV` on the coarse grid.
  - graphene high-symmetry path convention: **K = (2/3, 1/3, 0)**, **M = (1/2, 0, 0)**, Γ = (0,0,0).

## Module structure (`src/electron_defect_interaction/`)

- **io/qe_io.py** — the QE I/O backend. All compute functions take an `io=` module; pass `qe_io`.
  Functions take a `prefix.save/` dir: `get_C_nk`, `get_G_red`, `get_k_red`, `get_A_volume`,
  `get_B_volume`, `get_ecut`, `get_x_red`, `get_typat`, `get_eigenvalues`, `get_ngfft`, `get_pot`.
- **io/pseudo_io.py** — `read_upf` (QE UPF, the default `pseudo_reader`), `fq_from_fr` (Hankel
  transform of the radial projectors), and `read_psp8` (legacy ABINIT `.psp8` reader, kept but
  unused by default — safe to delete if ABINIT is fully dropped).
- **io/wannier_io.py** — Wannier90 readers: `read_w90_mat` (U / U_dis `.mat`, asserts
  unitarity/isometry), `read_w90_HR` (`_tb.dat`, returns `(HR, R, ndegen)`), `read_w90_hr`
  (`_hr.dat`, same return), `check_hermicity_HR`.
- **defects/local_R.py** — M^L in **real space**: `compute_ML_R` (serial, dense BLAS, fast for
  moderate supercells but O(D³)), `prep_realspace_inputs` + `compute_ML_R_mpi` (grid-distributed,
  Allreduce). Builds folded unit-cell Bloch parts via `wavefunctions/fold_wfk_to_sc`.
- **defects/local_G.py** — M^L in **reciprocal space**: `prep_reciprocal_inputs`, `compute_ML_G`
  (serial), `compute_ML_G_mpi` (distributed over (k',k) pairs, Gatherv), plus `*_interp` variants
  for non-commensurate grids. O(D²) but gather-bound; crossover with real-space ~D≈1000.
- **defects/non_local.py** — M^NL: `compute_M_NL` (serial), `compute_M_NL_mpi` (distributes the bra
  k'-index, Allreduce); helpers `build_K_vectors`, `compute_phase`, `compute_angular_part`. The
  pseudopotential path argument is `pseudo_path`; the reader is `pseudo_reader` (default `read_upf`).
- **wannier/wannier_hamiltonian.py** — `Hwr_to_Hwk(Hwr, Rw, k, ndegen=None)` → (Hwk, eigenvalues,
  eigenvectors); the eigenvectors are the per-k unitary to the smooth Bloch gauge.
- **wannier/wannier_interpolation.py** — the interpolation chain: `Mbk_to_Mwk` (V†MV, V=U_dis·U) →
  `Mwk_to_Mwr` (double FT) → `Mwr_to_Mwk` (inverse FT to fine grid) → `Mwk_to_Mbk` (back-rotate via
  H(k)). Top-level `wannier_interpolate(M, k_coarse, k_fine, wannier_tb, u_path, u_dis_path)` →
  `(M_bk_fine, E_fine)`. Helpers `_match_kpoint_order`, `_infer_mp_grid`.
- **wavefunctions/wfk.py** — `compute_psi_nk` (real-space ψ from C_nk on the FFT grid).
- **wavefunctions/fold_wfk_to_sc.py** — `compute_psi_nk_fold_sc` (unfold unit-cell ψ onto the
  supercell grid).
- **utils/** — `lattice` (`red_to_cart`, `monkhorst_pack_grid`), `planewaves` (`mask_invalid_G`),
  `fft_utils` (`map_G_to_fft_grid`, `fft_grid_from_G_red`), `interpolation` (periodic
  tri-linear / cubic-spline).
- **defects/many_body/single_defect.py** — orphan / work-in-progress (not imported anywhere).

Index convention everywhere: `M[bra_band, k', ket_band, k]`, shape `(nband, nk, nband, nk)`.

## Data / supercells (local mirror of the cluster runs)

Under `data/graphene/` (this whole dir is git-ignored):

- `unit_cell/qe/defect_unit_cell_5x5.save/` — unit cell on the 5×5 grid (25 k-points). Has
  `wfc*.hdf5`, `Vks_uc_5x5` (pp.x potential), `C.upf`, `data-file-schema.xml`. **Wannier files
  (`wannier_*.mat`, `wannier_tb.dat`, `wannier_hr.dat`) and `bands.dat` are not always present** —
  the Wannier tests need them. nb=16 disentanglement bands → nw=5 Wannier functions.
- `supercell/qe/defect_5x5_p.save/` (pristine, `Vks_5x5_p`) and `defect_5x5_d.save/` (defective,
  `Vks_5x5_d`) — the 5×5 pair used by the validation tests.
- `supercell/qe/defect_12x12_p.save/` — **pristine only**; no defective 12×12 yet.
- `M_ed.npy` at repo root — a previously computed matrix, shape `(16, 25, 16, 25)`.

The same `.save` directories live on `rorqual`; large supercells are computed there. Exact cluster
paths are not inspectable from this checkout — confirm names before launching.

## Validation tests & state

- **scripts/test_ks_reconstruction.py** — Test A reconstructs H_mn(k)=T+⟨ψ|V_loc|ψ⟩+V^NL and checks
  it equals diag(eps) from the XML (validates qe_io, kinetic term, get_pot, UPF projectors); Test B
  is the null-defect check (defect = pristine ⇒ M=0). **PASS** (reconstruction ~1e-8 Ha).
- **scripts/test_wannier.py** — 5 tests: parser properties (U unitary, U_dis isometry+projector),
  Wannier-gauge FT round-trip (exact), full pipeline gauge-invariant spectrum, fine-grid smoke
  test, and a round-trip on the **real** `M_ed.npy`. **PASS** when the Wannier `.mat`/`_tb.dat`
  files are present.
- **scripts/validate_wannier_bands.py** — Wannier-interpolated bands vs QE DFT eigenvalues on the
  coarse 5×5 grid (with/without ndegen contrast) + Γ–K–M–Γ Dirac-cone figure. **PASS**.
- **scripts/compare_bands_qe.py** — Wannier vs DFT along a continuous k-path from a `bands.x`
  `bands.dat`. Median agreement ~21 meV. Aligns each band structure on its own Dirac point.

Tests use hard-coded paths in `main()`/module constants (e.g. `DATA = ".../defect_5x5.save"`);
update them to the actual `.save` names (the unit cell is currently `defect_unit_cell_5x5.save`).

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
- **Do not commit `.npy` or `.save` data.** `data/`, `results/`, `notebooks/`, `jobs/` are
  git-ignored, but **`M_ed.npy` at the repo root is NOT covered by `.gitignore`** — never `git add`
  computed matrices, wavefunctions, potentials, or `.save` directories.
- **Do not reintroduce ABINIT** (netCDF4, `.psp8`, `WFK.nc`, `abinit_io`) into the main pipeline;
  the port to QE is complete.

## Typical workflow for a new supercell

1. **DFT (QE)**: relax/scf the defective and pristine supercells; run the unit cell scf/nscf on the
   full Γ-centered MP grid (`nosym`, `noinv`) matching the supercell folding.
2. **Potentials**: `pp.x` with `plot_num=1` (filplot) for both the pristine and defective
   supercells → `Vks_*` files.
3. **(For interpolation) Wannier90**: `pw2wannier90.x` + `wannier90.x` on the unit cell to produce
   `wannier_u.mat`, `wannier_u_dis.mat`, `wannier_tb.dat` (and `_hr.dat`).
4. **Matrix**: run `scripts/run.py` (serial) or `scripts/compute_M_cluster.py` (MPI on rorqual)
   with `io=qe_io`, `pseudo_reader=read_upf`, pointing at the `.save` dirs, the `Vks_*` potentials,
   and `C.upf`. Output is `M[bra_band, k', ket_band, k]`.
5. **Validate**: `test_ks_reconstruction.py` (sanity on the unit cell / null defect). For bands,
   `validate_wannier_bands.py` and, if a `bands.x` `bands.dat` is available, `compare_bands_qe.py`.
6. **Interpolate** M onto a fine grid with `wannier.wannier_interpolation.wannier_interpolate`.
