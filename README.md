# graphene-raman

First-principles pipeline for electron–defect and electron–phonon coupling in
graphene, developed for a master's thesis on resonant Raman spectroscopy of
defective graphene (Université de Montréal).

The central piece is an original implementation of the **electron–defect
coupling matrix** $M_{mn}(\mathbf{k'},\mathbf{k})$ from plane-wave DFT
supercell calculations, its **Wannier interpolation**, and the resulting
**defect-induced electron damping rates** in the full $T$-matrix formalism.
The repository also contains the post-processing of the electron–phonon
counterpart ($g$, $\Gamma^{\mathrm{ep}}$, phonon linewidths) computed with EPW,
and the electron–photon coupling (velocity matrix elements, optical
conductivity) from the same Wannier model.

## What it computes

- **Electron–defect matrix elements** $M = M^{\mathrm{L}} + M^{\mathrm{NL}}$
  for a vacancy in graphene. The local part is computed in real space on the
  supercell grid, with an exact Fourier refinement by zero padding for dense
  k-grids (no spline interpolation of the potential); the non-local part is
  reconstructed from the Kleinman–Bylander projectors. The defect potential is
  aligned far from the vacancy (Kumagai–Oba-type alignment).
- **Two-argument Wannier interpolation** of $M$ onto arbitrarily dense
  k-grids, with recentering on the defect site and real-space truncation.
- **Defect-induced damping rates** $\Gamma^{\mathrm{ed}}_{n\mathbf{k}}$ and
  spectral signatures (quasi-bound states, DOS changes) via an exact local
  $T$-matrix scheme valid beyond the Born approximation, which overestimates
  the vacancy rates by more than an order of magnitude (median factor ≈ 46),
  and disorder-averaged DOS and spectral functions.
- **Electron–phonon post-processing** (EPW): interpolation validation,
  $g$ vs. direct DFPT, Fan–Migdal damping rates $\Gamma^{\mathrm{ep}}$,
  phonon linewidths of the $E_{2g}$ and $A_1'$ modes.
- **Electron–photon coupling**: velocity operator in the Wannier basis
  (including the Berry-connection term) and Kubo optical conductivity
  $\sigma(\omega)$, checked against direct DFT and `postw90`.
- A systematic **verification suite**: Kohn–Sham eigenvalue reconstruction
  against Quantum ESPRESSO (µeV agreement), hermiticity, gauge and
  normalization checks, convergence studies.

## Method and references

The plane-wave formulation of $M$ follows Lu, Zhou & Bernardi
[PRM **3**, 033804 (2019); npj Comput. Mater. **6**, 17 (2020)], replacing
their B-spline interpolation of the defect potential by an exact zero-padding
refinement. The local $T$-matrix scheme follows Kaasbjerg
[PRB **101**, 045433 (2020)], formulated here in the Wannier basis.
Electron–phonon quantities use EPW [Poncé et al., CPC **209**, 116 (2016);
Lee et al., npj Comput. Mater. **9**, 156 (2023)].

## Stack

Quantum ESPRESSO (SCF/NSCF, `pp.x`, `ph.x`), Wannier90, EPW, Python
(NumPy, SciPy, mpi4py, matplotlib). Production runs are done on the Rorqual cluster
(Digital Research Alliance of Canada) via the SLURM scripts
in `scripts/`.

## Layout

    src/electron_defect_interaction/   core package (M, Wannier, T-matrix,
                                       electron–photon, EPW post-processing)
    scripts/                           production, analysis, figures, SLURM launchers
    tests/                             pytest suite
    config/production.json             frozen production parameters and paths
    wannier/                           Wannier90 outputs per k-grid
    results/M2_plateau/                final production data (npz/csv)
    results/M2/, results/M/            earlier, frozen versions (unaligned; v1)
    results/epw/                       electron–phonon post-processing (npz)
    figures/                           thesis figures (generated)
    campagnes/{R,EM,M}/                computation campaigns: inputs, drivers,
                                       reports, tables (mostly in French)
    article/                           plan of the follow-up article
    notes/                             working notes (T-matrix chain, EPW)
    admin/                             housekeeping logs
    data/                              raw QE data, not tracked (see below)

Install with `pip install -e ".[mpi,test]"` (`mpi4py` needs an MPI library and is
imported by the $M$ kernels), then
`python -m pytest tests`.

## Data availability

Raw DFT data (wavefunctions, potentials) and the $M$ matrices themselves are
**not** in this repository: on the cluster, `data/` holds symlinks to the
storage, rebuilt with `scripts/link_data.sh`. Committed npz/csv files under
`results/` contain the final production quantities. Everything else is
regenerable from the scripts given the QE outputs.

## Status

This is research code accompanying a thesis and not a packaged library. Paths and
conventions are tailored to the author's cluster environment. The physics
and the verification suite are documented in the thesis (in French,
forthcoming).

## Author

Grégoire Barrette, M.Sc. in computational condensed matter physics,
Université de Montréal. Thesis supervised by Michel Côté.

## Acknowledgments

Parts of this codebase (analysis scripts, figure generation, verification
tooling, and cleanup/automation workflows) were developed with the
assistance of AI tools (Claude and Claude Code, Anthropic). All physics,
derivations, and scientific decisions are the author's own. AI-assisted code
was reviewed by the author and validated through the verification suite
described above.

Computations were performed on the Rorqual cluster of Calcul Québec and the
Digital Research Alliance of Canada. This work was supported by NSERC.

## License

MIT
