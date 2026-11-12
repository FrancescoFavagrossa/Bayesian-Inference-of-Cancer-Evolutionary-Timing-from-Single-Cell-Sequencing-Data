# Bayesian Inference of Cancer Evolutionary Timing from Single-Cell Sequencing Data

Code for the thesis (University of Milano-Bicocca). The project develops and
validates a Bayesian hierarchical model that, from a **single snapshot** of
single-cell sequencing data, reconstructs **when** the genomic events of a tumour
happened (clone-birth and copy-number times) and **how fast** each clone grew.

The model (`base_rw`) ties four signals to one latent timeline: the somatic mutation
counts on each branch, the split of copy-number-region mutations into pre- and
post-event counts, the clone abundances, and the phylogenetic tree. Growth rates get
a hierarchical log-scale random walk along the tree; event times use stick-breaking
on pseudotime `[0,1]`.

Validation is staged: first the model is fitted to data from its **own** generative
process (well-specified: parameter recovery, simulation-based calibration,
robustness sweeps), then to synthetic data from the spatial simulator **ProCESS**
(misspecified: localized-biopsy scenarios).

## Repository layout

| path | content |
|---|---|
| `stan_models/` | Stan model sources. `base_rw.stan` (complete model: joint times + growth), `base_rw_nosub.stan` (ablation: NegBin abundance without the children-subtraction). |
| `robustness_checks/` | Well-specified simulation study: data-generating process, recovery, SBC, hyperparameter sweeps, mutation posterior-predictive checks, and the no-subtraction ablation. See `robustness_checks/README.md`. |
| `Clone 4 DLP/`, `Clone 6 DLP/` | The two ProCESS scenarios used in the thesis (four-clone linear chain, six-clone branching tree), each sampled by a localized biopsy (DLP+). Each holds its simulation/inference script (`*.R`), figures and `results/`. |
| `Helpers and Utils/` | Shared helper functions sourced by the Clone scenario scripts (tree setup, bias tables, plotting). |
| `results/` | Aggregated recovery/bias tables (`bias_summary.csv`, `bias_detailed.csv`, `bias_wide.csv`). |

The LaTeX thesis and large reference data (e.g. `GRCh38/`) are kept outside this
code repository.

## Requirements

- **R** with `cmdstanr` (and a working `CmdStan` install), `posterior`, `ggplot2`,
  `dplyr`, `tidyr`.
- **ProCESS** (R package) only to regenerate the Clone scenarios from scratch.

## Reproducing the results

Run every script **from the project root**, e.g.

```bash
Rscript robustness_checks/01_dgp_recovery_coverage.R   # parameter recovery
Rscript robustness_checks/06_sbc.R                      # simulation-based calibration
Rscript robustness_checks/05_threepanel_sweeps.R        # robustness sweeps (data)
Rscript robustness_checks/05b_threepanel_plots.R        # robustness sweep figures
Rscript robustness_checks/04_dgp_posterior_plots.R      # ideal-DGP posterior figures
Rscript robustness_checks/09b_dlp_mut_refresh.R         # mutation posterior-predictive figures
```

The long runs (`01`, `06`) are **checkpointed**: re-launching resumes from the last
checkpoint (seeds are deterministic, so the continuation is identical).

Figures are written into the thesis image folder (`Thesis_final/immagini/`); results
tables into `results/` or each scenario's `results/`.
