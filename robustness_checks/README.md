# Robustness checks

Simulation study for the `base_rw` inference model: data are generated from the
model's own (well-specified) generative process and the model is re-fitted, to
assess parameter recovery, calibration, and sensitivity to hyperparameters. A few
scripts also re-fit the ProCESS scenarios or run targeted ablations.

Run every script from the **project root**, e.g. `Rscript robustness_checks/01_dgp_recovery_coverage.R`.

Dependencies: `cmdstanr`, `posterior`, `ggplot2`, `dplyr`, `tidyr`;
`stan_models/base_rw.stan` (and `base_rw_nosub.stan` for the ablation).

Two helper modules are sourced at the top of the DGP scripts:

- **`_dgp.R`** — the real-time frame (`Tf`, `t0`, `mu`, `phi`), the priors, the
  example tree (4 clones, `A -> B -> {C, D}`, CNAs on C and D) and the generative
  process `simulate_dgp()`. One call returns a `base_rw` `stan_data` list (growth,
  times, mutations **and** NegBin abundances) plus the ground truth, in normalised
  pseudotime and in real time (`t_real = t0 + t*(Tf-t0)`). `prior_exact = TRUE`
  draws from the exact model prior for SBC.
- **`_plot_style.R`** — `make_posterior_ci_plot()`, `theme_thesis()`, palette.

Figures are written into the thesis image folder `Thesis_final/immagini/`.

## Scripts
| file | what it does | output |
|---|---|---|
| `01_dgp_recovery_coverage.R` | 100 DGP replicates (checkpointed/resumable); aggregate `\|bias\|`, 95% coverage and CI width of event times (real time) and growth | `results/dgp_recovery.rds` |
| `04_dgp_posterior_plots.R` | faceted posteriors (event times, growth) for the showcase DGP replicate (seed 305) | `Thesis_final/immagini/dgp_posterior_{times,growth}.png` |
| `05_threepanel_sweeps.R` | 4 hyperparameter sweeps under the DGP (ω sd, σ-rate, n. cells, φ dispersion), per-event bias/coverage/CI-width in pseudotime | `results/sweep3.rds` |
| `05b_threepanel_plots.R` | 3-panel figures (`\|bias\|`, coverage, CI width) from `05`, one per sweep | `Thesis_final/immagini/robust_sweep_*.png` |
| `06_sbc.R` | Simulation-Based Calibration (Talts et al. 2018), 1000 sims (checkpointed): truth from the **exact** model prior, fit, rank the truth among thinned posterior draws; rank histograms with the uniform band | `Thesis_final/immagini/sbc_rank_histograms.png` |
| `09_mutation_dist.R` / `09b_dlp_mut_refresh.R` | per-equation posterior-predictive distributions of the mutation counts (DGP replicate and the two Clone DLP scenarios), with `m_X` / `alpha_k` / `beta_k` labels | `Thesis_final/immagini/dgp_posterior_mutations.png`, `mut_dist_c{4,6}_dlp.png` |
| `08_nosub_compare.R` | ablation: `base_rw` (with children-subtraction) vs `base_rw_nosub` (without), on DGP replicates and the ProCESS scenarios | `results/nosub_compare.rds` |
| `08b_nosub_dgp_test.R` | focused ideal-DGP ablation of the children-subtraction (full vs nosub), recovery of times and growth | `results/nosub_dgp_test.rds` |

Exploratory / legacy (not used in the final thesis): `02_prior_sensitivity.R`,
`03_realsim_hyperparam_sweeps.R` + `03b_realsim_sweep_plots.R` (older real-sim
sweeps, superseded by `05`/`05b`), `07_ppc_mutations.R` (PPC scatter, figure dropped
from the thesis).

The long runs (`01`, `06`) are **checkpointed**: re-launching resumes from the last
checkpoint (deterministic seeds, so the continuation is identical). Outside SBC,
`simulate_dgp()` clamps `sigma_log_w` to `[0.12, 0.6]` (biologically plausible growth
heterogeneity), since the untruncated prior occasionally generates implausible trees.

## Main findings
- **Recovery / calibration:** event times are recovered accurately and well
  calibrated (mean absolute error ≈ 6% of the `[t0, Tf]` window; ≈ 97% coverage at
  the 95% level). Growth rates keep a residual identifiability limit (≈ 23% median
  relative error) intrinsic to the single snapshot, but stay calibrated (≈ 97%
  coverage): the model widens the intervals rather than committing to a point value.
- **SBC:** of 1000 simulations, 966 converged and the rank histograms are consistent
  with uniformity for every parameter.
- **Priors:** robust to a wider ω prior and to the σ-rate; an over-tight ω prior
  (sd 0.2) degrades coverage → `ω ~ LogNormal(log 1.25, 0.5)` and `λ = 1` are used.
- **Sampling depth:** more cells reduce bias and tighten the credible intervals.
- **Abundance dispersion (φ):** event-time recovery is robust to the assumed φ across
  two orders of magnitude; only very strong overdispersion (φ = 10) mildly widens the
  intervals.
- **Children-subtraction ablation:** removing it wrecks the ideal-DGP recovery
  (time coverage 0.97 → 0.24, growth relative error 0.35 → ~8), confirming the
  subtraction is required for the disjoint single-cell counts.
