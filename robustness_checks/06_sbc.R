## ----setup, include=FALSE-----------------------------------------------------
# Simulation-Based Calibration (Talts, Betancourt, Simpson, Vehtari, Gelman 2018,
# arXiv:1804.06788): draw the truth from the EXACT model prior, fit, and rank the
# truth among thinned posterior draws. Calibrated inference => uniform ranks.
knitr::opts_chunk$set(
  echo    = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.width  = 11,
  fig.height = 7.5
)

## ----libraries----------------------------------------------------------------
library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)
library(tidyr)

## ----dgp----------------------------------------------------------------------
# Real-time frame, tree, priors and the generative process (simulate_dgp()).
source("robustness_checks/_dgp.R")
source("robustness_checks/_plot_style.R")   # theme_thesis() + palette

## ----sbc-params---------------------------------------------------------------
N_sim  <- 1000     # SBC simulations (very deep run: ~50 per bin -> very sharp histograms)
L      <- 1000     # thinned posterior draws per simulation -> ranks in 0..L
N_bins <- 20       # bins for the rank histograms

## ----model--------------------------------------------------------------------
model <- cmdstan_model("stan_models/base_rw.stan")
tracked <- c("w", "sigma_log_w", "omega_out",
             paste0("tevents[",    seq_len(Nevents), "]"),
             paste0("w_node_out[", seq_len(Nnodes),  "]"))

## ----sbc-loop (checkpointed / resumable)--------------------------------------
# A checkpoint is written every 25 simulations; if the run is interrupted, simply
# launch the script again and it resumes from the last checkpoint (the per-sim seeds
# are deterministic, so the continuation is identical to an uninterrupted run).
dir.create("robustness_checks/results", showWarnings = FALSE, recursive = TRUE)
ckpt <- "robustness_checks/results/sbc_ckpt.rds"

ranks <- NULL
if (file.exists(ckpt)) {
  cp <- readRDS(ckpt)
  if (nrow(cp$ranks) == N_sim && identical(colnames(cp$ranks), tracked)) {
    ranks <- cp$ranks; diag <- cp$diag; L_thin <- cp$L; start <- cp$done + 1L
    cat(sprintf("resuming SBC from sim %d / %d\n", start, N_sim))
  }
}
if (is.null(ranks)) {
  ranks  <- matrix(NA_integer_, N_sim, length(tracked), dimnames = list(NULL, tracked))
  diag   <- data.frame(sim = 1:N_sim, max_rhat = NA_real_, divergences = NA_integer_)
  L_thin <- NA_integer_; start <- 1L
}

for (s in start:N_sim) {
  sim <- simulate_dgp(seed = 1000 + s, prior_exact = TRUE)   # exact-prior truth
  fit <- model$sample(data = sim$stan_data, chains = 4, parallel_chains = 4,
                      iter_warmup = 1500, iter_sampling = 2000, adapt_delta = 0.95,
                      max_treedepth = 13, seed = 1, refresh = 0, show_messages = FALSE)

  # true value for each tracked variable (tevents are in normalised [0,1])
  truth <- c(w = sim$true_w_bar, sigma_log_w = sim$true_sigma, omega_out = sim$true_omega)
  truth[paste0("tevents[",    seq_len(Nevents), "]")] <- sim$true_t
  truth[paste0("w_node_out[", seq_len(Nnodes),  "]")] <- sim$true_w

  dm  <- as_draws_matrix(fit$draws(tracked))
  idx <- unique(round(seq(1, nrow(dm), length.out = L)))   # thin to ~independent draws
  dm  <- dm[idx, , drop = FALSE]; L_thin <- length(idx)
  for (v in tracked) ranks[s, v] <- sum(dm[, v] < truth[[v]])   # rank in 0..L_thin

  ds <- fit$diagnostic_summary()
  diag$max_rhat[s]    <- max(fit$summary(tracked, "rhat")$rhat, na.rm = TRUE)
  diag$divergences[s] <- sum(ds$num_divergent)
  cat(sprintf("SBC sim %4d / %d  (max rhat %.2f, div %d)\n", s, N_sim, diag$max_rhat[s], diag$divergences[s]))

  if (s %% 25 == 0 || s == N_sim)               # checkpoint
    saveRDS(list(ranks = ranks, diag = diag, L = L_thin, done = s), ckpt)
}

saveRDS(list(ranks = ranks, diag = diag, L = L_thin), "robustness_checks/results/sbc_ranks.rds")
if (file.exists(ckpt)) file.remove(ckpt)

## ----rank-histograms----------------------------------------------------------
L_rank <- L_thin
ok     <- diag$max_rhat < 1.1                       # drop non-converged (invalid ranks)
N_ok   <- sum(ok)
long   <- as.data.frame(ranks[ok, , drop = FALSE]) |>
  pivot_longer(everything(), names_to = "param", values_to = "rank")
long$param <- factor(long$param, levels = tracked)
band <- data.frame(lo = qbinom(0.005, N_ok, 1 / N_bins), hi = qbinom(0.995, N_ok, 1 / N_bins))

p <- ggplot(long, aes(rank)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
            fill = thesis_pal$ci_band, alpha = .55) +
  geom_hline(yintercept = N_ok / N_bins, colour = "grey50", linetype = "dashed") +
  geom_histogram(bins = N_bins, boundary = 0,
                 fill = thesis_pal$hist_fill, colour = thesis_pal$hist_line, linewidth = .25) +
  facet_wrap(~ param, scales = "free_y") +
  labs(title = "Simulation-Based Calibration: rank histograms",
       subtitle = sprintf("%d simulations; shaded band = 99%% interval expected under uniformity (calibration)", N_ok),
       x = sprintf("Rank of true value among %d posterior draws", L_rank), y = "Count") +
  theme_thesis(13)

ggsave("Thesis_final/immagini/sbc_rank_histograms.png", p, width = 11, height = 7.5, dpi = 150, bg = "white")
cat("SBC_DONE  (", N_ok, "of", N_sim, "converged )\n")
