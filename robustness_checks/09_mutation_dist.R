## ----setup--------------------------------------------------------------------
# Posterior-predictive DISTRIBUTIONS of the mutation counts: per-equation
# histograms of mutations_rep with the observed count marked, for the ideal DGP
# and the ProCESS scenarios whose stan_data is available. Free per-facet scales,
# because mutation counts span 0 to ~2000 across equations.
suppressMessages({library(cmdstanr); library(posterior); library(ggplot2); library(dplyr); library(tidyr)})
source("robustness_checks/_plot_style.R")
source("robustness_checks/_dgp.R")
model <- cmdstan_model("stan_models/base_rw.stan")

## ----helpers------------------------------------------------------------------
fit_sd <- function(sd)
  model$sample(data = sd, chains = 4, parallel_chains = 4, iter_warmup = 1500,
               iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 13,
               seed = 1, refresh = 0, show_messages = FALSE)

# canonical base_rw priors on a ProCESS template (drops the old clock field)
fix_priors <- function(sd) {
  sd$sigma_nu_rate <- NULL
  sd$w_bar_prior_mean <- 0.06; sd$w_bar_prior_sd <- 0.35; sd$sigma_log_w_rate <- 1
  sd$omega_prior_mean <- 1.25; sd$omega_prior_sd <- 0.5; sd$phi <- 100
  sd
}

# long data frame (mutations_rep draws + observed) for plotting
mut_long <- function(fit, obs) {
  yr <- as_draws_matrix(fit$draws("mutations_rep"))
  yr <- yr[, order(as.integer(gsub("\\D", "", colnames(yr)))), drop = FALSE]
  lab <- factor(sprintf("Eq %02d", seq_along(obs)), levels = sprintf("Eq %02d", seq_along(obs)))
  data.frame(index = rep(seq_along(obs), each = nrow(yr)),
             label = rep(lab, each = nrow(yr)),
             value = as.vector(yr),
             true_value = rep(obs, each = nrow(yr)))
}

# per-equation predictive distribution plot (free scales; thesis palette)
mut_plot <- function(df, title) {
  ci <- df |> group_by(index, label) |>
    summarise(q025 = quantile(value, .025), q50 = median(value),
              q975 = quantile(value, .975), true_value = first(true_value), .groups = "drop")
  ggplot(df, aes(value)) +
    geom_rect(data = ci, aes(xmin = q025, xmax = q975, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = thesis_pal$ci_band, alpha = 0.45) +
    geom_histogram(bins = 30, fill = thesis_pal$hist_fill, color = thesis_pal$hist_line, linewidth = 0.2) +
    geom_vline(data = ci, aes(xintercept = true_value), inherit.aes = FALSE,
               color = thesis_pal$true, linetype = "dashed", linewidth = 0.9) +
    geom_vline(data = ci, aes(xintercept = q50), inherit.aes = FALSE,
               color = thesis_pal$median, linewidth = 0.8) +
    facet_wrap(~ label, scales = "free", ncol = 4) +
    labs(title = title,
         subtitle = "dashed red = observed count; solid black = posterior-predictive median; band = 95% interval",
         x = "Mutation count", y = "Posterior-predictive frequency") +
    theme_thesis(13)
}

save_mut <- function(sd, title, file, w = 11, h = 6.5) {
  fit <- fit_sd(sd)
  df  <- mut_long(fit, sd$mutations)
  ggsave(file, mut_plot(df, title), width = w, height = h, dpi = 150, bg = "white")
  cov <- df |> group_by(index) |>
    summarise(lo = quantile(value, .025), hi = quantile(value, .975), obs = first(true_value), .groups = "drop") |>
    summarise(c = mean(obs >= lo & obs <= hi)) |> pull(c)
  cat(sprintf("%-34s obs-in-95%%: %.2f  -> %s\n", title, cov, file))
}

## ----DGP (same replicate as the times/growth figures) ------------------------
sim <- simulate_dgp(seed = 27)
save_mut(sim$stan_data, "Posterior-predictive mutation counts (DGP replicate)",
         "Thesis_final/immagini/dgp_posterior_mutations.png", h = 5.5)

## ----ProCESS scenarios with available stan_data ------------------------------
c4p <- fix_priors(readRDS("Clone 4 proportional/results/Clone_4_prop_stan_diagnostics.rds")$stan_data)
c4d <- fix_priors(readRDS("robustness_checks/templates/c4d_sd.rds"))
c6d <- fix_priors(readRDS("robustness_checks/templates/c6d_sd.rds"))
save_mut(c4p, "Posterior-predictive mutation counts: Clone 4 (prop)", "Thesis_final/immagini/mut_dist_c4_prop.png", h = 5.5)
save_mut(c4d, "Posterior-predictive mutation counts: Clone 4 (DLP+)", "Thesis_final/immagini/mut_dist_c4_dlp.png", h = 5.5)
save_mut(c6d, "Posterior-predictive mutation counts: Clone 6 (DLP+)", "Thesis_final/immagini/mut_dist_c6_dlp.png", h = 7)
cat("MUT_DIST_DONE\n")
