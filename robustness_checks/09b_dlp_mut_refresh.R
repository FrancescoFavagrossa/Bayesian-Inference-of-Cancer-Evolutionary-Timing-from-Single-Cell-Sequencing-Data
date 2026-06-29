# Mutation posterior-predictive figures (DGP + Clone 4/6 DLP), with meaningful
# per-equation names: m_X for an MRCA (clonal) equation, alpha_k / beta_k for the
# pre- and post-CNA equations of the k-th copy-number event. Written into
# Thesis_final/immagini; the equation names are rendered with plotmath (subscripts).
suppressMessages({library(cmdstanr); library(posterior); library(ggplot2); library(dplyr); library(tidyr)})
source("robustness_checks/_plot_style.R")
source("robustness_checks/_dgp.R")
model <- cmdstan_model("stan_models/base_rw.stan")

fit_sd <- function(sd) model$sample(data = sd, chains = 4, parallel_chains = 4, iter_warmup = 1500,
  iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 13, seed = 1, refresh = 0, show_messages = FALSE)

# equation labels from the tree structure (plotmath strings: m[A], alpha[1], beta[1], ...)
build_eq_labels <- function(sd, clones = LETTERS[seq_len(sd$Nnodes)]) {
  labs <- character(sd$Nequations); k <- 0L
  for (node in seq_len(sd$Nnodes)) {
    i0 <- sd$idx_first_event[node]
    for (ev in seq_len(sd$number_of_events[node])) {
      idx <- i0 + ev - 1L
      es  <- sd$idx_first_mutation[idx]
      if (sd$karyotypes[idx] == 0) {
        labs[es] <- paste0("m[", clones[node], "]")
      } else {
        k <- k + 1L
        labs[es]     <- paste0("alpha[", k, "]")
        labs[es + 1] <- paste0("beta[",  k, "]")
      }
    }
  }
  labs
}

mut_long <- function(fit, obs, eqlab) {
  yr <- as_draws_matrix(fit$draws("mutations_rep"))
  yr <- yr[, order(as.integer(gsub("\\D", "", colnames(yr)))), drop = FALSE]
  data.frame(index = rep(seq_along(obs), each = nrow(yr)),
             label = factor(rep(eqlab, each = nrow(yr)), levels = eqlab),
             value = as.vector(yr), true_value = rep(obs, each = nrow(yr)))
}
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
    geom_vline(data = ci, aes(xintercept = q50), inherit.aes = FALSE, color = thesis_pal$median, linewidth = 0.8) +
    facet_wrap(~ label, scales = "free", ncol = 4, labeller = label_parsed) +
    labs(title = title,
         subtitle = "dashed red = observed count; solid black = posterior-predictive median; band = 95% interval",
         x = "Mutation count", y = "Posterior-predictive frequency") + theme_thesis(13)
}
save_mut <- function(sd, title, file, h = 6.5) {
  fit <- fit_sd(sd); df <- mut_long(fit, sd$mutations, build_eq_labels(sd))
  ggsave(file, mut_plot(df, title), width = 11, height = h, dpi = 150, bg = "white")
  cov <- df |> group_by(index) |>
    summarise(lo = quantile(value, .025), hi = quantile(value, .975), obs = first(true_value), .groups = "drop")
  cat(sprintf("%-34s in95: %d/%d\n", title, sum(cov$obs >= cov$lo & cov$obs <= cov$hi), nrow(cov)))
}

## ---- DGP (same replicate as the times/growth figures) ----
sim <- simulate_dgp(seed = 305)
save_mut(sim$stan_data, "Posterior-predictive mutation counts (DGP replicate)",
         "Thesis_final/immagini/dgp_posterior_mutations.png", 5.5)

## ---- ProCESS DLP scenarios (new focal-biopsy stan_data) ----
c4d <- readRDS("Clone 4 DLP/results/Clone_4_DLP_stan_data.rds")
c6d <- readRDS("Clone 6 DLP/results/Clone_6_DLP_stan_data.rds")
save_mut(c4d, "Posterior-predictive mutation counts: Linear chain (4 clones, DLP+)", "Thesis_final/immagini/mut_dist_c4_dlp.png", 5.5)
save_mut(c6d, "Posterior-predictive mutation counts: Branching tree (6 clones, DLP+)", "Thesis_final/immagini/mut_dist_c6_dlp.png", 7)
cat("REFRESH_DONE\n")
