## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.width  = 8,
  fig.height = 5
)

## ----libraries----------------------------------------------------------------
library(cmdstanr)
library(posterior)

## ----dgp----------------------------------------------------------------------
# Real-time frame, tree, priors and the generative process (simulate_dgp()).
source("robustness_checks/_dgp.R")

## ----model--------------------------------------------------------------------
model <- cmdstan_model("stan_models/base_rw.stan")

## ----fit-and-score------------------------------------------------------------
# Light fit; score per-event time recovery in pseudotime [0, 1] (the model frame).
fit_one <- function(stan_data)
  model$sample(data = stan_data, chains = 2, parallel_chains = 2,
               iter_warmup = 1000, iter_sampling = 1200, adapt_delta = 0.92,
               max_treedepth = 12, seed = 1, refresh = 0, show_messages = FALSE)

score_times <- function(fit, true_t) {
  dr  <- as_draws_matrix(fit$draws("tevents"))
  dr  <- dr[, order(as.integer(gsub("\\D", "", colnames(dr)))), drop = FALSE]   # pseudotime
  med <- apply(dr, 2, median); lo <- apply(dr, 2, quantile, 0.025); hi <- apply(dr, 2, quantile, 0.975)
  data.frame(event = seq_along(true_t),
             bias  = abs(med - true_t),
             in95  = true_t >= lo & true_t <= hi,
             ciw   = hi - lo,
             rhat  = max(fit$summary("tevents", "rhat")$rhat, na.rm = TRUE))
}

## ----sweep-grid---------------------------------------------------------------
N_rep <- 15
sweep_omega_sd   <- c(0.2, 0.35, 0.5, 0.75, 1.0)   # fit prior: omega_prior_sd
sweep_sigma_rate <- c(0.5, 1, 2, 3)                # fit prior: sigma_log_w_rate
sweep_ncells     <- c(300, 1000, 3000, 10000, 30000)
sweep_phi        <- c(10, 30, 100, 300, 1000)      # fit NegBin dispersion (data generated at phi = 100)
rows <- list(); k <- 1

## ----sweep-priors-------------------------------------------------------------
# Fixed replicates; refit each under the alternative fitting prior.
for (r in 1:N_rep) {
  sim <- simulate_dgp(seed = 300 + r)
  for (v in sweep_omega_sd) {
    sd <- sim$stan_data; sd$omega_prior_sd <- v
    d <- score_times(fit_one(sd), sim$true_t); d$sweep <- "omega_sd"; d$value <- v; d$rep <- r
    rows[[k]] <- d; k <- k + 1
  }
  for (v in sweep_sigma_rate) {
    sd <- sim$stan_data; sd$sigma_log_w_rate <- v
    d <- score_times(fit_one(sd), sim$true_t); d$sweep <- "sigma_rate"; d$value <- v; d$rep <- r
    rows[[k]] <- d; k <- k + 1
  }
  for (v in sweep_phi) {
    sd <- sim$stan_data; sd$phi <- v
    d <- score_times(fit_one(sd), sim$true_t); d$sweep <- "phi"; d$value <- v; d$rep <- r
    rows[[k]] <- d; k <- k + 1
  }
  cat(sprintf("priors: replicate %2d / %d\n", r, N_rep))
}

## ----sweep-ncells-------------------------------------------------------------
# Sampling depth: regenerate data with a different number of cells.
for (v in sweep_ncells) {
  for (r in 1:N_rep) {
    sim <- simulate_dgp(seed = 5000 + v + r, Ncells_target = v)
    d <- score_times(fit_one(sim$stan_data), sim$true_t)
    d$sweep <- "Ncells"; d$value <- v; d$rep <- r; rows[[k]] <- d; k <- k + 1
  }
  cat(sprintf("Ncells %5d done\n", v))
}

## ----save---------------------------------------------------------------------
df <- do.call(rbind, rows)
dir.create("robustness_checks/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(df = df, event_labels = event_labels), "robustness_checks/results/sweep3.rds")
cat("SWEEP3_DONE  rows =", nrow(df), "\n")
