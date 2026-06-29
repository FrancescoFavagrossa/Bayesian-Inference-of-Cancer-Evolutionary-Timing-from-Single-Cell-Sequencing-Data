## ----setup, include=FALSE-----------------------------------------------------
# Experiment: does removing the children-subtraction in the NegBin abundance
# (mass = exp growth only, no terminal-mass subtraction) help or hurt inference?
# Compares base_rw (full) vs base_rw_nosub on (a) well-specified DGP replicates and
# (b) the real ProCESS scenarios, scoring recovery of event times and growth rates.
suppressMessages({library(cmdstanr); library(posterior)})
source("robustness_checks/_dgp.R")

## ----models-------------------------------------------------------------------
m_full  <- cmdstan_model("stan_models/base_rw.stan")
m_nosub <- cmdstan_model("stan_models/base_rw_nosub.stan")

fit_it <- function(model, sd)
  model$sample(data = sd, chains = 4, parallel_chains = 4, iter_warmup = 1500,
               iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 13,
               seed = 1, refresh = 0, show_messages = FALSE)

score <- function(f, true_t, true_w) {
  g <- function(var, truth) {
    dr <- as_draws_matrix(f$draws(var))
    dr <- dr[, order(as.integer(gsub("\\D", "", colnames(dr)))), drop = FALSE]
    med <- apply(dr, 2, median); lo <- apply(dr, 2, quantile, .025); hi <- apply(dr, 2, quantile, .975)
    c(bias = mean(abs(med - truth)), cov = mean(truth >= lo & truth <= hi),
      relbias = mean(abs(med - truth) / pmax(abs(truth), 1e-9)))
  }
  ti <- g("tevents", true_t); wi <- g("w_node_out", true_w)
  data.frame(t_bias = unname(ti["bias"]), t_cov = unname(ti["cov"]),
             w_bias = unname(wi["bias"]), w_cov = unname(wi["cov"]),
             w_relbias = unname(wi["relbias"]),
             rhat = max(f$summary(c("tevents", "w_node_out"), "rhat")$rhat, na.rm = TRUE))
}

rows <- list(); k <- 1

## ----dgp (well-specified; data generated WITH subtraction) -------------------
NR <- 12
for (r in 1:NR) {
  S <- simulate_dgp(seed = 400 + r)
  for (mod in c("full", "nosub")) {
    f <- fit_it(if (mod == "full") m_full else m_nosub, S$stan_data)
    d <- score(f, S$true_t, S$true_w); d$model <- mod; d$rep <- r; d$data <- "DGP"
    rows[[k]] <- d; k <- k + 1
  }
  cat(sprintf("DGP rep %2d / %d\n", r, NR))
}

## ----real ProCESS scenarios (truth from saved bias.rds) ----------------------
scen <- list(
  list(d = "robustness_checks/templates/c4d_sd.rds", b = "Clone 4 DLP/results/Clone_4_DLP_bias.rds", nm = "Clone4DLP"),
  list(d = "robustness_checks/templates/c6d_sd.rds", b = "Clone 6 DLP/results/Clone_6_DLP_bias.rds", nm = "Clone6DLP"))
for (sc in scen) {
  sd <- readRDS(sc$d); sd$sigma_nu_rate <- NULL
  sd$w_bar_prior_mean <- 0.06; sd$w_bar_prior_sd <- 0.35; sd$sigma_log_w_rate <- 1
  sd$omega_prior_mean <- 1.25; sd$omega_prior_sd <- 0.5; sd$phi <- 100
  bias <- readRDS(sc$b)
  tt <- bias$true[bias$quantity == "event_time"]
  tw <- bias$true[bias$quantity == "growth_rate"]
  for (mod in c("full", "nosub")) {
    f <- fit_it(if (mod == "full") m_full else m_nosub, sd)
    d <- score(f, tt, tw); d$model <- mod; d$rep <- NA; d$data <- sc$nm
    rows[[k]] <- d; k <- k + 1
  }
  cat("scenario", sc$nm, "done\n")
}

res <- do.call(rbind, rows)
dir.create("robustness_checks/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "robustness_checks/results/nosub_compare.rds")

## ----summary------------------------------------------------------------------
cat("\n=== DGP recovery (mean over converged reps): full vs nosub ===\n")
dd <- subset(res, data == "DGP" & rhat < 1.1)
print(aggregate(cbind(t_bias, t_cov, w_bias, w_cov, w_relbias) ~ model, dd, mean), digits = 3)
cat("\n=== Real ProCESS scenarios: full vs nosub ===\n")
print(subset(res, data != "DGP", select = c(data, model, t_bias, t_cov, w_bias, w_cov, w_relbias, rhat)),
      digits = 3, row.names = FALSE)
cat("\nNOSUB_DONE\n")
