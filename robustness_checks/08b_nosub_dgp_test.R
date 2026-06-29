# Focused experiment: effect of REMOVING the children-subtraction in the NegBin
# abundance, on the ideal DGP only. Data are drawn from the well-specified process
# (which uses terminal/disjoint masses); we fit the same data with the full model
# (base_rw, with subtraction) and the nosub model (clade mass, no subtraction) and
# compare recovery of event times (real scale) and clone growth rates.
suppressMessages({library(cmdstanr); library(posterior)})
source("robustness_checks/_dgp.R")

m_full  <- cmdstan_model("stan_models/base_rw.stan")
m_nosub <- cmdstan_model("stan_models/base_rw_nosub.stan")

fit_it <- function(model, sd)
  model$sample(data = sd, chains = 4, parallel_chains = 4, iter_warmup = 1500,
               iter_sampling = 2000, adapt_delta = 0.95, max_treedepth = 13,
               seed = 1, refresh = 0, show_messages = FALSE)

score <- function(f) {
  g <- function(var, truth, tr = identity) {
    dr <- as_draws_matrix(f$draws(var))
    dr <- tr(dr[, order(as.integer(gsub("\\D", "", colnames(dr)))), drop = FALSE])
    med <- apply(dr, 2, median); lo <- apply(dr, 2, quantile, .025); hi <- apply(dr, 2, quantile, .975)
    list(med = med, lo = lo, hi = hi, truth = truth)
  }
  list(t = g("tevents", S$true_t_real, to_real_time), w = g("w_node_out", S$true_w))
}
metrics <- function(s) c(
  t_mab = mean(abs(s$t$med - s$t$truth)),
  t_cov = mean(s$t$truth >= s$t$lo & s$t$truth <= s$t$hi),
  w_mab = mean(abs(s$w$med - s$w$truth)),
  w_cov = mean(s$w$truth >= s$w$lo & s$w$truth <= s$w$hi),
  w_rel = median(abs(s$w$med - s$w$truth) / pmax(abs(s$w$truth), 1e-9)))

NR <- 12
rows <- list(); k <- 1L
for (r in 1:NR) {
  S <- simulate_dgp(seed = 500 + r)
  for (mod in c("full", "nosub")) {
    f  <- fit_it(if (mod == "full") m_full else m_nosub, S$stan_data)
    rh <- max(f$summary(c("tevents", "w_node_out"), "rhat")$rhat, na.rm = TRUE)
    mt <- metrics(score(f))
    rows[[k]] <- data.frame(rep = r, model = mod, t(mt), rhat = rh); k <- k + 1L
  }
  cat(sprintf("rep %2d / %d done\n", r, NR))
}
res <- do.call(rbind, rows)
dir.create("robustness_checks/results", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "robustness_checks/results/nosub_dgp_test.rds")

ok <- subset(res, rhat < 1.1)
cat("\n=== Ideal DGP: full (with subtraction) vs nosub (no subtraction) ===\n")
cat(sprintf("converged fits: full %d/%d, nosub %d/%d\n",
            sum(ok$model=="full"), NR, sum(ok$model=="nosub"), NR))
agg <- aggregate(cbind(t_mab, t_cov, w_mab, w_cov, w_rel) ~ model, ok, mean)
print(agg, digits = 3, row.names = FALSE)
cat("\n(t_mab = time mean abs bias in REAL time, window =", Tf - t0,
    "; t_cov/w_cov = 95% coverage; w_rel = median growth rel. error)\n")
cat("NOSUB_DGP_DONE\n")
