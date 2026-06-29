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

## ----recovery-fn--------------------------------------------------------------
# Fit one DGP replicate and score recovery of event times (in real time) and
# clone growth rates: |bias|, 95% credible-interval coverage and CI width.
recover_one <- function(seed) {
  sim <- simulate_dgp(seed = seed)
  fit <- model$sample(data = sim$stan_data, chains = 4, parallel_chains = 4,
                      iter_warmup = 1500, iter_sampling = 2000, adapt_delta = 0.95,
                      max_treedepth = 13, seed = 1, refresh = 0, show_messages = FALSE)

  score <- function(var, truth, transform = identity) {
    dr  <- as_draws_matrix(fit$draws(var))
    idx <- order(as.integer(gsub("\\D", "", colnames(dr))))
    dr  <- transform(dr[, idx, drop = FALSE])
    data.frame(med  = apply(dr, 2, median),
               lo   = apply(dr, 2, quantile, 0.025),
               hi   = apply(dr, 2, quantile, 0.975),
               true = truth)
  }
  # event times: scale posterior from [0,1] to real time before scoring
  ti <- score("tevents", sim$true_t_real, to_real_time); ti$quantity <- "time";   ti$label <- event_labels
  wi <- score("w_node_out", sim$true_w);                 wi$quantity <- "growth"; wi$label <- clone_labels

  out <- rbind(ti, wi)
  out$abs_bias <- abs(out$med - out$true)
  out$covered  <- out$true >= out$lo & out$true <= out$hi
  out$ci_width <- out$hi - out$lo
  out$rhat     <- max(fit$summary(c("tevents", "w_node_out"), "rhat")$rhat, na.rm = TRUE)
  out$seed     <- seed
  out
}

## ----replicates (checkpointed / resumable)------------------------------------
# Checkpoint every 20 replicates; re-launch to resume from the last checkpoint
# (per-replicate seeds are deterministic, so the continuation is identical).
N_rep <- 100
dir.create("robustness_checks/results", showWarnings = FALSE, recursive = TRUE)
ckpt <- "robustness_checks/results/dgp_recovery_ckpt.rds"
rows <- vector("list", N_rep); start <- 1L
if (file.exists(ckpt)) {
  cp <- readRDS(ckpt)
  if (isTRUE(cp$N_rep == N_rep)) {
    rows <- cp$rows; start <- cp$done + 1L
    cat(sprintf("resuming recovery from replicate %d / %d\n", start, N_rep))
  }
}
for (r in start:N_rep) {
  cat(sprintf("replicate %3d / %d\n", r, N_rep))
  rows[[r]] <- recover_one(300 + r)
  if (r %% 20 == 0 || r == N_rep) saveRDS(list(rows = rows, done = r, N_rep = N_rep), ckpt)
}
res <- do.call(rbind, rows)
res <- subset(res, rhat < 1.1)   # drop non-converged replicates

## ----recovery-summary---------------------------------------------------------
agg <- aggregate(cbind(abs_bias, covered, ci_width) ~ quantity, res, FUN = mean)
names(agg) <- c("quantity", "mean_abs_bias", "coverage_95", "mean_ci_width")
print(agg, digits = 3)
res$rel <- res$abs_bias / pmax(abs(res$true), 1e-9)
cat(sprintf("time MAE / window = %.3f\n", mean(res$abs_bias[res$quantity == "time"]) / (Tf - t0)))
cat(sprintf("growth median rel error = %.3f\n", median(res$rel[res$quantity == "growth"])))
cat(sprintf("time coverage = %.3f ; growth coverage = %.3f\n",
            mean(res$covered[res$quantity == "time"]), mean(res$covered[res$quantity == "growth"])))

saveRDS(res, "robustness_checks/results/dgp_recovery.rds")
if (file.exists(ckpt)) file.remove(ckpt)
cat("RECOVERY_DONE  (", length(unique(res$seed)), "converged replicates of", N_rep, ")\n")
