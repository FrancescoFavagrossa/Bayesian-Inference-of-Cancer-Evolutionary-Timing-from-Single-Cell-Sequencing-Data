# ----------------------------------------------------------------------------
# Bias tables for the clone-pseudotime simulations.
#
# Posterior estimates vs true values, for BOTH event times (pseudotime) and
# clone growth rates, for any scenario. The 4 scenario scripts each leave in the
# session the objects this needs: `fit`, `true_events`, `true_growth`.
#
# Workflow:
#   Each scenario script, at its end, sources this file and calls save_bias(),
#   which writes <scenario folder>/results/<scenario>_bias.rds. After running all
#   four scenarios, combine them:
#     source(".../Thesis/Helpers and Utils/bias_table.R")
#     res <- combine_bias()
#     print(res$detailed, n = Inf)   # per-parameter bias
#     print(res$summary)             # per-scenario / per-quantity summary
#     print(res$wide)                # one row per scenario (the table you want)
# ----------------------------------------------------------------------------

# Each scenario keeps its results in a "results/" subfolder INSIDE its own folder.
SCENARIO_RESULTS_DIRS <- c(
  "/Users/francescofavagrossa/Documents/Projects/Thesis/Clone 4 proportional/results",
  "/Users/francescofavagrossa/Documents/Projects/Thesis/Clone 4 DLP/results",
  "/Users/francescofavagrossa/Documents/Projects/Thesis/Clone 6 proportional/results",
  "/Users/francescofavagrossa/Documents/Projects/Thesis/Clone 6 DLP/results"
)

# Per-parameter posterior summary + bias for one fit / one scenario.
build_bias_table <- function(fit, true_events, true_growth, scenario) {

  one <- function(var, truth, quantity) {
    dm  <- as.matrix(fit$draws(var, format = "draws_matrix"))
    idx <- as.integer(stringr::str_extract(colnames(dm), "\\d+"))
    ord <- order(idx)
    dm  <- dm[, ord, drop = FALSE]
    idx <- idx[ord]

    tibble::tibble(
      scenario    = scenario,
      quantity    = quantity,
      index       = idx,
      label       = truth$label[idx],
      true        = truth$true_value[idx],
      post_mean   = apply(dm, 2, mean),
      post_median = apply(dm, 2, stats::median),
      q025        = apply(dm, 2, function(x) unname(stats::quantile(x, 0.025))),
      q975        = apply(dm, 2, function(x) unname(stats::quantile(x, 0.975)))
    )
  }

  dplyr::bind_rows(
    one("tevents",    true_events, "event_time"),
    one("w_node_out", true_growth, "growth_rate")
  ) |>
    dplyr::mutate(
      bias     = post_median - true,                  # signed bias (posterior median - true)
      abs_bias = abs(bias),
      rel_bias = bias / true,                         # relative bias
      ci_width = q975 - q025,
      covered  = true >= q025 & true <= q975          # is the truth inside the 95% CI?
    )
}

# Compute + save one scenario's table. `dir` defaults to "results" RELATIVE to
# the current working directory, i.e. a results/ folder inside the scenario's own
# folder (each scenario script setwd()s to its folder before calling this).
save_bias <- function(fit, true_events, true_growth, scenario, dir = "results") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  bt <- build_bias_table(fit, true_events, true_growth, scenario)
  file <- file.path(dir, paste0(gsub("[^A-Za-z0-9]+", "_", scenario), "_bias.rds"))
  saveRDS(bt, file)
  message("saved: ", normalizePath(file))
  invisible(bt)
}

# Per-scenario / per-quantity summary of the bias.
bias_summary <- function(detailed) {
  detailed |>
    dplyr::group_by(scenario, quantity) |>
    dplyr::summarise(
      n               = dplyr::n(),
      mean_bias       = mean(bias),
      mean_abs_bias   = mean(abs_bias),
      rmse            = sqrt(mean(bias^2)),
      median_rel_bias = stats::median(rel_bias),
      coverage_95     = mean(covered),
      .groups = "drop"
    )
}

# One row per scenario (event-time bias next to growth-rate bias): the table.
bias_summary_wide <- function(detailed) {
  bias_summary(detailed) |>
    tidyr::pivot_wider(
      id_cols     = scenario,
      names_from  = quantity,
      values_from = c(mean_abs_bias, rmse, median_rel_bias, coverage_95)
    )
}

# Load every saved *_bias.rds from all scenario folders, assemble the tables,
# PRINT summary + wide to the console, and write all three as CSV to `out_dir`.
combine_bias <- function(dirs = SCENARIO_RESULTS_DIRS,
                         out_dir = "/Users/francescofavagrossa/Documents/Projects/Thesis/results",
                         write_csv = TRUE) {
  files <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = "_bias\\.rds$", full.names = TRUE)))
  if (length(files) == 0) stop("No *_bias.rds found -- run save_bias() in each scenario first.")

  detailed <- dplyr::bind_rows(lapply(files, readRDS))
  out <- list(
    detailed = detailed,
    summary  = bias_summary(detailed),
    wide     = bias_summary_wide(detailed)
  )

  if (isTRUE(write_csv)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(out$detailed, file.path(out_dir, "bias_detailed.csv"), row.names = FALSE)
    utils::write.csv(out$summary,  file.path(out_dir, "bias_summary.csv"),  row.names = FALSE)
    utils::write.csv(out$wide,     file.path(out_dir, "bias_wide.csv"),     row.names = FALSE)
    message("CSV scritti in: ", normalizePath(out_dir))
  }

  cat("\n==== per-scenario / quantity summary ====\n")
  print(out$summary, n = Inf)
  cat("\n==== wide (one row per scenario) ====\n")
  print(out$wide, n = Inf)

  invisible(out)
}
