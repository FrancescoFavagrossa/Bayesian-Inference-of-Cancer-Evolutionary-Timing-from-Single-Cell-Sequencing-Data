# ----------------------------------------------------------------------------
# Pairwise comparison of clone growth-rate posteriors, per scenario.
#
# Uses the JOINT w_node_out draws saved by save_bias() as *_growth_draws.rds
# (a t-test is NOT valid here: those are MCMC draws, not i.i.d. data, and the
# two clones are correlated in the joint posterior). For every clone pair (i<j),
# from the PAIRED draws  Delta = w_i - w_j, it reports all three diagnostics:
#   1) mean_diff, ci_lo, ci_hi  -> 95% credible interval of Delta
#   2) p_i_gt_j                 -> P(w_i > w_j), the probability of direction
#   3) pct_in_rope              -> fraction of |Delta| <= rope_eps (equivalence)
# plus the TRUE difference, whether the pair is truly equal, and a verdict.
#
# NOTE on the model's monotonicity constraint (delta_log_w >= 0): for a
# parent->child pair Delta is sign-constrained, so read p_i_gt_j together with
# the CI/ROPE, not on its own.
#
# Workflow:
#   # make sure the draws exist (re-run scenarios with the updated save_bias, or
#   # call save_growth_draws(fit, true_growth, scenario) in each open session)
#   source(".../Thesis/Helpers and Utils/growth_compare.R")
#   res <- compare_growth(rope_eps = 0.005)
# ----------------------------------------------------------------------------

suppressMessages({library(dplyr); library(tibble)})

GC_SCENARIO_DIRS <- c(
  "/Users/francescofavagrossa/Desktop/RW SINGLE CELL/Thesis/Clone 4 proportional/results",
  "/Users/francescofavagrossa/Desktop/RW SINGLE CELL/Thesis/Clone 4 DLP/results",
  "/Users/francescofavagrossa/Desktop/RW SINGLE CELL/Thesis/Clone 6 proportional/results",
  "/Users/francescofavagrossa/Desktop/RW SINGLE CELL/Thesis/Clone 6 DLP/results"
)
GC_OUT_DIR <- "/Users/francescofavagrossa/Desktop/RW SINGLE CELL/Thesis/results"

# All pairwise comparisons for one scenario's growth-draws object.
growth_pairwise <- function(gd, rope_eps = 0.005, ci = 0.95) {
  draws  <- gd$draws
  labels <- colnames(draws)
  true   <- gd$true
  a      <- (1 - ci) / 2
  cb     <- utils::combn(seq_along(labels), 2)

  dplyr::bind_rows(lapply(seq_len(ncol(cb)), function(k) {
    i <- cb[1, k]; j <- cb[2, k]
    d  <- draws[, i] - draws[, j]          # paired difference, per draw
    td <- true[i] - true[j]
    tibble::tibble(
      scenario    = gd$scenario,
      pair        = paste0(labels[i], " - ", labels[j]),
      true_diff   = round(td, 4),
      truly_equal = abs(td) < 1e-9,
      mean_diff   = round(mean(d), 4),
      ci_lo       = round(unname(stats::quantile(d, a)), 4),
      ci_hi       = round(unname(stats::quantile(d, 1 - a)), 4),
      p_i_gt_j    = round(mean(d > 0), 3),
      pct_in_rope = round(mean(abs(d) <= rope_eps), 3)
    )
  })) |>
    dplyr::mutate(
      verdict = dplyr::case_when(
        pct_in_rope >= 0.95                   ~ "equivalent",   # ~all of Delta in ROPE
        ci_lo > rope_eps | ci_hi < -rope_eps  ~ "different",    # CI entirely outside ROPE
        TRUE                                  ~ "uncertain"
      )
    )
}

# Run the comparison for every scenario, print and (by default) write a CSV.
compare_growth <- function(dirs = GC_SCENARIO_DIRS, rope_eps = 0.005,
                           out_dir = GC_OUT_DIR, write_csv = TRUE) {
  files <- unlist(lapply(dirs, function(d)
    list.files(d, pattern = "_growth_draws\\.rds$", full.names = TRUE)))
  if (length(files) == 0)
    stop("No *_growth_draws.rds found. Re-run the scenarios with the updated ",
         "save_bias(), or call save_growth_draws(fit, true_growth, scenario) ",
         "in each open session.")

  res <- dplyr::bind_rows(lapply(files, function(f)
    growth_pairwise(readRDS(f), rope_eps = rope_eps)))

  if (isTRUE(write_csv)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(res, file.path(out_dir, "growth_pairwise.csv"), row.names = FALSE)
    message("CSV: ", normalizePath(file.path(out_dir, "growth_pairwise.csv")),
            "   (ROPE eps = ", rope_eps, ")")
  }

  print(res, n = Inf)
  invisible(res)
}
