# ============================================================================
# Shared plotting style for ALL robustness-check figures.
#
# Anchor: the ProCESS-style posterior credible-interval histogram
# (make_posterior_ci_plot). Every other figure inherits theme_thesis() and the
# same palette so the whole chapter has a single, uniform look.
#
# source("robustness_checks/_plot_style.R") at the top of every plotting script.
# ============================================================================

suppressMessages({library(ggplot2); library(dplyr)})

# ---- palette ---------------------------------------------------------------
thesis_pal <- list(
  ci_band   = "#FFF2A1",    # 95% credible-interval shading (and SBC uniformity band)
  hist_fill = "lightblue",  # histogram bars
  hist_line = "black",
  true      = "#D7301F",    # true value  (dashed red)
  median    = "#1F2937",    # posterior median (solid near-black)
  grid      = "grey88"
)
# categorical palette for per-event / per-clone lines (blue-first, no alarm-red)
thesis_event_cols <- c("#4C78C8", "#E45756", "#54A24B", "#B279A2",
                       "#FF9D45", "#17BECF", "#8C564B", "#E377C2")

# ---- shared theme (identical to the make_posterior_ci_plot theme block) -----
theme_thesis <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 12, color = "grey35", hjust = 0.5),
      plot.margin      = margin(8, 30, 8, 30),
      strip.background = element_rect(fill = "grey90", color = "grey35"),
      strip.text       = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = thesis_pal$grid),
      axis.title       = element_text(size = 14),
      axis.text        = element_text(color = "grey25")
    )
}

# ---- canonical posterior CI histogram --------------------------------------
# posterior_df must have columns: index, label, value, true_value
make_posterior_ci_plot <- function(posterior_df, plot_title, x_label, x_limits = NULL) {
  posterior_ci <- posterior_df %>%
    group_by(index, label) %>%
    summarise(q025 = quantile(value, 0.025),
              q50  = quantile(value, 0.50),
              q975 = quantile(value, 0.975),
              true_value = first(true_value), .groups = "drop")

  # Auto x-range from the 95% CIs and true values (not the raw posterior support):
  # heavy-tailed posteriors otherwise collapse every facet to a spike near zero.
  if (is.null(x_limits)) {
    x_hi <- max(c(posterior_ci$q975, posterior_ci$true_value), na.rm = TRUE) * 1.12
    x_lo <- min(c(posterior_ci$q025, posterior_ci$true_value, 0), na.rm = TRUE)
    x_limits <- c(x_lo, x_hi)
  }
  posterior_df <- posterior_df %>% filter(value >= x_limits[1], value <= x_limits[2])

  p <- ggplot(posterior_df, aes(x = value)) +
    geom_rect(data = posterior_ci,
              aes(xmin = q025, xmax = q975, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = thesis_pal$ci_band, alpha = 0.45) +
    geom_histogram(bins = 40, fill = thesis_pal$hist_fill,
                   color = thesis_pal$hist_line, linewidth = 0.25) +
    geom_vline(data = posterior_ci, aes(xintercept = true_value),
               inherit.aes = FALSE, color = thesis_pal$true,
               linetype = "dashed", linewidth = 1) +
    geom_vline(data = posterior_ci, aes(xintercept = q50),
               inherit.aes = FALSE, color = thesis_pal$median, linewidth = 0.9) +
    facet_wrap(~ label, scales = "free_y", ncol = 2) +
    labs(title = plot_title,
         subtitle = "Shaded regions represent 95% credible intervals; dashed red lines indicate true values; solid black lines denote posterior medians",
         x = x_label, y = "Posterior frequency") +
    theme_thesis(14)

  if (!is.null(x_limits)) p <- p + coord_cartesian(xlim = x_limits)
  p
}
