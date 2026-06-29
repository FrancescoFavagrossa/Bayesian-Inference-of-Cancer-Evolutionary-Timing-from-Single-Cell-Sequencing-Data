## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo    = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  fig.width  = 8.5,
  fig.height = 8.5
)

## ----libraries----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)

## ----style--------------------------------------------------------------------
source("robustness_checks/_plot_style.R")   # theme_thesis() + palette

## ----data---------------------------------------------------------------------
R   <- readRDS("robustness_checks/results/sweep3.rds")
df  <- R$df
df$Event <- factor(R$event_labels[df$event], levels = R$event_labels)
df  <- filter(df, rhat < 1.15)

## ----aggregate----------------------------------------------------------------
# Mean +/- standard error over replicates, per sweep value and event.
agg <- df %>% group_by(sweep, value, Event) %>%
  summarise(`|Bias|`          = mean(bias), `|Bias|_se`          = sd(bias) / sqrt(n()),
            `Coverage 95% CI` = mean(in95), `Coverage 95% CI_se` = sd(in95) / sqrt(n()),
            `CI width`        = mean(ciw),  `CI width_se`        = sd(ciw)  / sqrt(n()),
            .groups = "drop")
long <- agg %>%
  pivot_longer(c(`|Bias|`, `Coverage 95% CI`, `CI width`), names_to = "metric", values_to = "m") %>%
  rowwise() %>%
  mutate(se = agg[[paste0(metric, "_se")]][which(agg$sweep == sweep & agg$value == value & agg$Event == Event)][1]) %>%
  ungroup()
long$metric <- factor(long$metric, levels = c("|Bias|", "Coverage 95% CI", "CI width"))
hline   <- data.frame(metric = factor("Coverage 95% CI", levels = levels(long$metric)), y = 0.95)
ev_cols <- setNames(thesis_event_cols[seq_along(levels(long$Event))], levels(long$Event))

## ----plot-fn------------------------------------------------------------------
make_sweep_fig <- function(which_sweep, title, subtitle, xlab, log_x = FALSE, file) {
  d <- filter(long, sweep == which_sweep)
  p <- ggplot(d, aes(value, m, colour = Event)) +
    geom_hline(data = hline, aes(yintercept = y), linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0, alpha = .6) +
    geom_line(linewidth = .8) + geom_point(size = 2) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y", strip.position = "left") +
    scale_colour_manual(values = ev_cols) +
    labs(title = title, subtitle = paste0(subtitle, "   ·   |bias| and CI width in pseudotime"),
         x = xlab, y = NULL, colour = "Event") +
    theme_thesis(13) + theme(legend.position = "bottom", strip.placement = "outside")
  if (log_x) p <- p + scale_x_log10(breaks = unique(d$value))
  ggsave(file, p, width = 8.5, height = 8.5, dpi = 150, bg = "white")
}

## ----figures------------------------------------------------------------------
make_sweep_fig("omega_sd", "Sensitivity to the clock-factor prior width",
               "Mean over 15 replicates per value",
               expression("Prior sd of " * omega ~ "   ( smaller = tighter )"),
               FALSE, "Thesis_final/immagini/robust_sweep_omega_sd.png")
make_sweep_fig("sigma_rate", "Sensitivity to the growth random-walk step prior",
               "Mean over 15 replicates per value",
               expression(lambda ~ "   ( larger = tighter )"),
               FALSE, "Thesis_final/immagini/robust_sweep_sigma_rate.png")
make_sweep_fig("Ncells", "Sampling depth: number of cells",
               "Mean over 15 replicates per value",
               "Number of sampled cells   ( log scale )",
               TRUE, "Thesis_final/immagini/robust_sweep_ncells.png")
make_sweep_fig("phi", "Sensitivity to the abundance dispersion",
               "Mean over 15 replicates per value   ( data generated at phi = 100 )",
               expression(varphi ~ "   ( larger = less overdispersion )"),
               TRUE, "Thesis_final/immagini/robust_sweep_phi.png")
cat("PLOTS3_DONE\n")
