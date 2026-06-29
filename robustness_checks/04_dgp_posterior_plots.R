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
library(ggplot2)
library(dplyr)
library(tidyr)

## ----dgp----------------------------------------------------------------------
# Real-time frame, tree, priors and the generative process (simulate_dgp()).
source("robustness_checks/_dgp.R")
source("robustness_checks/_plot_style.R")   # make_posterior_ci_plot() + theme

## ----model--------------------------------------------------------------------
model <- cmdstan_model("stan_models/base_rw.stan")

## ----simulate-and-fit---------------------------------------------------------
sim <- simulate_dgp(seed = 305)  # representative replicate (full coverage, well-populated clones)
fit <- model$sample(data = sim$stan_data, chains = 4, parallel_chains = 4,
                    iter_warmup = 2000, iter_sampling = 3000, adapt_delta = 0.97,
                    max_treedepth = 13, seed = 1, refresh = 0, show_messages = FALSE)

## ----posterior-df-------------------------------------------------------------
# Long posterior_df (index, label, value, true_value) for make_posterior_ci_plot;
# event times stay in the model's native pseudotime [0, 1].
build_df <- function(var, labels, truth, transform = identity) {
  dr <- as_draws_df(fit$draws(var)) |> as.data.frame()
  cols <- grep(paste0("^", var, "\\["), names(dr), value = TRUE)
  long <- pivot_longer(dr[, cols], everything(), names_to = "v", values_to = "value")
  long$index      <- as.integer(gsub("\\D", "", long$v))
  long$value      <- transform(long$value)
  long$label      <- factor(labels[long$index], levels = labels)
  long$true_value <- truth[long$index]
  long
}
df_times  <- build_df("tevents",    event_labels, sim$true_t)   # pseudotime, no transform
df_growth <- build_df("w_node_out", clone_labels, sim$true_w)

## ----plots--------------------------------------------------------------------
p_times  <- make_posterior_ci_plot(df_times,  "Posterior event times (DGP replicate)", "Pseudotime")
p_growth <- make_posterior_ci_plot(df_growth, "Posterior growth rates (DGP replicate)", "Growth rate")
ggsave("Thesis_final/immagini/dgp_posterior_times.png",  p_times,  width = 12, height = 7, dpi = 150, bg = "white")
ggsave("Thesis_final/immagini/dgp_posterior_growth.png", p_growth, width = 12, height = 7, dpi = 150, bg = "white")
cat("POSTERIOR_DONE\n")
