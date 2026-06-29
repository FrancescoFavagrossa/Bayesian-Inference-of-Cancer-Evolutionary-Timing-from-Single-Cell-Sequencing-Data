## ----setup, include=FALSE-----------------------------------------------------
# Posterior predictive check on the mutation counts: refit base_rw on a ProCESS
# scenario and compare the observed per-equation mutation numbers with the
# posterior predictive replicates (mutations_rep, drawn in generated quantities).
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.align = "center", fig.width = 9, fig.height = 5)

## ----libraries----------------------------------------------------------------
library(cmdstanr)
library(posterior)
library(ggplot2)
library(scales)
source("robustness_checks/_plot_style.R")

## ----data---------------------------------------------------------------------
sd <- readRDS("robustness_checks/templates/c6d_sd.rds"); sd$sigma_nu_rate <- NULL
sd$w_bar_prior_mean <- 0.06; sd$w_bar_prior_sd <- 0.35; sd$sigma_log_w_rate <- 1
sd$omega_prior_mean <- 1.25; sd$omega_prior_sd <- 0.5; sd$phi <- 100

## ----fit----------------------------------------------------------------------
model <- cmdstan_model("stan_models/base_rw.stan")
fit <- model$sample(data = sd, chains = 4, parallel_chains = 4,
                    iter_warmup = 1500, iter_sampling = 2000, adapt_delta = 0.95,
                    max_treedepth = 13, seed = 1, refresh = 0, show_messages = FALSE)

## ----ppc-summary--------------------------------------------------------------
yrep <- as_draws_matrix(fit$draws("mutations_rep"))
yrep <- yrep[, order(as.integer(gsub("\\D", "", colnames(yrep)))), drop = FALSE]
y    <- sd$mutations
q    <- apply(yrep, 2, quantile, c(.025, .25, .5, .75, .975))
df   <- data.frame(eq = seq_along(y), obs = y,
                   lo95 = q[1, ], lo50 = q[2, ], med = q[3, ], hi50 = q[4, ], hi95 = q[5, ])
# fraction of observed counts inside the 95% predictive interval
cover <- mean(df$obs >= df$lo95 & df$obs <= df$hi95)
cat(sprintf("PPC coverage (obs in 95%% predictive): %d/%d\n", sum(df$obs >= df$lo95 & df$obs <= df$hi95), nrow(df)))

## ----plot---------------------------------------------------------------------
brk <- c(0, 1, 10, 100, 1000)
p <- ggplot(df, aes(x = obs, y = med)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(aes(ymin = lo95, ymax = hi95), width = 0, colour = "#A9C7E8", linewidth = 1.1) +
  geom_errorbar(aes(ymin = lo50, ymax = hi50), width = 0, colour = "#4C78C8", linewidth = 2.4) +
  geom_point(colour = thesis_pal$median, size = 2) +
  scale_x_continuous(trans = pseudo_log_trans(base = 10), breaks = brk) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), breaks = brk) +
  coord_equal() +
  labs(title = "Posterior predictive check: mutation counts",
       subtitle = "Clone 6 (DLP+): dashed = perfect fit, bands = 50/95% predictive intervals",
       x = "Observed mutation count", y = "Posterior predictive mutation count") +
  theme_thesis(13)

ggsave("Thesis_final/immagini/ppc_mutations.png", p, width = 8, height = 7, dpi = 150, bg = "white")
cat("PPC_DONE\n")
