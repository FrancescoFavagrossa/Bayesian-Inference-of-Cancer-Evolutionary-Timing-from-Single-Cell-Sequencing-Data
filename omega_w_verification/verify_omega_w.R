# ----------------------------------------------------------------------------
# Verification: is the original omega*w model (n growth fixed negative binomial)
# equivalent to the coupled clock (nu = kappa * w^gamma) once the per-clone growth
# is genuinely distinct? Both fitted on the SAME distinct-growth Clone 6 DLP data.
# ----------------------------------------------------------------------------
setwd("/Users/francescofavagrossa/Documents/Projects/Thesis")
suppressMessages({library(cmdstanr); library(posterior); library(tidyverse); library(stringr)})
source("Helpers and Utils/Helpers.R")
out <- "omega_w_verification"
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)

sd <- readRDS("/tmp/c6d_sd.rds")                 # distinct-growth Clone 6 DLP stan_data
elab <- c("MRCA_A","CNA_A_chr12_40Mb","MRCA_B","MRCA_C","CNA_C_hidden_chr5_20Mb",
          "MRCA_D","CNA_D_hidden_chr18_20Mb","MRCA_E","MRCA_F")
te <- c(0,0,0.291,0.602,0.602,0.715,0.715,0.742,0.760)
clab <- paste0("Clone_", c("A","B","C","D","E","F"))
tg <- c(0.0493,0.092,0.0946,0.0969,0.1879,0.1165)
sm <- function(f,v){s<-f$summary(v); s[order(as.integer(str_extract(s$variable,"[0-9]+"))),]}
fit_model <- function(stan, data) cmdstan_model(stan)$sample(
  data=data, chains=4, parallel_chains=4, iter_warmup=4000, iter_sampling=6000,
  adapt_delta=0.99, max_treedepth=14, seed=42, refresh=0, show_messages=FALSE)
score <- function(f, tag){
  s<-sm(f,"tevents"); g<-sm(f,"w_node_out")
  ae<-abs(s$median-te)[3:9]; ag<-abs(g$median-tg)
  cat(sprintf("%-22s | tempi err medio=%.3f max=%.3f | growth err=%.3f | rhat=%.3f\n",
      tag, mean(ae), max(ae), mean(ag), max(c(s$rhat,g$rhat),na.rm=TRUE)))
  tibble(model=tag, t_mean=mean(ae), t_max=max(ae), g_mean=mean(ag),
         rhat=max(c(s$rhat,g$rhat),na.rm=TRUE))
}

# original omega*w model: needs the standard data (no sigma_nu_rate)
sd_orig <- sd; sd_orig$sigma_nu_rate <- NULL
sd_coup <- sd; sd_coup$sigma_nu_rate <- 2

cat("=== Clone 6 (distinct growth) — same data, two models ===\n")
f_orig <- fit_model("stan_models/n growth fixed negative binomial.stan", sd_orig)
f_coup <- fit_model("stan_models/n growth fixed nb clock coupled.stan",   sd_coup)
res <- bind_rows(score(f_orig, "ORIGINALE omega*w"), score(f_coup, "coupled nu=k*w^g"))
write_csv(res, file.path(out, "results", "comparison.csv"))

# per-event detail for the original model
s<-sm(f_orig,"tevents"); g<-sm(f_orig,"w_node_out")
det <- tibble(event=elab, post_median=round(s$median,3), true=te, err=round(abs(s$median-te),3))
write_csv(det, file.path(out, "results", "omega_w_events.csv"))
cat("\n--- ORIGINALE omega*w, per evento ---\n"); print(det, n=9)
cat("\nomega stimato:", round(f_orig$summary("omega_out")$median,3), "\n")

# figures for the original model
getdf <- function(f,var,labs,trues) f$draws(var,format="draws_df") %>% as_tibble() %>%
  select(starts_with(paste0(var,"["))) %>% mutate(draw=row_number()) %>%
  pivot_longer(starts_with(paste0(var,"[")),names_to="parameter",values_to="value") %>%
  mutate(index=as.integer(str_extract(parameter,"[0-9]+")),
         label=factor(labs[index],levels=labs), true_value=trues[index])
ggsave(file.path(out,"figures","omega_w events time.png"),
  make_posterior_ci_plot(getdf(f_orig,"tevents",elab,te),"Posterior event times (original omega*w)","Pseudotime"),
  width=14,height=8,dpi=120,bg="white")
ggsave(file.path(out,"figures","omega_w growth rates.png"),
  make_posterior_ci_plot(getdf(f_orig,"w_node_out",clab,tg),"Posterior growth rates (original omega*w)","Growth rate"),
  width=12,height=7,dpi=120,bg="white")
cat("\nVERIFY_DONE — salvato in", out, "\n")
