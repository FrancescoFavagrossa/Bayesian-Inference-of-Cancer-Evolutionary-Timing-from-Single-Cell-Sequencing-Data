## ----dgp-frame----------------------------------------------------------------
# Real-time frame. The model works on normalised pseudotime t in [0, 1]; these
# constants map it back to real time:  t_real = t0 + t * (Tf - t0).
Tf  <- 30          # time of sequencing      (real time)
t0  <- 10          # tumour origin           (real time)
mu  <- 3           # background mutation rate (per unit length, per unit real time)
phi <- 100         # NegBin dispersion for clone abundances

## ----dgp-priors---------------------------------------------------------------
# base_rw priors. The ground truth is drawn from these too, so the DGP is
# well-specified (the model's own generative process).
w_bar_prior_mean <- 0.06
w_bar_prior_sd   <- 0.35
sigma_log_w_rate <- 1
omega_prior_mean <- 1.25
omega_prior_sd   <- 0.5

## ----dgp-tree-----------------------------------------------------------------
# Phylogeny: 4 clones, A -> B -> {C, D}; copy-number alterations on C and D.
Nnodes           <- 4
parents          <- c(0, 1, 2, 2)
number_of_events <- c(1, 1, 2, 2)
Nevents          <- 6
idx_first_event  <- c(1, 2, 3, 5)

karyotypes         <- c(0, 0, 0, 1, 0, 1)         # 0 = MRCA, 1 = CNA
lengths            <- c(10, 10, 10, 5, 10, 5)     # genomic length per event
Nequations         <- 8
coeff              <- c(1, 1, 1, 1, 1.5, 1, 1, 1.5)   # allelic-copy coefficients
idx_first_mutation <- c(1, 2, 3, 4, 6, 7)         # first equation of each event

## ----dgp-labels---------------------------------------------------------------
# Readable labels for events (t_N* = MRCA, t_CNA* = CNA) and clones (w_A, w_B...).
event_labels <- {
  lab <- character(Nevents); n <- 0L; k <- 0L
  for (i in seq_len(Nevents))
    if (karyotypes[i] == 0) { n <- n + 1L; lab[i] <- paste0("t_N", n) }
    else                    { k <- k + 1L; lab[i] <- paste0("t_CNA", k) }
  lab
}
clone_labels <- paste0("w_", LETTERS[seq_len(Nnodes)])

## ----dgp-simulator------------------------------------------------------------
# One draw from the base_rw generative process (the ideal DGP):
#   growth : log-random-walk on the tree,  w_node = exp(cumsum of N(0, sigma))
#   times  : stick-breaking  ->  normalised event times in [0, 1]
#   clock  : nu_node = omega * w_node
#   muts   : Poisson( nu * mu * (Tf-t0) * length * coeff * branch-duration );
#            CNA branches are split into pre (alpha) and post (beta) the event
#   counts : NegBin on the terminal-clone-mass simplex
# prior_exact = TRUE draws the truth from the EXACT model prior (truncations
# w_bar in (1e-9,1), sigma in (0,2)) for simulation-based calibration; otherwise
# sigma is clamped to a biologically plausible range to avoid pathological trees.
simulate_dgp <- function(seed = NA, Ncells_target = 3000,
                         prior_exact = FALSE, sigma_clamp = c(0.12, 0.6)) {
  if (!is.na(seed)) set.seed(seed)
  dt <- Tf - t0
  rtrunc <- function(r, lo, hi, ...) repeat { x <- r(1, ...); if (x > lo && x < hi) return(x) }

  ## ground-truth draws from the priors
  if (prior_exact) {
    w_bar <- rtrunc(rlnorm, 1e-9, 1, log(w_bar_prior_mean), w_bar_prior_sd)
    sigma <- rtrunc(rexp,   0,    2, sigma_log_w_rate)
  } else {
    w_bar <- rlnorm(1, log(w_bar_prior_mean), w_bar_prior_sd)
    sigma <- min(max(rexp(1, sigma_log_w_rate), sigma_clamp[1]), sigma_clamp[2])
  }
  omega <- rlnorm(1, log(omega_prior_mean), omega_prior_sd)
  u     <- runif(Nevents)

  ## growth: log-random-walk along the tree
  log_w <- numeric(Nnodes); log_w[1] <- log(w_bar)
  for (i in 2:Nnodes) log_w[i] <- log_w[parents[i]] + rnorm(1, 0, sigma)
  w  <- exp(log_w)
  nu <- omega * w

  ## times: stick-breaking -> normalised [0, 1]
  t <- numeric(Nevents)
  for (node in 1:Nnodes) {
    i0 <- idx_first_event[node]; p <- parents[node]
    tp <- if (p == 0) 0 else t[idx_first_event[p]]
    for (ev in 1:number_of_events[node]) {
      idx <- i0 + ev - 1
      if (p == 0) t[idx] <- if (ev == 1) u[idx] else u[idx] * t[i0]
      else        t[idx] <- if (ev == 1) tp + u[idx]*(1 - tp) else tp + u[idx]*(t[i0] - tp)
    }
  }

  ## mutations: Poisson, CNA events split into alpha (pre) and beta (post)
  m <- integer(Nequations)
  for (node in 1:Nnodes) {
    i0 <- idx_first_event[node]; p <- parents[node]
    tp <- if (p == 0) 0 else t[idx_first_event[p]]
    for (ev in 1:number_of_events[node]) {
      idx <- i0 + ev - 1
      if (ev == 1) {
        eq <- idx_first_mutation[i0]
        m[eq] <- rpois(1, max(nu[node]*mu*dt*lengths[i0]*coeff[eq]*(t[i0]-tp), 1e-9))
      } else {
        ea <- idx_first_mutation[idx]; eb <- ea + 1; tc <- t[idx]
        m[ea] <- rpois(1, max(nu[node]*mu*dt*lengths[idx]*coeff[ea]*(tc - tp),    1e-9))
        m[eb] <- rpois(1, max(nu[node]*mu*dt*lengths[idx]*coeff[eb]*(t[i0] - tc), 1e-9))
      }
    }
  }

  ## abundances: NegBin on the terminal-clone-mass simplex
  M <- numeric(Nnodes)
  for (node in 1:Nnodes) M[node] <- exp(min(w[node]*dt*(1 - t[idx_first_event[node]]), 200))
  for (ch in 2:Nnodes) {
    p <- parents[ch]
    M[p] <- M[p] - exp(min(w[p]*dt*(1 - t[idx_first_event[ch]]), 200))
  }
  M[M < 1e-12] <- 1e-12
  Nc <- rnbinom(Nnodes, size = phi, mu = Ncells_target * (M / sum(M)))

  ## assemble the stan_data list (every field base_rw expects)
  stan_data <- list(
    Tf_real = Tf, t0_real = t0, mu = mu, phi = phi,
    Nnodes = Nnodes, parents = parents, number_of_events = number_of_events,
    idx_first_event = idx_first_event,
    Nevents = Nevents, karyotypes = karyotypes, lengths = lengths,
    idx_first_mutation = idx_first_mutation,
    Nequations = Nequations, coeff = coeff, mutations = as.integer(m),
    Ncells = as.integer(sum(Nc)), Ncells_per_clone = as.integer(Nc),
    w_bar_prior_mean = w_bar_prior_mean, w_bar_prior_sd = w_bar_prior_sd,
    sigma_log_w_rate = sigma_log_w_rate,
    omega_prior_mean = omega_prior_mean, omega_prior_sd = omega_prior_sd)

  list(stan_data = stan_data,
       true_t      = t,                       # normalised event times
       true_t_real = t0 + t * dt,             # event times in real time
       true_w      = w, true_omega = omega, true_sigma = sigma, true_w_bar = w_bar)
}

## ----dgp-scale-helper---------------------------------------------------------
# Map a posterior of normalised times to real time: t_real = t0 + t * (Tf - t0).
to_real_time <- function(t_norm) t0 + t_norm * (Tf - t0)
