// ============================================================================
// BASE MODEL — random walk on growth, NO per-clone residual.
//   growth:     log w_i = log w_parent + sigma_log_w * z_i      (RW on the tree)
//   clock:      nu_i = omega * w_i                              (rigid, no residual)
//   mutations:  m ~ Poisson(nu_i * mu * dt * l * coeff * dt_interval)   (alpha uses
//               the clone's own w_node, NOT w_parent — alpha/beta both on the stem)
//   abundance:  N_i ~ NegBinomial(N * p_i, phi),  p_i from exp(w_i*dt*(1-t_i))
//   event times: nested-uniform stick-breaking (uniform prior).
//
// CONVERGENCE NOTE: with no residual the global clock factor omega is weakly
// identified on *chain* topologies (the telescoping abundance under-pins the
// w-scale), so a wide omega prior goes multimodal (all params rhat ~1.6). Fix:
// a moderately tight omega prior. Recommended: omega ~ LogNormal(log 1.25, 0.5)
// (omega_prior_mean=1.25, omega_prior_sd=0.5) -> both chain & branching converge
// (rhat <=1.003), full CI coverage, event-time |median-true| ~0.05. Tightening
// omega's *mean* trades chain vs branching accuracy (different sims want different
// omega); 1.25 is the balance point on Clone 4/6 DLP.
// ============================================================================
functions {
  real positive_rate(real rate) { return fmax(rate, 1e-9); }
}

data {
  real Tf_real;
  real t0_real;
  real<lower=0> mu;
  real<lower=0> phi;

  int<lower=1> Nnodes;
  array[Nnodes] int<lower=0, upper=Nnodes> parents;
  array[Nnodes] int<lower=1> number_of_events;
  array[Nnodes] int<lower=1> idx_first_event;

  int<lower=1> Nevents;
  array[Nevents] int<lower=0, upper=1> karyotypes;     // 0 = MRCA, 1 = CNA
  array[Nevents] real<lower=0> lengths;
  array[Nevents] int<lower=1> idx_first_mutation;

  int<lower=1> Nequations;
  array[Nequations] real<lower=0> coeff;               // allelic-copy coefficients
  array[Nequations] int<lower=0> mutations;            // observed mutation counts

  int<lower=0> Ncells;
  array[Nnodes] int<lower=0> Ncells_per_clone;

  real<lower=0> w_bar_prior_mean;
  real<lower=0> w_bar_prior_sd;
  real<lower=0> sigma_log_w_rate;
  real<lower=0> omega_prior_mean;   // Yule genealogical factor (~2)
  real<lower=0> omega_prior_sd;     // tightness; small -> nearly fixed
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;                 // stick-breaking auxiliaries

  real<lower=1e-9, upper=1> w_bar;                     // root growth rate
  real<lower=0, upper=2> sigma_log_w;                  // RW step sd on log growth
  vector[Nnodes - 1] delta_log_w;                      // RW increments (symmetric)

  real<lower=0> omega;                                 // global clock-to-growth factor
}

transformed parameters {
  real<lower=0> delta_t = Tf_real - t0_real;
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;
  vector<lower=0>[Nnodes] nu_node;

  vector<lower=0>[Nnodes] terminal_clone_mass;
  simplex[Nnodes] clone_prob_cells;

  // --- growth: random walk on log w along the tree ---
  log_w_node[1] = log(w_bar);
  for (node in 2:Nnodes) {
    log_w_node[node] = log_w_node[parents[node]] + delta_log_w[node - 1];
  }
  w_node = exp(log_w_node);

  // --- clock: rigid coupling nu = omega * w (no residual) ---
  // computed in log space + clamped to avoid Inf*0 = NaN during warmup
  {
    real log_omega = log(omega);
    for (node in 1:Nnodes) {
      nu_node[node] = exp(fmin(log_omega + log_w_node[node], 50));
    }
  }

  // --- event times: nested-uniform stick-breaking ---
  for (idx in 1:Nevents) tevents[idx] = 0;
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;
      if (p == 0) {
        tevents[idx] = (event == 1) ? u[idx] : u[idx] * tevents[idx0];
      } else {
        real t_parent = tevents[idx_first_event[p]];
        tevents[idx] = (event == 1)
          ? t_parent + u[idx] * (1.0 - t_parent)
          : t_parent + u[idx] * (tevents[idx0] - t_parent);
      }
    }
  }

  // --- abundance: exponential mass minus the children's contribution ---
  for (node in 1:Nnodes) {
    real t_node = tevents[idx_first_event[node]];
    terminal_clone_mass[node] = exp(fmin(w_node[node] * delta_t * (1.0 - t_node), 200));
  }
  for (child in 2:Nnodes) {
    int p = parents[child];
    real t_child = tevents[idx_first_event[child]];
    terminal_clone_mass[p] -= exp(fmin(w_node[p] * delta_t * (1.0 - t_child), 200));
  }
  for (node in 1:Nnodes) {
    if (terminal_clone_mass[node] < 1e-12) terminal_clone_mass[node] = 1e-12;
  }
  clone_prob_cells = terminal_clone_mass / sum(terminal_clone_mass);
}

model {
  // priors
  w_bar ~ lognormal(log(w_bar_prior_mean), w_bar_prior_sd);
  sigma_log_w ~ exponential(sigma_log_w_rate);
  delta_log_w ~ normal(0, sigma_log_w);
  omega ~ lognormal(log(omega_prior_mean), omega_prior_sd);
  // u ~ Uniform(0,1) implicit

  // mutation likelihood (Poisson), rate proportional to nu = omega * w
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_p = (p == 0) ? 0.0 : tevents[idx_first_event[p]];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;
      if (event == 1) {                                // MRCA -> diploid clonal m
        int idx_m = idx_first_mutation[idx0];
        real rate_m = nu_node[node] * mu * delta_t * lengths[idx0]
                      * coeff[idx_m] * (tevents[idx0] - t_p);
        target += poisson_lpmf(mutations[idx_m] | positive_rate(rate_m));
      } else {                                         // CNA -> alpha (pre) + beta (post)
        int idx_a = idx_first_mutation[idx];
        int idx_b = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];
        real rate_a = nu_node[node] * mu * delta_t * lengths[idx]
                      * coeff[idx_a] * (t_cna - t_p);
        real rate_b = nu_node[node] * mu * delta_t * lengths[idx]
                      * coeff[idx_b] * (tevents[idx0] - t_cna);
        target += poisson_lpmf(mutations[idx_a] | positive_rate(rate_a));
        target += poisson_lpmf(mutations[idx_b] | positive_rate(rate_b));
      }
    }
  }

  // abundance likelihood (Negative Binomial)
  for (node in 1:Nnodes) {
    target += neg_binomial_2_lpmf(
      Ncells_per_clone[node] | positive_rate(Ncells * clone_prob_cells[node]),
      positive_rate(phi)
    );
  }
}

generated quantities {
  real w = w_bar;
  real omega_out = omega;
  vector[Nnodes] w_node_out = w_node;
  vector[Nnodes] nu_node_out = nu_node;
  vector[Nnodes] clone_prob_cells_out = clone_prob_cells;

  // per-node RW increment (drop-in with the simple model's output interface)
  vector[Nnodes] delta_log_w_out;
  delta_log_w_out[1] = 0;
  for (node in 2:Nnodes) {
    delta_log_w_out[node] = log_w_node[node] - log_w_node[parents[node]];
  }
}
