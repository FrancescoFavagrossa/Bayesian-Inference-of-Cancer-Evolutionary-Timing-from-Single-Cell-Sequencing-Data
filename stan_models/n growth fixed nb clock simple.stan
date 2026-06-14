functions {
  real positive_rate(real rate) {
    return fmax(rate, 1e-9);
  }
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
  array[Nevents] int<lower=0, upper=1> karyotypes;
  array[Nevents] real<lower=0> lengths;
  array[Nevents] int<lower=1> idx_first_mutation;

  int<lower=1> Nequations;
  array[Nequations] real<lower=0> coeff;
  array[Nequations] int<lower=0> mutations;

  int<lower=0> Ncells;
  array[Nnodes] int<lower=0> Ncells_per_clone;

  real<lower=0> w_bar_prior_mean;
  real<lower=0> w_bar_prior_sd;
  real<lower=0> sigma_log_w_rate;
  real<lower=0> sigma_nu_rate;   // exponential rate on sigma_nu (clock spread); higher = tighter
}

transformed data {
  if (Tf_real <= t0_real) {
    reject("Tf_real must be greater than t0_real.");
  }

  if (parents[1] != 0) {
    reject("parents[1] must be 0.");
  }

  for (node in 2:Nnodes) {
    if (parents[node] < 1 || parents[node] >= node) {
      reject("parents must reference an earlier node; bad node: ", node);
    }
  }

  if (idx_first_event[1] != 1) {
    reject("idx_first_event[1] must be 1.");
  }

  for (node in 1:Nnodes) {
    if (idx_first_event[node] + number_of_events[node] - 1 > Nevents) {
      reject("Event range exceeds Nevents for node: ", node);
    }

    if (node < Nnodes) {
      int next_expected = idx_first_event[node] + number_of_events[node];
      if (idx_first_event[node + 1] != next_expected) {
        reject("Events must be contiguous and ordered by node; bad node: ", node + 1);
      }
    }
  }

  for (event in 1:Nevents) {
    int n_equations = karyotypes[event] == 1 ? 2 : 1;
    if (idx_first_mutation[event] + n_equations - 1 > Nequations) {
      reject("Mutation index range exceeds Nequations for event: ", event);
    }
  }

  {
    int Ncells_from_clones = 0;

    for (node in 1:Nnodes) {
      Ncells_from_clones += Ncells_per_clone[node];
    }

    if (Ncells_from_clones != Ncells) {
      reject("Ncells must equal sum(Ncells_per_clone).");
    }
  }
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;

  real<lower=1e-9, upper=1> w_bar;
  real<lower=0, upper=2> sigma_log_w;
  // Symmetric (unconstrained) random walk on log growth. Forcing the increments
  // >= 0 (monotone-increasing fitness) introduces a boundary/folding bias when a
  // true increment is ~0: the increment cannot be returned near zero, so it is
  // pushed positive and the error compounds down the lineage, badly over-
  // estimating deep clones. Allowing negative increments removes that bias.
  vector[Nnodes - 1] delta_log_w;

  // Per-branch molecular clock, coupled to growth: nu_node = kappa * w_node *
  // exp(sigma_nu * eps_node). The clock is taken simply proportional to growth
  // (the old omega*w idea); kappa is the global proportionality constant. The key
  // extra term is the per-clone residual exp(sigma_nu * eps): the genealogical
  // clock is only *approximately* proportional to w (the Yule/sampling depth
  // varies per branch), and rigid omega*w (no residual) pushes that deviation
  // onto the event times, biasing them. The residual absorbs it. (We do not fit
  // a free exponent gamma: gamma=1 + residual matches gamma free to within 0.004
  // on the event-time error, so it is fixed at 1 for simplicity.)
  real log_kappa;            // log of the clock-to-growth constant kappa
  real<lower=0> sigma_nu;    // spread of the per-clone clock residual
  vector[Nnodes] log_nu_raw;
}

transformed parameters {
  real<lower=0> delta_t = Tf_real - t0_real;
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;

  vector<lower=0>[Nnodes] nu_node;

  vector<lower=0>[Nnodes] terminal_clone_mass;
  simplex[Nnodes] clone_prob_cells;

  log_w_node[1] = log(w_bar);
  for (node in 2:Nnodes) {
    log_w_node[node] =
      log_w_node[parents[node]] + delta_log_w[node - 1];
  }
  w_node = exp(log_w_node);

  // Clock proportional to growth (kappa * w) times a per-clone residual. The
  // coupling to w lets the clonal mutation counts (which pin nu_node) inform
  // w_node, breaking the w-t degeneracy that pure abundance leaves; the residual
  // lets nu deviate from exactly kappa*w per branch. fmin guards exp() overflow.
  for (node in 1:Nnodes) {
    nu_node[node] =
      exp(fmin(log_kappa + log_w_node[node] + sigma_nu * log_nu_raw[node], 50));
  }

  for (idx in 1:Nevents) {
    tevents[idx] = 0;
  }

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (p == 0) {
        if (event == 1) {
          tevents[idx] = u[idx];
        } else {
          tevents[idx] = u[idx] * tevents[idx0];
        }
      } else {
        real t_parent = tevents[idx_first_event[p]];

        if (event == 1) {
          tevents[idx] = t_parent + u[idx] * (1.0 - t_parent);
        } else {
          tevents[idx] = t_parent + u[idx] * (tevents[idx0] - t_parent);
        }
      }
    }
  }

  for (node in 1:Nnodes) {
    real t_node = tevents[idx_first_event[node]];

    // Clamp the exponent: during warmup the random walk can push w_node high,
    // overflowing exp() to +inf, and inf - inf below produces NaN. The real
    // exponents are ~10-15, so a cap at 200 only trims runaway warmup draws.
    terminal_clone_mass[node] =
      exp(fmin(w_node[node] * delta_t * (1.0 - t_node), 200));
  }

  for (child in 2:Nnodes) {
    int p = parents[child];
    real t_child = tevents[idx_first_event[child]];

    terminal_clone_mass[p] -=
      exp(fmin(w_node[p] * delta_t * (1.0 - t_child), 200));
  }

  for (node in 1:Nnodes) {
    if (terminal_clone_mass[node] < 1e-12) {
      terminal_clone_mass[node] = 1e-12;
    }
  }

  clone_prob_cells = terminal_clone_mass / sum(terminal_clone_mass);
}

model {
  w_bar ~ lognormal(log(w_bar_prior_mean), w_bar_prior_sd);
  sigma_log_w ~ exponential(sigma_log_w_rate);
  delta_log_w ~ normal(0, sigma_log_w);

  log_kappa ~ normal(0, 3);
  sigma_nu ~ exponential(sigma_nu_rate);    // residual clock spread around the w-coupling
  log_nu_raw ~ std_normal();

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_p = (p == 0) ? 0.0 : tevents[idx_first_event[p]];
    real w_parent = (p == 0) ? w_node[node] : w_node[p];
    real w_child = w_node[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];
        real rate_m =
          nu_node[node] * mu * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents[idx0] - t_p);

        target += poisson_lpmf(mutations[idx_m] | positive_rate(rate_m));
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];

        real rate_alpha =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_alpha]
          * (t_cna - t_p);

        real rate_beta =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents[idx0] - t_cna);

        target += poisson_lpmf(mutations[idx_alpha] | positive_rate(rate_alpha));
        target += poisson_lpmf(mutations[idx_beta] | positive_rate(rate_beta));
      }
    }
  }

  for (node in 1:Nnodes) {
    target += neg_binomial_2_lpmf(
      Ncells_per_clone[node] |
      positive_rate(Ncells * clone_prob_cells[node]),
      positive_rate(phi)
    );
  }
}

generated quantities {
  real w = w_bar;
  real omega_out = exp(log_kappa);
  real kappa_out = exp(log_kappa);
  real sigma_nu_out = sigma_nu;
  vector[Nnodes] nu_node_out = nu_node;
  real clone_abundance_nb_phi_out = phi;
  vector[Nnodes] w_node_out = w_node;
  vector[Nnodes] log_w_node_out = log_w_node;
  vector[Nnodes] delta_log_w_out;

  array[Nevents] real tevents_prior;
  array[Nequations] int mutations_prior;
  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;
  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;
  vector[Nnodes] clone_prob_cells_out = clone_prob_cells;
  vector[Nnodes] expected_clone_cells_out = Ncells * clone_prob_cells;
  vector<lower=0, upper=1>[Nevents] u_prior;
  vector[Nnodes] terminal_clone_mass_out = terminal_clone_mass;

  for (i in 1:Nevents) {
    u_prior[i] = uniform_rng(0, 1);
    tevents_prior[i] = 0;
  }

  for (q in 1:Nequations) {
    mutations_prior[q] = 0;
    mutations_ppcheck[q] = 0;
    log_lik[q] = 0;
  }

  delta_log_w_out[1] = 0;
  for (node in 2:Nnodes) {
    delta_log_w_out[node] = log_w_node[node] - log_w_node[parents[node]];
  }

  for (node in 1:Nnodes) {
    real expected_clone_cells = positive_rate(Ncells * clone_prob_cells[node]);

    log_lik_cells[node] = neg_binomial_2_lpmf(
      Ncells_per_clone[node] | expected_clone_cells, positive_rate(phi)
    );
    Ncells_ppcheck[node] = neg_binomial_2_rng(
      expected_clone_cells, positive_rate(phi)
    );
  }

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (p == 0) {
        tevents_prior[idx] =
          event == 1 ? u_prior[idx] : u_prior[idx] * tevents_prior[idx0];
      } else {
        real tp = tevents_prior[idx_first_event[p]];
        tevents_prior[idx] =
          event == 1
            ? tp + u_prior[idx] * (1.0 - tp)
            : tp + u_prior[idx] * (tevents_prior[idx0] - tp);
      }
    }
  }

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_p_prior = (p == 0) ? 0.0 : tevents_prior[idx_first_event[p]];
    real t_p_post = (p == 0) ? 0.0 : tevents[idx_first_event[p]];
    real w_parent = (p == 0) ? w_node[node] : w_node[p];
    real w_child = w_node[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];
        real r_m_prior =
          nu_node[node] * mu * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents_prior[idx0] - t_p_prior);
        real r_m_post =
          nu_node[node] * mu * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents[idx0] - t_p_post);

        mutations_prior[idx_m] = poisson_rng(positive_rate(r_m_prior));
        mutations_ppcheck[idx_m] = poisson_rng(positive_rate(r_m_post));
        log_lik[idx_m] = poisson_lpmf(
          mutations[idx_m] | positive_rate(r_m_post)
        );
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real r_a_prior =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_alpha]
          * (tevents_prior[idx] - t_p_prior);
        real r_b_prior =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents_prior[idx0] - tevents_prior[idx]);
        real r_a_post =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_alpha]
          * (tevents[idx] - t_p_post);
        real r_b_post =
          nu_node[node] * mu * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents[idx0] - tevents[idx]);

        mutations_prior[idx_alpha] = poisson_rng(positive_rate(r_a_prior));
        mutations_prior[idx_beta] = poisson_rng(positive_rate(r_b_prior));
        mutations_ppcheck[idx_alpha] = poisson_rng(positive_rate(r_a_post));
        mutations_ppcheck[idx_beta] = poisson_rng(positive_rate(r_b_post));

        log_lik[idx_alpha] = poisson_lpmf(
          mutations[idx_alpha] | positive_rate(r_a_post)
        );
        log_lik[idx_beta] = poisson_lpmf(
          mutations[idx_beta] | positive_rate(r_b_post)
        );
      }
    }
  }
}
