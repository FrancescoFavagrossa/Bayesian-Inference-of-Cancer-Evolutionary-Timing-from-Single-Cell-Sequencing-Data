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
  // MRCA event times via a Dirichlet simplex on the gaps between consecutive
  // MRCAs (0 -> t_1 -> ... -> t_Nnodes -> 1). Dirichlet(1,...,1) = uniform prior
  // over ORDERED MRCA times in [0,1]. Replaces the nested stick-breaking, whose
  // uniform u compounds into a prior that pushes deep events toward 1 and biases
  // deep clones late. The cumulative construction t_node = sum(mrca_gap[1:node])
  // plus parents[node] < node guarantees t_node >= t_parent.
  simplex[Nnodes + 1] mrca_gap;

  // Within-branch placement fraction for CNA sub-events (MRCA slots unused).
  vector<lower=0, upper=1>[Nevents] u;

  real<lower=1e-9, upper=1> w_bar;
  real<lower=0, upper=2> sigma_log_w;
  // Symmetric (unconstrained) random walk on log growth.
  vector[Nnodes - 1] delta_log_w;

  real<lower=0> omega;
}

transformed parameters {
  real<lower=0> delta_t = Tf_real - t0_real;
  array[Nevents] real<lower=0, upper=1> tevents;
  vector[Nnodes] mrca_time;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;

  vector<lower=0>[Nnodes] terminal_clone_mass;
  simplex[Nnodes] clone_prob_cells;

  log_w_node[1] = log(w_bar);
  for (node in 2:Nnodes) {
    log_w_node[node] =
      log_w_node[parents[node]] + delta_log_w[node - 1];
  }
  w_node = exp(log_w_node);

  // Cumulative Dirichlet gaps -> ordered MRCA times.
  mrca_time[1] = mrca_gap[1];
  for (node in 2:Nnodes) {
    mrca_time[node] = mrca_time[node - 1] + mrca_gap[node];
  }

  for (idx in 1:Nevents) {
    tevents[idx] = 0;
  }

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_node = mrca_time[node];
    real t_parent = (p == 0) ? 0.0 : mrca_time[p];

    tevents[idx0] = t_node;

    for (event in 2:number_of_events[node]) {
      int idx = idx0 + event - 1;
      tevents[idx] = t_parent + u[idx] * (t_node - t_parent);
    }
  }

  for (node in 1:Nnodes) {
    real t_node = tevents[idx_first_event[node]];

    terminal_clone_mass[node] =
      exp(w_node[node] * delta_t * (1.0 - t_node));
  }

  for (child in 2:Nnodes) {
    int p = parents[child];
    real t_child = tevents[idx_first_event[child]];

    terminal_clone_mass[p] -=
      exp(w_node[p] * delta_t * (1.0 - t_child));
  }

  for (node in 1:Nnodes) {
    if (terminal_clone_mass[node] < 1e-12) {
      terminal_clone_mass[node] = 1e-12;
    }
  }

  clone_prob_cells = terminal_clone_mass / sum(terminal_clone_mass);
}

model {
  mrca_gap ~ dirichlet(rep_vector(1.0, Nnodes + 1));

  w_bar ~ lognormal(log(w_bar_prior_mean), w_bar_prior_sd);
  sigma_log_w ~ exponential(sigma_log_w_rate);
  delta_log_w ~ normal(0, sigma_log_w);

  omega ~ lognormal(0, 1);

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
          omega * mu * w_child * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents[idx0] - t_p);

        target += poisson_lpmf(mutations[idx_m] | positive_rate(rate_m));
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];

        real rate_alpha =
          omega * mu * w_parent * delta_t * lengths[idx] * coeff[idx_alpha]
          * (t_cna - t_p);

        real rate_beta =
          omega * mu * w_child * delta_t * lengths[idx] * coeff[idx_beta]
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
  real omega_out = omega;
  real clone_abundance_nb_phi_out = phi;
  vector[Nnodes] w_node_out = w_node;
  vector[Nnodes] log_w_node_out = log_w_node;
  vector[Nnodes] delta_log_w_out;

  array[Nevents] real tevents_prior;
  array[Nnodes] real log_lik_cells;
  vector[Nnodes] clone_prob_cells_out = clone_prob_cells;
  vector[Nnodes] expected_clone_cells_out = Ncells * clone_prob_cells;

  delta_log_w_out[1] = 0;
  for (node in 2:Nnodes) {
    delta_log_w_out[node] = log_w_node[node] - log_w_node[parents[node]];
  }

  for (node in 1:Nnodes) {
    real expected_clone_cells = positive_rate(Ncells * clone_prob_cells[node]);
    log_lik_cells[node] = neg_binomial_2_lpmf(
      Ncells_per_clone[node] | expected_clone_cells, positive_rate(phi)
    );
  }

  // Prior-predictive event times under the Dirichlet prior.
  {
    vector[Nnodes + 1] mrca_gap_prior = dirichlet_rng(rep_vector(1.0, Nnodes + 1));
    array[Nnodes] real mrca_time_prior;
    vector[Nevents] u_prior;

    for (i in 1:Nevents) {
      u_prior[i] = uniform_rng(0, 1);
      tevents_prior[i] = 0;
    }

    mrca_time_prior[1] = mrca_gap_prior[1];
    for (node in 2:Nnodes) {
      mrca_time_prior[node] = mrca_time_prior[node - 1] + mrca_gap_prior[node];
    }

    for (node in 1:Nnodes) {
      int idx0 = idx_first_event[node];
      int p = parents[node];
      real t_node = mrca_time_prior[node];
      real t_parent = (p == 0) ? 0.0 : mrca_time_prior[p];

      tevents_prior[idx0] = t_node;
      for (event in 2:number_of_events[node]) {
        int idx = idx0 + event - 1;
        tevents_prior[idx] = t_parent + u_prior[idx] * (t_node - t_parent);
      }
    }
  }
}
