functions {
  real masked_dirichlet_multinomial_lpmf(
    array[] int y,
    vector alpha,
    array[] int use_cell_count
  ) {
    real alpha_sum = 0;
    int y_sum = 0;
    real lp = 0;

    for (k in 1:size(y)) {
      if (use_cell_count[k] == 1) {
        alpha_sum += alpha[k];
        y_sum += y[k];
      }
    }

    if (y_sum == 0) {
      return 0;
    }

    lp += lgamma(y_sum + 1);
    lp += lgamma(alpha_sum) - lgamma(y_sum + alpha_sum);

    for (k in 1:size(y)) {
      if (use_cell_count[k] == 1) {
        lp += lgamma(y[k] + alpha[k]) - lgamma(alpha[k]) - lgamma(y[k] + 1);
      }
    }

    return lp;
  }
}

data {
  real<lower=0> mu;

  int<lower=1> Nnodes;
  array[Nnodes] int<lower=0, upper=Nnodes> parents;

  int<lower=1> Nevents;
  array[Nnodes] int<lower=1> number_of_events;
  array[Nnodes] int<lower=1, upper=Nevents> idx_first_event;

  array[Nevents] int<lower=0, upper=1> karyotypes;
  array[Nevents] real<lower=0> lengths;

  int<lower=1> Nequations;
  array[Nevents] int<lower=1, upper=Nequations> idx_first_mutation;
  array[Nequations] real<lower=0> coeff;
  array[Nequations] int<lower=0> mutations;

  // Uniform(0, 1) event priors are obtained with alpha = beta = 1.
  array[Nevents] real<lower=0> u_alpha;
  array[Nevents] real<lower=0> u_beta;

  array[Nnodes] int<lower=0> Ncells_per_clone;
  array[Nnodes] int<lower=0, upper=1> use_cell_count;
  real<lower=0> phi;

  // W is the growth rate on the same time scale used by tevents.
  // With normalized times, pass W = real_growth_rate * (Tf - t0).
  real log_W_mean;
  real<lower=0> log_W_sd;

  // Usually pass zeros here. Non-zero values encode a known expected drift.
  vector[Nnodes - 1] delta_log_w_mean;

  // Half-normal prior scale for the random-walk step sd.
  real<lower=0> tau_scale;

  // Extra clone sampling noise in the cell-count likelihood.
  real<lower=0> sampling_bias_sd;
}

transformed data {
  int n_used_cell_counts = 0;

  if (parents[1] != 0) {
    reject("parents[1] must be 0 for the root node.");
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
      reject("Event index range exceeds Nevents for node: ", node);
    }

    if (node < Nnodes) {
      int next_expected = idx_first_event[node] + number_of_events[node];
      if (idx_first_event[node + 1] != next_expected) {
        reject("Events must be contiguous and ordered by node; bad node: ", node + 1);
      }
    }

    n_used_cell_counts += use_cell_count[node];
  }

  if (n_used_cell_counts == 0) {
    reject("At least one clone must have use_cell_count = 1.");
  }

  for (event in 1:Nevents) {
    int n_equations = karyotypes[event] == 1 ? 2 : 1;
    if (idx_first_mutation[event] + n_equations - 1 > Nequations) {
      reject("Mutation index range exceeds Nequations for event: ", event);
    }
  }
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;

  real<lower=-10, upper=10> log_W0;
  real<lower=0, upper=2> tau_log_w;
  vector[Nnodes - 1] z_log_w;

  vector[Nnodes] sampling_bias_raw;
}

transformed parameters {
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;
  vector[Nnodes] log_clone_mass;
  vector[Nnodes] sampling_bias;
  vector<lower=0, upper=1>[Nnodes] clone_prob_raw;
  vector<lower=0, upper=1>[Nnodes] clone_prob_cells;

  log_w_node[1] = log_W0;
  for (node in 2:Nnodes) {
    int p = parents[node];
    log_w_node[node] =
      log_w_node[p] + delta_log_w_mean[node - 1] + tau_log_w * z_log_w[node - 1];
  }

  w_node = exp(log_w_node);

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    if (p == 0) {
      for (event in 1:number_of_events[node]) {
        int idx = idx0 + event - 1;

        if (event == 1) {
          tevents[idx] = u[idx];
        } else {
          tevents[idx] = u[idx] * tevents[idx0];
        }
      }
    } else {
      real t_parent = tevents[idx_first_event[p]];

      for (event in 1:number_of_events[node]) {
        int idx = idx0 + event - 1;

        if (event == 1) {
          tevents[idx] = t_parent + u[idx] * (1.0 - t_parent);
        } else {
          tevents[idx] = t_parent + u[idx] * (tevents[idx0] - t_parent);
        }
      }
    }
  }

  sampling_bias = sampling_bias_sd * (
    sampling_bias_raw - rep_vector(mean(sampling_bias_raw), Nnodes)
  );

  for (node in 1:Nnodes) {
    int idx_endpoint = idx_first_event[node];
    log_clone_mass[node] = w_node[node] * (1.0 - tevents[idx_endpoint])
      + sampling_bias[node];
  }

  clone_prob_raw = softmax(log_clone_mass);

  {
    real used_prob_sum = 0;

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        used_prob_sum += clone_prob_raw[node];
      }
    }

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        clone_prob_cells[node] = clone_prob_raw[node] / used_prob_sum;
      } else {
        clone_prob_cells[node] = 0;
      }
    }
  }
}

model {
  for (event in 1:Nevents) {
    u[event] ~ beta(u_alpha[event], u_beta[event]);
  }

  log_W0 ~ normal(log_W_mean, log_W_sd);
  tau_log_w ~ normal(0, tau_scale);
  z_log_w ~ std_normal();
  sampling_bias_raw ~ std_normal();

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    real t_parent = p == 0 ? 0.0 : tevents[idx_first_event[p]];
    real t_node = tevents[idx0];
    real w_branch = w_node[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];
        real rate_m = mu * w_branch * lengths[idx0] * coeff[idx_m]
          * (t_node - t_parent);

        target += poisson_lpmf(mutations[idx_m] | rate_m + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];

        real rate_alpha = mu * w_branch * lengths[idx] * coeff[idx_alpha]
          * (t_cna - t_parent);
        real rate_beta = mu * w_branch * lengths[idx] * coeff[idx_beta]
          * (t_node - t_cna);

        target += poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        target += poisson_lpmf(mutations[idx_beta] | rate_beta + 1e-9);
      }
    }
  }

  {
    vector[Nnodes] dirichlet_concentration =
      phi * clone_prob_cells + rep_vector(1e-9, Nnodes);

    target += masked_dirichlet_multinomial_lpmf(
      Ncells_per_clone | dirichlet_concentration, use_cell_count
    );
  }
}

generated quantities {
  vector[Nnodes] w_node_out = w_node;
  vector[Nnodes] log_w_node_out = log_w_node;
  vector[Nnodes] delta_log_w_out;
  real tau_log_w_out = tau_log_w;

  vector[Nnodes] sampling_bias_out = sampling_bias;
  vector[Nnodes] clone_prob_raw_out = clone_prob_raw;
  vector[Nnodes] clone_prob_cells_out = clone_prob_cells;

  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;

  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;

  delta_log_w_out[1] = 0;
  for (node in 2:Nnodes) {
    delta_log_w_out[node] = log_w_node[node] - log_w_node[parents[node]];
  }

  for (q in 1:Nequations) {
    mutations_ppcheck[q] = 0;
    log_lik[q] = 0;
  }

  {
    vector[Nnodes] dirichlet_concentration =
      phi * clone_prob_cells + rep_vector(1e-9, Nnodes);
    vector[Nnodes] sampled_prob = dirichlet_rng(dirichlet_concentration);
    int Ncells_total = 0;
    real lp_cells = masked_dirichlet_multinomial_lpmf(
      Ncells_per_clone | dirichlet_concentration, use_cell_count
    );

    for (node in 1:Nnodes) {
      Ncells_ppcheck[node] = 0;
      log_lik_cells[node] = 0;
      if (use_cell_count[node] == 1) {
        Ncells_total += Ncells_per_clone[node];
      }
    }

    Ncells_ppcheck = multinomial_rng(sampled_prob, Ncells_total);

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        log_lik_cells[node] = lp_cells;
      }
    }
  }

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    real t_parent = p == 0 ? 0.0 : tevents[idx_first_event[p]];
    real t_node = tevents[idx0];
    real w_branch = w_node[node];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];
        real rate_m = mu * w_branch * lengths[idx0] * coeff[idx_m]
          * (t_node - t_parent);

        mutations_ppcheck[idx_m] = poisson_rng(rate_m + 1e-9);
        log_lik[idx_m] = poisson_lpmf(mutations[idx_m] | rate_m + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];

        real rate_alpha = mu * w_branch * lengths[idx] * coeff[idx_alpha]
          * (t_cna - t_parent);
        real rate_beta = mu * w_branch * lengths[idx] * coeff[idx_beta]
          * (t_node - t_cna);

        mutations_ppcheck[idx_alpha] = poisson_rng(rate_alpha + 1e-9);
        mutations_ppcheck[idx_beta] = poisson_rng(rate_beta + 1e-9);

        log_lik[idx_alpha] = poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        log_lik[idx_beta] = poisson_lpmf(mutations[idx_beta] | rate_beta + 1e-9);
      }
    }
  }
}
