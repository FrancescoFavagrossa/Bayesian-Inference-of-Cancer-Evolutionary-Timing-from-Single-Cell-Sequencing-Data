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

  // Per-event Beta prior shape parameters for u.
  // Defaults should normally be u_alpha = u_beta = 1 (Uniform), matching v2.
  // Stronger event-time priors can be supplied from R for sensitivity checks.
  array[Nevents] real<lower=0> u_alpha;
  array[Nevents] real<lower=0> u_beta;

  array[Nnodes] int<lower=0> Ncells_per_clone;
  array[Nnodes] int<lower=0, upper=1> use_cell_count;
  real<lower=0> phi;

  real log_W_mean;
  real<lower=0> log_W_sd;
  vector[Nnodes - 1] delta_log_w_mean;
  real<lower=0> tau_scale;
  real<lower=0> sampling_bias_sd;
}

transformed data {
  if (parents[1] != 0) {
    reject("parents[1] must be 0 for the root node.");
  }

  for (node in 2:Nnodes) {
    if (parents[node] < 1 || parents[node] >= node) {
      reject("parents must reference an earlier non-root node; bad node: ", node);
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

  real<lower=-5, upper=5> log_W0;
  vector<lower=-2, upper=2>[Nnodes - 1] delta_log_w;
  vector[Nnodes] sampling_bias_raw;
}

transformed parameters {
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;
  vector[Nnodes] log_clone_mass;
  vector[Nnodes] sampling_bias;
  vector<lower=0, upper=1>[Nnodes] clone_prob;

  for (node in 1:Nnodes) {
    int p = parents[node];

    if (p == 0) {
      log_w_node[node] = log_W0;
    } else {
      log_w_node[node] = log_w_node[p] + delta_log_w[node - 1];
    }
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

  for (node in 1:Nnodes) {
    int idx_endpoint = idx_first_event[node];
    log_clone_mass[node] = w_node[node] * (1.0 - tevents[idx_endpoint]);
  }

  sampling_bias = sampling_bias_sd * (
    sampling_bias_raw - rep_vector(mean(sampling_bias_raw), Nnodes)
  );

  clone_prob = softmax(log_clone_mass + sampling_bias);
}

model {
  for (i in 1:Nevents) {
    u[i] ~ beta(u_alpha[i], u_beta[i]);
  }

  log_W0 ~ normal(log_W_mean, log_W_sd);
  delta_log_w ~ normal(delta_log_w_mean, tau_scale);
  sampling_bias_raw ~ std_normal();

  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    real t_p = (p == 0)
      ? 0.0 : tevents[idx_first_event[p]];

    real t_node = tevents[idx0];

    real w_child = w_node[node];

    real w_parent = (p == 0)
      ? w_node[node] : w_node[p];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];

        real rate_m = mu * w_child * lengths[idx0] * coeff[idx_m] * (t_node - t_p);

        target += poisson_lpmf(mutations[idx_m] | rate_m + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta  = idx_first_mutation[idx] + 1;

        real t_cna = tevents[idx];

        real rate_alpha = mu * w_parent * lengths[idx] * coeff[idx_alpha] * (t_cna - t_p);
        real rate_beta  = mu * w_child  * lengths[idx] * coeff[idx_beta]  * (t_node - t_cna);

        target += poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        target += poisson_lpmf(mutations[idx_beta]  | rate_beta  + 1e-9);
      }
    }
  }

  {
    vector[Nnodes] dirichlet_concentration = phi * clone_prob
      + rep_vector(1e-9, Nnodes);
    target += masked_dirichlet_multinomial_lpmf(
      Ncells_per_clone | dirichlet_concentration, use_cell_count
    );
  }
}

generated quantities {
  vector[Nnodes] w_node_out;
  vector[Nnodes] log_w_node_out;
  vector[Nnodes] delta_log_w_out;
  vector[Nnodes] sampling_bias_out;
  vector[Nnodes] clone_prob_out;

  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;

  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;

  w_node_out     = w_node;
  log_w_node_out = log_w_node;
  sampling_bias_out = sampling_bias;
  clone_prob_out = clone_prob;
  delta_log_w_out[1] = 0;
  for (node in 2:Nnodes) {
    delta_log_w_out[node] = delta_log_w[node - 1];
  }

  for (q in 1:Nequations) {
    mutations_ppcheck[q] = 0;
    log_lik[q] = 0;
  }

  for (node in 1:Nnodes) {
    Ncells_ppcheck[node] = 0;
    log_lik_cells[node]  = 0;
  }

  {
    int Ncells_total = 0;
    vector[Nnodes] dirichlet_concentration = phi * clone_prob
      + rep_vector(1e-9, Nnodes);
    vector[Nnodes] sampled_prob;
    real lp_cells = masked_dirichlet_multinomial_lpmf(
      Ncells_per_clone | dirichlet_concentration, use_cell_count
    );

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        Ncells_total += Ncells_per_clone[node];
      }
    }

    sampled_prob   = dirichlet_rng(dirichlet_concentration);
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

    real t_p = (p == 0)
      ? 0.0 : tevents[idx_first_event[p]];

    real t_node = tevents[idx0];

    real w_child = w_node[node];

    real w_parent = (p == 0)
      ? w_node[node] : w_node[p];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];

        real rate_m = mu * w_child * lengths[idx0] * coeff[idx_m] * (t_node - t_p);

        mutations_ppcheck[idx_m] = poisson_rng(rate_m + 1e-9);
        log_lik[idx_m] = poisson_lpmf(mutations[idx_m] | rate_m + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta  = idx_first_mutation[idx] + 1;

        real t_cna = tevents[idx];

        real rate_alpha = mu * w_parent * lengths[idx] * coeff[idx_alpha] * (t_cna - t_p);
        real rate_beta  = mu * w_child  * lengths[idx] * coeff[idx_beta]  * (t_node - t_cna);

        mutations_ppcheck[idx_alpha] = poisson_rng(rate_alpha + 1e-9);
        mutations_ppcheck[idx_beta]  = poisson_rng(rate_beta  + 1e-9);

        log_lik[idx_alpha] = poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        log_lik[idx_beta]  = poisson_lpmf(mutations[idx_beta]  | rate_beta  + 1e-9);
      }
    }
  }
}
