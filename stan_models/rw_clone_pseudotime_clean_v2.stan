data {
  real<lower=0> mu;
  int Nnodes;
  array[Nnodes] int<lower=0, upper=Nnodes> parents;

  int Nevents;
  array[Nnodes] int number_of_events;
  array[Nnodes] int idx_first_event;

  array[Nevents] int karyotypes;
  array[Nevents] real lengths;

  int Nequations;
  array[Nevents] int idx_first_mutation;
  array[Nequations] real coeff;
  array[Nequations] int mutations;

  array[Nnodes] int Ncells_per_clone;
  array[Nnodes] int<lower=0, upper=1> use_cell_count;
  real<lower=0> phi;

  real log_W_mean;
  real<lower=0> log_W_sd;
  real<lower=0> tau_scale;
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;

  real<lower=-5, upper=5> log_W0;
  vector<lower=-2, upper=2>[Nnodes - 1] delta_log_w;
}

transformed parameters {
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;
  vector[Nnodes] log_clone_mass;
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

  clone_prob = softmax(log_clone_mass);
}

model {
  log_W0 ~ normal(log_W_mean, log_W_sd);
  delta_log_w ~ normal(0, tau_scale);

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
        int idx_beta = idx_first_mutation[idx] + 1;

        real t_cna = tevents[idx];

        real rate_alpha = mu * w_parent * lengths[idx] * coeff[idx_alpha] * (t_cna - t_p);

        real rate_beta = mu * w_child * lengths[idx] * coeff[idx_beta] * (t_node - t_cna);

        target += poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        target += poisson_lpmf(mutations[idx_beta] | rate_beta + 1e-9);
      }
    }
  }

  {
    int Ncells_used = 0;
    int K_used = 0;
    array[Nnodes] int y_used;
    vector[Nnodes] log_prob_used;
    real log_norm = negative_infinity();

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        K_used += 1;
        y_used[K_used] = Ncells_per_clone[node];
        Ncells_used += Ncells_per_clone[node];
        log_prob_used[K_used] = log(clone_prob[node]);
        log_norm = log_sum_exp(log_norm, log_prob_used[K_used]);
      }
    }

    target += lgamma(Ncells_used + 1);
    for (k in 1:K_used) {
      target += y_used[k] * (log_prob_used[k] - log_norm) - lgamma(y_used[k] + 1);
    }
  }
}

generated quantities {
  vector[Nnodes] w_node_out;
  vector[Nnodes] log_w_node_out;
  vector[Nnodes] delta_log_w_out;
  vector[Nnodes] clone_prob_out;

  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;

  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;

  w_node_out = w_node;
  log_w_node_out = log_w_node;
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
    log_lik_cells[node] = 0;
  }

  {
    int Ncells_total = 0;
    vector[Nnodes] prob_used = rep_vector(0, Nnodes);
    real prob_norm = 0;
    real lp_cells = 0;

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        Ncells_total += Ncells_per_clone[node];
        prob_norm += clone_prob[node];
      }
    }

    lp_cells += lgamma(Ncells_total + 1);

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        prob_used[node] = clone_prob[node] / prob_norm;
        lp_cells += Ncells_per_clone[node] * log(prob_used[node])
          - lgamma(Ncells_per_clone[node] + 1);
        log_lik_cells[node] = lp_cells;
      }
    }

    Ncells_ppcheck = multinomial_rng(prob_used, Ncells_total);
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
        int idx_beta = idx_first_mutation[idx] + 1;

        real t_cna = tevents[idx];

        real rate_alpha = mu * w_parent * lengths[idx] * coeff[idx_alpha] * (t_cna - t_p);

        real rate_beta = mu * w_child * lengths[idx] * coeff[idx_beta] * (t_node - t_cna);

        mutations_ppcheck[idx_alpha] = poisson_rng(rate_alpha + 1e-9);
        mutations_ppcheck[idx_beta] = poisson_rng(rate_beta + 1e-9);

        log_lik[idx_alpha] = poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        log_lik[idx_beta] = poisson_lpmf(mutations[idx_beta] | rate_beta + 1e-9);
      }
    }
  }
}
