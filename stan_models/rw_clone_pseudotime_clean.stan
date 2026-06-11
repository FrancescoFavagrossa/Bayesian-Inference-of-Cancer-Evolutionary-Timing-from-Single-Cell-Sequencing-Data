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

  real log_W0;
  vector<lower=0>[Nnodes - 1] delta_log_w; // positive log-growth increments
}

transformed parameters {
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;

  for (node in 1:Nnodes) {
    int p = parents[node];

    if (p == 0) {
      log_w_node[node] = log_W0; // Se p == 0, root dell'albero --> growth rate base
    } else {
      log_w_node[node] = log_w_node[p] + delta_log_w[node - 1]; // child clones can only increase growth rate
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
} // Ciclo for non ho cambiato nulla 

model {
  // u has an implicit Uniform(0,1) prior from its bounds

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

  /*
    Clone cell-count likelihood.

    Growth dynamics define relative clone masses. The softmax converts
    them to expected clone fractions, then the negative binomial allows
    overdispersed counts around the observed sequencing depth.
  */
  {
    int Ncells_total = 0;
    vector[Nnodes] log_N_per_clone;
    real log_norm = negative_infinity();

    for (node in 1:Nnodes) {
      int idx_endpoint = idx_first_event[node];
      log_N_per_clone[node] = w_node[node] * (1.0 - tevents[idx_endpoint]);

      if (use_cell_count[node] == 1) {
        Ncells_total += Ncells_per_clone[node];
        log_norm = log_sum_exp(log_norm, log_N_per_clone[node]);
      }
    }

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        real log_expected_cells = log(Ncells_total)
          + log_N_per_clone[node]
          - log_norm;

        target += neg_binomial_2_log_lpmf(
          Ncells_per_clone[node] | log_expected_cells, phi
        );
      }
    }
  }
}

generated quantities {
  vector[Nnodes] w_node_out;
  vector[Nnodes] log_w_node_out;
  vector[Nnodes] delta_log_w_out;

  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;

  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;

  w_node_out = w_node;
  log_w_node_out = log_w_node;
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
    vector[Nnodes] log_N_per_clone;
    real log_norm = negative_infinity();

    for (node in 1:Nnodes) {
      int idx_endpoint = idx_first_event[node];
      log_N_per_clone[node] = w_node[node] * (1.0 - tevents[idx_endpoint]);

      if (use_cell_count[node] == 1) {
        Ncells_total += Ncells_per_clone[node];
        log_norm = log_sum_exp(log_norm, log_N_per_clone[node]);
      }
    }

    for (node in 1:Nnodes) {
      if (use_cell_count[node] == 1) {
        real log_expected_cells = log(Ncells_total)
          + log_N_per_clone[node]
          - log_norm;

        Ncells_ppcheck[node] = neg_binomial_2_log_rng(log_expected_cells, phi);
        log_lik_cells[node] = neg_binomial_2_log_lpmf(
          Ncells_per_clone[node] | log_expected_cells, phi
        );
      }
    }
  }

  // Posterior predictive mutation counts
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
