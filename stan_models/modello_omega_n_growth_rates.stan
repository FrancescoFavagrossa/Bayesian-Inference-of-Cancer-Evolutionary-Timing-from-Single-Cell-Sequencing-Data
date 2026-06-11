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
  vector<lower=0, upper=1>[Nevents] u;

  // This is the scalar w from the working model.
  real<lower=1e-9, upper=1> w_bar;

  // Random-walk step size for clone-specific growth rates on the log scale.
  real<lower=0, upper=2> sigma_log_w;
  vector[Nnodes - 1] z_w;
}

transformed parameters {
  real<lower=0> delta_t = Tf_real - t0_real;
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[Nnodes] log_w_node;
  vector<lower=0>[Nnodes] w_node;

  log_w_node[1] = log(w_bar);
  for (node in 2:Nnodes) {
    log_w_node[node] =
      log_w_node[parents[node]] + sigma_log_w * z_w[node - 1];
  }
  w_node = exp(log_w_node);

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
}

model {
  log(w_bar) ~ normal(log(w_bar_prior_mean), w_bar_prior_sd);

  // Keeps the random walk close to constant growth unless mutations demand otherwise.
  sigma_log_w ~ exponential(sigma_log_w_rate);
  z_w ~ std_normal();

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
          mu * w_child * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents[idx0] - t_p);

        target += poisson_lpmf(mutations[idx_m] | rate_m + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];

        real rate_alpha =
          mu * w_parent * delta_t * lengths[idx] * coeff[idx_alpha]
          * (t_cna - t_p);

        real rate_beta =
          mu * w_child * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents[idx0] - t_cna);

        target += poisson_lpmf(mutations[idx_alpha] | rate_alpha + 1e-9);
        target += poisson_lpmf(mutations[idx_beta] | rate_beta + 1e-9);
      }
    }
  }

  {
    vector[Nnodes] log_clone_mass;
    vector[Nnodes] log_clone_prob;
    real log_total_cells = log(Ncells + 1e-9);

    for (node in 1:Nnodes) {
      int idx_endpoint = idx_first_event[node];
      log_clone_mass[node] =
        w_node[node] * delta_t * (1.0 - tevents[idx_endpoint]);
    }

    log_clone_prob = log_softmax(log_clone_mass);

    for (node in 1:Nnodes) {
      target += neg_binomial_2_log_lpmf(
        Ncells_per_clone[node] | log_total_cells + log_clone_prob[node], phi
      );
    }
  }
}

generated quantities {
  real w = w_bar;
  vector[Nnodes] w_node_out = w_node;
  vector[Nnodes] log_w_node_out = log_w_node;
  vector[Nnodes] delta_log_w_out;

  array[Nevents] real tevents_prior;
  array[Nequations] int mutations_prior;
  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;
  array[Nnodes] int Ncells_ppcheck;
  array[Nnodes] real log_lik_cells;
  vector[Nnodes] clone_prob_cells_out;
  vector[Nnodes] expected_clone_cells_out;
  vector<lower=0, upper=1>[Nevents] u_prior;

  for (i in 1:Nevents) {
    u_prior[i] = uniform_rng(0, 1);
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
    Ncells_ppcheck[node] = 0;
    log_lik_cells[node] = 0;
    clone_prob_cells_out[node] = 0;
    expected_clone_cells_out[node] = 0;
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
          mu * w_child * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents_prior[idx0] - t_p_prior);
        real r_m_post =
          mu * w_child * delta_t * lengths[idx0] * coeff[idx_m]
          * (tevents[idx0] - t_p_post);

        mutations_prior[idx_m] = poisson_rng(r_m_prior + 1e-9);
        mutations_ppcheck[idx_m] = poisson_rng(r_m_post + 1e-9);
        log_lik[idx_m] = poisson_lpmf(mutations[idx_m] | r_m_post + 1e-9);
      } else {
        int idx_alpha = idx_first_mutation[idx];
        int idx_beta = idx_first_mutation[idx] + 1;
        real r_a_prior =
          mu * w_parent * delta_t * lengths[idx] * coeff[idx_alpha]
          * (tevents_prior[idx] - t_p_prior);
        real r_b_prior =
          mu * w_child * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents_prior[idx0] - tevents_prior[idx]);
        real r_a_post =
          mu * w_parent * delta_t * lengths[idx] * coeff[idx_alpha]
          * (tevents[idx] - t_p_post);
        real r_b_post =
          mu * w_child * delta_t * lengths[idx] * coeff[idx_beta]
          * (tevents[idx0] - tevents[idx]);

        mutations_prior[idx_alpha] = poisson_rng(r_a_prior + 1e-9);
        mutations_prior[idx_beta] = poisson_rng(r_b_prior + 1e-9);
        mutations_ppcheck[idx_alpha] = poisson_rng(r_a_post + 1e-9);
        mutations_ppcheck[idx_beta] = poisson_rng(r_b_post + 1e-9);

        log_lik[idx_alpha] =
          poisson_lpmf(mutations[idx_alpha] | r_a_post + 1e-9);
        log_lik[idx_beta] =
          poisson_lpmf(mutations[idx_beta] | r_b_post + 1e-9);
      }
    }
  }

  {
    vector[Nnodes] log_clone_mass;
    vector[Nnodes] log_clone_prob;
    real log_total_cells = log(Ncells + 1e-9);

    for (node in 1:Nnodes) {
      int idx_endpoint = idx_first_event[node];
      log_clone_mass[node] =
        w_node[node] * delta_t * (1.0 - tevents[idx_endpoint]);
    }

    log_clone_prob = log_softmax(log_clone_mass);
    clone_prob_cells_out = exp(log_clone_prob);

    for (node in 1:Nnodes) {
      real log_expected_cells = log_total_cells + log_clone_prob[node];
      expected_clone_cells_out[node] = exp(log_expected_cells);
      Ncells_ppcheck[node] = neg_binomial_2_log_rng(log_expected_cells, phi);
      log_lik_cells[node] =
        neg_binomial_2_log_lpmf(
          Ncells_per_clone[node] | log_expected_cells, phi
        );
    }
  }
}
