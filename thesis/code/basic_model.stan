// ============================================================================
// BASIC MODEL -- known clone growth rates, only event times inferred.
//   The per-clone growth rates w_node are supplied as data; the model infers the
//   event pseudotimes via nested-uniform stick-breaking.
//   mutations:  m ~ Poisson(w_i * mu * dt * l * coeff * dt_interval)
//   abundance:  N_i ~ NegBinomial(N * p_i, phi),  p_i from exp(w_i*dt*(1-t_i))
//   This is the complete model with the growth random walk and the clock factor
//   removed; with growth fixed the event times are fully identified.
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
  array[Nnodes] real<lower=0> w_node;                  // KNOWN clone growth rates

  int<lower=1> Nevents;
  array[Nevents] int<lower=0, upper=1> karyotypes;     // 0 = MRCA, 1 = CNA
  array[Nevents] real<lower=0> lengths;
  array[Nevents] int<lower=1> idx_first_mutation;

  int<lower=1> Nequations;
  array[Nequations] real<lower=0> coeff;               // allelic-copy coefficients
  array[Nequations] int<lower=0> mutations;            // observed mutation counts

  int<lower=0> Ncells;
  array[Nnodes] int<lower=0> Ncells_per_clone;
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;                 // stick-breaking auxiliaries
}

transformed parameters {
  real<lower=0> delta_t = Tf_real - t0_real;
  array[Nevents] real<lower=0, upper=1> tevents;
  vector<lower=0>[Nnodes] terminal_clone_mass;
  simplex[Nnodes] clone_prob_cells;

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
  // u ~ Uniform(0,1) implicit; the growth rates w_node are known (data)

  // mutation likelihood (Poisson), rate proportional to the known growth w
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_p = (p == 0) ? 0.0 : tevents[idx_first_event[p]];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;
      if (event == 1) {                                // MRCA -> diploid clonal m
        int idx_m = idx_first_mutation[idx0];
        real rate_m = w_node[node] * mu * delta_t * lengths[idx0]
                      * coeff[idx_m] * (tevents[idx0] - t_p);
        target += poisson_lpmf(mutations[idx_m] | positive_rate(rate_m));
      } else {                                         // CNA -> alpha (pre) + beta (post)
        int idx_a = idx_first_mutation[idx];
        int idx_b = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];
        real rate_a = w_node[node] * mu * delta_t * lengths[idx]
                      * coeff[idx_a] * (t_cna - t_p);
        real rate_b = w_node[node] * mu * delta_t * lengths[idx]
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
  // posterior predictive replicated mutation counts (PPC on mutation numbers)
  array[Nequations] int mutations_rep;
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];
    real t_p = (p == 0) ? 0.0 : tevents[idx_first_event[p]];
    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;
      if (event == 1) {
        int idx_m = idx_first_mutation[idx0];
        mutations_rep[idx_m] = poisson_rng(positive_rate(
          w_node[node] * mu * delta_t * lengths[idx0] * coeff[idx_m] * (tevents[idx0] - t_p)));
      } else {
        int idx_a = idx_first_mutation[idx];
        int idx_b = idx_first_mutation[idx] + 1;
        real t_cna = tevents[idx];
        mutations_rep[idx_a] = poisson_rng(positive_rate(
          w_node[node] * mu * delta_t * lengths[idx] * coeff[idx_a] * (t_cna - t_p)));
        mutations_rep[idx_b] = poisson_rng(positive_rate(
          w_node[node] * mu * delta_t * lengths[idx] * coeff[idx_b] * (tevents[idx0] - t_cna)));
      }
    }
  }
}
