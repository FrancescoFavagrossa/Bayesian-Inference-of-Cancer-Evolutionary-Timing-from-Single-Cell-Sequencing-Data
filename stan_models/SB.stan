functions {
  real integrated_w_piecewise(real a, real b, int K_w, vector w_grid) {
    real out = 0;

    if (b <= a) {
      return 0;
    }

    for (k in 1:K_w) {
      real grid_left = (k - 1.0) / K_w;
      real grid_right = k * 1.0 / K_w;

      real left = fmax(a, grid_left);
      real right = fmin(b, grid_right);

      if (right > left) {
        out += w_grid[k] * (right - left);
      }
    }

    return out;
  }

  matrix gp_exp_quad_cov_grid(int K_w, real alpha, real rho, real jitter) {
    matrix[K_w, K_w] K;

    for (i in 1:K_w) {
      real ti = (i - 0.5) / K_w;

      for (j in i:K_w) {
        real tj = (j - 0.5) / K_w;
        real sqdist = square(ti - tj);

        K[i, j] = square(alpha) * exp(-0.5 * sqdist / square(rho));
        K[j, i] = K[i, j];
      }

      K[i, i] += jitter;
    }

    return K;
  }
}

data {
  /*
    Pseudo-time model: no real-time scaling.
    Therefore Tf_real and t0_real are not needed here.
  */

  real<lower=0> mu;
  real<lower=0> phi;

  int<lower=1> Nnodes;
  array[Nnodes] int parents;
  array[Nnodes] int number_of_events;
  array[Nnodes] int idx_first_event;

  int<lower=1> Nevents;
  array[Nevents] int karyotypes;
  array[Nevents] real<lower=0> lengths;
  array[Nevents] int idx_first_mutation;

  int<lower=1> Nequations;
  array[Nequations] real coeff;
  array[Nequations] int mutations;

  int<lower=0> Ncells;

  int<lower=2> K_w;

  /*
    Priors on log(w).
    Here w is defined per unit of pseudo-time.
  */
  real log_w_mean_mu;
  real<lower=0> log_w_mean_sd;

  /*
    GP priors.
  */
  real<lower=0> gp_alpha_scale;
  real<lower=0> gp_rho_scale;

  /*
    Priors on latent relative positions.
  */
  array[Nevents] real<lower=0> u_alpha;
  array[Nevents] real<lower=0> u_beta;
}

parameters {
  vector<lower=0, upper=1>[Nevents] u;

  real log_w_mean;

  real<lower=0> gp_alpha;
  real<lower=0> gp_rho;

  vector[K_w] z_w;
}

transformed parameters {
  array[Nevents] real<lower=0, upper=1> tevents;

  vector[K_w] t_grid;
  matrix[K_w, K_w] K_gp;
  matrix[K_w, K_w] L_K_gp;

  vector[K_w] log_w_grid;
  vector<lower=0>[K_w] w_grid;

  for (k in 1:K_w) {
    t_grid[k] = (k - 0.5) / K_w;
  }

  K_gp = gp_exp_quad_cov_grid(K_w, gp_alpha, gp_rho, 1e-8);
  L_K_gp = cholesky_decompose(K_gp);

  log_w_grid = rep_vector(log_w_mean, K_w) + L_K_gp * z_w;
  w_grid = exp(log_w_grid);

  /*
    Construction of event times in pseudo-time.
  */
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
        if (event == 1) {
          tevents[idx] =
            tevents[idx_first_event[p]]
            + u[idx] * (1 - tevents[idx_first_event[p]]);
        } else {
          tevents[idx] =
            tevents[idx_first_event[p]]
            + u[idx] * (
              tevents[idx0] - tevents[idx_first_event[p]]
            );
        }
      }
    }
  }
}

model {
  /*
    Priors on event positions.
  */
  for (i in 1:Nevents) {
    u[i] ~ beta(u_alpha[i], u_beta[i]);
  }

  /*
    Non-centered GP prior for log(w).
  */
  z_w ~ normal(0, 1);

  log_w_mean ~ normal(log_w_mean_mu, log_w_mean_sd);

  gp_alpha ~ normal(0, gp_alpha_scale);
  gp_rho ~ inv_gamma(5, 5 * gp_rho_scale);

  /*
    Mutation likelihood in pseudo-time.

    Important:
    No multiplication by (Tf_real - t0_real).
    The integral of w is already in pseudo-time units.
  */
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    real t_p = (p == 0) ? 0.0 : tevents[idx_first_event[p]];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        real exposure_m =
          integrated_w_piecewise(t_p, tevents[idx0], K_w, w_grid);

        real rate_m =
          mu
          * lengths[idx0]
          * coeff[idx_first_mutation[idx0]]
          * exposure_m;

        target += poisson_lpmf(
          mutations[idx_first_mutation[idx0]]
          | rate_m + 1e-9
        );

      } else {
        real exposure_a =
          integrated_w_piecewise(t_p, tevents[idx], K_w, w_grid);

        real exposure_b =
          integrated_w_piecewise(tevents[idx], tevents[idx0], K_w, w_grid);

        real rate_a =
          mu
          * lengths[idx]
          * coeff[idx_first_mutation[idx]]
          * exposure_a;

        real rate_b =
          mu
          * lengths[idx]
          * coeff[idx_first_mutation[idx] + 1]
          * exposure_b;

        target += poisson_lpmf(
          mutations[idx_first_mutation[idx]]
          | rate_a + 1e-9
        );

        target += poisson_lpmf(
          mutations[idx_first_mutation[idx] + 1]
          | rate_b + 1e-9
        );
      }
    }
  }

  /*
    Cell number likelihood in pseudo-time.

    This corresponds to:

      E[Ncells] = exp( integral_{tevents[1]}^1 w(s) ds )

    because neg_binomial_2_log uses the log of the mean.
  */
  {
    real log_expected_val =
      integrated_w_piecewise(tevents[1], 1.0, K_w, w_grid);

    target += neg_binomial_2_log_lpmf(Ncells | log_expected_val, phi);
  }
}

generated quantities {
  array[Nequations] int mutations_ppcheck;
  array[Nequations] real log_lik;
  real log_lik_cells;

  /*
    Posterior predictive mutations.
  */
  for (node in 1:Nnodes) {
    int idx0 = idx_first_event[node];
    int p = parents[node];

    real t_p_post = (p == 0) ? 0.0 : tevents[idx_first_event[p]];

    for (event in 1:number_of_events[node]) {
      int idx = idx0 + event - 1;

      if (event == 1) {
        real exposure_m_post =
          integrated_w_piecewise(
            t_p_post,
            tevents[idx0],
            K_w,
            w_grid
          );

        real r_m_post =
          mu
          * lengths[idx0]
          * coeff[idx_first_mutation[idx0]]
          * exposure_m_post;

        mutations_ppcheck[idx_first_mutation[idx0]] =
          poisson_rng(r_m_post + 1e-9);

        log_lik[idx_first_mutation[idx0]] =
          poisson_lpmf(
            mutations[idx_first_mutation[idx0]]
            | r_m_post + 1e-9
          );

      } else {
        real exposure_a_post =
          integrated_w_piecewise(
            t_p_post,
            tevents[idx],
            K_w,
            w_grid
          );

        real exposure_b_post =
          integrated_w_piecewise(
            tevents[idx],
            tevents[idx0],
            K_w,
            w_grid
          );

        real r_a_post =
          mu
          * lengths[idx]
          * coeff[idx_first_mutation[idx]]
          * exposure_a_post;

        real r_b_post =
          mu
          * lengths[idx]
          * coeff[idx_first_mutation[idx] + 1]
          * exposure_b_post;

        mutations_ppcheck[idx_first_mutation[idx]] =
          poisson_rng(r_a_post + 1e-9);

        mutations_ppcheck[idx_first_mutation[idx] + 1] =
          poisson_rng(r_b_post + 1e-9);

        log_lik[idx_first_mutation[idx]] =
          poisson_lpmf(
            mutations[idx_first_mutation[idx]]
            | r_a_post + 1e-9
          );

        log_lik[idx_first_mutation[idx] + 1] =
          poisson_lpmf(
            mutations[idx_first_mutation[idx] + 1]
            | r_b_post + 1e-9
          );
      }
    }
  }

  /*
    Log-likelihood for Ncells.
  */
  {
    real log_expected_val =
      integrated_w_piecewise(tevents[1], 1.0, K_w, w_grid);

    log_lik_cells =
      neg_binomial_2_log_lpmf(Ncells | log_expected_val, phi);
  }
}