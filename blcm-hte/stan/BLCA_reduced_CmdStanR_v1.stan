data {
  int<lower=1> N;
  int<lower=1> K;                    // includes intercept as col 1
  int<lower=2> C;
  matrix[N, K] x;                    // original-scale X (with intercept column=1)
  vector[N] y;                       // original-scale Y
  array[N] int<lower=0, upper=1> t;  // treatment indicator
}

transformed data {
  // Standardize X (columns) and Y; keep intercept and any constant columns safe
  vector[K] mu_x;
  vector[K] sd_x;
  matrix[N, K] x_std;

  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;

  // Robust column means/SDs for X
  for (k in 1:K) {
    mu_x[k] = mean(col(x, k));
    sd_x[k] = sd(col(x, k));

    if (k == 1 || sd_x[k] <= 1e-12) {
      sd_x[k] = 1;
      x_std[, k] = col(x, k);
    } else {
      x_std[, k] = (col(x, k) - mu_x[k]) / sd_x[k];
    }
  }
}

parameters {
  // Reduced model: covariate-INDEPENDENT class proportions
  // pi_k replaces the full multinomial logistic membership model.
  // Class prevalence is estimated as a free simplex — constant across individuals.
  // This is the key structural difference from the full model.
  simplex[C] pi_k;

  // Shared prognostic regression (std scale) — unchanged
  vector[K] beta_std;

  // Gap parameterisation for treatment effects — unchanged
  real          tau_mid;
  real<lower=0> delta_harm;
  real<lower=0> delta_benefit;

  // Shared residual SD (std scale) — unchanged
  real<lower=1e-6> sigma_std;
}

transformed parameters {
  // Reconstruct ordered tau_std from gap parameterisation — unchanged
  vector[C] tau_std;
  tau_std[1] = tau_mid - delta_harm;
  tau_std[2] = tau_mid;
  tau_std[3] = tau_mid + delta_benefit;

  // Log class probabilities: constant across individuals (no covariates)
  vector[C] log_pi_k = log(pi_k);

  // Back-transform shared prognostic effects to original scale — unchanged
  vector[K] beta;
  vector[C] tau;
  real sigma;

  for (k in 2:K) {
    beta[k] = sd_y * beta_std[k] / sd_x[k];
  }
  {
    real adj = sd_y * beta_std[1] + mu_y;
    for (k in 2:K) adj -= beta[k] * mu_x[k];
    beta[1] = adj;
  }

  for (c in 1:C) tau[c] = sd_y * tau_std[c];
  sigma = sd_y * sigma_std;
}

model {
  // Prior on class proportions: symmetric Dirichlet(1) = uniform over simplex
  // No information about which class is larger — let data decide.
  pi_k ~ dirichlet(rep_vector(1.0, C));

  // Remaining priors — unchanged from full model
  beta_std      ~ normal(0, 0.5);
  tau_mid       ~ normal(0, 0.2);
  delta_harm    ~ normal(0, 1.0);
  delta_benefit ~ normal(0, 1.0);
  sigma_std     ~ normal(0, 1);

  // Mixture likelihood (std scale) — unchanged except log_pi_k is now scalar
  {
    vector[C] lps;
    for (n in 1:N) {
      real mu_base = dot_product(row(x_std, n), beta_std);
      for (c in 1:C) {
        real mu_std = mu_base + tau_std[c] * t[n];
        lps[c] = log_pi_k[c] + normal_lpdf(y_std[n] | mu_std, sigma_std);
      }
      target += log_sum_exp(lps);
    }
  }
}

generated quantities {
  // Structure identical to full model for direct comparison.
  // NOTE: class_prior[n] = pi_k for all n — constant across individuals
  // since there is no covariate-dependent membership model.
  // Comparing class_prior between full and reduced models directly shows
  // the contribution of the covariate membership component.
  matrix[N, C] p_class;
  matrix[N, C] class_prior;
  vector[N] log_lik;
  vector[N] y_rep;

  {
    vector[C] lps;
    real mu_base;
    for (n in 1:N) {
      mu_base = dot_product(row(x_std, n), beta_std);

      for (c in 1:C)
        lps[c] = log_pi_k[c] + normal_lpdf(y_std[n] | mu_base + tau_std[c] * t[n], sigma_std);

      p_class[n]     = softmax(lps)';
      log_lik[n]     = log_sum_exp(lps);
      class_prior[n] = pi_k';   // constant — same for all individuals

      int cls = categorical_rng(p_class[n]');
      y_rep[n] = normal_rng(mu_base + tau_std[cls] * t[n], sigma_std) * sd_y + mu_y;
    }
  }
}
