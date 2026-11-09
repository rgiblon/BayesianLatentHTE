// BLCM_k1_homogeneous_shared.stan
// Single treatment effect, no latent classes, shared residual SD across arms.
// Comparator for LOO. Identical to BLCM_k1_homogeneous.stan except that
// sigma_std is a single scalar rather than one per arm.
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] x;
  vector[N] y;
  array[N] int<lower=0, upper=1> t;
}
transformed data {
  real mean_y = mean(y);
  real sd_y   = sd(y);
  vector[N] y_std = (y - mean_y) / sd_y;

  matrix[N, K] x_std = x;              // col 1 = intercept, left alone
  for (k in 2:K) {
    x_std[, k] = (x[, k] - mean(x[, k])) / sd(x[, k]);
  }
}
parameters {
  vector[K] beta_std;
  real tau_std;
  real<lower=0> sigma_std;             // shared across arms
}
model {
  beta_std  ~ normal(0, 0.5);
  tau_std   ~ normal(0, 0.5);
  sigma_std ~ normal(0, 1);

  y_std ~ normal(x_std * beta_std + tau_std * to_vector(t), sigma_std);
}
generated quantities {
  vector[N] log_lik;
  real tau   = tau_std * sd_y;
  real sigma = sigma_std * sd_y;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y_std[n] |
                   x_std[n] * beta_std + tau_std * t[n],
                   sigma_std);
  }
}
