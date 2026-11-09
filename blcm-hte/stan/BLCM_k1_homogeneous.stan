// BLCM_k1_homogeneous.stan
// Single treatment effect, no latent classes. Comparator for LOO.
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
  vector<lower=0>[2] sigma_std;      // [1] = control, [2] = treated
}
model {
  beta_std  ~ normal(0, 0.5);
  tau_std   ~ normal(0, 0.5);
  sigma_std ~ normal(0, 1);

  for (n in 1:N)
    y_std[n] ~ normal(x_std[n] * beta_std + tau_std * t[n],
                      sigma_std[t[n] + 1]);
}
generated quantities {
  vector[N] log_lik;
  real tau = tau_std * sd_y;
  vector[2] sigma = sigma_std * sd_y;      // [1] = control, [2] = treated
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y_std[n] |
                   x_std[n] * beta_std + tau_std * t[n],
                   sigma_std[t[n] + 1]);
  }
}
