// ============================================================
// Bayesian Latent Class Model — Binary Outcome
// ============================================================
// Adapted from BLCA_fixedtrt_CmdStanR_v2.stan (continuous outcome).
//
// Key differences from the continuous model:
//   - y is binary (0/1); likelihood is Bernoulli logit, not Gaussian
//   - No sigma parameter — Bernoulli has no residual variance
//   - No Y standardisation — tau and beta are on the log-odds scale
//   - X standardisation is retained for HMC geometry
//   - tau back-transformed to odds ratios (tau_or) for interpretability
//   - No log-Jacobian correction needed for log_lik (no scale transform)
//   - Prior on tau_mid loosened slightly: normal(0, 0.2) on std scale
//     corresponded to ~±0.4 log-odds at 2 sigma; on the log-odds scale
//     directly, normal(0, 0.5) allows OR ~0.37 to 2.7 — more appropriate
//     for a clinical null.
//
// The sigma/tau_3 confound present in the continuous model is eliminated
// by construction: Bernoulli likelihood has no free variance parameter
// to compete with tau in the upper tail, so benefit class bias is
// expected to attenuate relative to the continuous specification.
// ============================================================

data {
  int<lower=1> N;
  int<lower=1> K;                             // includes intercept as col 1
  int<lower=2> C;
  matrix[N, K] x;                             // original-scale X (with intercept column=1)
  array[N] int<lower=0, upper=1> y;           // binary outcome
  array[N] int<lower=0, upper=1> t;           // treatment indicator
}

transformed data {
  // Standardize X columns only (no Y standardisation for binary outcome)
  // Intercept column (k=1) and constant columns are left unchanged
  vector[K] mu_x;
  vector[K] sd_x;
  matrix[N, K] x_std;

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
  // Mixing model — identical NCP structure to continuous model
  matrix[K-1, C-1] gamma_slopes_raw;
  vector[C-1]       gamma_intercept_raw;
  real<lower=1e-6>  gamma_sd;

  // Shared prognostic regression (log-odds scale)
  // beta[1] is the intercept on the log-odds scale
  // beta[2:K] are covariate effects on the log-odds scale
  // NOTE: no standardisation of y, so beta is directly on log-odds scale
  // and does not require back-transformation. Covariates are standardised
  // so beta[k] for k>1 represents the log-OR per 1 SD change in covariate k.
  vector[K] beta_std;

  // Gap parameterisation for treatment effects (log-odds scale):
  //   tau[1] = tau_mid - delta_harm    (harm class,    log-OR < 0)
  //   tau[2] = tau_mid                 (null class,    log-OR ~ 0)
  //   tau[3] = tau_mid + delta_benefit (benefit class, log-OR > 0)
  //
  // tau_mid is on the log-odds scale directly — no standardisation needed.
  // delta_harm and delta_benefit are gaps in log-odds units.
  // Ordering enforced by lower=0 on deltas, as in the continuous model.
  real          tau_mid;
  real<lower=0> delta_harm;
  real<lower=0> delta_benefit;
}

transformed parameters {
  // Reconstruct ordered tau vector from gap parameterisation (log-odds scale)
  vector[C] tau_std;
  tau_std[1] = tau_mid - delta_harm;
  tau_std[2] = tau_mid;
  tau_std[3] = tau_mid + delta_benefit;

  // Odds ratios for interpretability
  // tau_or[c] = exp(tau_std[c]) — the class-specific treatment odds ratio
  vector[C] tau_or;
  for (c in 1:C) tau_or[c] = exp(tau_std[c]);

  // Non-centred reconstruction of gamma_free (identical to continuous model)
  matrix[K, C-1] gamma_free;
  gamma_free[1, ]   = (gamma_sd * gamma_intercept_raw)';
  gamma_free[2:K, ] = gamma_sd * gamma_slopes_raw;

  matrix[K, C] gamma;
  matrix[N, C] log_pi;

  // Back-transform beta to original (unstandardised X) scale
  // beta[1] = intercept adjusted for X standardisation
  // beta[k] for k>1 = log-OR per original unit of covariate k
  // NOTE: unlike the continuous model, there is no sd_y factor here —
  // beta is already on the log-odds scale.
  vector[K] beta;
  for (k in 2:K) {
    beta[k] = beta_std[k] / sd_x[k];
  }
  {
    real adj = beta_std[1];
    for (k in 2:K) adj -= beta[k] * mu_x[k];
    beta[1] = adj;
  }

  // Build gamma and log_pi (identical to continuous model)
  gamma[, 1] = rep_vector(0, K);
  for (c in 2:C) {
    gamma[, c] = gamma_free[, c-1];
  }

  for (n in 1:N) {
    vector[C] eta;
    for (c in 1:C) {
      eta[c] = dot_product(row(x_std, n), gamma[, c]);
    }
    log_pi[n] = log_softmax(eta)';
  }
}

model {
  // Priors — membership model (identical to continuous model)
  gamma_sd ~ normal(0, 1);
  to_vector(gamma_slopes_raw) ~ normal(0, 1);
  gamma_intercept_raw ~ normal(0, 1);

  // Prognostic coefficients (log-odds scale)
  // normal(0, 0.5) on standardised X scale: covariate effects up to ~1 log-OR
  // at 2 sigma, corresponding to OR ~0.14 to 7.4 — weakly informative.
  beta_std ~ normal(0, 0.5);

  // Gap parameterisation priors (log-odds scale):
  //
  // tau_mid ~ normal(0, 0.5): null class log-OR regularised toward 0.
  // Allows ~±1.0 log-OR at 2 sigma (OR ~0.37 to 2.7) — appropriate for
  // a clinical null. Wider than the continuous model's normal(0, 0.2)
  // because there is no sd_y scaling — we are working directly in log-odds
  // units rather than standardised outcome units.
  //
  // delta_harm / delta_benefit ~ normal(0, 1) with lower=0 (half-normal):
  // Half-normal on log-odds scale; most mass on gaps of 0-2 log-OR units
  // (OR ratios up to ~7.4), allowing clinically meaningful separation while
  // regularising against implausibly large effects. Tighten if prior
  // predictive checks suggest unrealistic extreme probabilities.
  tau_mid       ~ normal(0, 0.5);
  delta_harm    ~ normal(0, 1.0);
  delta_benefit ~ normal(0, 1.0);

  // Bernoulli logit mixture likelihood
  // mu_base = X_i * beta_std (log-odds of outcome, prognostic component)
  // mu_base + tau_std[c] * t[n] = log-odds for class c, treatment t[n]
  // No sigma parameter — Bernoulli variance is determined by the mean.
  {
    vector[C] lps;
    for (n in 1:N) {
      real mu_base = dot_product(row(x_std, n), beta_std);
      for (c in 1:C) {
        lps[c] = log_pi[n, c] +
                 bernoulli_logit_lpmf(y[n] | mu_base + tau_std[c] * t[n]);
      }
      target += log_sum_exp(lps);
    }
  }
}

generated quantities {
  // NOTE: p_class[n] = p(class | y, x, t) — posterior class probability
  // class_prior[n] = p(class | x) — membership model only, no outcome
  // log_lik[n] = log p(y_n | x_n, t_n) — pointwise log-likelihood for LOO/WAIC
  // y_rep[n] = posterior predictive draw (binary, 0 or 1)
  //
  // NOTE on classification evaluation: treated-arm accuracy is the primary
  // metric. Control arm patients have t[n]=0, so tau does not enter their
  // likelihood and class assignment is driven entirely by gamma.
  //
  // NOTE on log_lik scale: for binary outcomes log_lik is already on the
  // natural scale (log-probability of a Bernoulli observation). No Jacobian
  // correction is needed — LOO/WAIC comparisons are directly interpretable.

  matrix[N, C] p_class;
  matrix[N, C] class_prior;
  vector[N] log_lik;
  array[N] int y_rep;          // binary posterior predictive draws

  {
    vector[C] lps;
    real mu_base;
    for (n in 1:N) {
      mu_base = dot_product(row(x_std, n), beta_std);

      for (c in 1:C)
        lps[c] = log_pi[n, c] +
                 bernoulli_logit_lpmf(y[n] | mu_base + tau_std[c] * t[n]);

      p_class[n]     = softmax(lps)';
      log_lik[n]     = log_sum_exp(lps);
      class_prior[n] = exp(log_pi[n]);

      // Posterior predictive: sample class then draw binary outcome
      int cls = categorical_rng(p_class[n]');
      y_rep[n] = bernoulli_logit_rng(mu_base + tau_std[cls] * t[n]);
    }
  }
}
