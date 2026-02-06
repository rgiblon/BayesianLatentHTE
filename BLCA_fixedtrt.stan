data {
  int<lower=1> N;                  // individuals
  int<lower=1> K;                  // covariates (incl intercept as col 1)
  int<lower=2> C;                  // classes
  matrix[N, K] x;                  // original-scale X
  vector[N] y;                     // original-scale Y
  array[N] int<lower=0, upper=1> t;  // treatment indicator (fixed 1:1)
}

transformed data {
  // Standardize X (columns) and Y; keep first column as intercept=1
  vector[K] mu_x;
  vector[K] sd_x;
  matrix[N, K] x_std;
  real mu_y = mean(y);
  real sd_y = sd(y);
  vector[N] y_std = (y - mu_y) / sd_y;

  // column means/SDs for X (protect intercept)
  for (k in 1:K) {
    if (k == 1) {
      mu_x[k] = 0;
      sd_x[k] = 1;
      x_std[, k] = rep_vector(1, N);
    } else {
      mu_x[k] = mean(col(x, k));
      sd_x[k] = sd(col(x, k));
      x_std[, k] = (col(x, k) - mu_x[k]) / sd_x[k];
    }
  }
}

parameters {
  simplex[C] pi;                       // class proportions

  // Standardized-scale regression per class for outcome model
  // y_std ~ Normal( x_std * beta_y_std[c] + tau_std[c] * t, sigma_std[c] )
  array[C] vector[K] beta_y_std;       // class-specific standardized betas
  ordered[C] tau_std;                  // class-specific standardized treatment effects
  vector<lower=0>[C] sigma_std;        // class-specific residual SDs on standardized Y
}

transformed parameters {
  // Map outcome parameters back to the original scale for interpretability
  array[C] vector[K] beta_y;   // original-scale betas
  vector[C] tau;               // original-scale treatment effects
  vector[C] sigma;             // original-scale residual SDs

  // For each class c:
  for (c in 1:C) {
    // slopes j>=2
    for (k in 2:K) {
      beta_y[c][k] = sd_y * beta_y_std[c][k] / sd_x[k];
    }
    // treatment effect
    tau[c] = sd_y * tau_std[c];

    // intercept (adjust for standardization shifts)
    real adj = sd_y * beta_y_std[c][1] + mu_y;   // baseline on std scale back to original mean
    for (k in 2:K) adj -= beta_y[c][k] * mu_x[k];
    beta_y[c][1] = adj;

    // residual SD
    sigma[c] = sd_y * sigma_std[c];
  }
}

model {
  // Priors (standardized scale)
  for (c in 1:C) {
    beta_y_std[c] ~ normal(0, 1);
  }
  tau_std ~ normal(0, 1);
  sigma_std ~ student_t(4, 0, 1); 

  // Mixture likelihood on STANDARDIZED scale for better geometry
  {
    vector[C] lps;
    for (n in 1:N) {
      for (c in 1:C) {
        real mu_std = dot_product(row(x_std, n), beta_y_std[c]) + tau_std[c] * t[n];
        lps[c] = log(pi[c]) + normal_lpdf(y_std[n] | mu_std, sigma_std[c]);
      }
      target += log_sum_exp(lps);
    }
  }
}

generated quantities {
  // Posterior class probabilities, pointwise log-lik (for loo/waic), and posterior predictive
  matrix[N, C] p_class;
  vector[N] log_lik;
  vector[N] y_rep;

  {
    vector[C] lps;
    for (n in 1:N) {
      for (c in 1:C) {
        // compute on standardized scale
        real mu_std = dot_product(row(x_std, n), beta_y_std[c]) + tau_std[c] * t[n];
        lps[c] = log(pi[c]) + normal_lpdf(y_std[n] | mu_std, sigma_std[c]);
      }
      // posterior class probabilities
      p_class[n] = softmax(lps)';

      // pointwise log-likelihood (mixture on std scale)
      log_lik[n] = log_sum_exp(lps);

      // posterior predictive draw (mixture sampling via class probs)
      int cls = categorical_rng(p_class[n]');
      real mu_std_cls = dot_product(row(x_std, n), beta_y_std[cls]) + tau_std[cls] * t[n];
      real y_std_rep = normal_rng(mu_std_cls, sigma_std[cls]);
      y_rep[n] = y_std_rep * sd_y + mu_y;  // back to original scale
    }
  }
}
