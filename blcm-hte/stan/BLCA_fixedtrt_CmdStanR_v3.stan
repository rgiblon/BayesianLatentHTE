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
  // Protect intercept and any constant columns from division by zero
  // NOTE: mu_x[1] is computed here but is not used in the membership model or slope
  // back-transform — the intercept column is handled separately in the beta back-transform
  // adjustment block below. This is intentional and not a bug.
  for (k in 1:K) {
    mu_x[k] = mean(col(x, k));
    sd_x[k] = sd(col(x, k));

    if (k == 1 || sd_x[k] <= 1e-12) {
      // Intercept or constant column: do NOT standardize
      sd_x[k] = 1;
      x_std[, k] = col(x, k);   // keeps intercept as all 1s
    } else {
      x_std[, k] = (col(x, k) - mu_x[k]) / sd_x[k];
    }
  }
}

parameters {
  // Mixing model: multinomial logit with baseline class 1 (gamma[,1] = 0)
  // Non-centred parameterisation: gamma_free = gamma_sd * gamma_free_raw
  // gamma_free_raw ~ normal(0,1) always — scale (gamma_sd) and shape are decoupled,
  // eliminating the funnel geometry that arises in the centred form when gamma_sd -> 0.
  //
  // Intercept and slope raw parameters are split to allow free per-class intercepts:
  // gamma_intercept_raw[c] captures unconditional class prevalence for non-baseline
  // classes, independent of covariates. This is important when covariates have weak
  // class-discriminating signal (R^2 ~ 0.15-0.20), as the model can otherwise only
  // express unconditional prevalence imbalance by distorting covariate coefficients.
  // gamma_slopes_raw captures covariate-driven membership (rows 2:K, non-baseline classes).
  
  matrix[K-1, C-1] gamma_slopes_raw;      // slopes only (rows 2:K, non-baseline classes)
  vector[C-1]       gamma_intercept_raw;  // free intercepts for non-baseline classes
  real<lower=1e-6>  gamma_sd;            // learned global shrinkage scale

  // Shared prognostic regression (std scale)
  vector[K] beta_std;                // includes intercept at beta_std[1]

  // Gap parameterization for treatment effects (std scale):
  //   tau_std[1] = tau_mid - delta_harm    (harm class)
  //   tau_std[2] = tau_mid                 (null class, strongly shrunk to ~0)
  //   tau_std[3] = tau_mid + delta_benefit (benefit class)
  //
  // Replaces ordered[C] tau_std. The ordered constraint conflated two separate
  // problems: (a) where the null class sits, and (b) how far harm/benefit deviate
  // from null. When tau[1] and tau[3] were uncertain, their uncertainty dragged
  // tau[2] wide — producing the wide null interval (width ~4.8 vs ~3.3 for tails).
  // Decoupling via tau_mid + deltas lets us strongly regularise (a) independently,
  // forcing null-class assignment onto the membership model (gamma) rather than
  // the outcome likelihood. Ordering is enforced by construction via lower=0 on deltas.
  real          tau_mid;          // null class effect — strongly shrunk toward 0
  real<lower=0> delta_harm;       // harm separation:    tau[1] = tau_mid - delta_harm
  real<lower=0> delta_benefit;    // benefit separation: tau[3] = tau_mid + delta_benefit

  // Shared residual SD (std scale)
  real<lower=1e-6> sigma_std;
}

transformed parameters {
  // Reconstruct ordered tau_std vector from gap parameterisation
  vector[C] tau_std;
  tau_std[1] = tau_mid - delta_harm;
  tau_std[2] = tau_mid;
  tau_std[3] = tau_mid + delta_benefit;

  // Non-centred reconstruction of gamma_free:
  // Row 1 (intercept): free per non-baseline class — absorbs unconditional prevalence.
  // Rows 2:K (slopes): covariate-driven class separation.
  // Posterior geometry is separable from gamma_sd for both, eliminating funnel issues.
  matrix[K, C-1] gamma_free;
  gamma_free[1, ]   = (gamma_sd * gamma_intercept_raw)';  // free intercepts, NCP
  gamma_free[2:K, ] = gamma_sd * gamma_slopes_raw;        // slopes, NCP

  // Full gamma including baseline class 1 = 0
  // Only the slope rows are pinned to zero for the baseline class; the intercept
  // row is implicitly zero for class 1 via the softmax baseline convention.
  matrix[K, C] gamma;
  matrix[N, C] log_pi;               // per-person log class probabilities

  // Original-scale quantities for interpretability
  vector[K] beta;   // original-scale betas
  vector[C] tau;    // original-scale taus
  real sigma;       // original-scale sigma

  // Build gamma with baseline class 1 set to zero for identifiability
  gamma[, 1] = rep_vector(0, K);
  for (c in 2:C) {
    gamma[, c] = gamma_free[, c-1];
  }

  // Compute log pi_i = log softmax(x_i * gamma)
  // Note: full rows of gamma (including intercept row 1) enter eta here,
  // so both unconditional prevalence (via gamma[1,c]) and covariate-driven
  // separation (via gamma[2:K,c]) contribute to class probabilities.
  for (n in 1:N) {
    vector[C] eta;
    for (c in 1:C) {
      eta[c] = dot_product(row(x_std, n), gamma[, c]);
    }
    log_pi[n] = log_softmax(eta)';
  }

  // Back-transform shared prognostic effects to original scale
  for (k in 2:K) {
    beta[k] = sd_y * beta_std[k] / sd_x[k];
  }
  {
    real adj = sd_y * beta_std[1] + mu_y;
    for (k in 2:K) adj -= beta[k] * mu_x[k];
    beta[1] = adj;
  }

  // Treatment effects to original scale
  for (c in 1:C) tau[c] = sd_y * tau_std[c];

  // Residual SD to original scale
  sigma = sd_y * sigma_std;
}

model {
  // Priors
  // Hierarchical prior on membership model coefficients:
  // gamma_sd is learned from data rather than fixed, allowing the model to decide
  // how much covariate separation between classes is supported by the evidence.
  // More flexible than fixed normal(0, 0.3) — if covariates don't separate classes
  // well, gamma_sd shrinks toward 0; if they do, it can grow.
  // Non-centred parameterisation: gamma_free_raw ~ normal(0,1), scale applied in
  // transformed parameters via gamma_free = gamma_sd * gamma_free_raw.
  // The prior on gamma_free_raw is always standard normal regardless of gamma_sd,
  // so the sampler geometry is uniform — no funnel.

  gamma_sd ~ normal(0, 1);
  to_vector(gamma_slopes_raw) ~ normal(0, 1);
  gamma_intercept_raw ~ normal(0, 1);   // same NCP prior as slopes

  // NOTE: If classes are not separating well, consider loosening the
  // gamma_sd hyperprior further (e.g. normal(0, 2.0)). The current normal(0, 1)
  // allows moderate-to-strong covariate separation while providing some regularisation.

  // NOTE: beta_std use normal(0, 0.5) — this assumes most
  // effects are < 0.5 SD on the standardized scale. Verify against domain knowledge
  // via prior predictive checks; loosen to normal(0, 1) if effects may be larger.

  // NOTE: beta_std is SHARED across all classes — only tau varies by
  // class. This means all class-specific prognosis is identical; covariate-outcome
  // heterogeneity is attributed entirely to the membership model (gamma), not the
  // outcome model. This is a strong assumption.

  beta_std ~ normal(0, 0.5);

  // Gap parameterization priors (std scale):
  //
  // Original tau_mid ~ normal(0, 0.2): strongly regularises the null class to ~0 effect.
  // Allows ~+/-0.4 SD at 2 sigma — generous for a true null. The tight prior
  // forces null-class membership to be expressed via covariates (gamma) rather
  // than via outcome wiggle room, which was the root cause of poor null-class
  // classification. Tighten to normal(0, 0.1) if the true null is very close to 0.
  // For real data example (HIV trial) loosened tau_mid ~ norm(0,0.5) to match data.
  //
  // delta_harm / delta_benefit ~ normal(0, 1) with lower=0 (half-normal):
  // Places most mass on separations of 0–2 SD while allowing larger gaps if
  // the data support them. The lower=0 constraint enforces ordering by construction.

  tau_mid       ~ normal(0, 0.5);
  delta_harm    ~ normal(0, 1.0);
  delta_benefit ~ normal(0, 1.0);

  sigma_std ~ normal(0, 1);   

  // Mixture likelihood (std scale)
  {
    vector[C] lps;
    for (n in 1:N) {
      real mu_base = dot_product(row(x_std, n), beta_std);
      for (c in 1:C) {
        real mu_std = mu_base + tau_std[c] * t[n];
        lps[c] = log_pi[n, c] + normal_lpdf(y_std[n] | mu_std, sigma_std);
      }
      target += log_sum_exp(lps);
    }
  }
}

generated quantities {
  // NOTE: p_class[n] is the POSTERIOR class probability p(class | y, x, t)
  // — it reflects both the membership model and how well each class's tau explains
  // the observed outcome y. This is correct for patient classification.
  // class_prior[n] = exp(log_pi[n]) is the PRIOR class probability p(class | x)
  // from the membership model alone (no y). Comparing the two is a useful diagnostic:
  // if class_prior is near-uniform but p_class is confident, outcome y is driving
  // class assignment more than covariates x.
  //
  // NOTE on classification evaluation: posterior class accuracy should be assessed
  // on the TREATED arm only. Control arm patients have t[n]=0, so tau does not enter
  // their likelihood — their class assignment is driven entirely by the membership
  // model (gamma). Since latent classes are defined by treatment response, control-arm
  // classification is extrapolation by design, not a model performance test.
  //
  // NOTE on redundancy: lps must be recomputed here because Stan does not expose
  // intermediate likelihood values from the model block to generated quantities.
  // log_pi is already available from transformed parameters (no redundancy there).
  // The only unavoidable duplication is the normal_lpdf evaluation per observation.
  matrix[N, C] p_class;
  matrix[N, C] class_prior;           // p(class | x) — membership model only, no y
  vector[N] log_lik;
  vector[N] y_rep;

  // NOTE: log_lik is computed on the standardized outcome scale (y_std), not the
  // original scale. This is internally consistent for model comparison via LOO/WAIC —
  // ELPD differences between models are valid and unaffected by the scale choice.
  // However, absolute ELPD values are not interpretable in original-scale units.
  // If original-scale ELPD is needed, apply the log-Jacobian correction:
  //   elpd_original = elpd_std + N * log(1 / sd_y)  (i.e. subtract N * log(sd_y))
  // y_rep is already back-transformed to original scale (see below) and is correct
  // for posterior predictive checks against raw observed data.

  {
    vector[C] lps;
    real mu_base;
    for (n in 1:N) {
      // mu_base is shared across classes for observation n
      mu_base = dot_product(row(x_std, n), beta_std);

      // lps[c] = log p(class=c) + log p(y | class=c) — needed for posterior class probs
      for (c in 1:C)
        lps[c] = log_pi[n, c] + normal_lpdf(y_std[n] | mu_base + tau_std[c] * t[n], sigma_std);

      p_class[n]    = softmax(lps)';
      log_lik[n]    = log_sum_exp(lps);
      class_prior[n] = exp(log_pi[n]);  // p(class | x) — already computed in transformed parameters

      // Posterior predictive draw: sample class then draw outcome
      // y_rep is on the original scale — multiply by sd_y and add mu_y to back-transform.
      int cls = categorical_rng(p_class[n]');
      y_rep[n] = normal_rng(mu_base + tau_std[cls] * t[n], sigma_std) * sd_y + mu_y;
    }
  }
}
