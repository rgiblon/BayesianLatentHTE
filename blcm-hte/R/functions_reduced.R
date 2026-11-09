#Install rstan and dependencies (if you haven't already)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores


#Set working directory to dissertation


# Load the cmdrstanr library & neededx packages
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(matrixStats)
library(MASS)
library(dplyr)

mod <- cmdstan_model("BLCA_fixedtrt_CmdStanR_v2.stan")

# ---------- DGP function ----------
simulate_dataset <- function(seed,
                             N = 800, M = 10, C = 3,
                             pk = c(0.3, 0.3, 0.4),
                             mu = list(c(-0.6, -0.6), c(0, 0.6), c(0.6, 0)),
                             Sigma_within = diag(2) * 0.25,
                             p_bio = c(0.05, 0.5, 0.95),
                             beta_common = c(1.4, 0.9, -0.7, 0.4, 1.0),
                             beta_intercept = 0,
                             tau_true = c(-5, 0, 5),
                             sigma = c(6, 6, 6),
                             p_treat = 0.5) {
  set.seed(seed)
  
  # -------------------------------------------------------------------
  # 1. Generate X unconditionally via mixture marginal
  #  X generated first, then class assigned via exact Bayes p(class|X)
  # -------------------------------------------------------------------
  
  # Draw x1 from mixture marginal: sum_k pk * N(mu_k, Sigma_within)
  comp  <- sample(1:C, N, replace = TRUE, prob = pk)  
  x1_level1 <- do.call(rbind, lapply(comp, function(k)
    MASS::mvrnorm(1, mu[[k]], Sigma = Sigma_within)
  ))
  
  # Site-level covariates
  x2_level2 <- MASS::mvrnorm(M, mu = c(0, 0),
                             Sigma = matrix(c(1, -1, -1, 1), nrow = 2))
  site <- sample(1:M, N, replace = TRUE)
  
  # Draw B from marginal: p(B=1) = sum_k pk * p_bio[k]
  p_bio_marginal <- sum(pk * p_bio)
  B <- rbinom(N, 1, prob = p_bio_marginal)
  
  # Full X matrix
  X <- cbind(1, x1_level1, x2_level2[site, ], B)  # N x K=6
  K <- ncol(X)
  
  # -------------------------------------------------------------------
  # 2. Compute p(class | x1, B) via Bayes' theorem
  #    p(class=k | x) proportional to pk * p(x1 | class=k) * p(B | class=k)
  #    Only x1 and B are class-informative; site covariates are not
  # -------------------------------------------------------------------
  log_p_class_given_x <- matrix(NA_real_, nrow = N, ncol = C)
  
  for (k in 1:C) {
    # log p(class=k)
    log_prior_k <- log(pk[k])
    
    # log p(x1 | class=k) — bivariate normal
    log_lik_x1 <- mvtnorm::dmvnorm(x1_level1,
                                   mean  = mu[[k]],
                                   sigma = Sigma_within,
                                   log   = TRUE)
    
    # log p(B | class=k) — Bernoulli
    log_lik_B <- dbinom(B, size = 1, prob = p_bio[k], log = TRUE)
    
    log_p_class_given_x[, k] <- log_prior_k + log_lik_x1 + log_lik_B
  }
  
  # Normalise to get proper probabilities (numerically stable)
  pi_mat <- t(apply(log_p_class_given_x, 1, function(lp) {
    lp <- lp - max(lp)
    exp(lp) / sum(exp(lp))
  }))  # N x C — this is the TRUE p(class | x)
  
  # Draw class from exact posterior
  class <- apply(pi_mat, 1, function(p) sample(seq_len(C), 1, prob = p))
  
  # -------------------------------------------------------------------
  # 3. Treatment (random binomial)
  # -------------------------------------------------------------------
  T <- rbinom(N, 1, p_treat)
  
  # -------------------------------------------------------------------
  # 4. Outcome Y (continous) 
  # -------------------------------------------------------------------
  beta_full <- c(beta_intercept, beta_common)  # length K
  
  Y <- numeric(N)
  for (n in 1:N) {
    mu_n <- sum(X[n, ] * beta_full) + tau_true[class[n]] * T[n]
    Y[n] <- rnorm(1, mean = mu_n, sd = sigma[class[n]])
  }
  
  # -------------------------------------------------------------------
  # 5. Stan data 
  # -------------------------------------------------------------------
  stan_data <- list(
    N = N,
    K = K,
    C = C,
    x = X,
    y = as.numeric(Y),
    t = as.integer(T)
  )
  
  # -------------------------------------------------------------------
  # 6. Truth — includes pi_mat to evaluate against true p(class|x)
  # -------------------------------------------------------------------
  truth <- list(
    class  = class,
    pi_mat = pi_mat,
    tau    = tau_true,
    sigma  = sigma,
    sd_y   = sd(Y),      
    N      = N,
    C      = C,
    t      = T
  )
  
  list(stan_data = stan_data, truth = truth)
}

# ---------- output function ----------

# Existing function — unchanged
get_pbar_from_fit <- function(fit, N, C) {
  pc_draws <- posterior::as_draws_matrix(fit$draws("p_class"))
  cn <- colnames(pc_draws)
  
  keep <- grepl("^p_class\\[", cn)
  pc_draws <- pc_draws[, keep, drop = FALSE]
  cn <- colnames(pc_draws)
  
  idx <- do.call(
    rbind,
    regmatches(cn, regexec("^p_class\\[([0-9]+),([0-9]+)\\]$", cn))
  )
  
  i_idx <- as.integer(idx[, 2])
  c_idx <- as.integer(idx[, 3])
  
  p_means <- colMeans(pc_draws)
  
  pbar <- matrix(NA_real_, nrow = N, ncol = C)
  for (j in seq_along(p_means)) {
    pbar[i_idx[j], c_idx[j]] <- p_means[j]
  }
  
  pbar
}

# NEW: identical structure but pulls class_prior instead of p_class
# class_prior = p(class | x) from membership model only, no outcome
get_prior_from_fit <- function(fit, N, C) {
  cp_draws <- posterior::as_draws_matrix(fit$draws("class_prior"))
  cn <- colnames(cp_draws)
  
  keep <- grepl("^class_prior\\[", cn)
  cp_draws <- cp_draws[, keep, drop = FALSE]
  cn <- colnames(cp_draws)
  
  idx <- do.call(
    rbind,
    regmatches(cn, regexec("^class_prior\\[([0-9]+),([0-9]+)\\]$", cn))
  )
  
  i_idx <- as.integer(idx[, 2])
  c_idx <- as.integer(idx[, 3])
  
  p_means <- colMeans(cp_draws)
  
  prior_mat <- matrix(NA_real_, nrow = N, ncol = C)
  for (j in seq_along(p_means)) {
    prior_mat[i_idx[j], c_idx[j]] <- p_means[j]
  }
  
  prior_mat
}

extract_metrics <- function(fit, truth, compute_class_acc = FALSE) {
  
  s_tau <- fit$summary("tau")
  C <- nrow(s_tau)
  
  # diagnostics
  s_key <- fit$summary(variables = c("tau", "sigma_std", "tau_std"))
  max_rhat     <- max(s_key$rhat, na.rm = TRUE)
  min_ess_bulk <- min(s_key$ess_bulk, na.rm = TRUE)
  n_div        <- sum(fit$sampler_diagnostics()[,, "divergent__"], na.rm = TRUE)
  runtime_sec  <- sum(fit$time()$total)
  
  out <- data.frame(
    max_rhat     = max_rhat,
    min_ess_bulk = min_ess_bulk,
    n_divergent  = n_div,
    runtime_sec  = runtime_sec,
    stringsAsFactors = FALSE
  )
  
  # tau performance
  for (c in 1:C) {
    est <- s_tau$median[c]
    lo  <- s_tau$q5[c]
    hi  <- s_tau$q95[c]
    
    out[[paste0("tau_", c, "_med")]]     <- est
    out[[paste0("tau_", c, "_q5")]]      <- lo
    out[[paste0("tau_", c, "_q95")]]     <- hi
    out[[paste0("tau_", c, "_bias")]]    <- est - truth$tau[c]
    out[[paste0("tau_", c, "_cover90")]] <- as.integer(lo <= truth$tau[c] && truth$tau[c] <= hi)
    out[[paste0("tau_", c, "_width")]]   <- hi - lo
  }
  
  # sigma (original scale) — for diagnosing tau_3/sigma confound
  # sigma is back-transformed in the transformed parameters block
  # truth$sigma should be stored in simulate_dataset() as the true residual SD
  # sigma_bias will be NA if truth$sigma is not present — safe to add later
  
  s_sigma        <- fit$summary("sigma")
  out$sigma_med  <- s_sigma$median
  out$sigma_q5   <- s_sigma$q5
  out$sigma_q95  <- s_sigma$q95
  out$sigma_bias <- if (!is.null(truth$sigma)) s_sigma$median - truth$sigma[1] else NA_real_
  
  # class recovery 
  if (compute_class_acc) {
    N  <- length(truth$class)
    C2 <- length(truth$tau)
    
    if (is.null(truth$t)) {
      stop("truth$t is missing. Store treatment assignment in simulate_dataset().")
    }
    
    # Posterior class probs p(class | x, y, t)
    pbar      <- get_pbar_from_fit(fit, N = N, C = C2)
    class_hat <- max.col(pbar, ties.method = "first")
    
    # Prior class probs p(class | x) — membership model only, no outcome
    prior_mat       <- get_prior_from_fit(fit, N = N, C = C2)
    class_hat_prior <- max.col(prior_mat, ties.method = "first")
    
    idx_trt <- truth$t == 1
    idx_ctl <- truth$t == 0
    
    # --- Posterior-based accuracy (p_class) ---
    out$class_acc     <- mean(class_hat == truth$class)
    out$class_acc_trt <- if (sum(idx_trt) > 0) mean(class_hat[idx_trt] == truth$class[idx_trt]) else NA_real_
    out$class_acc_ctl <- if (sum(idx_ctl) > 0) mean(class_hat[idx_ctl] == truth$class[idx_ctl]) else NA_real_
    
    for (cc in 1:C2) {
      idx_class     <- truth$class == cc
      idx_class_trt <- truth$class == cc & truth$t == 1
      idx_class_ctl <- truth$class == cc & truth$t == 0
      
      out[[paste0("class_", cc, "_acc")]] <-
        if (sum(idx_class) > 0) mean(class_hat[idx_class] == cc) else NA_real_
      out[[paste0("class_", cc, "_acc_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(class_hat[idx_class_trt] == cc) else NA_real_
      out[[paste0("class_", cc, "_acc_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(class_hat[idx_class_ctl] == cc) else NA_real_
    }
    
    # --- Prior-based accuracy (class_prior) — membership model only ---
    # This is the accuracy you get using covariates alone, no outcome.
    # Key diagnostic: if prior_acc << class_acc, outcome is doing the work.
    # If prior_acc ~ class_acc for treated, gamma is contributing meaningfully.
    out$prior_acc     <- mean(class_hat_prior == truth$class)
    out$prior_acc_trt <- if (sum(idx_trt) > 0) mean(class_hat_prior[idx_trt] == truth$class[idx_trt]) else NA_real_
    out$prior_acc_ctl <- if (sum(idx_ctl) > 0) mean(class_hat_prior[idx_ctl] == truth$class[idx_ctl]) else NA_real_
    
    for (cc in 1:C2) {
      idx_class     <- truth$class == cc
      idx_class_trt <- truth$class == cc & truth$t == 1
      idx_class_ctl <- truth$class == cc & truth$t == 0
      
      out[[paste0("prior_class_", cc, "_acc")]] <-
        if (sum(idx_class) > 0) mean(class_hat_prior[idx_class] == cc) else NA_real_
      out[[paste0("prior_class_", cc, "_acc_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(class_hat_prior[idx_class_trt] == cc) else NA_real_
      out[[paste0("prior_class_", cc, "_acc_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(class_hat_prior[idx_class_ctl] == cc) else NA_real_
    }
    
    # --- Mean class probabilities by true class and treatment ---
    # Tells you whether the model is uncertain vs confidently wrong.
    # For each true class, what is the mean posterior prob assigned to correct class?
    # Values near 1/C = model is uncertain; values > 0.5 = meaningful signal.
    for (cc in 1:C2) {
      idx_class     <- truth$class == cc
      idx_class_trt <- truth$class == cc & truth$t == 1
      idx_class_ctl <- truth$class == cc & truth$t == 0
      
      # Mean posterior prob on correct class
      out[[paste0("mean_pclass_", cc, "_correct")]] <-
        if (sum(idx_class) > 0) mean(pbar[idx_class, cc]) else NA_real_
      out[[paste0("mean_pclass_", cc, "_correct_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(pbar[idx_class_trt, cc]) else NA_real_
      out[[paste0("mean_pclass_", cc, "_correct_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(pbar[idx_class_ctl, cc]) else NA_real_
      
      # Same for prior
      out[[paste0("mean_prior_", cc, "_correct")]] <-
        if (sum(idx_class) > 0) mean(prior_mat[idx_class, cc]) else NA_real_
      out[[paste0("mean_prior_", cc, "_correct_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(prior_mat[idx_class_trt, cc]) else NA_real_
      out[[paste0("mean_prior_", cc, "_correct_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(prior_mat[idx_class_ctl, cc]) else NA_real_
    }

    # --- Classification confidence: entropy and max posterior probability ---
    # Entropy = 0 for certain assignment, log(C) for uniform (maximally uncertain)
    # Max prob = 1 for certain, 1/C for uniform
    # Stratified by true class and treatment arm to assess where uncertainty lives
    entropy    <- apply(pbar, 1, function(p) -sum(p * log(p + 1e-10)))
    max_pclass <- apply(pbar, 1, max)
    
    for (cc in 1:C2) {
      for (trt in c(0, 1)) {
        arm_label  <- ifelse(trt == 1, "trt", "ctl")
        idx_cell   <- truth$class == cc & truth$t == trt
        
        out[[paste0("entropy_class", cc, "_", arm_label)]] <-
          if (sum(idx_cell) > 0) mean(entropy[idx_cell],    na.rm = TRUE) else NA_real_
        out[[paste0("maxprob_class", cc, "_", arm_label)]] <-
          if (sum(idx_cell) > 0) mean(max_pclass[idx_cell], na.rm = TRUE) else NA_real_
      }
    }    
   
    # --- Off-diagonal misclassification rates (treated arm, primary metric) ---
    # Named as: misclass_<true>_to_<assigned>_trt
    # Catastrophic errors: harm<->benefit (1<->3)
    # Adjacent errors:     harm<->null (1<->2), benefit<->null (3<->2)
    misclass_pairs <- list(
      c(1, 2), # harm    -> null
      c(1, 3), # harm    -> benefit  (catastrophic)
      c(2, 1), # null    -> harm
      c(2, 3), # null    -> benefit
      c(3, 1), # benefit -> harm     (catastrophic)
      c(3, 2)  # benefit -> null
    )
    
    for (pair in misclass_pairs) {
      true_c     <- pair[1]
      assigned_c <- pair[2]
      col_name   <- paste0("misclass_", true_c, "_to_", assigned_c, "_trt")
      idx_true_trt <- truth$class == true_c & truth$t == 1
      out[[col_name]] <-
        if (sum(idx_true_trt) > 0)
          mean(class_hat[idx_true_trt] == assigned_c)
      else NA_real_
    }
    
  }  # end if (compute_class_acc)
  
  out
}
# ---------- init function ----------
make_init_fun <- function(stan_data, truth) {
  K <- stan_data$K
  C <- stan_data$C
  
  function() {
    tau_std_true <- truth$tau / truth$sd_y
    list(
      beta_std            = rnorm(K, 0, 0.1),
     #gamma_slopes_raw    = matrix(rnorm((K-1) * (C-1), 0, 0.1), K-1, C-1),
     #gamma_intercept_raw = rnorm(C-1, 0, 0.1),
     #gamma_sd            = 0.3,
      tau_mid             = tau_std_true[2],
      delta_harm          = tau_std_true[2] - tau_std_true[1],  # must be > 0
      delta_benefit       = tau_std_true[3] - tau_std_true[2],  # must be > 0
      sigma_std           = 1
    )
  }
}

#make_init_fun <- function(stan_data, truth) {
  #K <- stan_data$K
  #C <- stan_data$C
  
  #function() {
    #list(
      #beta_std            = rnorm(K, 0, 0.1),
      #gamma_slopes_raw    = matrix(rnorm((K-1) * (C-1), 0, 0.1), K-1, C-1),
      #gamma_intercept_raw = rnorm(C-1, 0, 0.1),
      #gamma_sd            = 0.3,
      #tau_std             = truth$tau / truth$sd_y,
      #sigma_std           = 1
#    )
#  }
#}

# ---------- sampling args ----------
fit_args <- list(
  chains          = 3,
  parallel_chains = 3,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.90,
  max_treedepth   = 12,
  refresh         = 0
)

# ---------- retry helper ----------
fit_with_retry <- function(stan_data, truth, model, fit_args, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    fit <- do.call(model$sample, c(
      list(data = stan_data,
           init = make_init_fun(stan_data, truth)),
      fit_args
    ))
    max_rhat <- max(fit$summary(variables = c("tau", "tau_std", "beta_std","sigma_std"))$rhat,
                    na.rm = TRUE)
    message(sprintf("Attempt %d: max_rhat = %.4f", attempt, max_rhat))  # always print
    if (max_rhat <= 1.01) return(fit)
    message("Rhat > 1.01 — retrying...")
  }
  message("Max retries reached")
  fit
}

# ---------- one replication ----------
run_once <- function(rep_id, mod, fit_args, compute_class_acc = TRUE) {
  seed <- 10000 + rep_id
  dat  <- simulate_dataset(seed = seed)
  
  fit <- fit_with_retry(          # <-- replaces the old do.call block
    stan_data  = dat$stan_data,
    truth      = dat$truth,
    model      = mod,
    fit_args   = fit_args
  )
  
  out <- extract_metrics(fit, dat$truth, compute_class_acc = compute_class_acc)
  out$rep  <- rep_id
  out$seed <- seed
  out
}

# ---------- many replications ----------
run_many <- function(R, mod, fit_args, compute_class_acc = TRUE) {
  res <- lapply(1:R, function(r) {
    tryCatch(
      run_once(r, mod, fit_args, compute_class_acc = compute_class_acc),
      error = function(e) data.frame(rep = r, seed = 10000 + r, error = conditionMessage(e))
    )
  })
  bind_rows(res)
}


