# functions_binary.R
# Binary outcome BLCA simulation functions.
# Mirrors functions.R exactly — differences noted with # BINARY comments.

options(mc.cores = parallel::detectCores())

# ---------- DGP function ----------
# BINARY: logistic outcome, beta on log-odds scale, no sigma
simulate_dataset_binary <- function(seed,
                                    N = 800, M = 10, C = 3,
                                    pk = c(0.3, 0.3, 0.4),
                                    mu = list(c(-0.6, -0.6), c(0, 0.6), c(0.6, 0)),
                                    Sigma_within = diag(2) * 0.25,
                                    p_bio = c(0.05, 0.5, 0.95),
                                    beta_common = c(0.3, 0.2, -0.1, 0.1, 0.2),
                                    beta_intercept = qlogis(0.375),
                                    tau_true = c(log(0.25), 0, log(3)),
                                    p_treat = 0.5) {
  set.seed(seed)

  # 1. Generate X unconditionally — identical to continuous model
  comp      <- sample(1:C, N, replace = TRUE, prob = pk)
  x1_level1 <- do.call(rbind, lapply(comp, function(k)
    MASS::mvrnorm(1, mu[[k]], Sigma = Sigma_within)
  ))

  x2_level2 <- MASS::mvrnorm(M, mu = c(0, 0),
                             Sigma = matrix(c(1, -1, -1, 1), nrow = 2))
  site <- sample(1:M, N, replace = TRUE)

  p_bio_marginal <- sum(pk * p_bio)
  B <- rbinom(N, 1, prob = p_bio_marginal)

  X <- cbind(1, x1_level1, x2_level2[site, ], B)
  K <- ncol(X)

  # 2. p(class | x) via Bayes — identical to continuous model
  log_p_class_given_x <- matrix(NA_real_, nrow = N, ncol = C)
  for (k in 1:C) {
    log_p_class_given_x[, k] <- log(pk[k]) +
      mvtnorm::dmvnorm(x1_level1, mean = mu[[k]], sigma = Sigma_within, log = TRUE) +
      dbinom(B, size = 1, prob = p_bio[k], log = TRUE)
  }
  pi_mat <- t(apply(log_p_class_given_x, 1, function(lp) {
    lp <- lp - max(lp); exp(lp) / sum(exp(lp))
  }))
  class <- apply(pi_mat, 1, function(p) sample(seq_len(C), 1, prob = p))

  # 3. Treatment — identical to continuous model
  T <- rbinom(N, 1, p_treat)

  # 4. Binary outcome via logistic model
  # log-odds = X_i * beta + tau[class] * T_i
  # Y_i ~ Bernoulli(plogis(log-odds))
  beta_full <- c(beta_intercept, beta_common)
  Y <- integer(N)
  for (n in 1:N) {
    log_odds_n <- sum(X[n, ] * beta_full) + tau_true[class[n]] * T[n]
    Y[n]       <- rbinom(1, 1, prob = plogis(log_odds_n))
  }

  # 5. Stan data
  stan_data <- list(
    N = N, K = K, C = C,
    x = X,
    y = as.integer(Y),   # BINARY: integer not numeric
    t = as.integer(T)
  )

  # 6. Truth
  truth <- list(
    class         = class,
    pi_mat        = pi_mat,
    tau           = tau_true,        # log-OR scale
    tau_or        = exp(tau_true),   # odds ratios
    sigma         = NULL,            # no sigma for binary
    N             = N,
    C             = C,
    t             = T,
    p_control     = plogis(beta_intercept),
    p_harm_trt    = plogis(beta_intercept + tau_true[1]),
    p_null_trt    = plogis(beta_intercept + tau_true[2]),
    p_benefit_trt = plogis(beta_intercept + tau_true[3]),
    event_rate    = mean(Y)
  )

  list(stan_data = stan_data, truth = truth)
}

# ---------- output function ----------

get_pbar_from_fit <- function(fit, N, C) {
  pc_draws <- posterior::as_draws_matrix(fit$draws("p_class"))
  cn       <- colnames(pc_draws)
  keep     <- grepl("^p_class\\[", cn)
  pc_draws <- pc_draws[, keep, drop = FALSE]
  cn       <- colnames(pc_draws)
  idx      <- do.call(rbind, regmatches(cn, regexec("^p_class\\[([0-9]+),([0-9]+)\\]$", cn)))
  i_idx    <- as.integer(idx[, 2])
  c_idx    <- as.integer(idx[, 3])
  p_means  <- colMeans(pc_draws)
  pbar     <- matrix(NA_real_, nrow = N, ncol = C)
  for (j in seq_along(p_means)) pbar[i_idx[j], c_idx[j]] <- p_means[j]
  pbar
}

get_prior_from_fit <- function(fit, N, C) {
  cp_draws <- posterior::as_draws_matrix(fit$draws("class_prior"))
  cn       <- colnames(cp_draws)
  keep     <- grepl("^class_prior\\[", cn)
  cp_draws <- cp_draws[, keep, drop = FALSE]
  cn       <- colnames(cp_draws)
  idx      <- do.call(rbind, regmatches(cn, regexec("^class_prior\\[([0-9]+),([0-9]+)\\]$", cn)))
  i_idx    <- as.integer(idx[, 2])
  c_idx    <- as.integer(idx[, 3])
  p_means  <- colMeans(cp_draws)
  prior_mat <- matrix(NA_real_, nrow = N, ncol = C)
  for (j in seq_along(p_means)) prior_mat[i_idx[j], c_idx[j]] <- p_means[j]
  prior_mat
}

extract_metrics <- function(fit, truth, compute_class_acc = FALSE) {

  # BINARY: pull tau_std (log-OR scale) and tau_or (odds ratio scale)
  # No tau back-transform and no sigma in the binary model
  s_tau    <- fit$summary("tau_std")
  s_tau_or <- fit$summary("tau_or")
  C        <- nrow(s_tau)

  # Diagnostics — no sigma_std or tau in binary model
  s_key <- fit$summary(variables = c("tau_std", "tau_mid", "delta_harm",
                                     "delta_benefit", "beta_std",
                                     "gamma_sd", "gamma_intercept_raw",
                                     "gamma_slopes_raw"))
  max_rhat     <- max(s_key$rhat,     na.rm = TRUE)
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

  # Tau performance (log-OR scale) + odds ratios
  for (cls in 1:C) {
    est <- s_tau$median[cls]
    lo  <- s_tau$q5[cls]
    hi  <- s_tau$q95[cls]

    out[[paste0("tau_", cls, "_med")]]     <- est
    out[[paste0("tau_", cls, "_q5")]]      <- lo
    out[[paste0("tau_", cls, "_q95")]]     <- hi
    out[[paste0("tau_", cls, "_bias")]]    <- est - truth$tau[cls]
    out[[paste0("tau_", cls, "_cover90")]] <- as.integer(lo <= truth$tau[cls] && truth$tau[cls] <= hi)
    out[[paste0("tau_", cls, "_width")]]   <- hi - lo
    out[[paste0("tau_or_", cls, "_med")]]  <- s_tau_or$median[cls]
  }

  # Class recovery — identical to continuous model
  if (compute_class_acc) {
    N  <- length(truth$class)
    C2 <- length(truth$tau)

    if (is.null(truth$t)) stop("truth$t is missing.")

    pbar            <- get_pbar_from_fit(fit, N = N, C = C2)
    class_hat       <- max.col(pbar, ties.method = "first")
    prior_mat       <- get_prior_from_fit(fit, N = N, C = C2)
    class_hat_prior <- max.col(prior_mat, ties.method = "first")

    idx_trt <- truth$t == 1
    idx_ctl <- truth$t == 0

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

    for (cc in 1:C2) {
      idx_class     <- truth$class == cc
      idx_class_trt <- truth$class == cc & truth$t == 1
      idx_class_ctl <- truth$class == cc & truth$t == 0

      out[[paste0("mean_pclass_", cc, "_correct")]] <-
        if (sum(idx_class) > 0) mean(pbar[idx_class, cc]) else NA_real_
      out[[paste0("mean_pclass_", cc, "_correct_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(pbar[idx_class_trt, cc]) else NA_real_
      out[[paste0("mean_pclass_", cc, "_correct_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(pbar[idx_class_ctl, cc]) else NA_real_

      out[[paste0("mean_prior_", cc, "_correct")]] <-
        if (sum(idx_class) > 0) mean(prior_mat[idx_class, cc]) else NA_real_
      out[[paste0("mean_prior_", cc, "_correct_trt")]] <-
        if (sum(idx_class_trt) > 0) mean(prior_mat[idx_class_trt, cc]) else NA_real_
      out[[paste0("mean_prior_", cc, "_correct_ctl")]] <-
        if (sum(idx_class_ctl) > 0) mean(prior_mat[idx_class_ctl, cc]) else NA_real_
    }

    entropy    <- apply(pbar, 1, function(p) -sum(p * log(p + 1e-10)))
    max_pclass <- apply(pbar, 1, max)

    for (cc in 1:C2) {
      for (trt in c(0, 1)) {
        arm_label <- ifelse(trt == 1, "trt", "ctl")
        idx_cell  <- truth$class == cc & truth$t == trt
        out[[paste0("entropy_class", cc, "_", arm_label)]] <-
          if (sum(idx_cell) > 0) mean(entropy[idx_cell],    na.rm = TRUE) else NA_real_
        out[[paste0("maxprob_class", cc, "_", arm_label)]] <-
          if (sum(idx_cell) > 0) mean(max_pclass[idx_cell], na.rm = TRUE) else NA_real_
      }
    }

    misclass_pairs <- list(c(1,2), c(1,3), c(2,1), c(2,3), c(3,1), c(3,2))
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
# BINARY: no sigma_std, tau on log-odds scale directly (no sd_y division)
make_init_fun_binary <- function(stan_data, truth) {
  K <- stan_data$K
  C <- stan_data$C

  function() {
    list(
      beta_std            = c(qlogis(0.375), rnorm(K-1, 0, 0.1)),
      gamma_slopes_raw    = matrix(rnorm((K-1) * (C-1), 0, 0.1), K-1, C-1),
      gamma_intercept_raw = rnorm(C-1, 0, 0.1),
      gamma_sd            = 0.3,
      tau_mid             = truth$tau[2],
      delta_harm          = truth$tau[2] - truth$tau[1],
      delta_benefit       = truth$tau[3] - truth$tau[2]
    )
  }
}

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
# BINARY: checks tau_std not tau; no sigma_std
fit_with_retry_bin <- function(stan_data, truth, model, fit_args, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    fit <- do.call(model$sample, c(
      list(data = stan_data,
           init = make_init_fun_binary(stan_data, truth)),
      fit_args
    ))
    max_rhat <- max(fit$summary(variables = c("tau_std", "beta_std",
                                              "gamma_sd", "gamma_intercept_raw",
                                              "gamma_slopes_raw"))$rhat,
                    na.rm = TRUE)
    message(sprintf("Attempt %d: max_rhat = %.4f", attempt, max_rhat))
    if (max_rhat <= 1.01) return(fit)
    message("Rhat > 1.01 — retrying...")
  }
  message("Max retries reached")
  fit
}

# ---------- summary functions ----------
summarise_tau <- function(data, cls) {
  med  <- data[[paste0("tau_", cls, "_med")]]
  lo   <- data[[paste0("tau_", cls, "_q5")]]
  hi   <- data[[paste0("tau_", cls, "_q95")]]
  bias <- data[[paste0("tau_", cls, "_bias")]]
  cov  <- data[[paste0("tau_", cls, "_cover90")]]

  tibble::tibble(
    class      = cls,
    mean_est   = mean(med,  na.rm = TRUE),
    emp_sd     = sd(med,    na.rm = TRUE),
    mean_bias  = mean(bias, na.rm = TRUE),
    rmse       = sqrt(mean(bias^2, na.rm = TRUE)),
    coverage90 = mean(cov,  na.rm = TRUE),
    mean_width = mean(hi - lo, na.rm = TRUE),
    mean_or    = mean(data[[paste0("tau_or_", cls, "_med")]], na.rm = TRUE)
  )
}

summarise_misclass <- function(data) {
  pairs <- list(
    c(1, 2, "Harm -> Null"),
    c(1, 3, "Harm -> Benefit (catastrophic)"),
    c(2, 1, "Null -> Harm"),
    c(2, 3, "Null -> Benefit"),
    c(3, 1, "Benefit -> Harm (catastrophic)"),
    c(3, 2, "Benefit -> Null")
  )
  tibble::tibble(
    misclassification = sapply(pairs, `[[`, 3),
    rate = sapply(pairs, function(p) {
      col <- paste0("misclass_", p[1], "_to_", p[2], "_trt")
      if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
    })
  ) %>%
    dplyr::mutate(rate_pct = scales::percent(rate, accuracy = 0.1))
}

summarise_confidence <- function(data) {
  expand.grid(class = 1:3, arm = c("trt", "ctl")) %>%
    dplyr::mutate(
      true_class    = factor(class, labels = c("Harm", "Null", "Benefit")),
      treated       = factor(arm, levels = c("trt", "ctl"),
                             labels = c("Treated", "Control")),
      mean_entropy  = mapply(function(c, a) {
        col <- paste0("entropy_class", c, "_", a)
        if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
      }, class, arm),
      mean_max_prob = mapply(function(c, a) {
        col <- paste0("maxprob_class", c, "_", a)
        if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
      }, class, arm)
    ) %>%
    dplyr::select(true_class, treated, mean_entropy, mean_max_prob) %>%
    dplyr::arrange(true_class, treated)
}

summarise_all <- function(data) {
  tau_summary   <- dplyr::bind_rows(lapply(1:3, function(cls) summarise_tau(data, cls)))

  diag_summary  <- data %>%
    dplyr::summarise(
      conv_rate          = mean(max_rhat < 1.01 & n_divergent == 0 & min_ess_bulk > 200, na.rm = TRUE),
      mean_max_rhat      = mean(max_rhat,     na.rm = TRUE),
      median_runtime_sec = median(runtime_sec, na.rm = TRUE)
    )

  class_summary <- data %>%
    dplyr::summarise(
      mean_class_acc     = mean(class_acc,     na.rm = TRUE),
      mean_class_acc_trt = mean(class_acc_trt, na.rm = TRUE),
      mean_class_acc_ctl = mean(class_acc_ctl, na.rm = TRUE)
    )

  per_class_summary <- data %>%
    dplyr::summarise(
      class_1_acc     = mean(class_1_acc,     na.rm = TRUE),
      class_2_acc     = mean(class_2_acc,     na.rm = TRUE),
      class_3_acc     = mean(class_3_acc,     na.rm = TRUE),
      class_1_acc_trt = mean(class_1_acc_trt, na.rm = TRUE),
      class_2_acc_trt = mean(class_2_acc_trt, na.rm = TRUE),
      class_3_acc_trt = mean(class_3_acc_trt, na.rm = TRUE)
    )

  prior_summary <- data %>%
    dplyr::summarise(
      mean_prior_acc     = mean(prior_acc,     na.rm = TRUE),
      mean_prior_acc_trt = mean(prior_acc_trt, na.rm = TRUE),
      mean_prior_acc_ctl = mean(prior_acc_ctl, na.rm = TRUE)
    )

  per_class_prior_summary <- data %>%
    dplyr::summarise(
      prior_class_1_acc     = mean(prior_class_1_acc,     na.rm = TRUE),
      prior_class_2_acc     = mean(prior_class_2_acc,     na.rm = TRUE),
      prior_class_3_acc     = mean(prior_class_3_acc,     na.rm = TRUE),
      prior_class_1_acc_trt = mean(prior_class_1_acc_trt, na.rm = TRUE),
      prior_class_2_acc_trt = mean(prior_class_2_acc_trt, na.rm = TRUE),
      prior_class_3_acc_trt = mean(prior_class_3_acc_trt, na.rm = TRUE)
    )

  mean_prob_summary <- data %>%
    dplyr::summarise(
      post_prob_class1  = mean(mean_pclass_1_correct, na.rm = TRUE),
      post_prob_class2  = mean(mean_pclass_2_correct, na.rm = TRUE),
      post_prob_class3  = mean(mean_pclass_3_correct, na.rm = TRUE),
      prior_prob_class1 = mean(mean_prior_1_correct,  na.rm = TRUE),
      prior_prob_class2 = mean(mean_prior_2_correct,  na.rm = TRUE),
      prior_prob_class3 = mean(mean_prior_3_correct,  na.rm = TRUE)
    )

  misclass_summary    <- summarise_misclass(data)
  confidence_summary  <- summarise_confidence(data)

  list(
    tau_summary             = tau_summary,
    diag_summary            = diag_summary,
    class_summary           = class_summary,
    per_class_summary       = per_class_summary,
    prior_summary           = prior_summary,
    per_class_prior_summary = per_class_prior_summary,
    mean_prob_summary       = mean_prob_summary,
    misclass_summary        = misclass_summary,
    confidence_summary      = confidence_summary
  )
}
