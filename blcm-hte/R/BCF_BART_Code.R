#Install rstan and dependencies (if you haven't already)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores
cmdstanr::set_cmdstan_path("/Users/yournamehere/.cmdstan/cmdstan-2.37.0")

#Set working directory to dissertation
setwd("~/LCM/Simulations_Code")

# Load the cmdrstanr library & needed packages
#library(cmdstanr) - not available for R 4.6.0
library(posterior)
library(bayesplot)
library(ggplot2)
library(matrixStats)
library(MASS)
library(dplyr)
library(bcf)
library(dbarts)

mod <- cmdstan_model("BLCA_fixedtrt_CmdStanR_v2.stan")
mod_reduced <- cmdstan_model("BLCA_reduced_CmdStanR_v1.stan")

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

# ---------- Output functions ----------
extract_metrics_bcf <- function(fit_bcf, truth, runtime_sec) {
  
  # Individual posterior treatment effects: matrix nsim x N
  tau_draws <- fit_bcf$tau   # rows = draws, cols = individuals
  
  out <- list()
  out$runtime_sec <- runtime_sec
  
  for (cc in 1:truth$C) {
    idx <- truth$class == cc
    
    # Posterior draws for this class: average across individuals within class
    tau_class_draws <- rowMeans(tau_draws[, idx, drop = FALSE])
    
    # Point estimate and uncertainty
    tau_med   <- median(tau_class_draws)
    tau_q5    <- quantile(tau_class_draws, 0.05)
    tau_q95   <- quantile(tau_class_draws, 0.95)
    tau_true  <- truth$tau[cc]
    
    out[[paste0("tau_", cc, "_med")]]     <- tau_med
    out[[paste0("tau_", cc, "_q5")]]      <- tau_q5
    out[[paste0("tau_", cc, "_q95")]]     <- tau_q95
    out[[paste0("tau_", cc, "_bias")]]    <- tau_med - tau_true
    out[[paste0("tau_", cc, "_cover90")]] <- as.integer(tau_q5 <= tau_true & tau_true <= tau_q95)
    out[[paste0("tau_", cc, "_width")]]   <- as.numeric(tau_q95 - tau_q5)
    out[[paste0("tau_", cc, "_rmse")]]    <- sqrt(mean((tau_class_draws - tau_true)^2))
  }
  
  as.data.frame(out)
}
extract_metrics_bart <- function(tau_draws_bart, truth, runtime_sec) {
  
  # Compute individual tau draws inside the function
  X_treated      <- cbind(X, 1)
  X_control      <- cbind(X, 0)
  tau_draws_bart <- predict(fit_bart, X_treated) - predict(fit_bart, X_control)
  
  out <- list()
  out$runtime_sec <- runtime_sec
  
  for (cc in 1:truth$C) {
    idx <- truth$class == cc
    
    # Average across individuals within true class for each draw
    tau_class_draws <- rowMeans(tau_draws_bart[, idx, drop = FALSE])
    
    tau_med  <- median(tau_class_draws)
    tau_q5   <- quantile(tau_class_draws, 0.05)
    tau_q95  <- quantile(tau_class_draws, 0.95)
    tau_true <- truth$tau[cc]
    
    out[[paste0("tau_", cc, "_med")]]     <- tau_med
    out[[paste0("tau_", cc, "_q5")]]      <- tau_q5
    out[[paste0("tau_", cc, "_q95")]]     <- tau_q95
    out[[paste0("tau_", cc, "_bias")]]    <- tau_med - tau_true
    out[[paste0("tau_", cc, "_cover90")]] <- as.integer(tau_q5 <= tau_true & tau_true <= tau_q95)
    out[[paste0("tau_", cc, "_width")]]   <- as.numeric(tau_q95 - tau_q5)
    out[[paste0("tau_", cc, "_rmse")]]    <- sqrt(mean((tau_class_draws - tau_true)^2))
  }
  
  as.data.frame(out)
}

#simulate data set for testing
dat <- simulate_dataset(10001)

#Setup for BCF
X      <- dat$stan_data$x[, -1]   # drop intercept column
y      <- dat$stan_data$y
t      <- dat$stan_data$t
pihat  <- rep(0.5, nrow(X))   # known randomization probability

#Fit BCF model 
fit_bcf <- bcf(y          = y,
               z          = t,
               x_control  = X,
               x_moderate = X,
               pihat      = pihat,
               nburn      = 1000,
               nsim       = 1000)

#Fit BART model 
fit_bart <- bart(x.train = cbind(X, t),   # covariates + treatment as predictors
                 y.train = y,
                 ntree   = 200,
                 ndpost  = 1000,
                 nskip   = 1000,
                 keeptrees = TRUE)

out_bcf  <- extract_metrics_bcf(fit_bcf, dat$truth, runtime_sec=96)
out_bart <- extract_metrics_bart(tau_draws_bart, dat$truth, runtime_sec=10)

# Rename to distinguish
names(out_bcf)  <- paste0(names(out_bcf),  "_bcf")
names(out_bart) <- paste0(names(out_bart), "_bart")

result <- cbind(out_bcf, out_bart)
print(result)

# Both BCF and BART recover marginal ATE reasonably but fail to identify 
# discrete subgroup-specific effects that are clinically actionable 

######################### RESULTS #################################

#downloaded aggregate results from HPC
results_comp <- readRDS("~/Desktop/PhD/Dissertation/Paper 2 - LCM/lcm_sim/results_combined_comp.rds")

#################### Summary Functions ################################
summarise_comp <- function(results_comp) {
  # For each condition and class, compute mean est, bias, RMSE for BCF and BART
  results_comp %>%
    group_by(cond_id, N, pk_str, tau_str) %>%
    summarise(
      # BCF
      tau_1_mean_bcf  = mean(tau_1_med_bcf,  na.rm = TRUE),
      tau_1_bias_bcf  = mean(tau_1_bias_bcf, na.rm = TRUE),
      tau_1_rmse_bcf  = sqrt(mean(tau_1_rmse_bcf^2, na.rm = TRUE)),
      
      tau_2_mean_bcf  = mean(tau_2_med_bcf,  na.rm = TRUE),
      tau_2_bias_bcf  = mean(tau_2_bias_bcf, na.rm = TRUE),
      tau_2_rmse_bcf  = sqrt(mean(tau_2_rmse_bcf^2, na.rm = TRUE)),
      
      tau_3_mean_bcf  = mean(tau_3_med_bcf,  na.rm = TRUE),
      tau_3_bias_bcf  = mean(tau_3_bias_bcf, na.rm = TRUE),
      tau_3_rmse_bcf  = sqrt(mean(tau_3_rmse_bcf^2, na.rm = TRUE)),
      
      # BART
      tau_1_mean_bart = mean(tau_1_med_bart,  na.rm = TRUE),
      tau_1_bias_bart = mean(tau_1_bias_bart, na.rm = TRUE),
      tau_1_rmse_bart = sqrt(mean(tau_1_rmse_bart^2, na.rm = TRUE)),
      
      tau_2_mean_bart = mean(tau_2_med_bart,  na.rm = TRUE),
      tau_2_bias_bart = mean(tau_2_bias_bart, na.rm = TRUE),
      tau_2_rmse_bart = sqrt(mean(tau_2_rmse_bart^2, na.rm = TRUE)),
      
      tau_3_mean_bart = mean(tau_3_med_bart,  na.rm = TRUE),
      tau_3_bias_bart = mean(tau_3_bias_bart, na.rm = TRUE),
      tau_3_rmse_bart = sqrt(mean(tau_3_rmse_bart^2, na.rm = TRUE)),
      
      n_reps = n(),
      .groups = "drop"
    )
}

comp_summary <- summarise_comp(results_comp)
print(comp_summary)





