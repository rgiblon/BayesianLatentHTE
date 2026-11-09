#Install rstan and dependencies (if you haven't already)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores
cmdstanr::set_cmdstan_path("/Users/rachelpim/.cmdstan/cmdstan-2.37.0")

#Set working directory to dissertation
setwd("~/Desktop/PhD/Dissertation/Paper 2 - LCM/Simulations_Code")

# Load the cmdrstanr library & needed packages
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(matrixStats)
library(MASS)
library(dplyr)
library(bcf)
library(BART)
library(speff2trial) #<--- Trial data is in this package
library(tidyr); 
library(loo)

## ------------------------------------------------------------------
## Exploring ACTG175 (speff2trial) ahead of the BLCM applied analysis
## ------------------------------------------------------------------
data(ACTG175)
df <- ACTG175

## 1. Basic structure -------------------------------------------------
str(df)
dim(df)
colnames(df)

## 2. Arm sizes --------------------------------------------------------
# "arms" is the 4-level treatment factor (0-3); "treat" is a collapsed
# binary indicator in some versions of the package. Check both exist:
table(df$arms)
if ("treat" %in% colnames(df)) table(df$treat)

## 3. Reduce to Yuki et al.'s two-arm comparison ------------------------
# Their split: zidovudine monotherapy (control) vs. zidovudine + didanosine (test)
# Confirm the arms coding from ?ACTG175 or the table() output above, then adjust
# the values below (this assumes 0 = ZDV alone, 1 = ZDV+ddI, consistent with
# the most common coding, but VERIFY against the codebook before trusting it):
df_2arm <- subset(df, arms %in% c(0, 1))
df_2arm$treat_binary <- ifelse(df_2arm$arms == 1, 1, 0)  # 1 = ZDV+ddI, 0 = ZDV monotherapy
table(df_2arm$treat_binary)

## 4. Outcome of interest: cd420 (CD4 count at 20±5 weeks) --------------
summary(df_2arm$cd420)
sum(is.na(df_2arm$cd420))
hist(df_2arm$cd420, main = "cd420 distribution (2-arm subset)", xlab = "CD4 count")

## 5. Candidate covariates (Yuki et al.'s 15-covariate set) -------------
covars <- c("age","wtkg","hemo","homo","drugs","karnof","oprior",
            "z30","zprior","race","gender","str2","symptom","cd40","cd80")
covars_present <- covars[covars %in% colnames(df_2arm)]
setdiff(covars, covars_present)  # flag any that don't exist under this name

summary(df_2arm[, covars_present])

## Drop zero-variance zprior, finalize covariate set, check cd420 shape
covars_final <- setdiff(covars_present, "zprior")
covars_final

# sanity check: confirm all remaining covariates actually vary
sapply(df_2arm[, covars_final], function(x) length(unique(x)))

## cd420 shape check
hist(df_2arm$cd420, breaks = 30, main = "cd420 distribution (n = 1054)",
     xlab = "CD4 count at 20±5 weeks")

# quick skewness check (base R, no extra package needed)
m <- mean(df_2arm$cd420); s <- sd(df_2arm$cd420)
skew <- mean(((df_2arm$cd420 - m) / s)^3)
skew

# also worth a quick look split by arm, in case the two treatment arms
# have differently-shaped distributions (would matter for your DGP assumptions)
par(mfrow = c(1, 2))
hist(df_2arm$cd420[df_2arm$treat_binary == 0], breaks = 20,
     main = "cd420, ZDV monotherapy", xlab = "CD4 count")
hist(df_2arm$cd420[df_2arm$treat_binary == 1], breaks = 20,
     main = "cd420, ZDV+ddI", xlab = "CD4 count")
par(mfrow = c(1, 1))

## 6. Missingness check --------------------------------------------------
sapply(df_2arm[, c("cd420", covars_final)], function(x) sum(is.na(x)))

## 7. Quick covariate x treatment balance check ------------------------
# sanity check that randomization looks balanced on baseline covariates
for (v in covars_present) {
  cat("\n---", v, "---\n")
  print(tapply(df_2arm[[v]], df_2arm$arms, function(x) round(mean(x, na.rm = TRUE), 2)))
}

## =====================================================================
## BLCM applied analysis — ACTG175 (continuous outcome, K = 3)
##
## Replaces the simulation harness. No ground truth here: no true class
## labels, no true tau, so no bias / coverage / classification-accuracy.
## What we keep: convergence diagnostics, posterior summaries, and
## class-membership probabilities.
## =====================================================================

mod <- cmdstan_model("BLCA_fixedtrt_CmdStanR_v3.stan") #tau_mid loosened to N(0,0.5)
#mod <- cmdstan_model("BLCA_fixedtrt_CmdStanR_v2.stan") #tau_mid tighten to N(0,0.2)

## ---- 1. Data prep (replaces simulate_dataset) ------------------------
prep_actg175 <- function(outcome = "cd420",
                         covars = c("age","wtkg","hemo","homo","drugs",
                                    "karnof","oprior","z30","race","gender",
                                    "str2","symptom","cd40","cd80"),
                         C = 3) {
  data(ACTG175)
  df <- subset(ACTG175, arms %in% c(0, 1))
  df$treat_binary <- ifelse(df$arms == 1, 1, 0)   # 1 = ZDV+ddI, 0 = ZDV mono
  # zprior excluded — constant (== 1) within this two-arm subset
  
  stopifnot(all(covars %in% names(df)))
  if (any(is.na(df[, c(outcome, covars)]))) stop("Unexpected missingness.")
  
  # intercept in col 1; Stan standardizes internally and protects col 1,
  # so pass raw scale
  X <- cbind(1, as.matrix(df[, covars]))
  colnames(X) <- c("(Intercept)", covars)
  
  stan_data <- list(N = nrow(df), K = ncol(X), C = C, x = X,
                    y = as.numeric(df[[outcome]]),
                    t = as.integer(df$treat_binary))
  list(stan_data = stan_data, df = df, covars = covars, outcome = outcome)
}

## ---- 2. Init (replaces make_init_fun) --------------------------------
## Sim version initialised tau from truth$tau / truth$sd_y. With no truth,
## anchor tau_mid at 0 and size the gaps from the observed naive effect —
## neutral, assumes classes straddle the marginal effect.
make_init_fun_real <- function(stan_data) {
  K <- stan_data$K; C <- stan_data$C
  y <- stan_data$y; t <- stan_data$t; sd_y <- sd(y)
  naive_tau_std <- (mean(y[t == 1]) - mean(y[t == 0])) / sd_y
  gap0 <- max(abs(naive_tau_std), 0.3)   # floor so gaps don't start at ~0
  
  function() list(
    beta_std            = rnorm(K, 0, 0.1),
    gamma_slopes_raw    = matrix(rnorm((K-1)*(C-1), 0, 0.1), K-1, C-1),
    gamma_intercept_raw = rnorm(C-1, 0, 0.1),
    gamma_sd            = 0.3,
    tau_mid             = 0,
    delta_harm          = gap0 * runif(1, 0.5, 1.5),
    delta_benefit       = gap0 * runif(1, 0.5, 1.5),
    sigma_std           = 1
  )
}

## ---- 3. Sampling args ------------------------------------------------
fit_args <- list(chains = 4, parallel_chains = 4,
                 iter_warmup = 2000, iter_sampling = 2000,
                 adapt_delta = 0.90, max_treedepth = 12, refresh = 200)

## ---- 4. Fit with retry ------------------------------------------------
fit_real <- function(stan_data, model, fit_args, seed = 1, max_retries = 3) {
  for (attempt in seq_len(max_retries)) {
    fit <- do.call(model$sample, c(
      list(data = stan_data, init = make_init_fun_real(stan_data),
           seed = seed + attempt - 1), fit_args))
    s <- fit$summary(c("tau","tau_std","beta_std","gamma_sd","sigma_std",
                       "gamma_intercept_raw","gamma_slopes_raw"))
    max_rhat <- max(s$rhat, na.rm = TRUE)
    min_ess  <- min(s$ess_bulk, na.rm = TRUE)
    n_div    <- sum(fit$sampler_diagnostics()[,, "divergent__"], na.rm = TRUE)
    message(sprintf("Attempt %d: max_rhat = %.4f, min_ess = %.0f, n_div = %d",
                    attempt, max_rhat, min_ess, n_div))
    if (max_rhat < 1.01 && min_ess > 200 && n_div == 0) return(fit)
    message("Convergence criteria not met, retrying with a new seed...")
  }
  message("Max retries reached, returning last fit.")
  fit
}

## ---- 5. Diagnostics + posterior summary (replaces extract_metrics) ----
summarise_real_fit <- function(fit, dat) {
  key <- c("tau_std","tau_mid","delta_harm","delta_benefit",
           "sigma_std","beta_std","gamma_sd")
  s <- fit$summary(key)
  diag <- list(
    max_rhat     = max(s$rhat, na.rm = TRUE),
    min_ess_bulk = min(s$ess_bulk, na.rm = TRUE),
    n_divergent  = sum(fit$sampler_diagnostics()[,, "divergent__"], na.rm = TRUE),
    runtime_sec  = sum(fit$time()$total))
  
  tau_tab <- fit$summary("tau") %>%
    mutate(class = c("Harm","Null","Benefit")) %>%
    select(class, median, q5, q95, rhat, ess_bulk)
  
  beta_tab <- fit$summary("beta") %>%
    mutate(term = colnames(dat$stan_data$x)) %>%
    select(term, median, q5, q95, rhat)
  
  sigma_tab <- fit$summary("sigma") %>% select(variable, median, q5, q95)
  
  list(diagnostics = diag, tau = tau_tab, beta = beta_tab, sigma = sigma_tab)
}

## ---- 6. Class membership (these two carry over unchanged) ------------
get_pbar_from_fit <- function(fit, N, C) {
  pc <- posterior::as_draws_matrix(fit$draws("p_class"))
  pc <- pc[, grepl("^p_class\\[", colnames(pc)), drop = FALSE]
  idx <- do.call(rbind, regmatches(colnames(pc),
                                   regexec("^p_class\\[([0-9]+),([0-9]+)\\]$", colnames(pc))))
  m <- matrix(NA_real_, N, C)
  m[cbind(as.integer(idx[,2]), as.integer(idx[,3]))] <- colMeans(pc); m
}

get_prior_from_fit <- function(fit, N, C) {
  cp <- posterior::as_draws_matrix(fit$draws("class_prior"))
  cp <- cp[, grepl("^class_prior\\[", colnames(cp)), drop = FALSE]
  idx <- do.call(rbind, regmatches(colnames(cp),
                                   regexec("^class_prior\\[([0-9]+),([0-9]+)\\]$", colnames(cp))))
  m <- matrix(NA_real_, N, C)
  m[cbind(as.integer(idx[,2]), as.integer(idx[,3]))] <- colMeans(cp); m
}

## Descriptive class summary — replaces the accuracy metrics.
## Modal posterior class. Per your Stan comment, assignment is only really
## informative in the TREATED arm (tau doesn't enter the likelihood at t=0).
describe_classes <- function(fit, dat) {
  N <- dat$stan_data$N; C <- dat$stan_data$C
  pbar <- get_pbar_from_fit(fit, N, C)
  prior_mat <- get_prior_from_fit(fit, N, C)
  class_hat <- max.col(pbar, ties.method = "first")
  labels <- c("Harm","Null","Benefit")
  
  df <- dat$df
  df$class_hat <- factor(labels[class_hat], levels = labels)
  df$max_pprob <- apply(pbar, 1, max)
  df$entropy   <- apply(pbar, 1, function(p) -sum(p * log(p + 1e-10)))
  
  sizes <- df %>% count(class_hat, treat_binary) %>%
    pivot_wider(names_from = treat_binary, values_from = n,
                names_prefix = "arm_", values_fill = 0)
  
  confidence <- df %>% group_by(class_hat) %>%
    summarise(n = n(), mean_max_pprob = mean(max_pprob),
              mean_entropy = mean(entropy), .groups = "drop")
  
  # this becomes the "class characterization" table in the paper
  profile_trt <- df %>% filter(treat_binary == 1) %>%
    group_by(class_hat) %>%
    summarise(across(all_of(dat$covars), ~ round(mean(.x), 2)),
              n = n(), .groups = "drop")
  
  # membership model alone vs. full posterior — low agreement means y is
  # driving assignment more than covariates are
  agreement <- mean(max.col(prior_mat, ties.method = "first") == class_hat)
  
  list(sizes = sizes, confidence = confidence, profile_treated = profile_trt,
       prior_post_agree = agreement, df = df, pbar = pbar, prior_mat = prior_mat)
}

## ---- 7. Seed stability — the real-data analogue of "run 100 reps" -----
## One dataset, a few refits. Confirms the posterior is stable under
## sampler stochasticity; a handful of seeds is plenty.
check_seed_stability <- function(dat, model, fit_args, seeds = 1:5) {
  res <- lapply(seeds, function(s) {
    set.seed(s)
    fit <- fit_real(dat$stan_data, model, fit_args, seed = s, max_retries = 1)
    cls <- describe_classes(fit, dat)
    list(tau = fit$summary("tau")$median,
         class_hat = as.integer(cls$df$class_hat))
  })
  ref <- res[[1]]$class_hat
  data.frame(seed = seeds,
             tau_harm = sapply(res, function(r) r$tau[1]),
             tau_null = sapply(res, function(r) r$tau[2]),
             tau_ben  = sapply(res, function(r) r$tau[3]),
             n_harm   = sapply(res, function(r) sum(r$class_hat == 1)),
             n_moved_vs_seed1 = sapply(res, function(r) sum(r$class_hat != ref)))
}

## ---- 8. Run -----------------------------------------------------------
set.seed(20261005)   # R-side inits; Stan's seed is set inside fit_real()
dat <- prep_actg175()
fit <- fit_real(dat$stan_data, mod, fit_args, seed = 1)
fit$save_object("actg175_blcm_k3_fit.rds")

out <- summarise_real_fit(fit, dat)
print(out$diagnostics); print(out$tau); print(out$beta); print(out$sigma)

cls <- describe_classes(fit, dat)
print(cls$sizes); print(cls$confidence); print(cls$profile_treated)
cat("prior/posterior class agreement:", round(cls$prior_post_agree, 3), "\n")

# LOO — log_lik already in generated quantities. Standardized scale, so
# fine for comparison, not interpretable in CD4 units.
loo_k3 <- fit$loo()
print(loo_k3)

saveRDS(list(fit_summary = out, classes = cls, loo = loo_k3),
        "actg175_blcm_k3.rds")

## ---- 9. Membership model: sigma_gamma and gamma coefficients ----------
fit <- readRDS("actg175_blcm_k3_fit.rds") 
dat <- prep_actg175()                        

# 0. List the gamma-related variables in the Stan file
grep("gamma", fit$metadata()$stan_variables, value = TRUE)

# 1. sigma_gamma: median, 90% CrI and posterior SD
fit$summary("gamma_sd", "median", "sd", ~quantile(.x, c(.05, .95)))

# 2. Scaled membership slopes (gamma = gamma_sd * gamma_slopes_raw), draw by draw
sd_draws <- as.vector(as_draws_matrix(fit$draws("gamma_sd")))
raw      <- as_draws_matrix(fit$draws("gamma_slopes_raw"))
K1 <- length(dat$covars); C1 <- dat$stan_data$C - 1
contrast <- c("Null vs harm", "Benefit vs harm")   # harm is the reference class

rows <- list()
for (c in seq_len(C1)) for (k in seq_len(K1)) {
  d <- sd_draws * as.vector(raw[, sprintf("gamma_slopes_raw[%d,%d]", k, c)])
  rows[[length(rows) + 1]] <- data.frame(
    contrast = contrast[c], covariate = dat$covars[k],
    median = median(d), q5 = quantile(d, .05), q95 = quantile(d, .95))
}
gamma_tab <- do.call(rbind, rows); rownames(gamma_tab) <- NULL
gamma_tab$excludes_zero <- gamma_tab$q5 > 0 | gamma_tab$q95 < 0

print(gamma_tab, digits = 3)
sum(gamma_tab$excludes_zero)                            # expect 0
gamma_tab[order(-abs(gamma_tab$median)), ][1:6, ]       # largest medians for the text

# 3. LaTeX rows for the appendix table
lab <- c(age = "Age", wtkg = "Weight", hemo = "Haemophilia",
         homo = "Homosexual activity", drugs = "Intravenous drug use",
         karnof = "Karnofsky score", oprior = "Prior non-zidovudine ART",
         z30 = "Zidovudine in prior 30 days", race = "Race", gender = "Sex",
         str2 = "Antiretroviral history", symptom = "Symptomatic status",
         cd40 = "Baseline CD4", cd80 = "Baseline CD8")
fmt <- function(x) sprintf("$%.2f$", x)
ci  <- function(a, b) sprintf("$(%.2f,\\ %.2f)$", a, b)
nh <- gamma_tab[gamma_tab$contrast == "Null vs harm", ]
bh <- gamma_tab[gamma_tab$contrast == "Benefit vs harm", ]
cat(sprintf("%s & %s & %s & %s & %s \\\\",
            lab[nh$covariate], fmt(nh$median), ci(nh$q5, nh$q95),
            fmt(bh$median), ci(bh$q5, bh$q95)), sep = "\n")

## ------------------------------------------------------------------
## K=1 Comparison 
## ------------------------------------------------------------------

mod_k1 <- cmdstan_model("BLCM_k1_homogeneous.stan")
mod_k1_shared <- cmdstan_model("BLCM_k1_homogeneous_shared.stan")

d1 <- dat$stan_data[c("N","K","x","y","t")]   # drop C

fit_k1 <- mod_k1$sample(
  data = d1, chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  adapt_delta = 0.9, seed = 1, refresh = 500
)

fit_k1$summary(c("tau","sigma"))
loo_k1 <- fit_k1$loo()
print(loo_k1) #check Pareto k estimates <0.7

fit_k1_shared <- mod_k1_shared$sample(
  data = d1, chains = 4, parallel_chains = 4,
  iter_warmup = 2000, iter_sampling = 2000,
  adapt_delta = 0.9, seed = 1, refresh = 500
)

fit_k1_shared$summary(c("tau", "sigma"))
loo_k1_shared <- fit_k1_shared$loo()
print(loo_k1_shared)

loo_k3 <- fit$loo()   # or: loo_k3 <- readRDS("actg175_blcm_k3.rds")$loo
loo::loo_compare(list(K3 = loo_k3, K1_arm = loo_k1))

loo::loo_compare(list(K3 = loo_k3, K1_arm = loo_k1, K1_shared = loo_k1_shared))
arm_split(loo_k3, loo_k1_shared, dat$stan_data$t)

arm_split <- function(loo_a, loo_b, t) {
  d <- loo_a$pointwise[, "elpd_loo"] - loo_b$pointwise[, "elpd_loo"]
  c(total = sum(d), control = sum(d[t == 0]), treated = sum(d[t == 1]))
}

arm_split(loo_k3, loo_k1, dat$stan_data$t)          # vs arm-specific sigma
# once the shared-sigma model is fitted:
loo_k1_shared <- fit_k1_shared$loo()
loo::loo_compare(list(K3 = loo_k3, K1_arm = loo_k1, K1_shared = loo_k1_shared))
arm_split(loo_k3, loo_k1_shared, dat$stan_data$t)  

## ------------------------------------------------------------------
## BCF and BART Comparison
## ------------------------------------------------------------------

df <- dat$df
covars <- dat$covars
X  <- as.matrix(df[, covars])          # no intercept column
y  <- df$cd420
z  <- df$treat_binary

## ---- BCF -------------------------------------------------------------
## Randomized trial, so the true propensity is known and constant
pihat <- rep(mean(z), nrow(df))

set.seed(20261005)
t_bcf <- system.time(
  fit_bcf <- bcf(y = y, z = z, x_control = X, x_moderate = X,
                 pihat = pihat, nburn = 2000, nsim = 2000)
)

round(c(bcf = t_bcf[["elapsed"]], bart = t_bart[["elapsed"]]) / 60, 1)   # minutes

tau_bcf <- colMeans(fit_bcf$tau)       # posterior mean CATE per subject

## ---- BART ------------------------------------------------------------
Xtr <- cbind(trt = z, X)
set.seed(20261005)
t_bart <- system.time(
  bf <- wbart(x.train = Xtr, y.train = y, ndpost = 2000, nskip = 1000)
)

p1 <- predict(bf, cbind(trt = 1, X))
p0 <- predict(bf, cbind(trt = 0, X))
tau_bart <- colMeans(p1 - p0)

## ---- The number that matters: spread of individual CATEs -------------
data.frame(
  method = c("BCF", "BART"),
  mean   = c(mean(tau_bcf),  mean(tau_bart)),
  sd     = c(sd(tau_bcf),    sd(tau_bart)),
  q05    = c(quantile(tau_bcf, .05),  quantile(tau_bart, .05)),
  q95    = c(quantile(tau_bcf, .95),  quantile(tau_bart, .95))
)

tau_draws <- unclass(posterior::as_draws_matrix(fit$draws("tau")))

expected_tau <- function(fit, var, tau_draws, N, C) {
  P <- unclass(posterior::as_draws_matrix(fit$draws(var)))
  out <- numeric(N)
  for (c in seq_len(C)) {
    Pc <- P[, sprintf("%s[%d,%d]", var, seq_len(N), c), drop = FALSE]
    out <- out + colMeans(Pc * as.numeric(tau_draws[, c]))
  }
  out
}

N <- dat$stan_data$N; C <- dat$stan_data$C
tau_prosp <- expected_tau(fit, "class_prior", tau_draws, N, C)
tau_post  <- expected_tau(fit, "p_class",     tau_draws, N, C)

summ <- function(x) c(mean = mean(x), sd = sd(x),
                      q05 = unname(quantile(x, .05)), q95 = unname(quantile(x, .95)))
round(rbind(BCF = summ(tau_bcf), BART = summ(tau_bart),
            BLCM_prospective = summ(tau_prosp), BLCM_posterior = summ(tau_post)), 2)

round(c(bcf_bart   = cor(tau_bcf, tau_bart),
        prosp_bcf  = cor(tau_prosp, tau_bcf),
        prosp_bart = cor(tau_prosp, tau_bart),
        post_bcf   = cor(tau_post, tau_bcf)), 3)

## ------------------------------------------------------------------
## Seed Stability
## ------------------------------------------------------------------
conv_pars <- c("tau", "tau_std", "beta_std", "gamma_sd", "sigma_std",
               "gamma_intercept_raw", "gamma_slopes_raw")

main_class <- describe_classes(fit, dat)$df$class_hat
stab_list  <- list()

for (s in 1:5) {
  set.seed(s)
  f   <- fit_real(dat$stan_data, mod, fit_args, seed = 100 * s)
  cls <- describe_classes(f, dat)
  tt  <- f$summary("tau")
  stab_list[[s]] <- data.frame(
    seed = s,
    max_rhat = max(f$summary(conv_pars)$rhat, na.rm = TRUE),
    tau_harm = tt$median[1], tau_null = tt$median[2], tau_ben = tt$median[3],
    n_harm = sum(cls$df$class_hat == "Harm"),
    n_ben  = sum(cls$df$class_hat == "Benefit"),
    n_moved_vs_main = sum(cls$df$class_hat != main_class))
  saveRDS(stab_list, "seed_stability_partial.rds")
  message(sprintf("Seed %d done", s))
}

stab <- do.call(rbind, stab_list)
stab

