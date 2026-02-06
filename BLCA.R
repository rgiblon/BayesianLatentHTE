#Install rstan and dependencies (if you haven't already)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores
cmdstanr::set_cmdstan_path("/Users/rachelpim/.cmdstan/cmdstan-2.37.0")

#Set working directory to dissertation
setwd("~/Desktop/PhD/Dissertation/Paper 2 - LCM/Simulations_Code")

# Load the rstan library
library(rstan)
library(cmdstanr)
library(MASS)
library(bayesplot)
library(ggplot2)
library(label.switching)
library(posterior)
library(clue)

#Set random seed for reproducibility
set.seed(6417)

###################### DATA SIMULATION #############################
#Parameters
N <- 1000   # Number of individuals
M <- 10     # Number of sites
C <- 3      # Number of latent classes
K <- 6      # Number of covariates

# Step 1: Assign individuals to classes
pk <- c(0.3,0.3,0.4)  # Class probabilities (larger groups = easier for model)
C_n <- sample(1:C, N, replace = TRUE, prob = pk)

# Define mean vectors for each class 
mu <- list(
  c(-6,-6), # Mean vector for class 1
  c(0, 6),  # Mean vector for class 2
  c(6, 0))  # Mean vector for class 3

#Can increase or decrease within-class covariate variance (cluster tighter)
Sigma_within <- diag(2) * 0.25

# Generate Level 1 covariates (at the individual level)
x1_level1 <- do.call(rbind, lapply(C_n, function(k) mvrnorm(1, mu[[k]], Sigma = Sigma_within)))
                                   
# Step 2: Generate Level 2 covariates (at the site level)
x2_level2 <- matrix(NA, M, 2)  # Assuming two level-2 covariates
x2_level2[1:M, ] <- mvrnorm(M, mu = c(0, 0), Sigma = matrix(c(1,-1,-1,1), nrow = 2))

#Biomarker predicted by class
p_bio <- c(0.05, 0.5, 0.95)
B <- sapply(C_n, function(k) rbinom(1, 1, p_bio[k]))

# Step 3: Assign individuals to sites
cluster_assignments <- sample(1:M, N, replace = TRUE)  # simpler assign
X_ij_raw <- cbind(1, x1_level1, x2_level2[cluster_assignments, ], B)  #now K=6

# Standardize the non-intercept covariates (center & scale)
#X_nonint <- scale(X_ij_raw[, -1])    # N x 4
#X_ij <- cbind(1, X_nonint)           # N x 5 (intercept + standardized covs)
X_ij <- X_ij_raw                      # Try unstandardized for better seperation

#Step 4: Generate treatment indicators - random treatment assignment
T_ij <- rbinom(N,1,0.5)

# True taus unchanged
tau_true <- c(-5, 0, 5)

#Beta coefficients
beta1 <- c(5, 7, -5, 1, 3)   # include biomarker coefficient last (strong)
beta2 <- c(8, -5, 0, -1, 0.5)
beta3 <- c(-1, 2, 3, 1, -2)

# Option: reduce residual sds to improve signal/noise
sigma <- c(1.5, 1.5, 1.5)  # smaller than original c(4,3,5)

# Step 5: Generate outcomes via treatment & covars
Y_ij <- numeric(N)
for (n in 1:N) {
  c_idx <- C_n[n]
  Xn <- X_ij[n, ]
  if (c_idx == 1) {
    mean_val <- tau_true[1] * T_ij[n] + sum(Xn * c(0, beta1))  # intercept handled
    Y_ij[n] <- rnorm(1, mean = mean_val, sd = sigma[1])
  } else if (c_idx == 2) {
    mean_val <- tau_true[2] * T_ij[n] + sum(Xn * c(0, beta2))
    Y_ij[n] <- rnorm(1, mean = mean_val, sd = sigma[2])
  } else {
    mean_val <- tau_true[3] * T_ij[n] + sum(Xn * c(0, beta3))
    Y_ij[n] <- rnorm(1, mean = mean_val, sd = sigma[3])
  }
}

simulated_data <- data.frame(
  id = 1:N,
  class = C_n,
  treatment = T_ij,
  outcome = Y_ij,
  X_ij[, -1]   # standardized covariates
)

# Print the first few rows of the simulated dataset
head(simulated_data)

#################### DATA CHECKS - CLASS SEPERATION ###################

check_separation <- function(X, classes, reg_eps = 1e-6, K = NULL) {
  # X: N x K matrix of covariates (incl. intercept in col 1, will be dropped)
  # classes: vector of true class labels (length N)
  # reg_eps: ridge term for regularization
  # K: number of clusters (default = number of unique classes)
  
  if (is.null(K)) K <- length(unique(classes))
  
  # Drop intercept
  X_cov <- X[, -1, drop = FALSE]
  
  # Regularized pooled covariance
  S_pool <- cov(X_cov)
  S_reg  <- S_pool + diag(reg_eps, ncol(S_pool))
  
  # Class means
  mu_list <- lapply(split(data.frame(X_cov), classes), colMeans)
  
  # Pairwise Mahalanobis distances
  dists <- matrix(NA, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      if (i < j) {
        dists[i, j] <- mahalanobis(mu_list[[i]], mu_list[[j]], cov = S_reg)
        dists[j, i] <- dists[i, j]
      }
    }
  }
  rownames(dists) <- colnames(dists) <- paste0("Class", 1:K)
  
  # Quick clustering check (k-means)
  set.seed(123)
  km <- kmeans(X_cov, centers = K, nstart = 20)
  tab <- table(km$cluster, classes)
  
  return(list(
    maha_distances = dists,
    confusion_matrix = tab,
    kmeans_accuracy = sum(diag(tab)) / length(classes)
  ))
}

results <- check_separation(X_ij, C_n)

print(results$maha_distances)   # pairwise class separation
print(results$confusion_matrix) # how k-means clusters vs truth
print(results$kmeans_accuracy)  # rough "recoverability"

#Relabel - kmeans also label switches
cluster_accuracy <- function(tab) {
  assignment <- solve_LSAP(tab, maximum = TRUE)
  sum(tab[cbind(seq_along(assignment), assignment)]) / sum(tab)}

cluster_accuracy(results$confusion_matrix)

####################### RUN MODEL - HMC/NUTS #######################

#Prepare data list for Stan
stan_data <- list(
  N = N,
  M = M,
  K = ncol(X_ij),
  C = C,
  x = X_ij,
  y = Y_ij,
  t = T_ij,
  s = cluster_assignments
)

#Compile Stan model
mod <-cmdstan_model("BLCA_fixedtrt_CmdStanR.stan") 
mod2 <-cmdstan_model("BLCA_fixedtrt_CmdStanR_2.stan") 

#Pilot first for quick fixes
fit_pilot_cmd <- mod$sample(data=stan_data, chains=3, parallel_chains=3,
                            iter_warmup=3000, iter_sampling=3000,
                            adapt_delta=0.995, max_treedepth=15)

#Pilot first for quick fixes
fit_pilot_cmd2 <- mod2$sample(data=stan_data, chains=3, parallel_chains=3,
                            iter_warmup=1000, iter_sampling=1000,
                            adapt_delta=0.995, max_treedepth=15)


#1. Get model stable on one dataset  ← you are here
#2. Build simulation wrapper
#3. Scale to many datasets (500?)

################### TAU ESTIMATES - NO ECR ##################################

#CMDSTAN
tau_samples <- as_draws_matrix(fit_pilot_cmd$draws("tau"))
pi_samples  <- as_draws_matrix(fit_pilot_cmd$draws("pi"))
colnames(tau_samples) <- c("tau_1", "tau_2", "tau_3")
tau_summary <- fit_pilot_cmd$summary("tau")
print(tau_summary)


mcmc_trace(tau_samples) #chain diagnostics
mcmc_areas(pi_samples)
mcmc_trace(pi_samples)

#Plot posterior distribution of Tau
mcmc_areas(as.matrix(tau_samples), 
           pars = colnames(tau_samples)) + 
  labs(title = "Posterior Distribution of Tau",
       x = "Tau",
       y = "Density") +
  theme_minimal()

#Diagnostics: distribution of R-hat across parameters, effective sample sizes
mcmc_rhat(rhat(fit_pilot)) 
mcmc_neff(neff_ratio(fit_pilot))
mcmc_trace(as.array(fit_pilot), pars = c("tau[1]", "tau[2]", "tau[3]"))

# Posterior predictive checks
y_rep <- extract(fit_pilot, pars = "y_rep")$y_rep        # draws x N
bayesplot::ppc_dens_overlay(y = Y_ij, yrep = y_rep[1:100, ])   # quick PPC

#Fit checks
fit_summary <- summary(fit_pilot)$summary
# Identify parameters with low ESS (e.g., less than 400)
low_ess_params <- rownames(fit_summary)[fit_summary[, "n_eff"] < 400]
low_ess_params  # List of parameter names with low ESS





