#Install rstan and dependencies (if you haven't already)
#install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
options(mc.cores = parallel::detectCores())  # Allows Stan to use multiple cores

#Set working directory to dissertation
setwd("~/Desktop/PhD/Dissertation")

# Load the rstan library
library(rstan)
library(MASS)
library(bayesplot)
library(ggplot2)

# Set random seed for reproducibility
set.seed(6417)

###################### DATA SIMULATION #############################
# Parameters
N <- 1000   # Number of individuals
M <- 10     # Number of sites
C <- 3      # Number of latent classes
K <- 5      # Number of covariates

# Step 1: Assign individuals to classes
pk <- c(0.2,0.3,0.5)  # Class probabilities
C_n <- sample(1:C, N, replace = TRUE, prob = pk)

# Define mean vectors for each class #Make means more different? and reduce variance
#mean_X1 = I(c=1)*0 + I(C=2)*5 + I(C=2)*10
#mean_X1 = I(c=1)*0 + I(C=2)*5 + I(C=3)*10, x1 ~N(mean_X1, 1)
mu <- list(
  c(-4,-5),  # Mean vector for class 1
  c(0, 0),  # Mean vector for class 2
  c(5, 4))   # Mean vector for class 3

# Generate Level 1 covariates (at the individual level)
x1_level1 <- do.call(rbind, lapply(C_n, function(k) mvrnorm(1, mu[[k]], Sigma=matrix(c(1,0,0,1), nrow = 2))))

# Step 2: Generate Level 2 covariates (at the site level)
x2_level2 <- matrix(NA, M, 2)  # Assuming two level-2 covariates
x2_level2[1:M, ] <- mvrnorm(M, mu = c(0, 0), Sigma = matrix(c(1,-1,-1,1), nrow = 2))

# Step 3: Assign individuals to sites
# Create an empty vector to store site assignments
cluster_assignments <- numeric(N)
# Define equal probabilities for each cluster
cluster_probs <- rep(1/M, M)  # 1/M probability for each cluster
X_ij <- matrix(NA, N, ncol = 5)  # 5 covariates (1, x1, 0, x2)

  for (n in 1:N) {
    # Assign each participant to a cluster randomly with equal probability
    cluster_assignments[n] <- sample(1:M, 1, prob = cluster_probs)
    X_ij[n, ] <- c(1, x1_level1[n, ], x2_level2[cluster_assignments[n], ])
  }

#Step 4: Generate treatment indicators - random treatment assignment
T_ij <- rbinom(1000,1,0.5)

# Step 5: Generate outcomes via treatment & covars
Y_ij <- numeric(N)  # Initialize the outcome vector
tau <- c(-5,0,5)  # True treatment effects (t1, t2, t3)
sigma <- c(4,3,5)  # Standard deviations (s1, s2, s3)

for (n in 1:N) {
  j <- cluster_assignments[n]  # Get the cluster assignment for the observation
  c <- C_n[n]  # Get the categorical variable for the observation
  Xij <- X_ij[n, ]  # Get the covariates for the observation
  Tij <- T_ij[n]  # Get the treatment assignment for the observation
  
  # Calculate the mean based on the treatment effect, and covariates
  if (c == 1) {
    mean_value <- tau[1] * Tij + sum(c(5,7,-5,1,-2) * Xij)
    Y_ij[n] <- rnorm(1, mean = mean_value, sd = sigma[1])  # Generate outcome
  } else if (c == 2) {
    mean_value <- tau[2] * Tij + sum(c(8,-5,0,-1,1) * Xij)
    Y_ij[n] <- rnorm(1, mean = mean_value, sd = sigma[2])  # Generate outcome
  } else if (c == 3) {
    mean_value <- tau[3] * Tij + sum(c(-1,2,3,1,-1) * Xij)
    Y_ij[n] <- rnorm(1, mean = mean_value, sd = sigma[3])  # Generate outcome
  }
}

# Combine all data into a data frame
simulated_data <- data.frame(
  individual_id = 1:N,
  #cluster = cluster_assignments,
  class = C_n,
  treatment = T_ij,
  outcome = Y_ij,
  x1 = X_ij[, 2],
  x2 = X_ij[, 3],
  x3 = X_ij[, 4],
  x4 = X_ij[, 5]
)

# Print the first few rows of the simulated dataset
head(simulated_data)

####################### END OF DATA SIMULATION #######################

# Prepare data list for Stan
stan_data <- list(
  N = N,
  M = M,
  K = K,
  C = C,
  x = X_ij,
  y = Y_ij,
  t = T_ij,
  s = cluster_assignments
)

#Compile Stan model - modified 
model.blca <-stan_model("BLCA_fixedtrt.stan")

#Fit the model with the data - no random effects or clusters
fit.blca <- sampling(model.blca, data=stan_data, chains=3, iter=10000, warmup=2000, 
                     control = list(adapt_delta = 0.99, max_treedepth = 15))
#Check number of parameters

#Print results - estimates for treatment effects
tau_samples <- extract(fit.blca, pars = "tau")$tau
colnames(tau_samples) <- c("tau_1", "tau_2", "tau_3")
tau_summary <- summary(fit.blca, pars = "tau")$summary
print(tau_summary)

#Plot posterior distribution of Tau
mcmc_areas(as.matrix(tau_samples), 
           pars = colnames(tau_samples)) + 
  labs(title = "Posterior Distribution of Tau",
       x = "Tau",
       y = "Density") +
  theme_minimal()

#Fit checks
fit_summary <- summary(fit.blca)$summary
# Identify parameters with low ESS (e.g., less than 400)
low_ess_params <- rownames(fit_summary)[fit_summary[, "n_eff"] < 400]
low_ess_params  # List of parameter names with low ESS

#Label switching - need to check proportion of individuals correctly assigned to each group 
#Look at cluster label for each output




