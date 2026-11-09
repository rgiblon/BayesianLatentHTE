# run_sim_comparators.R
# Called by SLURM as: Rscript run_sim_comparators.R <job_id>
# Fits BCF and BART on each simulated dataset and saves results.

# ── Setup ────────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
job_id  <- as.integer(args[1])

#PROJECT_DIR <- "~/Desktop/PhD/Dissertation/Paper 2 - LCM/lcm_sim"
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
CONDITIONS  <- file.path(PROJECT_DIR, "conditions.csv")
RESULTS_DIR <- file.path(PROJECT_DIR, "results_comparators")

# ── Libraries ────────────────────────────────────────────────────────────────
library(bcf)
library(dbarts)
library(MASS)
library(dplyr)
library(mvtnorm)
library(matrixStats)

# ── Source functions ──────────────────────────────────────────────────────────
source(file.path(PROJECT_DIR, "functions.R"))

# ── Extraction functions ──────────────────────────────────────────────────────
extract_metrics_bcf <- function(fit_bcf, truth, runtime_sec) {
  tau_draws <- fit_bcf$tau
  out <- list()
  out$runtime_sec <- runtime_sec
  for (cc in 1:truth$C) {
    idx             <- truth$class == cc
    tau_class_draws <- rowMeans(tau_draws[, idx, drop = FALSE])
    tau_med         <- median(tau_class_draws)
    tau_q5          <- quantile(tau_class_draws, 0.05)
    tau_q95         <- quantile(tau_class_draws, 0.95)
    tau_true        <- truth$tau[cc]
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
  out <- list()
  out$runtime_sec <- runtime_sec
  for (cc in 1:truth$C) {
    idx             <- truth$class == cc
    tau_class_draws <- rowMeans(tau_draws_bart[, idx, drop = FALSE])
    tau_med         <- median(tau_class_draws)
    tau_q5          <- quantile(tau_class_draws, 0.05)
    tau_q95         <- quantile(tau_class_draws, 0.95)
    tau_true        <- truth$tau[cc]
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

# ── Read condition ────────────────────────────────────────────────────────────
conditions <- read.csv(CONDITIONS, stringsAsFactors = FALSE)
row        <- conditions[conditions$job_id == job_id, ]

if (nrow(row) == 0) stop("job_id not found in conditions.csv: ", job_id)

pk  <- as.numeric(strsplit(row$pk_str,  "_")[[1]])
tau <- as.numeric(strsplit(row$tau_str, "_")[[1]])
N   <- row$N
rep <- row$rep

cat(sprintf("[job %d] cond=%d rep=%d N=%d pk=%s tau=%s\n",
            job_id, row$cond_id, rep, N, row$pk_str, row$tau_str))

seed <- row$cond_id * 100000 + rep

# ── Run ───────────────────────────────────────────────────────────────────────
result <- tryCatch({
  
  # Simulate data
  dat <- simulate_dataset(seed = seed, N = N, pk = pk, tau_true = tau)
  
  X     <- dat$stan_data$x[, -1]
  y     <- dat$stan_data$y
  t     <- dat$stan_data$t
  pihat <- rep(0.5, nrow(X))
  
  # ── BCF ──────────────────────────────────────────────────────────────────
  t_bcf <- proc.time()
  fit_bcf <- bcf(y          = y,
                 z          = t,
                 x_control  = X,
                 x_moderate = X,
                 pihat      = pihat,
                 nburn      = 1000,
                 nsim       = 1000,
                 verbose    = FALSE)
  runtime_bcf <- (proc.time() - t_bcf)["elapsed"]
  out_bcf <- extract_metrics_bcf(fit_bcf, dat$truth, runtime_bcf)
  
  # ── BART ─────────────────────────────────────────────────────────────────
  t_bart <- proc.time()
  fit_bart <- bart(x.train   = cbind(X, t),
                   y.train   = y,
                   ntree     = 200,
                   ndpost    = 1000,
                   nskip     = 1000,
                   keeptrees = TRUE,
                   verbose   = FALSE)
  runtime_bart <- (proc.time() - t_bart)["elapsed"]
  
  X_treated        <- cbind(X, 1)
  X_control        <- cbind(X, 0)
  tau_draws_bart   <- predict(fit_bart, X_treated) - predict(fit_bart, X_control)
  out_bart         <- extract_metrics_bart(tau_draws_bart, dat$truth, runtime_bart)
  
  # ── Combine ───────────────────────────────────────────────────────────────
  names(out_bcf)  <- paste0(names(out_bcf),  "_bcf")
  names(out_bart) <- paste0(names(out_bart), "_bart")
  
  result <- cbind(
    data.frame(job_id  = job_id,
               cond_id = row$cond_id,
               rep     = rep,
               seed    = seed,
               N       = N,
               pk_str  = row$pk_str,
               tau_str = row$tau_str,
               factor  = row$factor),
    out_bcf,
    out_bart
  )
  result
  
}, error = function(e) {
  data.frame(
    job_id  = job_id,
    cond_id = row$cond_id,
    rep     = rep,
    seed    = seed,
    N       = N,
    pk_str  = row$pk_str,
    tau_str = row$tau_str,
    factor  = row$factor,
    error   = conditionMessage(e),
    stringsAsFactors = FALSE
  )
})

# ── Save ──────────────────────────────────────────────────────────────────────
out_file <- file.path(RESULTS_DIR, sprintf("sim_%06d.rds", job_id))
saveRDS(result, out_file)
cat(sprintf("[job %d] saved -> %s\n", job_id, out_file))
