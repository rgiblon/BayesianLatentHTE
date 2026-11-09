# run_sim_binary.R
# Called by SLURM as: Rscript run_sim_binary.R <job_id>
# Reads one row from conditions_binary.csv, runs one rep, saves result to results_binary/.
# Mirrors run_sim.R exactly — only differences are:
#   - Stan model: BLCA_binary_CmdStanR.stan
#   - Conditions file: conditions_binary.csv
#   - Results dir: results_binary/
#   - DGP: simulate_dataset_binary() instead of simulate_dataset()
#   - Retry: fit_with_retry_bin() instead of fit_with_retry()
#   - tau values stored as log-OR (e.g. log(0.25)_0_log(4))
#
# ── Setup ────────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
job_id  <- as.integer(args[1])

# Paths — adjust if your project folder name differs
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")
cat("PROJECT_DIR:", PROJECT_DIR, "\n")
cat("exe exists:", file.exists(file.path(PROJECT_DIR, "BLCA_binary_CmdStanR")), "\n")
STAN_MODEL  <- file.path(PROJECT_DIR, "BLCA_binary_CmdStanR.stan")
CONDITIONS  <- file.path(PROJECT_DIR, "conditions_binary.csv")
RESULTS_DIR <- file.path(PROJECT_DIR, "results_binary")

# ── Libraries ────────────────────────────────────────────────────────────────
library(cmdstanr)
library(posterior)
library(bayesplot)
library(matrixStats)
library(MASS)
library(dplyr)
library(mvtnorm)

# ── Source functions ──────────────────────────────────────────────────────────
# functions_binary.R should contain:
#   simulate_dataset_binary(), extract_metrics(), get_pbar_from_fit(),
#   get_prior_from_fit(), make_init_fun_binary(), fit_with_retry_bin(),
#   summarise_tau(), summarise_misclass(), summarise_confidence(),
#   summarise_all(), fit_args
source(file.path(PROJECT_DIR, "functions_binary.R"))

# ── Load Stan model ───────────────────────────────────────────────────────────
# cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")  # uncomment if needed
mod_binary <- cmdstan_model(
  STAN_MODEL,
  exe_file = file.path(PROJECT_DIR, "BLCA_binary_CmdStanR")
)

# ── Read this job's condition ─────────────────────────────────────────────────
conditions <- read.csv(CONDITIONS, stringsAsFactors = FALSE)
row        <- conditions[conditions$job_id == job_id, ]

if (nrow(row) == 0) stop("job_id not found in conditions_binary.csv: ", job_id)

# Parse pk and tau back from strings
# tau values are stored as log-OR: e.g. "-1.386_0_1.386" for OR (0.25, 1, 4)
pk  <- as.numeric(strsplit(row$pk_str,  "_")[[1]])
tau <- as.numeric(strsplit(row$tau_str, "_")[[1]])
N   <- row$N
rep <- row$rep

cat(sprintf("[job %d] cond=%d rep=%d N=%d pk=%s tau=%s\n",
            job_id, row$cond_id, rep, N, row$pk_str, row$tau_str))

# ── Unique seed ───────────────────────────────────────────────────────────────
seed <- row$cond_id * 100000 + rep

# ── Run ───────────────────────────────────────────────────────────────────────
result <- tryCatch({
  dat <- simulate_dataset_binary(
    seed     = seed,
    N        = N,
    pk       = pk,
    tau_true = tau
    # all other DGM params stay at their defaults (baseline values)
  )

  fit <- fit_with_retry_bin(
    stan_data = dat$stan_data,
    truth     = dat$truth,
    model     = mod_binary,
    fit_args  = fit_args
  )

  out           <- extract_metrics(fit, dat$truth, compute_class_acc = TRUE)
  out$job_id    <- job_id
  out$cond_id   <- row$cond_id
  out$rep       <- rep
  out$seed      <- seed
  out$N         <- N
  out$pk_str    <- row$pk_str
  out$tau_str   <- row$tau_str
  out$factor    <- row$factor
  out

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
out_file <- file.path(RESULTS_DIR, sprintf("sim_binary_%06d.rds", job_id))
saveRDS(result, out_file)
cat(sprintf("[job %d] saved -> %s\n", job_id, out_file))
