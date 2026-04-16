# run_sim.R
# Called by SLURM as: Rscript run_sim.R <job_id>
# Reads one row from conditions.csv, runs one rep, saves result to results/.
#
# ── Setup ────────────────────────────────────────────────────────────────────
args    <- commandArgs(trailingOnly = TRUE)
job_id  <- as.integer(args[1])

# Paths — adjust if your project folder name differs
PROJECT_DIR <- Sys.getenv("PROJECT_DIR", unset = ".")   # set in submit.sh
STAN_MODEL  <- file.path(PROJECT_DIR, "BLCA_fixedtrt_CmdStanR_v2.stan")
CONDITIONS  <- file.path(PROJECT_DIR, "conditions.csv")
RESULTS_DIR <- file.path(PROJECT_DIR, "results")

# ── Libraries ────────────────────────────────────────────────────────────────
library(cmdstanr)
library(posterior)
library(bayesplot)
library(matrixStats)
library(MASS)
library(dplyr)
library(mvtnorm)

# ── Source your existing functions (everything above run_once in RepsCode.R) ─
# We keep all your existing functions in a separate sourced file so nothing
# needs to change in them. Just split RepsCode.R at the "trial replications"
# section: keep simulate_dataset(), extract_metrics(), helpers, fit_args, etc.
# in "functions.R" and source it here.
source(file.path(PROJECT_DIR, "functions.R"))

# ── Load Stan model ──────────────────────────────────────────────────────────
# CmdStan path on Alliance — module load r sets this automatically if configured;
# otherwise set explicitly. Check with: which cmdstan
# cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.37.0")  # uncomment if needed

mod <- cmdstan_model(STAN_MODEL)

# ── Read this job's condition ────────────────────────────────────────────────
conditions <- read.csv(CONDITIONS, stringsAsFactors = FALSE)
row        <- conditions[conditions$job_id == job_id, ]

if (nrow(row) == 0) stop("job_id not found in conditions.csv: ", job_id)

# Parse pk and tau back from strings (e.g. "0.3_0.3_0.4" -> c(0.3, 0.3, 0.4))
pk  <- as.numeric(strsplit(row$pk_str,  "_")[[1]])
tau <- as.numeric(strsplit(row$tau_str, "_")[[1]])
N   <- row$N
rep <- row$rep

cat(sprintf("[job %d] cond=%d rep=%d N=%d pk=%s tau=%s\n",
            job_id, row$cond_id, rep, N, row$pk_str, row$tau_str))

# ── Unique seed: combines condition and rep so seeds never collide ────────────
seed <- row$cond_id * 100000 + rep

# ── Run ──────────────────────────────────────────────────────────────────────
result <- tryCatch({
  dat <- simulate_dataset(
    seed     = seed,
    N        = N,
    pk       = pk,
    tau_true = tau
    # all other DGM params stay at their defaults (baseline values)
  )

  fit <- fit_with_retry(
    stan_data = dat$stan_data,
    truth     = dat$truth,
    model     = mod,
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

# ── Save ─────────────────────────────────────────────────────────────────────
out_file <- file.path(RESULTS_DIR, sprintf("sim_%06d.rds", job_id))
saveRDS(result, out_file)
cat(sprintf("[job %d] saved -> %s\n", job_id, out_file))
