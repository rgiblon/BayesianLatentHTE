# make_conditions.R
# Run once locally to generate the condition grid.
# Based on ADEMP preregistration: partially factorial design, ~12-16 conditions.
# Each row = one simulation condition. Reps are handled by SLURM array offset.

# ── Baseline values (held fixed unless varied) ──────────────────────────────
BASELINE <- list(
  N      = 800,
  pk     = list(c(0.30, 0.30, 0.40)),   # stored as list so expand.grid works
  tau    = list(c(-5, 0, 5))
)

# ── Factor levels (varied one/two at a time around baseline) ────────────────

# Factor 1: Sample size (vary N, hold pk and tau at baseline)
cond_N <- expand.grid(
  N      = c(200, 400, 800, 1000),
  pk     = list(c(0.30, 0.30, 0.40)),
  tau    = list(c(-5, 0, 5)),
  factor = "sample_size"
)

# Factor 2: Subgroup prevalence (vary pk, hold N and tau at baseline)
cond_pk <- expand.grid(
  N      = 800,
  pk     = list(c(0.30, 0.30, 0.40),   # baseline (already in cond_N, kept for clarity)
                c(0.10, 0.20, 0.70)),   # imbalanced: small harm subgroup
  tau    = list(c(-5, 0, 5)),
  factor = "prevalence"
)

# Factor 3: Treatment-effect separation (vary tau, hold N and pk at baseline)
cond_tau <- expand.grid(
  N      = 800,
  pk     = list(c(0.30, 0.30, 0.40)),
  tau    = list(c(-7,  0,  7),    # very strong separation
                c(-6,  0,  6),    # strong separation
                c(-5,  0,  5),    # baseline (already covered)
                c(-3,  0,  3),    # moderate separation
                c(-2,  0,  2),    # weak separation
                c( 0,  2,  5)),   # asymmetric/exploratory
  factor = "tau_separation"
)

# ── Combine and deduplicate ──────────────────────────────────────────────────
conditions_raw <- rbind(cond_N, cond_pk, cond_tau)

# Convert list columns to strings for deduplication
conditions_raw$pk_str  <- sapply(conditions_raw$pk,  paste, collapse = "_")
conditions_raw$tau_str <- sapply(conditions_raw$tau, paste, collapse = "_")

# Deduplicate on the actual parameter values (baseline appears in multiple factors)
conditions_unique <- conditions_raw[
  !duplicated(conditions_raw[, c("N", "pk_str", "tau_str")]), 
]

# Assign condition ID and keep only what's needed
conditions_unique$cond_id <- seq_len(nrow(conditions_unique))

# ── Add reps ─────────────────────────────────────────────────────────────────
REPS_PER_COND <- 1000   # adjust after pilot; preregistration minimum

conditions_full <- conditions_unique[
  rep(seq_len(nrow(conditions_unique)), each = REPS_PER_COND), 
]
conditions_full$rep <- rep(seq_len(REPS_PER_COND), times = nrow(conditions_unique))
conditions_full$job_id <- seq_len(nrow(conditions_full))

# ── Save ─────────────────────────────────────────────────────────────────────
# Keep only columns needed by run_sim.R
out <- data.frame(
  job_id  = conditions_full$job_id,
  cond_id = conditions_full$cond_id,
  rep     = conditions_full$rep,
  N       = conditions_full$N,
  pk_str  = conditions_full$pk_str,
  tau_str = conditions_full$tau_str,
  factor  = conditions_full$factor,
  stringsAsFactors = FALSE
)

write.csv(out, "conditions.csv", row.names = FALSE)

cat("Conditions summary:\n")
cat("  Unique conditions:", nrow(conditions_unique), "\n")
cat("  Reps per condition:", REPS_PER_COND, "\n")
cat("  Total jobs:", nrow(out), "\n\n")
print(conditions_unique[, c("cond_id", "N", "pk_str", "tau_str", "factor")])
