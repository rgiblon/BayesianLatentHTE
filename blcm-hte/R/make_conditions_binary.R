# make_conditions_binary.R
# Run once locally to generate the condition grid for the binary outcome model.
# Mirrors make_conditions.R exactly — tau values are log-OR (e.g. log(0.25), 0, log(3)).
# Partially factorial design varying N, pk, and tau separation around baseline.

# ── Baseline values ──────────────────────────────────────────────────────────
# Baseline tau: OR = (0.25, 1, 3) on log scale
# Corresponds to 60-day survival probabilities of ~13%, 37.5%, 64% on treated arm
BASELINE <- list(
  N   = 800,
  pk  = list(c(0.30, 0.30, 0.40)),
  tau = list(c(log(0.25), 0, log(3)))
)

# ── Factor levels ─────────────────────────────────────────────────────────────

# Factor 1: Sample size (vary N, hold pk and tau at baseline)
cond_N <- expand.grid(
  N      = c(200, 400, 800, 1000),
  pk     = list(c(0.30, 0.30, 0.40)),
  tau    = list(c(log(0.25), 0, log(3))),
  factor = "sample_size"
)

# Factor 2: Subgroup prevalence (vary pk, hold N and tau at baseline)
cond_pk <- expand.grid(
  N      = 800,
  pk     = list(c(0.30, 0.30, 0.40),
                c(0.10, 0.20, 0.70)),
  tau    = list(c(log(0.25), 0, log(3))),
  factor = "prevalence"
)

# Factor 3: Treatment-effect separation (vary tau, hold N and pk at baseline)
# ORs chosen to mirror continuous model separation conditions:
#   Strong:        OR (0.1,  1, 4.0)  -- very strong, mirrors tau=(-7,0,7)
#   Strong-mod:    OR (0.20, 1, 3.5) -- mirrors tau=(-6,0,6)
#   Baseline:      OR (0.25, 1, 3.0)   -- baseline
#   Moderate:      OR (0.4,  1, 2.5) -- mirrors tau=(-3,0,3)
#   Weak:          OR (0.6,  1, 1.7) -- mirrors tau=(-2,0,2)
#   Landmark (Ewan):OR(0.85, 1, 1.3) -- exploratory, small real-world effect
cond_tau <- expand.grid(
  N      = 800,
  pk     = list(c(0.30, 0.30, 0.40)),
  tau    = list(
    c(log(0.15),  0, log(4.0)),   # very strong separation
    c(log(0.20),  0, log(3.5)),    # strong separation
    c(log(0.25),  0, log(3.0)),    # baseline
    c(log(0.40),  0, log(2.5)),    # moderate separation
    c(log(0.60),  0, log(1.70)),   # weak separation
    c(log(0.85),  0, log(1.30))    # exploratory: Ewan landmark OR setting
  ),
  factor = "tau_separation"
)

# ── Combine and deduplicate ───────────────────────────────────────────────────
conditions_raw <- rbind(cond_N, cond_pk, cond_tau)

# Round tau values before string conversion to avoid floating point mismatches
conditions_raw$pk_str  <- sapply(conditions_raw$pk,
                                  paste, collapse = "_")
conditions_raw$tau_str <- sapply(conditions_raw$tau, function(x)
  paste(round(x, 6), collapse = "_"))

conditions_unique <- conditions_raw[
  !duplicated(conditions_raw[, c("N", "pk_str", "tau_str")]),
]

conditions_unique$cond_id <- seq_len(nrow(conditions_unique))

# ── Add reps ──────────────────────────────────────────────────────────────────
REPS_PER_COND <- 1000

conditions_full <- conditions_unique[
  rep(seq_len(nrow(conditions_unique)), each = REPS_PER_COND),
]
conditions_full$rep    <- rep(seq_len(REPS_PER_COND), times = nrow(conditions_unique))
conditions_full$job_id <- seq_len(nrow(conditions_full))

# ── Save ──────────────────────────────────────────────────────────────────────
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

write.csv(out, "conditions_binary.csv", row.names = FALSE)

cat("Binary conditions summary:\n")
cat("  Unique conditions:", nrow(conditions_unique), "\n")
cat("  Reps per condition:", REPS_PER_COND, "\n")
cat("  Total jobs:", nrow(out), "\n\n")
print(conditions_unique[, c("cond_id", "N", "pk_str", "tau_str", "factor")])
