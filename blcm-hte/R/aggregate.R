# aggregate.R
# Run after all jobs complete (or periodically to check progress).
# Collects all sim_*.rds files into a single dataframe.

library(dplyr)

RESULTS_DIR  <- "results"
CONDITIONS   <- read.csv("conditions.csv", stringsAsFactors = FALSE)
TOTAL_JOBS   <- nrow(CONDITIONS)

files <- list.files(RESULTS_DIR, pattern = "^sim_.*\\.rds$", full.names = TRUE)
cat(sprintf("Found %d / %d result files (%.1f%%)\n",
            length(files), TOTAL_JOBS, 100 * length(files) / TOTAL_JOBS))

# ── Load all results ──────────────────────────────────────────────────────────
results_list <- lapply(files, function(f) {
  tryCatch(readRDS(f), error = function(e) {
    warning("Could not read: ", f)
    NULL
  })
})
results_list <- Filter(Negate(is.null), results_list)

# ── Separate successful runs from errors ──────────────────────────────────────
is_error <- sapply(results_list, function(x) "error" %in% names(x))

results_ok  <- bind_rows(results_list[!is_error])
results_err <- bind_rows(results_list[is_error])

cat(sprintf("Successful: %d | Errors: %d\n",
            nrow(results_ok), nrow(results_err)))

# ── Find missing job IDs ───────────────────────────────────────────────────────
completed_ids <- as.integer(gsub(".*sim_(\\d+)\\.rds$", "\\1",
                                 basename(files)))
missing_ids   <- setdiff(seq_len(TOTAL_JOBS), completed_ids)

if (length(missing_ids) > 0) {
  cat(sprintf("\nMissing %d jobs. Resubmit with:\n", length(missing_ids)))
  cat(sprintf("  sbatch --array=%s submit.sh\n\n",
              paste(missing_ids, collapse = ",")))
  write(missing_ids, "missing_jobs.txt", ncolumns = 1)
} else {
  cat("\nAll jobs completed.\n")
}

# ── Save combined results ─────────────────────────────────────────────────────
saveRDS(results_ok,  "results_combined.rds")
if (nrow(results_err) > 0) saveRDS(results_err, "results_errors.rds")

cat("Saved: results_combined.rds\n")

# ── Quick sanity check ────────────────────────────────────────────────────────
if (nrow(results_ok) > 0) {
  cat("\nReps per condition:\n")
  print(table(results_ok$cond_id))

  cat("\nConvergence rate by condition:\n")
  results_ok %>%
    group_by(cond_id, N, pk_str, tau_str) %>%
    summarise(
      n_reps    = n(),
      conv_rate = mean(max_rhat < 1.01 & n_divergent == 0, na.rm = TRUE),
      .groups   = "drop"
    ) %>%
    print()
}
