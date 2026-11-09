# aggregate_comp.R
# Collects all comparator (BCF + BART) simulation results into one dataframe.

library(dplyr)

RESULTS_DIR  <- "results_comparators"
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

# ── Separate successes from errors ────────────────────────────────────────────
is_error <- sapply(results_list, function(x) {
  "error" %in% names(x) && is.character(x$error)
})

results_ok  <- bind_rows(results_list[!is_error])
results_err <- bind_rows(results_list[is_error])

cat(sprintf("Successful: %d | Errors: %d\n", nrow(results_ok), nrow(results_err)))

# ── Find missing jobs ─────────────────────────────────────────────────────────
completed_ids <- as.integer(gsub(".*sim_(\\d+)\\.rds$", "\\1", basename(files)))
missing_ids   <- setdiff(seq_len(TOTAL_JOBS), completed_ids)

if (length(missing_ids) > 0) {
  cat(sprintf("\nMissing %d jobs.\n", length(missing_ids)))
  write(missing_ids, "missing_jobs_comp.txt", ncolumns = 1)
} else {
  cat("\nAll jobs completed.\n")
}

# ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(results_ok, "results_combined_comp.rds")
if (nrow(results_err) > 0) saveRDS(results_err, "results_errors_comp.rds")
cat("Saved: results_combined_comp.rds\n")

# ── Quick summary ─────────────────────────────────────────────────────────────
if (nrow(results_ok) > 0) {
  cat("\nReps per condition:\n")
  print(table(results_ok$cond_id))
  
  cat("\nMean BCF runtime by condition (seconds):\n")
  results_ok %>%
    group_by(cond_id, N, tau_str) %>%
    summarise(mean_runtime_bcf  = mean(runtime_sec_bcf,  na.rm = TRUE),
              mean_runtime_bart = mean(runtime_sec_bart, na.rm = TRUE),
              .groups = "drop") %>%
    print()
}