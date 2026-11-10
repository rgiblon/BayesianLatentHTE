# tau_tilde.R
# Within-class targets of the CATE methods (BCF, BART) for the continuous
# simulation, tau_tilde_c = E{ tau(X) | C = c } = sum_k W[c, k] * tau_k,
# where W[c, k] = E{ pi_k(X) | C = c } and pi_k(x) = P(C = k | X = x).
#
# Covariates and class membership are generated exactly as in functions.R
# (simulate_dataset): x1 from the class mixture with weights pk, the biomarker
# B from its marginal, and class drawn from the exact p(C | x1, B).
# Site covariates carry no class information and do not enter W.
#
# Note: pk is the mixture weight used to generate covariates. Because B is drawn
# from its marginal and class is then reassigned by Bayes' theorem, the realised
# class prevalences P(C = c) differ from pk. The identity
#   sum_c P(C = c) * tau_tilde_c = sum_c P(C = c) * tau_c   ( = E{tau(X)} )
# holds with the REALISED prevalences, and is checked below.
#
# Optional: if the aggregated comparator results are available, the script also
# reports BCF/BART bias and point-estimate RMSE relative to tau_c and tau_tilde_c.

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
})

set.seed(20261109)
n_mc <- 2e6   # Monte Carlo draws per prevalence condition

# ---------- DGM settings (as in functions.R) ----------
mu           <- list(c(-0.6, -0.6), c(0, 0.6), c(0.6, 0))
Sigma_within <- diag(2) * 0.25
p_bio        <- c(0.05, 0.50, 0.95)
C            <- 3

# ---------- Class probabilities pi_k(x) under the DGM ----------
class_probs <- function(pk, n) {
  comp <- sample(seq_len(C), n, replace = TRUE, prob = pk)
  x1   <- matrix(NA_real_, n, 2)
  for (k in seq_len(C)) {
    idx <- which(comp == k)
    if (length(idx) > 0)
      x1[idx, ] <- MASS::mvrnorm(length(idx), mu[[k]], Sigma_within)
  }
  B <- rbinom(n, 1, sum(pk * p_bio))
  lp <- sapply(seq_len(C), function(k)
    log(pk[k]) +
      mvtnorm::dmvnorm(x1, mean = mu[[k]], sigma = Sigma_within, log = TRUE) +
      dbinom(B, 1, p_bio[k], log = TRUE))
  lp <- lp - apply(lp, 1, max)
  P  <- exp(lp)
  P / rowSums(P)                       # n x C matrix of pi_k(x)
}

# W[c, k] = E{ pi_k(X) | C = c } = E{ pi_c(X) pi_k(X) } / E{ pi_c(X) }
# (class is drawn from pi(x), so conditioning on C = c weights by pi_c(x))
w_matrix <- function(P) crossprod(P) / colSums(P)

# ---------- Prevalence conditions ----------
pk_list <- list(baseline   = c(0.30, 0.30, 0.40),
                imbalanced = c(0.10, 0.20, 0.70))

W_list    <- list()
prev_list <- list()
for (nm in names(pk_list)) {
  P <- class_probs(pk_list[[nm]], n_mc)
  W_list[[nm]]    <- w_matrix(P)
  prev_list[[nm]] <- colMeans(P)       # realised class prevalence
  rm(P); invisible(gc())
}

cat("Realised class prevalences P(C = c):\n")
print(round(do.call(rbind, prev_list), 4))
cat("\nW = E{ pi_k(X) | C = c }, rows c, columns k:\n")
for (nm in names(W_list)) { cat(nm, "\n"); print(round(W_list[[nm]], 4)) }

# ---------- Simulation conditions (continuous outcome) ----------
conds <- data.frame(
  label  = c("N=200", "N=400", "N=800", "N=1000",
             "pi=(0.10,0.20,0.70)",
             "tau=(-7,0,7)", "tau=(-6,0,6)", "tau=(-3,0,3)", "tau=(-2,0,2)",
             "tau=(0,2,5)"),
  pk     = c(rep("baseline", 4), "imbalanced", rep("baseline", 5)),
  stringsAsFactors = FALSE
)
conds$tau <- list(c(-5, 0, 5), c(-5, 0, 5), c(-5, 0, 5), c(-5, 0, 5),
                  c(-5, 0, 5),
                  c(-7, 0, 7), c(-6, 0, 6), c(-3, 0, 3), c(-2, 0, 2),
                  c(0, 2, 5))

tt_rows <- lapply(seq_len(nrow(conds)), function(i) {
  tau  <- conds$tau[[i]]
  W    <- W_list[[conds$pk[i]]]
  prev <- prev_list[[conds$pk[i]]]
  tt   <- as.vector(W %*% tau)
  data.frame(condition = conds$label[i], class = c("Harm", "Null", "Benefit"),
             tau = tau, tau_tilde = tt, gap = tt - tau,
             ate_check_tau = sum(prev * tau), ate_check_tilde = sum(prev * tt))
})
tau_tilde_tab <- do.call(rbind, tt_rows)

cat("\nWithin-class targets tau_tilde_c:\n")
print(transform(tau_tilde_tab,
                tau_tilde = round(tau_tilde, 3), gap = round(gap, 3),
                ate_check_tau = round(ate_check_tau, 4),
                ate_check_tilde = round(ate_check_tilde, 4)),
      row.names = FALSE)

# Identity check: population-weighted tau_tilde must recover the ATE
max_dev <- max(abs(tau_tilde_tab$ate_check_tau - tau_tilde_tab$ate_check_tilde))
cat(sprintf("\nMax |sum P(C=c) tau_tilde_c - sum P(C=c) tau_c| over conditions: %.2e\n",
            max_dev))

write.csv(tau_tilde_tab, "tau_tilde_targets.csv", row.names = FALSE)

# ---------- Optional: BCF/BART performance against both targets ----------
res_file <- "results_combined_comp.rds" # stored results from HPC 
if (file.exists(res_file)) {
  rc  <- readRDS(res_file)
  key <- function(tau) paste(round(tau, 6), collapse = "_")
  pk_key <- c(baseline = "0.3_0.3_0.4", imbalanced = "0.1_0.2_0.7")

  perf <- list()
  for (i in seq_len(nrow(conds))) {
    tau <- conds$tau[[i]]
    tt  <- as.vector(W_list[[conds$pk[i]]] %*% tau)
    d   <- rc[rc$pk_str == pk_key[[conds$pk[i]]] & rc$tau_str == key(tau), ]
    N_i <- switch(conds$label[i], "N=200" = 200, "N=400" = 400,
                  "N=1000" = 1000, 800)
    d   <- d[d$N == N_i, ]
    if (nrow(d) == 0) next
    for (m in c("bcf", "bart")) for (k in 1:3) {
      b  <- d[[paste0("tau_", k, "_bias_", m)]]     # posterior median - tau
      bt <- b - (tt[k] - tau[k])                    # posterior median - tau_tilde
      perf[[length(perf) + 1]] <- data.frame(
        condition = conds$label[i], method = toupper(m),
        class = c("Harm", "Null", "Benefit")[k],
        mean_est   = mean(d[[paste0("tau_", k, "_med_", m)]], na.rm = TRUE),
        bias       = mean(b,  na.rm = TRUE),
        rmse       = sqrt(mean(b^2,  na.rm = TRUE)),   # point-estimate RMSE
        tau_tilde  = tt[k],
        bias_tilde = mean(bt, na.rm = TRUE),
        rmse_tilde = sqrt(mean(bt^2, na.rm = TRUE)),
        n_reps     = sum(!is.na(b)))
    }
  }
  perf_tab <- do.call(rbind, perf)
  cat("\nBCF/BART performance against tau_c and tau_tilde_c:\n")
  print(perf_tab, digits = 3, row.names = FALSE)
  write.csv(perf_tab, "comparator_perf_tau_tilde.csv", row.names = FALSE)
} else {
  cat("\n", res_file, "not found: skipping comparator performance.\n")
}

sessionInfo()
