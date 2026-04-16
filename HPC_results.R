#Set working directory to dissertation
setwd("~/Desktop/PhD/Dissertation/Paper 2 - LCM/Simulations_Code")

# Load the cmdrstanr library & neededx packages
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(matrixStats)
library(MASS)
library(dplyr)

#################### Summary Functions ################################
summarise_tau <- function(data, c) {
  med  <- data[[paste0("tau_", c, "_med")]]
  lo   <- data[[paste0("tau_", c, "_q5")]]
  hi   <- data[[paste0("tau_", c, "_q95")]]
  bias <- data[[paste0("tau_", c, "_bias")]]
  cov  <- data[[paste0("tau_", c, "_cover90")]]
  
  tibble(
    class      = c,
    mean_est   = mean(med,  na.rm = TRUE),
    emp_sd     = sd(med,    na.rm = TRUE),
    mean_bias  = mean(bias, na.rm = TRUE),
    rmse       = sqrt(mean(bias^2, na.rm = TRUE)),
    coverage90 = mean(cov,  na.rm = TRUE),
    mean_width = mean(hi - lo, na.rm = TRUE)
  )
}
summarise_misclass <- function(data) {
  pairs <- list(
    c(1, 2, "Harm -> Null"),
    c(1, 3, "Harm -> Benefit (catastrophic)"),
    c(2, 1, "Null -> Harm"),
    c(2, 3, "Null -> Benefit"),
    c(3, 1, "Benefit -> Harm (catastrophic)"),
    c(3, 2, "Benefit -> Null")
  )
  
  tibble(
    misclassification = sapply(pairs, `[[`, 3),
    rate = sapply(pairs, function(p) {
      col <- paste0("misclass_", p[1], "_to_", p[2], "_trt")
      if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
    })
  ) %>%
    mutate(rate_pct = scales::percent(rate, accuracy = 0.1))
}
summarise_confidence <- function(data) {
  expand.grid(
    class = 1:3,
    arm   = c("trt", "ctl")
  ) %>%
    mutate(
      true_class    = factor(class, labels = c("Harm", "Null", "Benefit")),
      treated       = factor(arm,   levels = c("trt", "ctl"),
                             labels = c("Treated", "Control")),   # <-- fixed order
      mean_entropy  = mapply(function(c, a) {
        col <- paste0("entropy_class", c, "_", a)
        if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
      }, class, arm),
      mean_max_prob = mapply(function(c, a) {
        col <- paste0("maxprob_class", c, "_", a)
        if (col %in% names(data)) mean(data[[col]], na.rm = TRUE) else NA_real_
      }, class, arm)
    ) %>%
    select(true_class, treated, mean_entropy, mean_max_prob) %>%
    arrange(true_class, treated)
}
summarise_all <- function(data) {
  
  tau_summary <- bind_rows(lapply(1:3, function(c) summarise_tau(data, c)))
  
  diag_summary <- data %>%
    summarise(
      conv_rate          = mean(max_rhat < 1.01 & n_divergent == 0 & min_ess_bulk > 200, na.rm = TRUE),
      mean_max_rhat      = mean(max_rhat,    na.rm = TRUE),
      median_runtime_sec = median(runtime_sec, na.rm = TRUE)
    )
  
  class_summary <- data %>%
    summarise(
      mean_class_acc     = mean(class_acc,     na.rm = TRUE),
      mean_class_acc_trt = mean(class_acc_trt, na.rm = TRUE),
      mean_class_acc_ctl = mean(class_acc_ctl, na.rm = TRUE)
    )
  
  per_class_summary <- data %>%
    summarise(
      class_1_acc     = mean(class_1_acc,     na.rm = TRUE),
      class_2_acc     = mean(class_2_acc,     na.rm = TRUE),
      class_3_acc     = mean(class_3_acc,     na.rm = TRUE),
      class_1_acc_trt = mean(class_1_acc_trt, na.rm = TRUE),
      class_2_acc_trt = mean(class_2_acc_trt, na.rm = TRUE),
      class_3_acc_trt = mean(class_3_acc_trt, na.rm = TRUE)
    )
  
  prior_summary <- data %>%
    summarise(
      mean_prior_acc     = mean(prior_acc,     na.rm = TRUE),
      mean_prior_acc_trt = mean(prior_acc_trt, na.rm = TRUE),
      mean_prior_acc_ctl = mean(prior_acc_ctl, na.rm = TRUE)
    )
  
  per_class_prior_summary <- data %>%
    summarise(
      prior_class_1_acc     = mean(prior_class_1_acc,     na.rm = TRUE),
      prior_class_2_acc     = mean(prior_class_2_acc,     na.rm = TRUE),
      prior_class_3_acc     = mean(prior_class_3_acc,     na.rm = TRUE),
      prior_class_1_acc_trt = mean(prior_class_1_acc_trt, na.rm = TRUE),
      prior_class_2_acc_trt = mean(prior_class_2_acc_trt, na.rm = TRUE),
      prior_class_3_acc_trt = mean(prior_class_3_acc_trt, na.rm = TRUE)
    )
  
  mean_prob_summary <- data %>%
    summarise(
      post_prob_class1  = mean(mean_pclass_1_correct, na.rm = TRUE),
      post_prob_class2  = mean(mean_pclass_2_correct, na.rm = TRUE),
      post_prob_class3  = mean(mean_pclass_3_correct, na.rm = TRUE),
      prior_prob_class1 = mean(mean_prior_1_correct,  na.rm = TRUE),
      prior_prob_class2 = mean(mean_prior_2_correct,  na.rm = TRUE),
      prior_prob_class3 = mean(mean_prior_3_correct,  na.rm = TRUE)
    )
  
  misclass_summary <- summarise_misclass(data)
  confidence_summary = summarise_confidence(data)
  
  list(
    tau_summary             = tau_summary,
    diag_summary            = diag_summary,
    class_summary           = class_summary,
    per_class_summary       = per_class_summary,
    prior_summary           = prior_summary,
    per_class_prior_summary = per_class_prior_summary,
    mean_prob_summary       = mean_prob_summary,
    misclass_summary        = misclass_summary,
    confidence_summary      = confidence_summary   # <-- new
  )
}

#################### Check HPC Results ################################
results <- readRDS("results_combined.rds")
by_cond <- split(results, results$cond_id)

invisible(lapply(names(by_cond), function(id) {
  cond <- by_cond[[id]]
  cat("\n========================================\n")
  cat(sprintf("Condition %s | N=%d | pk=%s | tau=%s\n", 
              id, cond$N[1], cond$pk_str[1], cond$tau_str[1]))
  cat("========================================\n")
  out <- summarise_all(cond)
  invisible(lapply(names(out), function(nm) { 
    cat("\n>", nm, "\n")
    print(out[[nm]]) 
  }))
}))