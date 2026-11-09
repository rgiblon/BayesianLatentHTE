#!/bin/bash
# submit_comp.sh — SLURM job array for Bayesian LCM simulation study comparators - BART & BCF
#
# Usage:
#   sbatch submit_comp.sh                      # full run
#   sbatch --array=45,112,334 submit_comp.sh   # resubmit specific failed jobs
#
# ── SLURM directives ─────────────────────────────────────────────────────────
#SBATCH --job-name=blcm_comp
#SBATCH --account=def-heathann
#SBATCH --array=1-9999                    # ← set to nrow(conditions.csv)
#SBATCH --time=0:45:00
#SBATCH --mem=3G                          # cmdstanr with 3 chains needs ~5-6G
#SBATCH --cpus-per-task=3                 # 1 per chain (parallel_chains=3)
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err
#SBATCH --mail-type=FAIL                  # email only on failure (optional)
#SBATCH --mail-user=rachel.giblon@mail.utoronto.ca

# ── Environment ───────────────────────────────────────────────────────────────
module load gcc/12.3    # required for Stan C++ compilation
module load r/4.4.0

# Project directory — edit to match your path on the cluster
export PROJECT_DIR=$HOME/scratch/lcm_sim

# ── Run ───────────────────────────────────────────────────────────────────────
cd $PROJECT_DIR
mkdir -p results_comp logs

Rscript run_sim_comp.R $SLURM_ARRAY_TASK_ID

# ── Notes ─────────────────────────────────────────────────────────────────────
# BEFORE full submission:
#   1. Run pilot: sbatch --array=1-10 submit.sh
#   2. Check timing: grep "Elapsed" logs/sim_*_1.out
#      OR: sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
#   3. Set --time generously (1.5x your observed runtime) to avoid timeouts
#   4. Set --array=1-$(wc -l < conditions.csv) to match your grid exactly
#
# INSTALLING R PACKAGES (do once in an interactive job, not here):
#   salloc --time=1:00:00 --mem=8G --cpus-per-task=4 --account=def-heathann
#   module load gcc/12.3 r/4.4.0
#   R
#   > install.packages(c("cmdstanr","posterior","bayesplot","matrixStats",
#                        "MASS","dplyr","mvtnorm"),
#                      lib = "~/R/library",
#                      repos = "https://cloud.r-project.org")
#   > cmdstanr::install_cmdstan()
