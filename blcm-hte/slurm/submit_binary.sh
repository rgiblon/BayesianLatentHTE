#!/bin/bash
# submit_binary.sh — SLURM job array for binary outcome BLCA simulation study
# Mirrors submit.sh exactly — points to run_sim_binary.R and results_binary/
#
# Usage:
#   sbatch submit_binary.sh                      # full run
#   sbatch --array=45,112,334 submit_binary.sh   # resubmit specific failed jobs
#
# ── SLURM directives ─────────────────────────────────────────────────────────
#SBATCH --job-name=blcm_binary
#SBATCH --account=def-heathann
#SBATCH --array=1-10000                   # ← set to nrow(conditions_binary.csv)
#SBATCH --time=0:20:00                    # binary ~62s/rep, 20min is generous
#SBATCH --mem=3G
#SBATCH --cpus-per-task=3                 # 1 per chain (parallel_chains=3)
#SBATCH --output=logs_binary/sim_%A_%a.out
#SBATCH --error=logs_binary/sim_%A_%a.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rachel.giblon@mail.utoronto.ca

# ── Environment ───────────────────────────────────────────────────────────────
module load gcc/12.3
module load r/4.4.0

# Project directory — edit to match your path on the cluster
export PROJECT_DIR=$HOME/scratch/lcm_sim

# ── Run ───────────────────────────────────────────────────────────────────────
cd $PROJECT_DIR
mkdir -p results_binary logs_binary

Rscript run_sim_binary.R $SLURM_ARRAY_TASK_ID

# ── Notes ─────────────────────────────────────────────────────────────────────
# BEFORE full submission:
#   1. Run make_conditions_binary.R locally to generate conditions_binary.csv
#   2. Update --array=1-N where N = nrow(conditions_binary.csv)
#   3. Run pilot: sbatch --array=1-10 submit_binary.sh
#   4. Check timing: grep "Elapsed" logs_binary/sim_*_1.out
#      OR: sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
#   5. Binary model is faster than continuous (~62s vs ~90s) so 20min is safe
#
# RESUBMITTING FAILURES:
#   After running aggregate_binary.R, check missing_jobs_binary.txt then:
#   sbatch --array=$(cat missing_jobs_binary.txt | tr '\n' ',') submit_binary.sh
