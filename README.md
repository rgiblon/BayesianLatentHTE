# Bayesian Latent Class Model for Heterogeneity of Treatment Effect

Code accompanying:

> Giblon R, Liu K, Heath A. *A Bayesian latent class model for heterogeneity of treatment effect: ordered harm, null and benefit subgroups in randomized trials.* [Journal, year, DOI to be added.]

The Bayesian latent class model (BLCM) estimates class-specific treatment effects in randomized trials jointly with each participant's probability of belonging to a harm, no-effect or benefit class. Prognostic effects are shared across classes, a gap parametrization orders the class-specific effects, and a multinomial membership model gives each participant a class probability from baseline covariates. Because treatment is randomized, the treatment-assignment mechanism is not modelled.

This repository contains the Stan programs, the simulation study (continuous and binary outcomes), the comparator analyses (reduced latent class model, BCF, BART) and the ACTG175 application.

## Repository structure

```
stan/
  BLCA_fixedtrt_CmdStanR_v2.stan      BLCM, continuous outcome (simulation; s_mid = 0.2)
  BLCA_fixedtrt_CmdStanR_v3.stan      BLCM, continuous outcome (ACTG175; s_mid = 0.5)
  BLCA_binary_CmdStanR.stan           BLCM, binary outcome (s_mid = 0.5)
  BLCA_reduced_CmdStanR_v1.stan       Reduced model with fixed mixing proportions
  BLCM_k1_homogeneous_shared.stan     One-class model, shared residual SD (ACTG175)
  BLCM_k1_homogeneous.stan            One-class model, arm-specific residual SDs (ACTG175)
R/         Data-generating mechanisms, simulation drivers, aggregation, tables, ACTG175 analysis, tau tilde
slurm/     SLURM array scripts used on the Digital Research Alliance of Canada Nibi cluster
results/   Simulation condition grids
```

## Requirements

- R (simulations run with R 4.4.0 on Nibi; aggregation in R 4.6.0)
- CmdStan 2.38.0 via `cmdstanr`
- R packages: `cmdstanr`, `posterior`, `bcf`, `dbarts`, `MASS`, `mvtnorm`, `dplyr`, `loo`, `speff2trial`

```r
install.packages(c("posterior", "bcf", "dbarts", "MASS", "mvtnorm", "dplyr", "loo", "speff2trial"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Reproducing the simulation study

Each model has its own pipeline: full BLCM (continuous), reduced model, comparators (BCF and BART) and BLCM (binary).

1. **Compile the Stan models:** `R/compile.R`
2. **Generate the condition grids:** `R/make_conditions.R` (continuous) and `R/make_conditions_binary.R` (binary). These write `conditions.csv` and `conditions_binary.csv`, included in `results/`. There are 10 conditions with 1,000 replicates each.
3. **Run the replicates.** Each replicate is an independent SLURM array job:

   | Model | Driver | Functions | SLURM script |
   |---|---|---|---|
   | BLCM (continuous) | `run_sim.R` | `functions.R` | `submit.sh` |
   | Reduced model | `run_sim_reduced.R` | `functions_reduced.R` | `submit_reduced.sh` |
   | BCF and BART | `run_sim_comp.R` | `functions.R` | `submit_comp.sh` |
   | BLCM (binary) | `run_sim_binary.R` | `functions_binary.R` | `submit_binary.sh` |

   Random seeds are set deterministically from the condition index and replicate number.
4. **Aggregate replicate output:** `aggregate.R`, `aggregate_reduced.R`, `aggregate_comp.R` and `aggregate_binary.R`
5. **Produce the manuscript tables:** `results_tables.R`

Simulation output is not included in this repository. Because seeds are set deterministically from the condition index and replicate number, running the pipeline reproduces the results reported in the paper. Aggregated results are available from the corresponding author on request.

The SLURM scripts are written for the Nibi cluster (three cores, 3 GB memory, 30-minute limit per replicate). Adjust the account and module settings for other systems.

## ACTG175 application

`R/blcm_actg175.R` fits the BLCM, the one-class comparison models, BCF and BART to the ACTG175 trial. The data are publicly available in the `speff2trial` R package.

## Licence

MIT. See `LICENSE`.

## Contact

Rachel Giblon, Dalla Lana School of Public Health, University of Toronto
