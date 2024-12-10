# Bayesian Latent Class Model for Treatment Heterogeneity

This repository contains R code for simulating and implementing a Bayesian latent class model to detect heterogeneity of treatment effect (HTE) in critical care patients. The model uses a mixed-mode latent class approach and leverages Hamiltonian Monte Carlo (HMC) through RStan for robust Bayesian inference.

---

## Table of Contents
1. [Project Overview](#project-overview)
2. [Features](#features)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Methodology](#methodology)
6. [Contributing](#contributing)
7. [License](#license)

---

## Project Overview
Understanding treatment heterogeneity is essential for tailoring medical interventions to specific patient subgroups. This Bayesian latent class model:
- Identifies unobserved subgroups (latent classes) in clinical trial populations.
- Jointly models treatment allocation and outcomes.
- Uses Bayesian inference and Hamiltonian Monte Carlo for precise estimation.

The approach supports adaptive trial designs and precision medicine by uncovering latent treatment-response patterns.

---

## Features
- Bayesian latent class modeling implemented in R.
- Efficient sampling via Hamiltonian Monte Carlo (HMC) using **RStan**.
- Flexible model specification with customizable priors and covariates.
- Simulation-based validation and performance metrics.
- Application-ready framework for real-world critical care datasets.

---

## Installation
To set up the project locally, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/username/bayesian-latent-class-model.git
   cd bayesian-latent-class-model

install.packages(c("rstan", "ggplot2", "dplyr"))
