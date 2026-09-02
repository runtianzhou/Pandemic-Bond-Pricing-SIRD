# Pandemic Bond Pricing under Stochastic SIRD Modelling

This repository contains the R code supporting the preprint

**Pandemic Bond Pricing under Stochastic SIRD Modelling**

by **Runtian Zhou, Gail S. K. Wolkowicz, and Traian A. Pirvu**.

The manuscript is currently a preprint and has been submitted for journal publication.

> **Software implementation and repository:** Runtian Zhou

## Overview

This project develops a pandemic bond pricing framework based on a stochastic SIRD epidemic model in which the disease transmission rate is modelled by a mean-reverting Ornstein–Uhlenbeck (OU) process.

The repository contains the computational implementation used to generate the numerical results, tables, figures, and sensitivity analyses reported in the manuscript.

The main components of the analysis include:

- reproduction of benchmark results based on the stochastic logistic growth epidemic model;
- implementation of the stochastic SIRD–OU epidemic framework;
- recursive valuation of pandemic bond cash flows;
- numerical simulation of stochastic epidemic trajectories;
- Monte Carlo estimation of pandemic bond prices and trigger probabilities;
- sensitivity analysis, including partial rank correlation coefficient (PRCC) analysis; and
- generation of the numerical tables and figures reported in the manuscript.

## Repository Structure

The main computational files are contained in the `R/` directory.

```text
Pandemic-Bond-Pricing-SIRD/
├── R/                  # R scripts, data inputs, and intermediate computational objects
├── renv/               # renv activation and project settings
├── .Rprofile           # activates the project-specific R environment
├── .gitignore
├── CITATION.cff        # citation metadata
├── LICENSE             # MIT License
├── README.md
└── renv.lock           # R package versions used in the project
