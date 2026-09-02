# Pandemic Bond Pricing under Stochastic SIRD Modelling

This repository contains the computational code supporting the preprint:

**Pandemic Bond Pricing under Stochastic SIRD Modelling**

by **Runtian Zhou, Gail S. K. Wolkowicz, and Traian A. Pirvu**.

The manuscript is currently a preprint and has been submitted for journal publication.

> **Software implementation and repository maintenance:** Runtian Zhou

## Overview

This repository provides the R implementation used for the numerical analysis in the associated manuscript.

The study develops a pandemic bond pricing framework based on a stochastic SIRD epidemic model in which the disease transmission rate evolves according to a mean-reverting Ornstein–Uhlenbeck (OU) process.

The computational work includes:

- benchmark calculations based on the stochastic logistic growth model;
- simulation of the stochastic SIRD–OU epidemic model;
- epidemic-trigger calculations;
- recursive pandemic bond valuation;
- Monte Carlo estimation of bond prices and trigger probabilities;
- generation of numerical tables and figures; and
- parameter sensitivity analysis using partial rank correlation coefficients (PRCC).

## Repository Structure

```text
Pandemic-Bond-Pricing-SIRD/
├── R/                  # R scripts, data inputs, and computational objects
├── renv/               # renv activation and project settings
├── .Rprofile           # activates the project-specific R environment
├── .gitignore
├── CITATION.cff        # software citation metadata
├── LICENSE             # MIT License
├── README.md
└── renv.lock           # recorded R package versions
```

The main computational files are located in the [`R/`](R/) directory.

## Reproducibility

This repository uses [`renv`](https://rstudio.github.io/renv/) to record the R package environment used for the analysis.

### 1. Clone the repository

```bash
git clone https://github.com/runtianzhou/Pandemic-Bond-Pricing-SIRD.git
cd Pandemic-Bond-Pricing-SIRD
```

### 2. Restore the R environment

Open R or RStudio in the repository directory and run:

```r
install.packages("renv")
renv::restore()
```

This restores the package versions recorded in `renv.lock`.

### 3. Run the analysis

The scripts in `R/` correspond to the main computational components of the manuscript, including:

- stochastic logistic growth benchmark calculations;
- stochastic SIRD–OU simulations;
- pandemic-trigger calculations;
- recursive pandemic bond pricing;
- numerical tables and figures; and
- PRCC sensitivity analysis.

A script-to-manuscript output map will be maintained in this repository to make individual results easier to reproduce.

## Methodological Components

### Stochastic SIRD–OU Model

The epidemic component is based on a stochastic SIRD framework with a mean-reverting Ornstein–Uhlenbeck process for the disease transmission rate.

### Pandemic Bond Pricing

The simulated epidemic dynamics are incorporated into a recursive pandemic bond pricing framework in which epidemic-trigger events affect future bond cash flows.

### Numerical Simulation

The stochastic model is evaluated numerically using Euler–Maruyama-type simulation methods.

### Monte Carlo Analysis

Monte Carlo simulation is used to estimate epidemic-trigger probabilities and pandemic bond values under stochastic epidemic dynamics.

### Sensitivity Analysis

Partial rank correlation coefficient (PRCC) analysis is used to investigate the sensitivity of model outputs to epidemiological and stochastic-model parameters.

## Associated Manuscript

**Zhou, R., Wolkowicz, G. S. K., and Pirvu, T. A.**  
*Pandemic Bond Pricing under Stochastic SIRD Modelling.*  
Preprint, 2026. Submitted for journal publication.

The manuscript contains the mathematical formulation, theoretical development, pricing framework, numerical methodology, and interpretation of the results implemented in this repository.

## Software Authorship

The R implementation, computational organization, and GitHub repository were developed and are maintained by:

**Runtian Zhou**

Software authorship is distinct from authorship of the associated manuscript.

## Citation

Software citation metadata are provided in [`CITATION.cff`](CITATION.cff).

GitHub also provides a **Cite this repository** option in the repository sidebar.

If you use this code in academic work, please cite:

1. the archived software release; and
2. the associated manuscript.

## Archived Version

An archived version of this research software is available through Zenodo:

https://zenodo.org/records/19685421

The Zenodo record provides a persistent archival version of the software and its metadata.

## Data

Some data files included in or referenced by this repository originate from third-party public data sources.

Such data remain subject to the terms and conditions of their original providers and are not necessarily covered by the MIT License that applies to the software in this repository.

## License

The software in this repository is released under the **MIT License**.

See [`LICENSE`](LICENSE) for the full license terms.

## Keywords

Pandemic bond · Catastrophe bond · Recursive pandemic bond pricing · Stochastic SIRD model · Ornstein–Uhlenbeck process · Stochastic differential equations · Euler–Maruyama method · Monte Carlo simulation · PRCC sensitivity analysis
