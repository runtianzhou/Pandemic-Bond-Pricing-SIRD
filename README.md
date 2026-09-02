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
```

The `R/` directory contains the computational material for the benchmark stochastic logistic growth model, the stochastic SIRD–OU framework, pandemic-trigger calculations, pandemic bond valuation, numerical experiments, figures, tables, and sensitivity analysis.

## Reproducibility

This repository uses [`renv`](https://rstudio.github.io/renv/) to record the R package environment used for the computational analysis.

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

The analysis files are located in:

```text
R/
```

The scripts cover the principal computational components of the manuscript, including:

- benchmark stochastic logistic growth calculations;
- stochastic SIRD–OU simulations;
- epidemic-trigger calculations;
- recursive pandemic bond pricing;
- numerical tables and figures; and
- PRCC sensitivity analysis.

Because several scripts correspond to specific numerical experiments and manuscript outputs, users should consult the script names and comments within the `R/` directory when reproducing individual results.

## Methodological Components

### Stochastic SIRD–OU Epidemic Model

The proposed epidemic framework is based on a stochastic SIRD model in which the disease transmission rate evolves according to a mean-reverting Ornstein–Uhlenbeck process.

The stochastic specification is intended to incorporate time-varying uncertainty in disease transmission while retaining the compartmental structure of the SIRD model.

### Pandemic Bond Pricing

Simulated epidemic trajectories are incorporated into a recursive pandemic bond pricing framework.

Pandemic-dependent trigger conditions affect future bond cash flows, allowing epidemic dynamics and financial valuation to be considered jointly.

### Numerical Simulation

The stochastic epidemic system is simulated numerically using Euler–Maruyama-type methods, with particular attention to the positivity requirements of epidemiological state variables.

### Monte Carlo Analysis

Monte Carlo simulation is used to evaluate stochastic epidemic paths, trigger events, and pandemic bond values under parameter uncertainty.

### Sensitivity Analysis

Parameter sensitivity is investigated using partial rank correlation coefficients (PRCC), allowing the relative influence of epidemiological and stochastic parameters on model outputs and pandemic bond values to be assessed.

## Associated Manuscript

**Zhou, R., Wolkowicz, G. S. K., and Pirvu, T. A.**  
*Pandemic Bond Pricing under Stochastic SIRD Modelling.*  
Preprint, 2026. Submitted for journal publication.

The manuscript provides the mathematical formulation, theoretical development, pandemic bond pricing framework, numerical methodology, and interpretation of the computational results implemented in this repository.

## Software Authorship

The R implementation, computational repository, and software organization were developed and are maintained by:

**Runtian Zhou**

The authorship of the associated manuscript is separate from the authorship of the software implementation.

## Citation

Citation metadata for the software are provided in [`CITATION.cff`](CITATION.cff).

GitHub also provides a **Cite this repository** option in the repository sidebar.

If you use the code in academic work, please cite:

1. the archived software release; and
2. the associated manuscript.

## Archived Version

An archived version of this research software is available through Zenodo:

https://zenodo.org/records/19685421

The Zenodo archive provides a persistent record of the software and associated metadata.

## Data

Some data files included or referenced in this repository originate from third-party public data sources. Such data remain subject to the terms and conditions of their original providers and are not necessarily covered by the MIT License applicable to the software in this repository.

## License

The code in this repository is released under the **MIT License**.

See [`LICENSE`](LICENSE) for the full license terms.

## Keywords

Pandemic bond · Catastrophe bond · Recursive pandemic bond pricing · Stochastic SIRD model · Ornstein–Uhlenbeck process · Stochastic differential equations · Euler–Maruyama method · Monte Carlo simulation · PRCC sensitivity analysis
