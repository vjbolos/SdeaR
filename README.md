# SdeaR

**Stochastic Data Envelopment Analysis with R**

[![CRAN version](https://www.r-pkg.org/badges/version/SdeaR)](https://cran.r-project.org/package=SdeaR)
[![Downloads](https://cranlogs.r-pkg.org/badges/SdeaR)](https://cran.r-project.org/package=SdeaR)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

`SdeaR` is an open-source R package for **Stochastic Data Envelopment Analysis (SDEA)** based on **chance-constrained programming**.

Unlike conventional DEA approaches that assume deterministic input/output data, `SdeaR` explicitly incorporates **uncertainty** by modelling inputs and outputs as **multivariate normally distributed random variables**.

The package provides a comprehensive framework for estimating stochastic efficiency under alternative model specifications, including:

- Radial chance-constrained DEA models
- Radial super-efficiency chance-constrained models
- Directional chance-constrained DEA models
- Additive chance-constrained DEA models (E-model and P-model)
- Additive super-efficiency chance-constrained models

with support for:

- Constant Returns to Scale (CRS)
- Variable Returns to Scale (VRS)
- Non-increasing Returns to Scale (NIRS)
- Non-decreasing Returns to Scale (NDRS)
- Generalized Returns to Scale (GRS)

---

# Motivation

Data Envelopment Analysis (DEA) is a widely used nonparametric methodology for measuring the relative efficiency of homogeneous decision-making units (DMUs).

Classical DEA assumes exact observations. However, many real-world datasets contain **measurement uncertainty, noise, or stochastic variability**.

`SdeaR` fills an important gap in the DEA software ecosystem by providing a robust and user-friendly implementation of **chance-constrained DEA models for stochastic data**.

---

# Installation

Install from CRAN:

```r
install.packages("SdeaR")
```

Development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("vjbolos/SdeaR")
```

Load package:

```r
library(SdeaR)
library(deaR)
```

---

# Workflow

The standard workflow consists of:

1. Build a deterministic DEA dataset with `deaR`
2. Construct stochastic data with `make_deadata_stoch()`
3. Select a stochastic DEA model
4. Estimate efficiency scores
5. Extract and analyse results

---

# Minimal Example

Using the classical **Program Follow Through** dataset:

```r
library(SdeaR)
library(deaR)

data("PFT1981")

# Select Program Follow Through sites
PFT <- PFT1981[1:49, ]

# Create deterministic DEA data
PFT <- make_deadata(PFT, ni = 5, no = 3)

# Define stochastic output variances
c <- 0.5
var_output <- matrix(c^2, nrow = 3, ncol = 49)

# Create stochastic dataset
PFT_stoch <- make_deadata_stoch(
  datadea = PFT,
  var_output = var_output
)

# Run radial stochastic DEA model
results <- modelstoch_radial(PFT_stoch)

# Extract efficiencies
efficiencies(results)
```

---

# Main Functions

## Data preparation

- `make_deadata_stoch()`

Constructs stochastic DEA datasets from deterministic `deaR` objects and covariance specifications.

Supports:

- Full covariance matrix
- Input-input covariance arrays
- Output-output covariance arrays
- Input-output covariance arrays
- Diagonal variance structures

---

## Radial models

- `modelstoch_radial()`
- `modelstoch_radial_supereff()`

Implements chance-constrained radial DEA models and super-efficiency extensions.

---

## Directional models

- `modelstoch_dir()`
- `modelstoch_dir_dd()`

Supports:

- Stochastic directions
- Deterministic directions

---

## Additive models

- `modelstoch_additive()`
- `modelstoch_additive_p()`
- `modelstoch_addsupereff()`

Includes:

- Expected value models (E-models)
- Probability models (P-models)
- Super-efficiency additive models

---

# Documentation

Full reference manual:

https://cran.r-project.org/package=SdeaR

Package PDF manual:

https://cran.r-project.org/web/packages/SdeaR/SdeaR.pdf

---

# Mathematical Framework

All models assume:

- Inputs and outputs follow a **multivariate normal distribution**
- Chance constraints are transformed into deterministic equivalents
- Constraints are satisfied with probability at least **1 − α**

This allows explicit efficiency analysis under uncertainty while preserving tractable optimization formulations.

---

# Software Requirements

- R >= 3.5
- Depends on:

  - `deaR`

Platform-independent and tested on:

- Windows
- Linux
- macOS

---

# Citation

If you use `SdeaR` in academic work, please cite:

Bolós, V.J., Coll-Serrano, V., & Benítez, R. (2026)

**SdeaR: Stochastic Data Envelopment Analysis. R package version 1.0.2.**

---

# Related Work

The package extends previous work on:

- Chance-constrained DEA
- Stochastic efficiency analysis
- Directional stochastic DEA models
- Additive stochastic DEA formulations

and complements the deterministic DEA package:

- `deaR`

---

# License

GPL-3

---

# Authors

**Vicente J. Bolós**  
Department of Business Mathematics  
University of Valencia

**Vicente Coll-Serrano**  
Department of Applied Economics  
University of Valencia

**Rafael Benítez**  
Department of Business Mathematics  
University of Valencia

---

# Support

Questions and bug reports:

vicente.bolos@uv.es

GitHub issues:

https://github.com/vjbolos/SdeaR/issues
