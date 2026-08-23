# Federated ADMM-LASSO
This repository contains the implementation and mathematical derivation of a federated ADMM-LASSO framework for high-dimensional sparse linear regression across multiple data centers.

Each participating center keeps its patient-level data locally and computes a center-specific coefficient update using only its own feature matrix and outcome vector. The local model-level updates are coordinated through the Alternating Direction Method of Multipliers (ADMM) to obtain a common global sparse coefficient vector. The central aggregation step applies soft-thresholding to perform LASSO-based variable selection, while the dual update enforces consensus across centers.

The framework enables collaborative coefficient estimation, variable selection, and prediction without pooling or directly sharing patient-level data. The current implementation provides data locality but does not include formal differential privacy, secure aggregation, or cryptographic privacy guarantees.

## Mathematical Derivation

The complete Federated ADMM-LASSO formulation, update derivations, convergence diagnostics, and algorithm are available on the project webpage:

[View the Federated ADMM-LASSO](https://atikur616.github.io/Fedearted-ADMM-for-LASSO-Regression/)

## Reproducibility

The repository includes reproducible R code for the primary simulation setting with p = 50 and n = 100 observations per center.

The simulation code includes:

- LASSO-1
- LASSO-2
- LASSO-Pool
- Federated ADMM-LASSO
- Regularization parameter tuning using minimum validation MSE
- ADMM convergence diagnostics
- Statistical comparisons between methods
- Rho-sensitivity analysis

The default code uses 100 Monte Carlo replications. 

## Repository Files

- `index.html` — mathematical formulation, derivation, convergence diagnostics, and algorithm
- `Federated_ADMM_LASSO_Simulation_Reproducibility.R` — reproducible R code for the simulation analysis
- `README.md` — project overview and reproducibility information
