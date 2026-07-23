# Federated ADMM-LASSO
This repository contains the implementation and mathematical derivation of a federated ADMM-LASSO framework for high-dimensional sparse linear regression across multiple data centers.

Each participating center keeps its patient-level data locally and computes a center-specific coefficient update using only its own feature matrix and outcome vector. The local model-level updates are then coordinated through the Alternating Direction Method of Multipliers (ADMM) to obtain a common global sparse coefficient vector. The central aggregation step applies soft-thresholding to perform LASSO-based variable selection, while the dual update enforces consensus across centers.

The framework enables collaborative coefficient estimation, variable selection, and prediction without pooling or directly sharing patient-level data. Its privacy-preserving property is based on federated data locality rather than formal cryptographic privacy guarantees.

## Mathematical Derivation

The complete Federated ADMM-LASSO formulation, update derivations, convergence diagnostics, and algorithm are available on the project webpage:

[View the Federated ADMM-LASSO webpage](https://atikur616.github.io/Federated-ADMM-for-LASSO-Regression/)

## Repository Files

- `index.html` — mathematical formulation, derivation, convergence diagnostics, and algorithm
- `ADMM-LASSO_FL.txt` — implementation of the Federated ADMM-LASSO method
- `README.md` — project overview
