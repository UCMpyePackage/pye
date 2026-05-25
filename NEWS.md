# pye 0.1.0

- **Initial Release**: First official release of the `pye` package, providing a unified toolkit for high-dimensional binary classification, feature selection, and covariate adjustment.

## New Implemented Methods

- **Penalized Youden index Estimator (PYE)**: Introduced an embedded feature selection method for low- and high-dimensional binary classification ($p \gg n$) that directly maximizes a differentiable, Kernel-Smoothed (KS) version of the Youden Index using a standard normal CDF kernel.
- **Covariate-adjusted Youden Index (covYI)**: Implemented an adaptive extension to incorporate covariates, allowing for observation-specific thresholding ($t_i = c_i^\top \gamma$) and automated covariate selection.

## Optimization & Penalties

- **Accelerated Proximal Gradient (APG)**: Implemented two efficient optimization algorithms tailored for non-convex and non-smooth objective functions: `mmAPG` (modified monotone variant) and `mnmAPG` (non-monotone variant).
- **Sparsity-Inducing Penalties**: Integrated closed-form proximal operators for a wide range of penalty functions, including $L_{1/2}$ norm, $L_1$ (Lasso), Elastic-Net, SCAD, and MCP.

## Core Functionality & Benchmarking Suite

- **Model Estimation**: Added core routines `pye_KS_estimation` and `covYI_KS_estimation` to perform simultaneous feature selection and coefficient estimation.
- **Unified Benchmarking**: Included wrapper functions to estimate and compare established high-dimensional binary decision engines under identical data handling: Penalized Logistic Regression (`plr_estimation`), Penalized Support Vector Machines (`psvm_estimation`), and AUC-based methods (`AucPR_estimation`).

## Tuning, Utilities & Simulations

- **Hyperparameter Selection**: Added automated $k$-fold cross-validation routines (`pye_KS_compute_cv`, `plr_compute_cv`, `psvm_compute_cv`, `AucPR_compute_cv`) to optimize tuning parameters ($\lambda$ and $\tau$) across grid searches.
- **Data Generation & Validation**: Included `create_sample_with_covariates` to generate synthetic high-dimensional datasets with controlled correlation structures.
- **Simulation Wrappers**: Added `pye_simulation_study` and `model_simulation_study` to automate repeated train-test splits for evaluating selection stability and performance metrics under varying sparsity constraints.
