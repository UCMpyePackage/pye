# pye: Penalized Youden Index Estimator R Package

The **pye** package implements the penalized variable selection and classification methodologies developed during my PhD research at the Complutense University of Madrid (UCM).

The package features the **Penalized Youden Index Estimator (PYE)** and the **Covariate-adjusted Youden Index (covYI)** estimator, designed for simultaneous feature selection, coefficient estimation, and threshold optimization, both in low- and high-dimensional settings.

## 📌 Repository Structure & Versions

* **`main` (Version 0.1.0 - 2026):** Clean production-ready codebase optimized according to CRAN guidelines. Prepared for submission.

## ⚙️ Core Features Included

* **Primary Estimators:** Smooth kernel-based estimation for standard/weighted PYE and covYI.
* **Sparsity Penalties:** Full implementation of $L_{1/2}$, $L_1$ (Lasso), Elastic-Net, SCAD, and MCP regularization.
* **Optimization Engines:** Modified monotone (mmAPG) and non-monotone (mnmAPG) Accelerated Proximal Gradient solvers.
* **Cross-Validation & Simulations:** Integrated automated grid search for hyperparameter tuning ($\lambda$ and $\tau$) and parallelized train-test validation splits.
* **Benchmark Modules:** Standardized interfaces for penalized logistic regression, penalized Support Vector Machines (SVM), and penalized AUC-based regression.

## 🚀 Installation (for Advisors & Collaborators)

To install the package directly from this private repository, make sure you have `devtools` installed and run the following command in R/RStudio:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("UCMpyePackage/pye")
