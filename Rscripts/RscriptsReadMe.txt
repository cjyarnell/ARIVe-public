This folder contains the R and Stan scripts necessary to perform the ARIVe analysis.

1) CreateCohortARIVe.R is a script that preprocesses, combines, transforms, centers, and scales the data from the two cohorts.
2) ARIVe_analysis.R is a script that calls the Stan program to fit the Bayesian multistate model, calls the Stan program to perform the mediation analysis, and extracts the primary and secondary outcomes from the results, including figures.
3) ARIVe_4state_msm.stan is a Stan program that fits the parametric Bayesian multistate model
4) ARIVe_stan_piecewise_functions_joint.R is an R script containing functions used by ARIVe_analysis.R
5) ARIVe_sensitivity_analyses.R is an R script that conducts the database-specific sensitivity analyses