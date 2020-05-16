#!/usr/bin/env Rscript

###################################
## filename: well_titer_estimates_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model used to infer titers non
## hierarchically from the well data
####################################

hyperparam_list <- list(
    titer_prior_mean = 3,
    titer_prior_sd = 3,
    debug = FALSE)

inits_seed <- 32
fixed_seed <- 2542
