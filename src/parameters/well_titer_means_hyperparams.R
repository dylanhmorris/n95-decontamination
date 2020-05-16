#!/usr/bin/env Rscript

###################################
## filename: well_titer_means_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model used to infer mean replicate
## titers directly from the raw well data
####################################

hyperparam_list <- list(
    titer_mean_prior_mean = 3,
    titer_mean_prior_sd = 5,
    titer_sd_prior_mean = 0,
    titer_sd_prior_sd = 3,
    debug = FALSE)

fixed_seed <- 235902
inits_seed <- 876543
