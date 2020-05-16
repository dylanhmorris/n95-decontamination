#!/usr/bin/env Rscript

###################################
## filename: decon_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of titer decay
## for decontamination experiments
####################################

hyperparam_list <- list(
    intercept_prior_mean = 4.5, # official dose was 10^5 TCID50
    intercept_prior_sd = 3,
    decay_rate_prior_mean = 0,
    decay_rate_prior_sd = 20,
    lower_lim_decay_rate = 0,
    first_dilution = 0,
    last_dilution = 7,
    n_wells = 4,
    debug = FALSE)

fixed_seed <- 23032
inits_seed <- 32
adapt_delta <- 0.9
