#!/usr/bin/env Rscript

###################################
## filename: mask_hyperparams.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: contains hyperparameters for
## the Bayesian model of mask fit/
## filtration decay with repeated
## use/decontamination
####################################

hyperparam_list <- list(
    upper_ff_limit = 200,  # data censored above 200
    ff_interval = 1, # data are rounded to integers
    fit_decay_mean = 0,
    fit_decay_sd = 0.33,
    treatment_effect_mean = 0.25, # treatments expected to hurt
    treatment_effect_sd = 0.5,
    mask_intercept_filt_mean_mean = .995,
    mask_intercept_filt_mean_sd = 0.0025,
    logit_mask_intercept_filt_sd_mean = 0,
    logit_mask_intercept_filt_sd_sd = 0.5, # 3 sigma --> uniform on [0, 1]
    sd_log_decay_mean = 0,
    sd_log_decay_sd = log(1.5) / 3, # 3 sigma --> uniform on [0, 1]
    sd_logit_obs_mean = 0,
    sd_logit_obs_sd = 0.5, # 3 sigma --> uniform on [0, 1]
    debug = debug)
