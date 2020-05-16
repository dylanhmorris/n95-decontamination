#!/usr/bin/env Rscript

###################################
## filename: fit_mask_model.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
## description: fit the specified stan model
## for mask fit factor decay
## to the specified clean dataset,
## and save the result as an .Rds
## file
##
####################################

script_packages <- c(
    'rstan',     # stan interface
    'parallel',  # parallelized MCMC
    'readr',     # csv read-in
    'magrittr',  # for pipe operator %>%
    'dplyr',     # for filter()
    'tidyr'      # for drop_na()
)


## set up hyperparameters for models

debug <- FALSE ## set to TRUE to diagnose sampler problems


## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <- args[1]
mask_data_path <- args[2]
hyperparam_path <- args[3]
mcmc_output_path <- args[4]


## read data
cat('reading in mask data from file ', mask_data_path, ' ...\n')
dat <- read_csv(mask_data_path,
                col_types = cols())
cat('data loaded successfully!\n')

cat('reading in hyperparameters from file ', hyperparam_path, ' ...\n')
source(hyperparam_path)
cat('hyperparameters loaded successfully!\n')

## set stan options
n_cores <- parallel::detectCores()
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)
niter <- 2000
nchains <- n_cores
adapt_d <- 0.95
max_tree <- 15
fixed_seed <- 2332
inits_seed <- 323122
set.seed(inits_seed) ## R's rng set for random inits

if (debug) {
    nchains <- 1
    niter <- 500
}



mask_t_ids <- dat %>%
    distinct(mask_id, treatment_id) %>%
    arrange(mask_id) %>%
    select(mask_id, treatment_id)

print(mask_t_ids, n = 100)

mask_treatment_id <- mask_t_ids$treatment_id

## make observation lists and set initial values
n_treatments <- max(dat$treatment_id)
n_subjects <- max(dat$subject_id)
n_masks <- max(dat$mask_id)
max_n_decons <- max(dat$decon_number)
n_observations <- length(dat$fit_factor)



observation_data_list <- list(
    n_subjects = n_subjects,
    n_treatments = n_treatments,
    n_masks = n_masks,
    max_n_decons = max_n_decons,
    n_observations = n_observations,
    subject_id = dat$subject_id,
    treatment_id = dat$treatment_id,
    mask_id = dat$mask_id,
    mask_treatment_id = mask_treatment_id,
    decon_number = dat$decon_number,
    obs_fit = dat$fit_factor,
    n_datapoints = length(dat$fit_factor))

init_val <- "random"

init_val <- lapply(
    1:nchains,
    function(x){ list(treatment_effect = matrix(
                          runif(n_treatments *
                                max_n_decons,
                                0, 0.1), nrow = n_treatments),
                      sd_logit_obs = runif(1, 20, 30)) })


###############################
## Compile, fit, and save model
###############################
## pass stan data and hyperparams
stan_data <- c(
    observation_data_list,
    hyperparam_list)

fit <- stan(
    model_src_path,
    data = stan_data,
    iter = niter,
    seed = fixed_seed,
    chains = nchains,
    init = init_val,
    control = list(max_treedepth = max_tree,
                   adapt_delta = adapt_d))

cat('\nsaving mcmc samples to', mcmc_output_path, '\n')

saveRDS(fit, mcmc_output_path)

warnings()
