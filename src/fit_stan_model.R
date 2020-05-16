#!/usr/bin/env Rscript

###################################
## fit the specified stan model
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

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


## read command line args
args <- commandArgs(trailingOnly=TRUE)
model_src_path <- args[1]
well_data_path <- args[2]
hyperparam_path <- args[3]
mcmc_output_path <- args[4]



## read data
cat('reading in well data from file ', well_data_path, ' ...\n')
dat <- read_csv(well_data_path,
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
max_tree <- 15
set.seed(inits_seed) ## R's rng set for random inits

debug <- TRUE

if("debug" %in% names(hyperparam_list) ) {
    debug <- hyperparam_list$debug
}

if (debug) {
    nchains <- 1
    niter <- 500
}

adapt_d <- 0.8

## calculate replicate counts

## observations are either status
## (+/-) of wells or log10 estimated
## titers

## calculate transformed data
replicates <- dat %>%
    group_by(trial_unique_id) %>%
    summarise(replicates = max(replicate),
              condition_name = factor(paste(virus, material))[1])

replicates$condition_id = as.numeric(replicates$condition_name)

replicates <- replicates %>%
    group_by(condition_id) %>%
    add_tally()

n_experiments <- max(dat$trial_unique_id)
n_experiments_check <- length(unique(dat$trial_unique_id))


if(abs(n_experiments - n_experiments_check) > 0.5){
    cat('warning: different number of unique ids from max unique id:',
        'max:', n_experiments, 'number unique', n_experiments_check,
        '\n')
}

n_conditions <- max(replicates$condition_id)
n_experiments_in_condition <- replicates$n


## make observation lists and set initial values,
#depending on which kind of data is used

## some models use raw well data

titer_replicate_set_id <- dat %>%
    distinct(titer_id, replicate_set_id) %>%
    arrange(titer_id) 
titer_replicate_set_id <- titer_replicate_set_id$replicate_set_id

titer_experiment_id <- dat %>%
    distinct(titer_id, trial_unique_id) %>%
    arrange(titer_id) 
titer_experiment_id <- titer_experiment_id$trial_unique_id

titer_replicate_id <- dat %>%
    distinct(titer_id, replicate) %>%
    arrange(titer_id) 
titer_replicate_id <- titer_replicate_id$replicate

titer_time <- dat %>%
    distinct(titer_id, time) %>%
    arrange(titer_id) 
titer_time <- titer_time$time

if("control_index" %in% names(dat)){
    titer_control_index <- dat %>%
        distinct(titer_id, control_index) %>%
        arrange(titer_id) 
    titer_control_index <- titer_control_index$control_index

    control_index <- dat$control_index
    
} else {
    titer_control_index = -99
    control_index = -99
}

if("treatment_end_time" %in% names(dat)){
    titer_treatment_end_time <- dat %>%
        distinct(titer_id, treatment_end_time) %>%
        arrange(titer_id) 
    titer_treatment_end_time <- titer_treatment_end_time$treatment_end_time

    treatment_end_time <- dat$treatment_end_time
    
} else {
    titer_treatment_end_time <- -99
    treatment_end_time <- -99
}

observation_data_list <- list(
    n_total_datapoints = length(dat$virus_detect),
    well_status = dat$virus_detect,
    dilution = dat$dilution,
    titer_id = dat$titer_id,
    titer_replicate_set_id = titer_replicate_set_id,
    titer_experiment_id = titer_experiment_id,
    titer_time = titer_time,
    titer_replicate_id = titer_replicate_id,
    titer_control_index = titer_control_index,
    titer_treatment_end_time = titer_treatment_end_time,
    n_replicate_sets = max(titer_replicate_set_id),
    n_experiments = n_experiments,
    n_conditions = n_conditions,
    n_replicates = replicates$replicates[1],
    n_experiments_in_condition = n_experiments_in_condition,
    condition_id = replicates$condition_id,
    experiment_id = dat$trial_unique_id,
    replicate_id = dat$replicate,
    control_index = control_index,
    treatment_end_time = treatment_end_time,
    times = dat$time,
    n_titers = length(titer_time))

init_val <- "random"

init_val <- lapply(
    1:nchains,
    function(x){ list(decay_rate = runif(n_experiments, 0, 0.2))})

## pass stan data and hyperparams
stan_data <- c(
    observation_data_list,
    hyperparam_list)


###############################
## Compile, fit, and save model
###############################
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
