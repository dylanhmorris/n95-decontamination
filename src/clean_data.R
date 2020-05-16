#!/usr/bin/env Rscript

#####################################
## name: clean_data.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## process raw data and save cleaned
## data for use in model fitting
##
####################################


script_packages <- c(
    'readr',     # csv read-in
    'dplyr',     # grouping and combining
    'tidyr',     # for replace_na
    'magrittr'   # for %>% 
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}

## functions

round_nearest <- function(x, stepsize){
    return( round(x / stepsize) * stepsize )
}

## specify needed cleaning
clean_titer_data <- function(dat){

    dat$n_wells <- ifelse((grepl("MERS", dat$virus) &
                           dat$material == "Copper"),
                          3,
                          4)

    dat$detection_limit_log10_titer <- ifelse(
        dat$material == "Copper",
        1.5,
        0.5)
    
    log_titer_stepsize = 1 / dat$n_wells
    
    if (!"log10_titer" %in% names(dat)){

        ## undetectable titers are always
        ## reported as 0.5 log10 units
        ## (3.16... unlogged), but others
        ## are in intervals of either 1/3s
        ## or 1/4s
        dat$log10_titer = ifelse(
            dat$titer > 3.17,
            0.5 + round_nearest(
                      log10(dat$titer) - 0.5,
                      log_titer_stepsize),
            0.5)
    }

    if (!"strain" %in% names(dat)) {
        dat$strain <- dat$virus

    }

    dat <- dat %>%
        filter(virus != "MERS-CoV")

    ## for multi-virus fitting,
    ## include both the current
    ## MERS-CoV EMC/12 strain experiments
    ## and the prior ones (but as different trials)
    if ("xp" %in% names(dat)) {
        dat <- dat %>%
            filter(xp == "new" |
                   strain == "MERS-CoV_EMC/12") %>%
            mutate(trial_unique_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH ",
                       sub("MERS_", "", xp)),
                   grouping_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"))
        
        dat <- dat %>%
            arrange(virus,
                    strain,
                    temperature,
                    humidity,
                    material,
                    xp) %>%
            group_by(virus,
                     strain,
                     temperature,
                     humidity,
                     material,
                     xp) %>%
            mutate(trial_unique_id = group_indices())
        
        dat <- dat %>%
            arrange(virus,
                    material,
                    temperature,
                    humidity,
                    replicate,
                    time,
                    xp) %>%
            group_by(virus,
                     material,
                     temperature,
                     humidity,
                     replicate,
                     time,
                     xp) %>%
            mutate(titer_id = group_indices())


    } else {
        dat <- dat %>%
            mutate(trial_unique_name = paste0(
                       virus,
                       " ",
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"),
                   grouping_name = paste0(
                       material,
                       " ",
                       temperature,
                       "°C ",
                       humidity,
                       "% RH"))

        dat <- dat %>%
            arrange(virus,
                    strain,
                    temperature,
                    humidity,
                    material) %>%
            group_by(virus,
                     strain,
                     temperature,
                     humidity,
                     material) %>%
            mutate(trial_unique_id =
                       group_indices())

        dat <- dat %>%
            arrange(virus,
                    material,
                    temperature,
                    humidity,
                    replicate,
                    time) %>%
            group_by(virus,
                     material,
                     temperature,
                     humidity,
                     replicate,
                     time) %>%
            mutate(titer_id = group_indices())
    }

    dat <- dat %>% arrange(titer_id)
    return(dat)
}


clean_well_data <- function(dat) {

    if( grepl("env", raw_data_path) ){
        
        dat <- dat %>%
            filter(material %in% c("Plastic", "Aerosols"))
    }


    ## for now, do not estimate dose deposited
    ## loss rate -- maybe add if everything
    ## else works
    dat <- dat %>% filter(time >= min_time,
                          !grepl("cytotoxic", comment) &
                          virus == "SARS-CoV-2" &
                          !(material %in% c("N95 mask",
                                           "DMEM")))

    ## add naming and trial ids
    dat <- dat %>%
        mutate(trial_unique_name = paste0(
                   virus,
                   " ",
                   material,
                   " ",
                   temperature,
                   "°C ",
                   humidity,
                   "% RH"),
               grouping_name = paste0(
                   material,
                   " ",
                   temperature,
                   "°C ",
                   humidity,
                   "% RH"),
               )
    dat$treatment <- dat$treatment %>%
        replace_na("Untreated")

    dat <- dat %>%
        filter(treatment == "Untreated")
    dat <- dat %>%
        group_by(virus,
                 temperature,
                 humidity,
                 material) %>%
        mutate(trial_unique_id = group_indices())

    dat <- dat %>%
        group_by(virus,
                 material,
                 temperature,
                 humidity,
                 replicate,
                 time) %>%
        mutate(titer_id = group_indices())

    dat <- dat %>%
        group_by(virus,
                 material,
                 temperature,
                 humidity,
                 time) %>%
        mutate(replicate_set_id = group_indices())


    dat <- dat %>% arrange(titer_id)
    return (dat)
}

clean_decon_data <- function(dat) {

    decon_materials <- c(
        "DMEM uncovered plate oven",
        "DMEM closed vial heat block",
        "DMEM closed vial oven",
        "DMEM covered plate oven",
        "DMEM",
        "Steel",
        "N95 mask")
    
    ## for now, do not estimate dose deposited
    ## loss rate -- maybe add if everything
    ## else works
    
    dat <- dat %>% filter(time >= min_time &
                          !grepl("cytotoxic", comment) &
                          material %in% decon_materials &
                          virus == "SARS-CoV-2")

    dat$treatment <- dat$treatment %>%
        replace_na("Control")


    ## VHP treatment ended at 7 minutes,
    ## then measured at 97
    dat$treatment_end_time <- ifelse(
        dat$treatment == "VHP",
        7 / 60,
        99999)

    print(dat %>% distinct(treatment))

    ## add naming and trial ids
    dat <- dat %>%
        mutate(trial_unique_name = paste0(
                   virus,
                   " ",
                   material,
                   " ",
                   treatment)
               )

    dat <- dat %>%
        group_by(virus, material, treatment) %>%
        mutate(trial_unique_id = group_indices())

    dat <- dat %>%
        group_by(virus,
                 material,
                 treatment,
                 replicate,
                 time) %>%
        mutate(titer_id = group_indices())

    dat <- dat %>%
        group_by(virus,
                 material,
                 treatment,
                 time) %>%
        mutate(replicate_set_id = group_indices())
    
    dat <- dat %>%
        group_by(virus,
                 material) %>%
        mutate(control_index =
                   ifelse(
                       length(trial_unique_id[treatment == "Control"]) > 0,
                       first(trial_unique_id[treatment == "Control"]),
                       trial_unique_id))
    
    dat <- dat %>% arrange(titer_id) %>%
        select(virus,
               material,
               treatment,
               time,
               treatment_end_time,
               dilution,
               replicate,
               trial_unique_id,
               control_index,
               replicate_set_id,
               titer_id,
               virus_detect)

    return (dat)
}


## load and clean data
args <- commandArgs(trailingOnly=TRUE)
raw_data_path <- args[1]
outpath <- args[2]
min_time <- args[3]

if(is.na(min_time)){
    min_time <- 0
}

delim <- ";"

col_types <- cols()

if ( grepl("decon", raw_data_path) |
     grepl("env", raw_data_path) ) {
    col_types <- cols_only(
        virus = col_factor(),
        material = col_factor(),
        treatment = col_character(),
        time = col_double(),
        dilution = col_integer(),
        replicate = col_integer(),
        virus_detect = col_integer(),
        temperature = col_character(),
        humidity = col_double(),
        comment = col_character())
}

dat <- read_delim(raw_data_path,
                  delim = delim,
                  col_types = col_types)
if ( grepl("decon", raw_data_path))  {

    cleaned <- clean_decon_data(dat)

} else if ( any(grepl("dilution", names(dat)))) {
    cat("cleaning well data...\n")
    cleaned <- clean_well_data(
        dat)

} else {

    cleaned <- clean_titer_data(dat)
}



write_csv(cleaned,
          outpath)

warnings()
