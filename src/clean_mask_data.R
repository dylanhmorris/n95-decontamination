#!/usr/bin/env Rscript

#####################################
## name: clean_mask_data.R
## author: Dylan H. Morris <dhmorris@princeton.edu>
##
## description: process raw mask data
## and save cleaned data for use in
## model fitting
##
####################################


####################################
## dependency loading
####################################

script_packages <- c(
    'readr',     # csv read-in
    'dplyr'      # grouping and combining
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}

####################################
## functions
####################################

## specify needed cleaning
clean_mask_data <- function(dat){
    ## calculate ids
    dat <- dat %>%
        ungroup() %>%
        group_by(treatment) %>%
        mutate(treatment_id =
                   group_indices()) %>%
        ungroup() %>%
        group_by(treatment,
                 decon_number) %>%
        ungroup() %>%
        rename(subject_id = subject,
               mask_id = item_id)

    dat$obs_id = as.numeric(rownames(dat))
    return(dat)
}

####################################
## script body
####################################

## load data
args <- commandArgs(trailingOnly=TRUE)
raw_data_path <- args[1]

delim <- ";"
col_types <- cols()

dat <- read_delim(raw_data_path,
                  delim = delim,
                  col_types = col_types)

## clean and write data

cleaned <- clean_mask_data(dat)

print(cleaned)

write_csv(cleaned,
          args[2])

warnings()
