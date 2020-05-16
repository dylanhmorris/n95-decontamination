#!/usr/bin/env Rscript

###################################
## filename: table_decon.R
## description: programmatically
## generate a table showing important
## info from decon experiments
####################################

script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidyr',      # for spread()
    'tidybayes',  # for spread_draws()
    'huxtable',   # for programmatic results output
    'officer',    # for outputing tables to docx
    'flextable'   # for table output
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly=TRUE)
data_path <- args[1]
results_path <- args[2]
plotting_params_path <- args[3]
outpath <- args[4]


## read data
chains <- readRDS(results_path)

source(plotting_params_path) # shared project style

dat <- read_data_for_plotting(data_path)
sorting_virus <- "SARS-CoV-1"

time_to_decay <- chains %>%
    spread_draws(decay_rate[trial_unique_id]) %>%
    mutate(
        time_to_half = 60 * log10(2) / decay_rate,
        time_to_thousand = 60 * 3 / decay_rate,
        time_to_million = 60 * 6 / decay_rate) %>%
    ungroup()

time_to_decay <- dat %>%
    distinct(trial_unique_id, .keep_all = TRUE) %>%
    select(trial_unique_id,
           virus,
           material,
           treatment) %>%
    inner_join(time_to_decay, by = "trial_unique_id") %>%
    group_by(trial_unique_id)

ordering <- time_to_decay %>%
    group_by(treatment) %>%
    summarise(med = median(-decay_rate)) %>%
    arrange(med)

full_summary <- time_to_decay %>%
    ungroup() %>%
    group_by(treatment, material) %>%
    summarise(
        hl_median = quantile(time_to_half, 0.5),
        hl_q025 = quantile(time_to_half, 0.025),
        hl_q975 = quantile(time_to_half, 0.925),
        thous_median = quantile(time_to_thousand, 0.5),
        thous_q025 = quantile(time_to_thousand, 0.025),
        thous_q975 = quantile(time_to_thousand, 0.925),
        mil_median = quantile(time_to_million, 0.5),
        mil_q025 = quantile(time_to_million, 0.025),
        mil_q975 = quantile(time_to_million, 0.925))

#################################
## Construct and format table
#################################

output_table <- as_hux(full_summary)

print(output_table)

header_row <- c("Treatment", "Material", rep(c("median", "2.5%", "97.5%"), 3))
unit_row <- c("", "",
              "half-life (min)", "", "",
              "time to one thousanth (min)", "", "",
              "time to one millionth (min)", "", "")

output_table <- output_table %>%
    insert_row(header_row,
               after = 0) %>%
    insert_row(unit_row,
               after = 0) %>%
    merge_cells(1, 1:2) %>%
    merge_cells(1, 3:5) %>%
    merge_cells(1, 6:8) %>%
    merge_cells(1, 9:11) %>%
    merge_cells(3:4, 1) %>%
    merge_cells(5:6, 1) %>%
    merge_cells(7:8, 1) %>%
    merge_cells(9:10, 1) %>%
    merge_cells(11:12, 1)

    


## style the table a bit

bottom_border(output_table)[2, ] <- 1



cat('Saving table to ', outpath, '...\n')
quick_docx(output_table,
           file = outpath)
warnings()
