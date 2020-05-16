#!/usr/bin/env Rscript

########################################
## filename: figure_mask_pp_check.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot supplementary figure of individual
## mask fitted trajectories with prior
## and posterior predictive checks
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'cowplot'     # publication ready ggplot
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
integrity_data_path <- args[1]
mask_results_path <- args[2]
plotting_params_path <- args[3]
outpath <- args[4]

## read data / style files
source(plotting_params_path) # shared project style
text_size <- 35

integrity_chains <- readRDS(mask_results_path)
integrity_dat <- read_csv(integrity_data_path,
                          col_types = cols())

#################################
## overall plot styling
#################################
set.seed(989327) # reproducible! (since we use random draws)



##################################################
## calculate predicted integrity (needed
## for panel order)
##################################################
intercept_predictions <- integrity_chains %>%
    spread_draws(predicted_intercept[treatment_id]) %>%
    mutate(fit_factor_pred = 1 / (1 - predicted_intercept),
           decon_number = 0) %>%
    rename(predicted_integrity = predicted_intercept)

decay_predictions <- integrity_chains %>%
    spread_draws(predicted_integrity[treatment_id, decon_number]) %>%
    ungroup() %>%
    mutate(fit_factor_pred = 1 / (1 - predicted_integrity))


integrity_predictions <- bind_rows(intercept_predictions,
                         decay_predictions)

integrity_predictions <- integrity_dat %>%
    distinct(treatment_id,
             decon_number,
             .keep_all = TRUE) %>%
    inner_join(integrity_predictions,
               by = c("treatment_id", "decon_number"))


## sort facets by the posterior median decay
## of the first virus in our virus order
## in the plotted data
cat("sorting panels by treatment integrity...\n")

ordering <- integrity_predictions %>%
    group_by(treatment) %>%
    filter(decon_number == 3) %>%
    summarise(med = median(predicted_integrity)) %>%
    arrange(med)

integrity_dat$treatment <- factor(
    integrity_dat$treatment,
    levels = ordering$treatment)


###################################
## plot mask trajectories
###################################
cat('plotting individual mask trajectories...\n')


obs_predictions <- integrity_chains %>%
    spread_draws(post_pred_obs_integrity[obs_id]) %>%
    ungroup() %>%
    mutate(fit_factor_uncensored = 1 / (1 - post_pred_obs_integrity),
           below_limit = fit_factor_uncensored < 200,
           fit_factor_censored = ifelse(
               below_limit,
               round(fit_factor_uncensored),
               200)
           )

mask_order <- integrity_dat %>%
    filter(decon_number == 3) %>%
    group_by(mask_id) %>%
    summarise(ff_rank = fit_factor) %>%
    ungroup()


obs_predictions <- obs_predictions %>%
    inner_join(integrity_dat,
               by = "obs_id") %>%
    ungroup()

integrity_labels <- integrity_dat %>%
    distinct(mask_id,
             .keep_all = TRUE) %>%
    inner_join(mask_order,
               by = "mask_id") %>%
    arrange(desc(ff_rank)) %>%
    group_by(treatment) %>%
    mutate(treatment_mask_id =
               row_number()) %>%
    ungroup()

integrity_dat <- integrity_dat %>%
    inner_join(integrity_labels %>%
               select(mask_id,
                      treatment_mask_id),
               by = "mask_id") %>%
    mutate(below_limit = fit_factor < 200) %>%
    arrange(treatment_mask_id) %>%
    ungroup()

random_draws <- sample(max(obs_predictions$.draw), 25)

print(random_draws)

obs_predictions <- obs_predictions %>%
    inner_join(integrity_labels %>%
               select(mask_id,
                      treatment_mask_id),
               by = "mask_id") %>%
    filter(.draw %in% random_draws) %>%
    group_by(.draw, mask_id) %>%
    mutate(line_id = group_indices()) %>%
    ungroup()

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 24,
                         "TRUE" = 21)))

fig <- obs_predictions %>% ggplot(aes(
                               x = factor(decon_number),
                               y = fit_factor_censored,
                               shape = below_limit,
                               fill = treatment,
                               group = line_id)) +
    geom_line(data = integrity_dat,
              mapping = aes(
                  x = factor(decon_number),
                  y = fit_factor,
                  group = mask_id),
              color = "black",
              size = 2,
              alpha = 0.75) +
    geom_point(data = integrity_dat,
               mapping = aes(
                   x = factor(decon_number),
                   y = fit_factor,
                   fill = treatment,
                   shape = below_limit,
                   group = mask_id),
               size = 8,
               stroke = 1,
               alpha = 0.8) +
    geom_hline(aes(yintercept = 100),
               color = "red",
               linetype = "longdash",
               size = 1.5) +
    geom_line(color = interval_grey,
              size = 0.5,
              alpha = 0.1) +
    geom_point(size = 2,
               stroke = 1,
               alpha = 0.1) +
    shape_scale +
    scale_fill_manual(values = unlist(treatment_colors)) +
    scale_x_discrete() +
    scale_y_continuous(trans = "log2",
                       breaks = c(25, 50, 100, 200, 400)) + 
    coord_cartesian(ylim = c(25, 600)) +
    facet_grid(vars(treatment),
               vars(treatment_mask_id)) +
    theme_project(base_size = text_size) +
    xlab("Number of decontaminations") + 
    ylab("Mask filtration (fit factor)") +
    theme(legend.position = "none")



####################################
## compose full figure from panels
####################################

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          fig,
          base_height = 15,
          base_asp = 1)
warnings()
