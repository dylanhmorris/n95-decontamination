#!/usr/bin/env Rscript

########################################
## filename: figure_decon_controls.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot figure showing control decay all
## the way to undetectability
#######################################


script_packages <- c(
    'rstan',      # stan interface
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'modelr',     # for data_grid()
    'tidyr',      # for crossing()
    'cowplot',    # publication ready ggplot
    'scales',     # for trans_breaks(), etc.
    'forcats'     # for fct_relevel()
)

## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}

#################################
# functions
#################################

logit <- function(x) {
    return( log(x / (1 - x)) )
}

inv_logit <- function(x) {
    return( 1 / (1 + exp(-x)) )
}


#################################
# read in needed data
#################################

## read command line args
args <- commandArgs(trailingOnly=TRUE)
decay_data_path <- args[1]
integrity_data_path <- args[2]
decay_results_path <- args[3]
titer_means_path <- args[4]
mask_results_path <- args[5]
plotting_params_path <- args[6]
outpath <- args[7]

## read data / style files
source(plotting_params_path) # shared project style

decay_chains <- readRDS(decay_results_path)
titer_means_chains <- readRDS(titer_means_path)
integrity_chains <- readRDS(mask_results_path)
dat <- read_csv(decay_data_path,
                col_types = cols())
integrity_dat <- read_csv(integrity_data_path,
                          col_types = cols())

#################################
## overall plot styling
#################################
set.seed(989327) # reproducible! (since we use random draws)
text_size = 30
detection_linesize = 0.75
titer_ylab <- expression("Virus titer (TCID"[50] * "/mL media)")
bar_width <- 0.75 # width for bars in barplot
bar_cap_width <- 0.5 # cap width for error bar caps in barplot
bar_err_size <- 0.75 # line thickness for error bars in barplot
scatter_err_size <- 1 # line thickness for error bars in scatterplot
dodge_width <- 0.5 # how much to dodge dots
convert_mL_media_to_L_air <- 10 / 3



##################################################
## calculate predicted integrity (needed
## for panel order)
##################################################
intercept_predictions <- integrity_chains %>%
    spread_draws(predicted_intercept[treatment_id]) %>%
    mutate(fit_factor_pred = 1 / (1 - predicted_intercept),
           decon_number = 0) %>%
    rename(predicted_integrity = predicted_intercept)
print(intercept_predictions)

decay_predictions <- integrity_chains %>%
    spread_draws(predicted_integrity[treatment_id, decon_number]) %>%
    ungroup() %>%
    mutate(fit_factor_pred = 1 / (1 - predicted_integrity))
print(decay_predictions)


integrity_predictions <- bind_rows(intercept_predictions,
                         decay_predictions)

## integrity_predictions <- integrity_chains %>%
##     spread_draws(filtration_means[treatment_id, decon_id]) %>%
##     rename(predicted_integrity = filtration_means) %>%
##     mutate(fit_factor_pred = 1 / (1 - predicted_integrity)) %>%
##     mutate(decon_number = decon_id - 1)

integrity_predictions <- integrity_dat %>%
    distinct(treatment_id, decon_number,
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

print(ordering)

integrity_predictions$treatment <- factor(
    integrity_predictions$treatment,
    levels = ordering$treatment)

dat$treatment <- factor(
    dat$treatment,
    levels = ordering$treatment)
integrity_dat$treatment <- factor(
    integrity_dat$treatment,
    levels = ordering$treatment)
dat$material <- factor(
    dat$material,
    levels = c("Steel",
               "N95 mask"))


##################################################
## calculate posterior draws for regression lines
##################################################

decon_numbers <- tibble(decon_number = 0:3)

plot_times <- dat %>%
    data_grid(time = seq_range(c(0, 60),
                               n = fineness))

## get needed draws and add human readable names
tidy_draws <- decay_chains %>%
    spread_draws(decay_rate[trial_unique_id],
                 intercept[trial_unique_id][replicate])

pos_wells <- dat %>%
    group_by(replicate_set_id) %>%
    summarise(
        n_wells = n(),
        n_pos = sum(virus_detect))

tidy_draws <- dat %>%
    distinct(trial_unique_id,
             .keep_all = TRUE) %>%
    select(trial_unique_id,
           virus,
           material,
           treatment) %>% 
    inner_join(tidy_draws,
               by = "trial_unique_id")

titer_means_draws <- titer_means_chains %>%
    spread_draws(set_mean_log10_tcid50[replicate_set_id])


## get human readable names and detectability
titer_means_draws <- dat %>%
    distinct(time,
             replicate_set_id,
             .keep_all = TRUE) %>%
    select(replicate_set_id,
           virus,
           time,
           material,
           treatment) %>% 
    inner_join(titer_means_draws,
               by = "replicate_set_id") %>%
    inner_join(pos_wells,
               by = "replicate_set_id") %>%
    mutate(detectable = n_pos > 1)

## check that small and undetectable titers have
## higher variance if there are fewer well
print(titer_means_draws %>%
      filter(n_pos <= 1) %>%
      group_by(replicate_set_id) %>%
      mutate(med = median(set_mean_log10_tcid50),
             q025 = quantile(set_mean_log10_tcid50, 0.025),
             q975 = quantile(set_mean_log10_tcid50, 0.975)) %>%
      ungroup() %>%
      distinct(replicate_set_id, .keep_all = TRUE) %>%
      select(replicate_set_id,
             treatment,
             time,
             n_pos,
             n_wells,
             med,
             q025,
             q975),
      n = 50)

## convert from TCID50/(0.1mL) to TCID50/mL
## and visualize 0 positive well titers at
## the traditional LOD
LOD_log10 = 0.5
LOD = 10^LOD_log10
titer_means_draws <- titer_means_draws %>%
    mutate(set_mean_log10_tcid50 = ifelse(
               detectable,
               set_mean_log10_tcid50 + 1,
               LOD_log10)) %>%
    ungroup() %>%
    arrange(desc(time), desc(material)) ## so masks are on top 

print(titer_means_draws)

ytrans <- 'log10'
ylim <- c(LOD/10, 2 * 10^5)

ybreaks <- trans_breaks('log10', function(x) 10^x)
yformat <- trans_format(ytrans,
                        math_format(10^.x))

###################################
## plot panel showing raw surface
## data
###################################
cat('plotting raw data...\n')

###################################
## plot panel showing fit of
## regression lines to real data
###################################
cat('plotting regression lines...\n')
## draw n_lines random regression lines
func_samples <- tidy_draws %>%
    group_by(trial_unique_id,
             replicate) %>%
    sample_n(n_lines) %>%
    ungroup()

## annotate lines so that each
## has a unique id for ggplot overplotting
## (else two lines from the same draw but
## different replicates can get confused
## with each other)
func_samples <- func_samples %>%
    mutate(line_id = as.numeric(rownames(func_samples)))

## cross product decay_rates with x (time) values
## and calculate y (titer) values
cat('setting up x values...\n')

to_plot <- func_samples %>%
    filter(treatment == "Control") %>%
    crossing(plot_times)

## adding one to convert to per mL from per mL/10
to_plot <- to_plot %>%
    filter(treatment == "Control") %>%
    mutate(predicted_titer = 10^(1 + intercept - decay_rate * time))

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

lines_steel <- to_plot %>%
    filter(material == "Steel")

lines_mask <- to_plot %>%
    filter(material == "N95 mask")

## nudge steel points at the LOD slightly
## to avoid overplotting
steel_sub <- titer_means_draws %>%
    filter(material == "Steel") %>%
    filter(treatment == "Control")

mask_sub <- titer_means_draws %>%
    filter(material == "N95 mask") %>%
    filter(treatment == "Control")

fit_panel <- lines_mask %>%
    ggplot(aes(x = time,
               y = predicted_titer,
               color = material,
               group = line_id)) +
    geom_hline(aes(yintercept = LOD),
               size = 2,
               linetype = "dotted") +
    geom_line(data = lines_steel,
              alpha = line_alpha) +
    geom_line(alpha = line_alpha) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(x = time,
                      y = 10^set_mean_log10_tcid50,
                      shape = detectable,
                      fill = material,
                      group = replicate_set_id),
        data = steel_sub,
        point_size = 6,
        size = 8,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    stat_pointinterval(
        .width = 0.6827,
        mapping = aes(x = time,
                      y = 10^set_mean_log10_tcid50,
                      shape = detectable,
                      fill = material,
                      group = replicate_set_id),
        data = steel_sub,
        point_size = 6,
        size = 14,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    stat_pointinterval(
        .width = 0.95,
        mapping = aes(x = time,
                      y = 10^set_mean_log10_tcid50,
                      shape = detectable,
                      fill = material,
                      group = replicate_set_id),
        data = mask_sub,
        point_size = 6,
        size = 8,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    stat_pointinterval(
        .width = 0.6827,
        mapping = aes(x = time,
                      y = 10^set_mean_log10_tcid50,
                      shape = detectable,
                      fill = material,
                      group = replicate_set_id),
        data = mask_sub,
        point_size = 6,
        size = 14,
        stroke = 2,
        interval_color = interval_grey,
        interval_alpha = 1,
        color = "black",
        alpha = 0.9) +
    scale_fill_manual(values = unlist(material_colors)) +
    scale_fill_manual(values = unlist(material_colors),
                       aesthetics = "point_fill") +
    scale_color_manual(values = unlist(material_colors)) +
    shape_scale + 
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       label = yformat) +
    coord_cartesian(ylim = ylim,
                    xlim = c(0, 55)) +
    facet_wrap(vars(treatment), nrow = 1)
    
## styling: no facet labels because is background plot
fit_panel <- fit_panel +
    guides(fill = guide_legend(override.aes = list(alpha = 1,
                                                   stroke = 1,
                                                   point_size = 15,
                                                   shape = 21)),
           shape = "none") +
    theme_project(base_size = text_size) +
    theme(legend.key.size = unit(5, "lines"),
          legend.position = c(0.75, 0.7),
          legend.margin = margin(t = -10,
                                 l = 0,
                                 r = 5,
                                 b = 0),
          legend.title = element_blank(),
          legend.box.background = element_rect(color = "black",
                                               size = 2))+
    xlab("Time (hrs)") +
    ylab(titer_ylab)


####################################
## compose full figure from panels
####################################

labeled_panel_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("0.5", "in")) 

int_int_theme <- theme(
    strip.background = element_blank(),
    strip.text.x = element_text(size = text_size),
    strip.placement = "outside",
    strip.switch.pad.grid = unit("-0.5", "in"))

cat('making full figure...\n')

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          fit_panel + labeled_panel_theme,
          base_width = 10,
          base_height = 8)
warnings()
