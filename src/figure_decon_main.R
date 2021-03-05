#!/usr/bin/env Rscript

########################################
## filename: figure_decon_main.R
## author: Dylan Morris <dhmorris@princeton.edu>
## plot main figure, which contains raw
## data, fitted kill rates, and mask
## integrity
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
text_size <- 35

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

decay_predictions <- integrity_chains %>%
    spread_draws(predicted_integrity[treatment_id, decon_number]) %>%
    ungroup() %>%
    mutate(fit_factor_pred = 1 / (1 - predicted_integrity))


integrity_predictions <- bind_rows(intercept_predictions,
                         decay_predictions)

integrity_predictions <- integrity_dat %>%
    distinct(treatment_id, decon_number,
             .keep_all = TRUE) %>%
    inner_join(integrity_predictions,
               by = c("treatment_id", "decon_number"))


## sort facets by the posterior median decay
## of the first virus in our virus order
## in the plotted data
cat("sorting panels by treatment integrity...\n")

dmem_name <- "DMEM closed vial heat block"

material_order <- c(
    "N95 mask",
    "Steel",
    "DMEM uncovered plate oven",
    "DMEM closed vial heat block",
    "DMEM closed vial oven",
    "DMEM covered plate oven",
    "DMEM")


dat$material[dat$material == dmem_name] = "DMEM"

materials_to_include <- c(
    "N95 mask",
    "Steel",
    "DMEM")

dat <- dat %>%
    filter(material %in% materials_to_include)


ordering <- integrity_predictions %>%
    group_by(treatment) %>%
    filter(decon_number == 3) %>%
    summarise(med = median(predicted_integrity)) %>%
    arrange(med)

integrity_dat$treatment <- factor(
    integrity_dat$treatment,
    levels = ordering$treatment)

integrity_predictions$treatment <- factor(
    integrity_predictions$treatment,
    levels = ordering$treatment)

dat$treatment <- factor(
    dat$treatment,
    levels = ordering$treatment)
dat$material <- factor(
    dat$material,
    levels = material_order)

##################################################
## calculate posterior draws for regression lines
##################################################

decon_numbers <- tibble(decon_number = 0:3)

plot_times <- dat %>%
    data_grid(time_hrs = seq_range(c(0, 2),
                               n = fineness)) %>%
    mutate(time_min = time_hrs * 60)


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
           treatment,
           control_index,
           treatment_end_time) %>% 
    inner_join(tidy_draws,
               by = "trial_unique_id") %>%
    group_by(.draw, material) %>%
    mutate(control_decay = first(
               decay_rate[trial_unique_id == control_index])) %>%
    ungroup() %>%
    filter(material %in% materials_to_include)

##################################
## plot mean estimated titers to
## check goodness of fit
##################################
titer_means_draws <- titer_means_chains %>%
    spread_draws(set_mean_log10_tcid50[replicate_set_id])


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
    mutate(time_min = time * 60,
           detectable = n_pos > 0) %>%
    ungroup() %>%
    filter(material %in% materials_to_include)

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
    ungroup()

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
    crossing(plot_times)

## adding one to convert to per mL from per mL/10
to_plot <- to_plot %>%
    mutate(log_predicted_titer = ifelse(
               time_hrs < treatment_end_time,
               intercept - decay_rate * time_hrs,
               intercept - decay_rate * treatment_end_time -
               control_decay * (time_hrs - treatment_end_time)),
           predicted_titer = 10^(1 + log_predicted_titer))

shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 25,
                         "TRUE" = 21)))

lines_steel <- to_plot %>%
    filter(material == "Steel")

lines_mask <- to_plot %>%
    filter(material == "N95 mask")

lines_dmem <- to_plot %>%
    filter(material == "DMEM")

## nudge steel points at the LOD slightly
## to avoid overplotting
steel_sub <- titer_means_draws %>%
    filter(material == "Steel") %>%
    mutate(set_mean_log10_tcid50 = ifelse(
               detectable,
               set_mean_log10_tcid50,
               LOD_log10 + 0.1),
           time_min = ifelse(
               detectable,
               time_min,
               time_min - 3)) %>%
    mutate(
        set_mean_log10_tcid50 =
            ifelse(
            (treatment %in% c("Ethanol", "Heat") &
             time_min == 0),
            set_mean_log10_tcid50 + 0.2,
            set_mean_log10_tcid50),
        time_min = ifelse(
        (treatment %in% c("Ethanol", "Heat") &
         time_min == 0),
        time_min - 1.5,
        time_min))

dmem_sub <- titer_means_draws %>%
    filter(material == "DMEM")

mask_sub <- titer_means_draws %>%
    filter(material == "N95 mask")

fit_panel <- lines_mask %>%
    ggplot(aes(x = time_min,
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
        mapping = aes(x = time_min,
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
        mapping = aes(x = time_min,
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
        mapping = aes(x = time_min,
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
        mapping = aes(x = time_min,
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
    scale_fill_manual(values = unlist(material_colors),
                      breaks = levels(dat$material)) +
    scale_fill_manual(values = unlist(material_colors),
                      breaks = levels(dat$material),
                      aesthetics = "point_fill") +
    scale_color_manual(values = unlist(material_colors),
                       breaks = levels(dat$material)) +
    shape_scale + 
    scale_y_continuous(trans = ytrans,
                       breaks = ybreaks,
                       label = yformat) +
    coord_cartesian(ylim = ylim,
                    xlim = c(0, 96)) +
    facet_wrap(vars(treatment), nrow = 1)
    
## styling: no facet labels because is background plot
fit_panel <- fit_panel +
    guides(fill = guide_legend(override.aes = list(alpha = 1,
                                                   stroke = 1,
                                                   point_size = 10,
                                                   shape = 21)),
           shape = "none") +
    theme_project(base_size = text_size) +
    theme(legend.key.size = unit(3, "lines"),
          legend.position = c(0.135, 0.8),
          legend.margin = margin(t = -10,
                                 l = 0,
                                 r = 5,
                                 b = 0),
          legend.title = element_blank(),
          legend.box.background = element_rect(color = "black",
                                               size = 2))+
    xlab("Time (min)") +
    ylab(titer_ylab)


cat('plotting mask integrity predictions...\n')

## sort facets by the posterior median integrity
cat('Plotting integrity...\n')

big_point_size <- 25 # size for median points
big_point_stroke <- 2

integrity_dat <- integrity_dat %>%
    mutate(below_limit = fit_factor < 200)

mask_shape_scale = scale_shape_manual(
    values = unlist(list("FALSE" = 24,
                         "TRUE" = 21)))

integrity_panel <- integrity_predictions %>%
    ggplot(aes(
        x = factor(decon_number),
        y = fit_factor_pred,
        point_fill = treatment,
        fill = treatment)) +
    geom_hline(aes(yintercept = 100),
               linetype = "longdash",
               size = 3,
               color = "red") +
    stat_pointinterval(
        .width = c(0.6827, 0.95),
        point_size = big_point_size,
        interval_size_range = c(2.5, 5),
        interval_color = interval_grey,
        stroke = big_point_stroke,
        shape = 21,
        color = "black") +
    geom_line(data = integrity_dat,
              mapping = aes(
                  y = fit_factor,
                  group = subject_id),
              color = "black",
              size = 2,
              alpha = 0.75,
              position = position_dodge(width = 0.3)) +
    geom_point(data = integrity_dat,
               mapping = aes(
                   y = fit_factor,
                   fill = treatment,
                   group = subject_id,
                   shape = below_limit),
               size = 8,
               stroke = 1,
               alpha = 0.8,
               position = position_dodge(width = 0.3)) +  
    ylab("Mask integrity (fit factor)") +
    xlab("Number of decontaminations") +
    scale_fill_manual(values = unlist(treatment_colors),
                      aesthetics = "point_fill") +
    mask_shape_scale +
    scale_fill_manual(values = unlist(treatment_colors)) +
    scale_color_manual(values = unlist(treatment_colors)) +
    facet_grid(cols = vars(treatment)) +
    coord_cartesian(ylim = c(25, 400)) +
    scale_x_discrete() +
    scale_y_continuous(expand = c(0, 0),
                       trans = "log2",
                       breaks = c(25, 50, 100, 200, 400)) +
    theme_project(base_size = text_size) +
    theme(legend.position = "none")

####################################
## create inactivation vs integrity
## plot
####################################


cat('extracting posterior draws...\n')

inactivation <- decay_chains %>%
    spread_draws(decay_rate[trial_unique_id]) %>%
    mutate(kill_rate_min = decay_rate * 60) %>%
    crossing(decon_numbers)

inact_to_plot <- dat %>%
    distinct(trial_unique_id, .keep_all = TRUE) %>%
    inner_join(inactivation,
               by = "trial_unique_id") %>%
    filter(material == "N95 mask")

to_plot <- inact_to_plot %>%
    inner_join(integrity_predictions,
               by = c("treatment", "decon_number", ".draw")) %>%
    group_by(treatment) %>%
    mutate(kill_rate_med = median(kill_rate_min)) %>%
    ungroup() %>%
    group_by(treatment, decon_number) %>%
    mutate(fit_factor_med = median(fit_factor_pred)) %>%
    ungroup()


quant_lower = 0.025
quant_upper = 0.975
quant_mid_lower = (1 - 0.6827) / 2
quant_mid_upper = 1 - quant_mid_lower
ebar_size <- 2
ci_to_plot <- to_plot %>% group_by(treatment, decon_number) %>%
    summarise(med_ff = median(fit_factor_pred),
              med_kr = median(kill_rate_min),
              lb_ff = quantile(fit_factor_pred, quant_lower),
              uq_ff = quantile(fit_factor_pred, quant_mid_upper),
              lq_ff = quantile(fit_factor_pred, quant_mid_lower),
              ub_ff = quantile(fit_factor_pred, quant_upper),
              lb_kr = quantile(kill_rate_min, quant_lower),
              ub_kr = quantile(kill_rate_min, quant_upper)) %>%
    ungroup()

text_justs <- tibble(
    treatment = factor(
        c("Control",
          "Ethanol",
          "Heat",
          "VHP",
          "UV"),
        levels = ordering$treatment),
    hjusts = c(1.3, 1.3, 1.5, 1.3, 1.3),
    vjusts = c(2.6, -1.35, -1.5, 2.6, 2.9))
           

text_data <- ci_to_plot %>%
    filter(decon_number == 1) %>%
    distinct(treatment, .keep_all=TRUE) %>%
    inner_join(text_justs,
               by = "treatment")

cat('Plotting integrity vs kill rate...\n')

int_int_panel <- to_plot %>%
    filter(decon_number > 0) %>%
    arrange(desc(treatment)) %>%
    ggplot(aes(
        x = fit_factor_pred,
        y = kill_rate_med,
        group = treatment,
        fill = treatment)) +
    geom_vline(aes(xintercept = 100),
               color = "red",
               linetype = "longdash",
               size = 3) + 
    stat_pointintervalh(
        .width = 0.6827,
        interval_size_range = c(0, 6),
        interval_color = interval_grey,
        color = "black",
        point_size = big_point_size,
        stroke = big_point_stroke, 
        shape = 21) +
    stat_pointinterval(
        mapping = aes(x = fit_factor_med,
                      y = kill_rate_min),
        .width = 0.6827,
        interval_size_range = c(0, 6),
        interval_color = interval_grey,
        color = "black",
        point_size = big_point_size,
        stroke = big_point_stroke, 
        shape = 21) +
    geom_text(
        aes(x = med_ff,
            y = med_kr,
            color = treatment,
            label = treatment,
            hjust = hjusts,
            vjust = vjusts),
        data = text_data,
        size = 8) +
    xlab("Integrity after decontamination (fit factor)") +
    ylab("Kill rate (per virus per min)") +
    scale_x_continuous(trans = "log2",
                       breaks = c(50, 100, 200),
                       expand = c(0, 0)) +
    scale_y_continuous(trans = "log2",
                       breaks = c(10, 100, 1000, 10000)) +
        coord_cartesian(xlim = c(50, 400),
                        ylim = c(5, 10000)) +
    theme_project(base_size = text_size) +
    facet_wrap(vars(factor(decon_number)), nrow = 1,
               strip.position = "bottom") +
    scale_fill_manual(values = unlist(treatment_colors)) +
    scale_color_manual(values = unlist(treatment_colors)) +
    theme(legend.position = "none")



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

## no facet labels for lower panel because 
## can just read the upper ones
unlabeled_panel_theme <- theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = text_size))

cat('making full figure...\n')
full_fig <- plot_grid(
    fit_panel + labeled_panel_theme,
    integrity_panel + unlabeled_panel_theme,
    int_int_panel + int_int_theme,
    label_size = text_size,
    align = "v",
    axis = "lr",
    rel_heights = c(1, 1, 1.3),
    ncol = 1)

## save the plot to outpath
cat('saving figure to ', outpath, '...\n')
save_plot(outpath,
          full_fig,
          base_height = 22,
          base_asp = 1.2)
warnings()
