#!/usr/bin/env Rscript

###################################
## shared styling for all plots
##
####################################
script_packages <- c("ggplot2",    # for theme_set()
                     "scales",     # for scientific_format()
                     "readr",      # for read_csv()
                     "cowplot",    # for background_grid()
                     "extrafont")   # for nice fonts 


## load in packages without messages
for (package in script_packages){
    suppressPackageStartupMessages(
        library(package,
                character.only=TRUE))
}


## load in fonts

## styling for text
ncov_name <- "SARS-CoV-2"
sars_name <- "SARS-CoV-1"
mers_name <- "MERS-CoV"
flu_name <- "Influenza"
    
virus_order <- c(ncov_name,
                 sars_name,
                 mers_name,
                 flu_name)

strain_order <- c(
    ncov_name,
    sars_name,
    mers_name,
    "MERS-CoV_EMC/12",
    "MERS-CoV_U/14",
    "MERS-CoV_KSA/15",
    "MERS-CoV_SK/15",
    "MERS-CoV_KSA/18",
    "MERS-CoV_C/KSA/13",
    "MERS-CoV_C/E/13",
    "MERS-CoV_C/BF/15")

## color palette
interval_grey <- "#525252"
virus_colors <- list(
    "SARS-CoV-2" = "#cc4747",
    "SARS-CoV-1" = "#7570b3",
    "MERS-CoV" = "#518591")

material_colors <- list(
    "N95 mask" = "#554fa8",
    "Plastic" = "#c4ae8f",
    "Aerosols" = "#54e3ff",
    "Cardboard" = "#615242",
    "Copper" = "#da8a67",
    "Steel" = "#adbab8",
    "DMEM" = "#609c70",
    "Uncovered plate oven" = "#609c70",
    "Closed vial heat block" = "#609c70",
    "Closed vial oven" =  "#609c70",
    "Covered plate oven" =  "#609c70")

treatment_colors <- list(
    "Control" = "grey",
    "Ethanol" = "#d2db1f",
    "Heat" = "#ff8b0f",
    "VHP" = "#0fc7ff",
    "UV" = "#a503ab")


pp_color <- "#5564c2"

sci10_formatter <- function(x){
    return (ifelse(x == 0,
                   "0",
                   parse(text = gsub(
                             "[+]", "",
                             gsub("e", " %*% 10^",
                                  scientific_format()(x))))
                   ))
}


## function to read in cleaned data and format
## it as a tibble that can be used for plotting
read_data_for_plotting <- function(data_path) {
   dat <- read_csv(data_path,
                   col_types = cols())

   ## add unlogged titers for plotting
   if ("log10_titer" %in% names(dat) ) {
       dat$titer = 10^dat$log10_titer
   }
   
   ## make sure the viruses/strains display
   ## in our order or interest (set below)
   dat$virus <- factor(dat$virus,
                       virus_order)
   
   if ("strain" %in% names(dat) ) {
       dat$strain <- factor(dat$strain,
                            strain_order)
   }
   
   return (dat)
}


vector_formats = c(".pdf")
save_plot <- function(outpath, fig, ...){
    cowplot::save_plot(outpath,
                       fig,
                       ...)
    if(any(endsWith(outpath, vector_formats)))
        embed_fonts(file = outpath)
}

theme_project <- function(base_size = 30,
                          ...){
    x_axis_margin <- margin(t = base_size / 2)
    y_axis_margin <- margin(r = base_size / 2)

    theme_classic(
        base_size = base_size,
        ...) +
        background_grid(major = "xy",
                        minor = "none",
                        size.major = 0.5) +
        theme(axis.title.x = element_text(
                  margin = x_axis_margin),
              axis.title.y = element_text(
                  margin = y_axis_margin),
              strip.background = element_blank(),
              panel.border = element_rect(color = "black",
                                          fill = NA,
                                          size = 2),
              text = element_text(family = "CM Sans"))

}

theme_set(theme_project())

## styling for figures

panel_fig_height <- 20
panel_fig_aspect <- 2.5

## styling for lines
fineness <- 500  ## fineness of line plotting
n_lines <- 50    ## number of random posterior lines to plot
line_alpha = 5 / n_lines  ## transparency of overplotted lines
line_size <- 0.1 * panel_fig_height
lowest_log_titer <- 0.5 ## sets the lower y limit
detection_linestyle <- "dashed" ## sets the style of the detection line
n_breaks <- 4

## styling for points
pointfill <- "#56c9ff"
pointborder <- 'black'
pointsize <- 0.25 * panel_fig_height
pointstroke <- 0.5
pointalpha <- 0.5
jitwid <- 0.25
jith <- 0.0


cat("Loaded plotting style successfully!\n")
