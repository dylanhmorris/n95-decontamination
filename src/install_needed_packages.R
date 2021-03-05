#!/usr/bin/env Rscript

#####################################
## name: install_needed_packages.R
## author: Dylan Morris <dhmorris@princeton.edu>
##
## installs needed packages for
## reproducing SARS-CoV-2 decay
## study
##
####################################

install_if_absent <- function(package_name){
    if (!suppressPackageStartupMessages(
             require(package_name, character.only=TRUE))){
      install.packages(pkgs=package_name,
                       repos="http://cloud.r-project.org")
  }
  else
      cat(sprintf("Package %s already installed\n", package_name))
}

needed_packages <- c(
    'rstan',      # stan interface
    'parallel',   # parallelized MCMC
    'readr',      # csv read-in
    'magrittr',   # for pipe operator %>%
    'dplyr',      # for filter()
    'tidybayes',  # for for spread_draws(), etc.
    'ggplot2',    # for plotting
    'modelr',     # for data_grid()
    'tidyr',      # for crossing(), drop_na()
    'forcats',    # for fct_rev()
    'cowplot',    # publication ready ggplot
    'bayesplot',  # for rhat, neff_ratio stats
    'huxtable',   # for programmatic results output
    'officer',    # for outputing tables to docx
    'flextable',  # for table output
    'stats',      # for stable numerics
    'Rcpp',       # for *fast* numerics
    'extrafont',  # for nice plot fonts
    'viridis')    # for magma color palette

args <- commandArgs(trailingOnly=TRUE)

for (package in needed_packages)
    install_if_absent(package)


r <- getOption("repos")
r["CRAN"] <- "http://cloud.r-project.org"
options(repos = r)
library(extrafont)
font_install('fontcm', prompt = FALSE)
