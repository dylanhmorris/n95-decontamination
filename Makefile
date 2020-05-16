#####################################
# name: Makefile
# author: Dylan Morris <dhmorris@princeton.edu>
#
# Makefile to generate analyses
# for Fischer et al study of N95
# mask decontamination for SARS-CoV-2
####################################

#####################################
# Directory structure
####################################

default: all

SRC := src
OUTPUT := out
DATA := dat
RAW := $(DATA)/raw
CLEANED := $(DATA)/cleaned
PARAMS := $(SRC)/parameters

MCMC_CHAINS := $(OUTPUT)/mcmc_chains
FIGURE_DIR = $(OUTPUT)/figures
TABLE_DIR = $(OUTPUT)/tables

#####################################
# File extensions and the like
####################################
CHAINS_SUFFIX = _chains.Rds
PRIOR_CHECK_NAME = _prior_check

#####################################
# Expected bash settings
#
# Check these vs your local
# machine setup if you are having
# difficulty reproducing the
# analysis
#####################################

CXX := g++-8
MKDIR := @mkdir -p
RM := rm -rf

R_OPTIONS = --vanilla
R_COMMAND := Rscript $(R_OPTIONS)

# R creates a blank Rplots.pdf when run
# from the command line. This removes it.
FIG_CLEANUP = @$(RM) Rplots.pdf

#####################################
# Installation / dependencies
#
# Rules for prepping analysis
#####################################

.PHONY: depend

depend:
	$(R_COMMAND) $(SRC)/install_needed_packages.R

#####################################
# data locations
#####################################
DECON_DATAFILE = decon_data.csv
MASK_DATAFILE = mask_fit_data.csv

RAW_DECON_DATA = $(RAW)/$(DECON_DATAFILE)
CLEANED_DECON_DATA = $(CLEANED)/$(DECON_DATAFILE)

RAW_MASK_DATA = $(RAW)/$(MASK_DATAFILE)
CLEANED_MASK_DATA = $(CLEANED)/$(MASK_DATAFILE)

CLEANED_DATA = $(CLEANED_DECON_DATA) $(CLEANED_MASK_DATA) 

## use corrected titer data
DECON_DATA = $(CLEANED_DECON_DATA)
MASK_DATA = $(CLEANED_MASK_DATA)

#####################################
# code locations
#####################################
CLEANING_SCRIPT = $(SRC)/clean_data.R
MASK_CLEANING_SCRIPT = $(SRC)/clean_mask_data.R
FITTING_SCRIPT = $(SRC)/fit_stan_model.R
MASK_FITTING_SCRIPT = $(SRC)/fit_mask_model.R
MASK_PRIOR_CHECK_SCRIPT = $(SRC)/fit_mask_model_prior_check.R
DIAGNOSTIC_SCRIPT = $(SRC)/chain_diagnostics.R

SPLINE_MODEL_SRC = decon_decay_spline.stan
MASK_MODEL_SRC = mask_fit_nonlinear.stan
INFER_TITER_MODEL_SRC = well_titer_estimates.stan
TITER_MEAN_MODEL_SRC = well_titer_means.stan
#####################################
# parameter locations
#####################################
PLOT_PARAMS = $(PARAMS)/plotting_style.R

## paths to hyperparameters
MASK_HYPERS = $(PARAMS)/mask_hyperparams.R
DECON_HYPERS = $(PARAMS)/decon_hyperparams.R
INFER_TITER_HYPERS = $(PARAMS)/well_titer_estimates_hyperparams.R
INFER_MEAN_HYPERS = $(PARAMS)/well_titer_means_hyperparams.R
#####################################
# model names and output locations
#####################################

WELL_MODEL_NAME = well

DECON_MODEL_NAME = decon
DECON_CHAINS =  $(MCMC_CHAINS)/$(DECON_MODEL_NAME)$(CHAINS_SUFFIX)

DECON_TITER_NAME = inferred_decon_titers
DECON_TITER_CHAINS =  $(MCMC_CHAINS)/$(DECON_TITER_NAME)$(CHAINS_SUFFIX)

DECON_MEANS_NAME = inferred_decon_means
DECON_MEANS_CHAINS =  $(MCMC_CHAINS)/$(DECON_MEANS_NAME)$(CHAINS_SUFFIX)

MASK_MODEL_NAME = mask
MASK_CHAINS = $(MCMC_CHAINS)/$(MASK_MODEL_NAME)$(CHAINS_SUFFIX)
MASK_PRIOR_CHECK_CHAINS =  $(MCMC_CHAINS)/$(MASK_MODEL_NAME)$(PRIOR_CHECK_NAME)$(CHAINS_SUFFIX)

MODELS = $(WELL_MODEL_NAME) $(DECON_MODEL_NAME) $(MASK_MODEL_NAME)

CHAIN_PATHS = $(DECON_CHAINS) $(MASK_CHAINS) $(DECON_MEANS_CHAINS) $(MASK_PRIOR_CHECK_CHAINS)

CHAIN_DIAGNOSTICS = $(OUTPUT)/chain_diagnostics.csv


DECON_FIGURES = figure_decon_main.pdf figure_decon_controls.pdf figure_mask_posterior_check.pdf figure_mask_prior_check.pdf figure_decon_masks.pdf figure_decon_dmem.pdf 

FIGURES = $(DECON_FIGURES)

FIGURE_PATHS = $(addprefix $(FIGURE_DIR)/, $(FIGURES))


TABLES = table_decon.docx
TABLE_PATHS = $(addprefix $(TABLE_DIR)/, $(TABLES))

#####################################
# Rules
#
# definition of dependency
# tree and specification of
# rules for doing stuff
#####################################

##########################
# rules for data cleaning
##########################

$(CLEANED)/%.csv: $(RAW)/%.csv $(CLEANING_SCRIPT)
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $(CLEANING_SCRIPT) $< $@

$(CLEANED_MASK_DATA): $(MASK_CLEANING_SCRIPT) $(RAW_MASK_DATA)
	$(MKDIR) $(CLEANED)
	$(R_COMMAND) $^ $@


#####################################
# rules for model fitting and post-processing
#####################################


$(DECON_TITER_CHAINS): $(FITTING_SCRIPT) $(SRC)/$(INFER_TITER_MODEL_SRC) $(DECON_DATA) $(INFER_TITER_HYPERS)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(DECON_MEANS_CHAINS): $(FITTING_SCRIPT) $(SRC)/$(TITER_MEAN_MODEL_SRC) $(DECON_DATA) $(INFER_MEAN_HYPERS)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(DECON_CHAINS): $(FITTING_SCRIPT) $(SRC)/$(SPLINE_MODEL_SRC) $(DECON_DATA) $(DECON_HYPERS)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(MASK_CHAINS): $(MASK_FITTING_SCRIPT) $(SRC)/$(MASK_MODEL_SRC) $(MASK_DATA) $(MASK_HYPERS)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(MASK_PRIOR_CHECK_CHAINS): $(MASK_PRIOR_CHECK_SCRIPT) $(SRC)/$(MASK_MODEL_SRC) $(MASK_DATA) $(MASK_HYPERS)
	$(MKDIR) $(MCMC_CHAINS)
	$(R_COMMAND) $^ $@

$(CHAIN_DIAGNOSTICS): $(DIAGNOSTIC_SCRIPT) $(CHAIN_PATHS)
	$(MKDIR) $(OUTPUT)
	$(R_COMMAND) $^ $@


#####################################
# rules for table generation
#####################################

$(TABLE_DIR)/table_decon.docx: $(SRC)/table_decon.R $(DECON_DATA) $(DECON_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(TABLE_DIR)
	$(R_COMMAND) $^ $@


#####################################
# rules for figure generation
#####################################

$(FIGURE_DIR)/figure_decon_main.pdf: $(SRC)/figure_decon_main.R $(DECON_DATA) $(MASK_DATA) $(DECON_CHAINS) $(DECON_MEANS_CHAINS) $(MASK_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_decon_controls.pdf: $(SRC)/figure_decon_controls.R $(DECON_DATA) $(MASK_DATA) $(DECON_CHAINS) $(DECON_MEANS_CHAINS) $(MASK_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_decon_dmem.pdf: $(SRC)/figure_decon_dmem.R $(DECON_DATA) $(DECON_CHAINS) $(DECON_MEANS_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_decon_masks.pdf: $(SRC)/figure_decon_masks.R $(MASK_DATA) $(MASK_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_mask_posterior_check.pdf: $(SRC)/figure_mask_pp_check.R $(MASK_DATA) $(MASK_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)

$(FIGURE_DIR)/figure_mask_prior_check.pdf: $(SRC)/figure_mask_pp_check.R $(MASK_DATA) $(MASK_PRIOR_CHECK_CHAINS) $(PLOT_PARAMS)
	$(MKDIR) $(FIGURE_DIR)
	$(R_COMMAND) $^ $@
	$(FIG_CLEANUP)



#####################################
# convenience rules for making
# various quantities
#####################################
.PHONY: data
data: $(CLEANED_DATA)

.PHONY: chains
chains: $(CHAIN_PATHS)

.PHONY: diagnostics
diagnostics: $(CHAIN_DIAGNOSTICS)

.PHONY: figures
figures: $(FIGURE_PATHS) 

.PHONY: tables
tables: $(TABLE_PATHS)

.PHONY: echo_figures echo_chains

echo_figures:
	echo $(FIGURE_PATHS)
echo_chains:
	echo $(CHAIN_PATHS)

## remove emacs tempfiles, etc.
.PHONY deltemp:
	$(RM) $(SRC)/*~*
	$(RM) $(SRC)/*#*
	$(RM) $(PARAMS)/*~*
	$(RM) $(PARAMS)/*#*

.PHONY: clean
clean: deltemp
	$(RM) $(OUTPUT)
	$(RM) $(CLEANED)

all: depend data chains diagnostics figures tables
