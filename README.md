# Effectiveness of N95 Respirator Decontamination and Reuse against SARS-CoV-2 Virus
Robert J. Fischer(1\*), Dylan H. Morris(2), Neeltje van Doremalen(1), Shanda Sarchette(1), M. Jeremiah Matson(1), Claude Kwe Yinda(1), Stephanie N. Seifert(1), Amandine Gamble(3), Brandi N. Williamson(1), Seth D. Judson(4), Emmie de Wit1, James O. Lloyd-Smith3, Vincent J. Munster1

1. Laboratory of Virology, Division of Intramural Research, National Institute of Allergy and Infectious Diseases, National Institutes of Health, Hamilton, MT, USA
2. Dept. of Ecology \& Evolutionary Biology, Princeton University, Princeton, NJ, USA
3. Dept. of Ecology \& Evolutionary Biology, University of California, Los Angeles, Los Angeles, CA, USA
4. Dept. of Medecine, University of Washington, Seattle, WA, USA.

## Repository information
This repository accompanies the article "Effectiveness of N95 Respirator Decontamination and Reuse against SARS-CoV-2 Virus" (R.J. Fischer et al.). It provides code and data for reproducing all data analysis from the paper and recreating all associated display figures.

## Changelog
This repository was last updated March 2021 to add bibliographic information for the published article and to ensure that plotting code exactly matches published figure captions. This results in a cosmetic change to [Appendix Figure 1](out/figures/figure_decon_controls.pdf) relative to the published version. There are no other changes, and no conclusions change.

## License and citation information
If you use the code or data provided here, please make sure to do so in light of the project [license](LICENSE.txt) and please cite our work as below:

Fischer RJ, Morris DH, van Doremalen N, Sarchette S, Matson M, Bushmaker T, et al. Effectiveness of N95 Respirator Decontamination and Reuse against SARS-CoV-2 Virus. Emerg Infect Dis. 2020;26(9):2253-2255. https://dx.doi.org/10.3201/eid2609.201524

Bibtex record:
```
@electronic{fischer2020n95,
  title={Effectiveness of N95 respirator decontamination and reuse against SARS-CoV-2 virus},
  author={Fischer, Robert J and 
  Morris, Dylan H and 
  van Doremalen, Neeltje and 
  Sarchette, Shanda and 
  Matson, M Jeremiah and 
  Bushmaker, Trenton and 
  Yinda, Claude Kwe and 
  Seifert, Stephanie N and 
  Gamble, Amandine and 
  Williamson, Brandi N and
  Judson, Seth D and
  de Wit, Emmie and 
  Lloyd-Smith, James O and 
  Munster, Vincent J},
  journal={Emerging infectious diseases},
  volume={26},
  number={9},
  pages={2253},
  year={2020},
  DOI = {10.3201/eid2609.201524}
  URL = {https://github.com/dylanhmorris/n95-decontamination}
}
```

## Article abstract 
The coronavirus pandemic has created worldwide shortages of N95 respirators. We analyzed 4 decontamination methods for effectiveness in deactivating severe acute respiratory syndrome coronavirus 2 virus and effect on respirator function. Our results indicate that N95 respirators can be decontaminated and reused, but the integrity of respirator fit and seal must be maintained.


## Directories
- ``src``: all code, including data preprocessing, Bayesian model definition and fitting, and results post-processing and figure generation:
    - ``src/parameters``: parameters specification for models and style specification for plots
- ``dat``: data files in comma-separated values (``.csv``) formats
    - ``dat/raw``: raw data files (semicolon-separated)
    - ``dat/cleaned``: data files processed and prepared for model fitting
    - ``dat/fonts``: freely available fonts for figures
- ``out``: output files
    - ``out/mcmc_chains``: Markov Chain Monte Carlo (MCMC) output, as serialized R data (``.Rds``) files. 
    - ``out/figures``: figures generated from results
    - ``out/tables``: tables generated from results
    - ``out/chain_diagnostics.csv``: diagnostic tests for MCMC convergence.

## Reproducing analysis

A guide to reproducing the analysis from the paper follows.

### Getting the code
First download this repository. The recommended way is to ``git clone`` it from the command line:

    git clone https://github.com/dylanhmorris/n95-decontamination.git

Downloading it manually via Github's download button or from OSF should also work.

### Dependency installation
The analysis can be auto-run from the project ``Makefile``, but you may need to install some external dependencies first. See the **Dependency installation guide** below for a complete walkthrough. In the first instance, you'll need a working installation of the statistical programming language R, a working C++ compiler, and a working installation of Gnu Make or similar. A few external R packages can then be installed from the command line by typing.

    make depend

from within the project directory.

### Running the analysis

The simplest approach is simply to type ``make`` at the command line, which should produce a full set of figures and MCMC output (saved as R Dataset ``.Rds`` files in the ``out/mcmc-chains/`` directory as ``<model_name>_chains.Rds``). These can be loaded in any working R installation, as long as the package ``rstan`` is also installed.

If you want to do things piecewise, typing ``make <filename>`` for any of the files listed in the ``dat/cleaned`` or ``out`` directories below should run the steps needed to produce that file.

Some shortcuts are available:

- ``make data`` produces cleaned data files.
- ``make chains`` produces all MCMC output
- ``make diagnostics`` extracts MCMC diagnostic statistics
- ``make figures`` produces all figures
- ``make tables`` produces all tables
- ``make clean`` removes all generated files, leaving only source code (though it does not uninstall packages)

### Examining code

Examining the raw Stan code is the place to start to understand how models have been specified. But note that parameters for the prior distributions are set at runtime rather than hard-coded into the ``.stan`` files, so that recompilation is not required when parameter choices are changed (this makes it easier to try the models using different priors, for sensitivity analysis).

Prior parameter choices are specified in files found in the directory, ``src/parameters``.

## Project structure when complete

Once the full analysis has been run, you should be able to find a full set of figures in ``out/figures`` and a table of regression results in ``out/tables``.

## Dependency installation guide
You will need a working R installation with the command line interpreter ``Rscript`` (macOS and Linux) or ``Rscript.exe`` (Windows). On mac and Linux, you can check that you have an accessible ``Rscript`` by typing ``which Rscript``at the command line and seeing if one is found.

If you do not have an R installation, you can install it from [the R project website](https://www.r-project.org/) or from the command line using a package manager such as [Homebrew](https://brew.sh/) on macOS or ``apt-get`` on Linux. macOS users may also need to install the macOS "command line tools" by typing ``xcode-select --install`` at a command prompt.

Once R is installed, you can automatically install all other dependencies (including the Hamiltonian Monte Carlo software Stan and its R interface rstan) on most systems using ``make``. In the top level project directory, type the following at the command line:

    make depend

Alternatively, you can run the script ``src/install_needed_packages.R`` manually. 

Note that installing Stan and RStan can be time-consuming, Stan is a large program that must be compiled from source. Some of the packages in the very valuable [tidyverse](https://www.tidyverse.org/) may also take some time to install.
