# `fdcoexist` – a model of community coexistence with trits

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip) [![CircleCI](https://circleci.com/gh/Rekyt/fdcoexist.svg?style=svg&circle-token=fbf1bc61e0b42230c72d8e000cc728ec3843ebcf)](https://circleci.com/gh/Rekyt/fdcoexist) [![codecov](https://codecov.io/gh/Rekyt/fdcoexist/branch/master/graph/badge.svg?token=8IIrOJuF0u)](https://codecov.io/gh/Rekyt/fdcoexist)

The goal of `fdcoexist` is to provide for Denelle et al. *in prep.*.

## Authors

* Pierre Denelle
* Matthias Grenié
* Cyrille Violle
* François Munoz
* Caroline M. Tucker

## Content

`fdcoexist` is a population model of coexisting plants across patches that lets the user explore the consequences of trait contributing to different processes on the CWM <-> environment relationship as well as the relation between various estimators of trait-environment relationship.


## R package

The R functions contains the bare minimum in the package to run simulations while most of the actual analysis code is in the `meta_simulations.Rmd` vignette.
The vignettes contains most of the analysis code for the moment:
* `vignettes/equation_number_individuals_equilibrium_noncompetition.Rmd` contains the equation that demonstrates the number of individuals of each species expected in the absence of limiting similarity as well as hierarchical competition,
* `vignettes/equation_only.Rmd` contains the equation of the model and explanation of its development,
* `vignettes/meta_simulations.Rmd` contains most of the analysis code to run analysis and plot corresponding figures
* `vignettes/solving_A.Rmd` tries different approach to explore good values of A


Other interesting code is in `inst/test.R` and consists of just a scratchpad of different figures and analyses.
