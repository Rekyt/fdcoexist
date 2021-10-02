# `fdcoexist` – Multi-species Trait-based Coexistence Model in discrete time

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/Rekyt/fdcoexist/workflows/R-CMD-check/badge.svg)](https://github.com/Rekyt/fdcoexist/actions)
[![DOI](https://zenodo.org/badge/341563824.svg)](https://zenodo.org/badge/latestdoi/341563824)
[![codecov](https://codecov.io/gh/Rekyt/fdcoexist/branch/master/graph/badge.svg?token=8IIrOJuF0u)](https://codecov.io/gh/Rekyt/fdcoexist)
<!-- badges: end -->

`fdcoexist` is a population model of coexisting plants across patches that lets the user explore the consequences of trait contributing to different processes on the CWM <-> environment relationship as well as the relation between various estimators of trait-environment relationship.

## Reference

> How competition shape trait-environment relationships? Insights from a population-dynamic model. Denelle\* P., Grenié\* M.*, Violle C., Munoz F., Tucker C.M., _in prep_.

## Running the analyses

To reproduce the figures shown in the manuscript copy this repository and with the `fdcoexist` directory selected as working directory run the code available in `inst/00-all_paper_figures.R`.

The `R/` folder contains all the functions to run the model while the `inst/` folder contains the code that actually run it.

