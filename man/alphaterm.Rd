% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multigen_function.R
\name{alphaterm}
\alias{alphaterm}
\title{Function definition for deterministic model run with global dispersal
Compute alpha term in Beverton-Holt function}
\usage{
alphaterm(distance, Nts)
}
\arguments{
\item{distance}{competiton matrix of species}

\item{Nts}{vector of abundances of species at time t}
}
\description{
From a competition matrix (for the moment the distance between species
traits) and a vector of abundances by species, return the alpha term in the
Beverton-Holt equation. The order of species between the two should be the
same as no checks are done.
}