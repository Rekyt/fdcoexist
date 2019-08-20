# Script with support functions

#' Simpler data.frame extracted from population
#'
#' @param given_composition a result of a call to `multigen()`
#' @param time              the timestep at which population abundances should
#'                          be extracted
#'
#' @return a data.frame with the following columns
#'         * `patch` the label of the patch
#'         * `species` the name of the species
#'         * `N` the local abundance of species in given patch
#'         * `patch_optim` environmental value of given patch
#' @export
landscape_df = function(given_composition, time = 150) {

    compo_df = as.data.frame(given_composition$compo[,,time])

    # Add patch names as column
    compo_df$patch = rownames(compo_df)


        # Transform to tidy data.frame (define name of columns and column name
        # for values)
    tidy_compo = tidyr::gather(compo_df, "species", "N", -"patch")
        # Add patch number as a new column and compute distance to optimal value
    tidy_compo$patch_optim = gsub("\\w+(\\d+)", "\\1", tidy_compo$patch)

    return(tidy_compo)
}

wtd_mean = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    weighted.mean(x, w)
}

wtd_var = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.var(x, w)
}

#' Weighted Skewness with na.rm
#'
#' Compute weighted Skewness using [`Weighted.Desc.Stat::w.skewness`] but add an
#' option to remove `NA` values
#'
#' @param x the values to be weighted
#' @param w the weighted
#' @param na.rm should NA be removed?
#'
#' @export
wtd_skewness = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.skewness(x, w)
}

#' Weighted Kurtosis with na.rm
#'
#' Compute weighted kurtosis using [`Weighted.Desc.Stat::w.kurtosis`] but add an
#' option to remove `NA` values
#'
#' @param x the values to be weighted
#' @param w the weighted
#' @param na.rm should NA be removed?
#'
#' @export
wtd_kurtosis = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.kurtosis(x, w)
}



#' Extract Contribution from Scenario Name
#'
#' All scenarios names are formatted the same way 'RXAYHZ' with X, Y and Z being
#' integers between 0 and 100
#' @param scenario a character of type `"R0A0H0"`
#'
#' @return a list with three elements:
#' * `R` an integer giving the contribution of trait 2 to growth
#' * `A` an integer giving the contribution of trait 2 to competition (limiting
#'       similarity)
#' * `H` an integer giving the contribution of trait 2 to hierarchical
#'       competition
#' @export
extract_trait_contrib_from_scenar = function(scenario) {
    R_scenar <- as.integer(sub("R(\\d+).*", "\\1", scenario, perl = TRUE))
    A_scenar <- as.integer(sub(".*A(\\d+).*", "\\1", scenario, perl = TRUE))
    H_scenar <- as.integer(sub(".*H(\\d+)", "\\1", scenario, perl = TRUE))

    list(R = R_scenar,
         A = A_scenar,
         H = H_scenar)
}
