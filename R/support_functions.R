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
