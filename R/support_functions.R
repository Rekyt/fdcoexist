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
    given_composition$compo[,,time] %>%
        as.data.frame() %>%
        rownames_to_column("patch") %>%  # Add patch names as column
        # Transform to tidy data.frame (define name of columns and column name
        # for values)
        gather("species", "N", -patch) %>%
        # Add patch number as a new column and compute distance to optimal value
        mutate(patch_optim = gsub("\\w+(\\d+)", "\\1", patch) %>%
                   as.numeric())
}
