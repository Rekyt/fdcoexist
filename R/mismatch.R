# Graphical function to growth of species due to hierarchical competition

#' Plot mismatch per species between environmental optimum and max. abundance
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param time an integer indicating the number of time steps
#' @param plot a boolean determining whether to plot or not the mismatch
#'
#' @import ggplot2
#' @export
#'
#'
mismatch <- function(simul, n_patches = 25, time = 50, plot = TRUE){
    # data.frame with environmental growth per species
    env_growth <- fdcoexist::r_env(simul, n_patches = n_patches, plot = FALSE)

    # Extract environment where species has its maximal R
    max_r <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        as.data.frame() %>%
        rename(max_r_env = env) %>%
        select(sp, max_r_env)

    # Remove species absent from any patch
    dat <- simul$compo[, , time]
    dat <- dat[, colSums(dat) > 0]

    # data.frame with observed maximal abundances
    max_ab <- apply(dat, 2, which.max)
    max_ab <- data.frame(
        sp =  gsub("species", "sp", as.character(names(max_ab))),
        max_env = as.integer(max_ab))

    # Merge environmental growth with maximal abundances
    mismatch <- left_join(max_r, max_ab, by = "sp") %>%
        mutate(mismatch = max_r_env - max_env) %>%
        mutate(sp_lab = as.integer(gsub("sp", "", sp)))

    if(plot == TRUE){
        return(list(
            mismatch,
            ggplot(mismatch, aes(mismatch, sp_lab)) +
                geom_vline(xintercept = 0,linetype = 1, color = "grey",
                           size = 1) +
                geom_point() +
                labs(x = "Mismatch", y = "Species",
                     title = "Optimal Environment vs Maximal abundance") +
                theme_classic() +
                theme(panel.border = element_rect(fill = NA))))
    } else{
        return(mismatch)
    }
}
