# Graphical function to plot environmental response curves of species

#' Plot environmental response curves of species in all the patches
#'
#' Automatically assigns one color per species.
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param plot a boolean determining whether to plot or not the response curves
#'
#' @import ggplot2
#' @export
#'
#'
r_env <- function(simul, n_patches, plot = TRUE){
    env_growth <- data.frame(simul[["rmatrix"]])
    env_growth$env <- 1:25
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    if(plot == TRUE){
        env_plot <- ggplot(env_growth, aes(env, r_env)) +
            geom_line(aes(color = as.factor(sp))) +
            geom_label(aes(label = ifelse(r_env == k, sp, NA)), size = 3) +
            labs(x = "Patch", y = expression("r"["env"]),
                 title = "Environmental response curve of species") +
            theme_classic() +
            guides(color = FALSE) +
            theme(panel.border = element_rect(fill = NA))

        return(list(env_growth, env_plot))
    }else {
        return(env_growth)
    }
}
