# Graphical function to plot environmental response curves of species

#' Plot environmental response curves of species in all the patches
#'
#' Automatically assigns one color per species.
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param sp integer for the number of species to consider
#' @param plot a boolean determining whether to plot or not the response curves
#'
#' @import ggplot2
#' @export
#'
#'
r_env <- function(simul, n_patches, sp, plot = TRUE){
    sp <- 1:sp

    env_growth <- data.frame(simul[["rmatrix"]][, sp])
    env_growth$env <- 1:n_patches
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        rename(max_r_env = env) %>%
        as.data.frame()

    env_growth <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")],
                            by = "sp")

    if(plot == TRUE){
        env_plot <- ggplot(env_growth, aes(env, r_env)) +
            geom_line(aes(color = as.factor(sp))) +
            geom_label(aes(label = ifelse(env == max_r_env, sp, NA)),
                       size = 3) +
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
