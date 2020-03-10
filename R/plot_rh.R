# Graphical function to growth of species due to hierarchical competition

#' Plot hierarchical growth of species in all the patches
#'
#' Automatically assigns one color per species.
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param time an integer indicating the number of time steps
#'
#' @import ggplot2
#' @export
#'
#'
plot_rh <- function(simul, n_patches, time){

    # Plot hierarchical coefficients in a given patch over time
    rh_data <- list()
    for(i in 1:(time-1)){
        rh_data[[i]] <- simul[[time-1]][n_patches, ]
    }
    # Convert arrays into data frames
    rh_data <- lapply(rh_data,
                      function(x) data.frame("sp" = names(x),
                                             "rh" = as.numeric(x))
    )
    # Add column to identify time
    rh_data <- mapply(cbind, rh_data, "time" = 2:time, SIMPLIFY = FALSE)
    # Only one data frame
    rh_data <- do.call("rbind", rh_data)

    # Plot
    time_end <- time
    rh_data$sp_label <- gsub("species", "sp", rh_data$sp)

    ggplot(rh_data, aes_string("time", "rh")) +
        geom_line(aes_string(color = "as.factor(sp_label)"), size = 1) +
        scale_color_discrete(name = "Species") +
        ylab("Hierarchical coefficient") +
        xlab("Time") +
        ggtitle(paste0("Patch ", n_patches)) +
        theme_classic() +
        geom_label(aes(label = ifelse(time == time_end, as.character(sp_label),
                                      NA))) +
        theme_classic() +
        guides(color = FALSE)
}
