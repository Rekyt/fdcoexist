# Graphical function to plot species dynamics through time in a focal patch

#' Plot patch dynamics
#'
#' Plot the dynamics of species in a single patch over time. Automatically
#' assigns one color per species.
#'
#' @param results results array of a given simulation (a [multigen()] output)
#' @param patch an integer indicating the patch number
#' @param time an integer indicating the maximum timestep to look at
#'
#' @import ggplot2
#' @export
plot_patch <- function(results, patch, time){

  patch_data <- list()
  for (i in 1:time) {
    patch_data[[i]] <- results[patch,, i]
  }

  # Convert arrays into data frames
  patch_data <- lapply(patch_data,
                       function(x) data.frame("sp" = names(x),
                                              "n_ind" = as.numeric(x))
                       )

  # Add column to identify time
  patch_data <- mapply(cbind, patch_data, "time" = 1:time, SIMPLIFY = FALSE)

  # Only one data frame
  patch_data <- do.call("rbind", patch_data)

  # Plot
  ggplot(patch_data, aes_string("time", "n_ind")) +
    geom_line(aes_string(color = "as.factor(sp)"), size = 1) +
    scale_color_discrete(name = "Species") +
    ylab("Number of individuals") +
    xlab("Time") +
    ggtitle(paste0("Patch ", patch)) +
    theme_classic()
}
