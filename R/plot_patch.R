# Graphical function to plot species dynamics through time in a focal patch

#' Plot patch dynamics
#'
#' Plot the dynamics of species in a single patch over time. Automatically
#' assigns one color per species.
#'
#' @param results results array of a given simulation (a [multigen()] output)
#' @param patch an integer indicating the patch number
#' @param time an integer indicating the maximum timestep to look at
#' @param equilibrium a boolean, if TRUE, displays a vertical line when
#' equilibrium is reached
#'
#' @import ggplot2
#' @export
plot_patch <- function(results, patch, time, equilibrium = FALSE){

  patch_data <- list()

  for (i in 1:time) {
    patch_data[[i]] <- results[patch,, i]
  }

  # Convert arrays into data frames
  patch_data <- lapply(
      patch_data,
      function(x) data.frame("sp" = names(x), "n_ind" = as.numeric(x))
  )

  # Add column to identify time
  patch_data <- mapply(cbind, patch_data, "time" = 1:time, SIMPLIFY = FALSE)

  # Only one data frame
  patch_data <- do.call("rbind", patch_data)

  # Reorder levels of factors so that species are ordered
  asc_order <- order(as.integer(gsub("species", "", unique(patch_data$sp))))
  patch_data$sp <- factor(patch_data$sp, unique(patch_data$sp)[asc_order])

  if (!isTRUE(equilibrium)) {
    # Plot
    ggplot(patch_data, aes_string("time", "n_ind")) +
      geom_line(aes_string(color = "as.factor(sp)"), linewidth = 1) +
      scale_color_discrete(name = "Species") +
      labs(
          x = "Time", y = "Number of individuals",
          title = paste0("Patch ", patch)
      ) +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))

  } else {
    nb_ind_t0 <- round(sum(results[patch, , 1]), 0)

    # data.frame of number of individuals per time step
    nb_ind_df <- data.frame(time_step = 1,
                            nb_ind = round(sum(results[patch, , 1]), 0),
                            diff = NA)
    for (i in 2:time) {
      nb_ind_df <- rbind(
        nb_ind_df,
        data.frame(
            time_step = i,
            nb_ind = round(sum(results[patch, , i]), 0),
            diff = sum(abs(round(results[patch, , i], 0) -
                               round(results[patch, , (i - 1)], 0)))
        )
      )
    }

    first_0 <- nb_ind_df[which(nb_ind_df$diff == 0), "time_step"][1]

    # Control for potential later outbreak of a species
    if (!is.na(first_0)) {
      # If more than 10 individuals difference during one time step,
      # search for another equilibrium after that 10 diff
      if (sum(abs(nb_ind_df[first_0:nrow(nb_ind_df), "diff"]) > 10) > 0) {
        last_10 <- nb_ind_df[which(abs(nb_ind_df$diff) > 10),
                             "time_step"]
        last_10 <- last_10[length(last_10)]

        nb_ind_df <- nb_ind_df[last_10:nrow(nb_ind_df), ]

        first_0 <- nb_ind_df[which(nb_ind_df$diff == 0), "time_step"][1]
      }
    }

    # Plot
    ggplot(patch_data, aes_string("time", "n_ind")) +
      geom_line(aes_string(color = "as.factor(sp)"), linewidth = 1) +
      geom_vline(xintercept = first_0, linetype = 2) +
      scale_color_discrete(name = "Species") +
      labs(x = "Time", y = "Number of individuals") +
      theme_classic() +
      theme(panel.border = element_rect(fill = NA))
  }
}
