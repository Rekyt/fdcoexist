
# Graphical function to plot species dynamics through time in a focal patch

plot_patch <- function(patch, time){
  require(ggplot2)
  
  patch_data <- list()
  for (i in 1:time){
    patch_data[[i]] <- results[patch,, i]
  }
  # Convert arrays into data frames
  patch_data <- lapply(patch_data, function(x) data.frame("sp" = names(x), "n_ind" = as.numeric(x)))
  # Add column to identify time
  patch_data <- mapply(cbind, patch_data, "time" = 1:time, SIMPLIFY = FALSE)
  # Only one data frame
  patch_data <- do.call("rbind", patch_data)
  # Plot
  ggplot(patch_data, aes(time, n_ind)) +
    geom_line(aes(color = as.factor(sp)), size = 1) +
    scale_color_discrete(name = "Species") +
    ylab("Number of individuals") +
    xlab("Time") +
    ggtitle(paste0("Patch ", patch)) +
    theme_classic()
}
