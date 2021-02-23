

#' Generate random traits
#' Compared to generate_cor_traits() introduce a little of variability in first
#' trait as instead of being directly determined by the species number it adds
#' little white noise to it and scale it to a minimum of 0 if negative and
#' maximum of 25 if maximum value is over 25.
#'
#' @export
generate_cor_traits_rand <- function (number_patches, number_species,
                                      number_other = 9, cor_coef = 0.7,
                                      min_value = 1)
 {
     x1 = seq(min_value, number_patches, length.out = number_species)
     # Add small random noise
     x1 = x1 + rnorm(length(x1), mean = 0, sd = 1)

     # If negative trait value replace negative trait values with random number
     # between 0 and first non-negative number
     x1 = sort(x1)
     first_non_neg = x1[which(x1 > 0)[1]]
     x1 = sort(ifelse(x1 < 0, runif(1, 0, first_non_neg), x1))
     x1 = ifelse(x1 < number_patches, x1, number_patches)
     cov_mat = matrix(c(1, cor_coef, 0, sqrt(1 - cor_coef^2)), nrow = 2)
     n_other <- number_other
     res <- lapply(seq(n_other), function(x) {
         x2 <- stats::runif(number_species, 1, number_patches)
         cov_mat %*% matrix(c(x1, x2), nrow = 2, byrow = TRUE)
     })
     res <- do.call(rbind, lapply(res, function(x) x[2, ]))
     res <- rbind(x1, res)
     res <- t(res)
     trait_mat <- cbind(res[, 1], sapply(2:ncol(res), function(x) {
         val <- res[, x]
         (val - min(val))/(max(val) - min(val)) * (diff(range(x1))) +
             min(x1)
     }))
     colnames(trait_mat) <- paste0("trait", seq(number_other +
         1))
     rownames(trait_mat) <- paste0("species", seq(number_species))
     return(trait_mat)
}

