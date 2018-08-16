#' Generate Traits data.frame
#'
#' Generates a random traits data.frame by default the traits are uniformly
#' & independently distributed with specifiex minimum and maximum values
#'
#' @param n_species the number of species to generate
#' @param min_val   minimum trait value (for all traits)
#' @param max_val   maximum trait value (for all traits)
#'
#' @importFrom stats runif
#' @export
generate_traits <- function(n_species, min_val, max_val) {

    traits_matrix <- mapply(seq, from = min_val, to = max_val,
                            length.out = n_species)
    row.names(traits_matrix) <- paste0("species", seq_len(n_species))
    colnames(traits_matrix) <- paste0("trait", seq_along(min_val))

    return(traits_matrix)

}

#' Generate correlated traits
#'
#' This function generates a matrix of traits based on the given number of
#' patches, species, the number of "additional" traits with a given correlation
#' coefficient. The first trait is always uniform between 1 and the number of
#' patches provided. Then the other traits are generated correlated to this
#' first one. In the end all the traits are scaled between 1 and the number of
#' patches.
#' @param number_patches a numeric value giving the total number of pathces
#' @param number_species a numeric value giving the total number of species to
#'                       simulate (number of rows in the trait tables)
#' @param number_other   a numeric value giving the number of additional traits
#'                       to generate in addition to the trait correlated to the
#'                       number of patches
#' @param cor_coef       a numeric value giving the correlation coefficient
#'                       between the first trait and the other ones
#'
#' @return a matrix of traits with species in rows and traits in columns
#' @export
#'
#' @examples
#' traits <- generate_cor_traits(25, 100, 3, 0.3)
generate_cor_traits = function(number_patches, number_species, number_other = 9,
                               cor_coef = 0.7) {

    # Initial Trait
    x1 = seq(1, number_patches, length.out = number_species)


    # Variance-covariance matrix between first and any other traits
    cov_mat = matrix(c(1, cor_coef, 0, sqrt(1 - cor_coef^2)), nrow = 2)

    # Number of additional trait to generate
    n_other <- number_other

    res <- lapply(seq(n_other), function(x) {

        # Independent trait
        x2 <- runif(number_species, 1, number_patches)

        # Correlated trait in second row
        cov_mat %*% matrix(c(x1, x2), nrow = 2, byrow = TRUE)

    })

    # Bind all the additional traits
    res <- do.call(rbind, lapply(res, function(x) x[2,]))
    res <- rbind(x1, res)

    res <- t(res)

    # Repaste first column
    trait_mat <- cbind(res[,1],

          # Scale all the traits between min and max value of x1
          sapply(2:ncol(res), function(x) {
              val <- res[,x]

              (val - min(val))/(max(val) - min(val)) * (diff(range(x1))) +
                  min(x1)
          }))

    colnames(trait_mat) <- paste0("trait", seq(number_other + 1))
    rownames(trait_mat) <- paste0("species", seq(number_species))

    return(trait_mat)
}
