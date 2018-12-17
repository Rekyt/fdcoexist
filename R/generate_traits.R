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
        x2 <- stats::runif(number_species, 1, number_patches)

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


#' Generates a data.frame of trait weights
#'
#' Only in the special case of 3 traits with 1 always driving growth and
#' another one only driving competition
#' @param R an integer value giving the contribution of trait in growth
#' @param A an integer value giving the contribution of trait in limiting
#'          similarity
#' @param H an integer value giving the contribution of trait in hierarchical
#'          competition
#'
#' @export
#'
#' @examples create_trait_weights(50, 50, 0)
create_trait_weights = function(R, A, H) {
    if (any(!is.numeric(R), !is.numeric(A), !is.numeric(H))) {
        stop("All provided weights should be integer")
    }

    if (any(R > 100, R < 0, A > 100, A < 0, H > 100, H < 0)) {
        stop("All weights should be in [0; 100]")
    }

    g_weights = c(100 - R, R,       0,       0)
    c_weights = c(      0, A, 100 - A,       0)
    h_weights = c(      0, H,       0, 100 - H)

    data.frame(
        trait            = paste0("trait", 1:4),
        growth_weight    = g_weights/100,
        compet_weight    = c_weights/100,
        hierarchy_weight = h_weights/100
    )
}
