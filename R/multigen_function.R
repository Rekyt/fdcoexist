## Define algorithms and fcts for communtiy assembly
#' Function definition for deterministic model run with global dispersal

#' Compute alpha term in Beverton-Holt function
#'
#' From a competition matrix (for the moment the distance between species
#' traits) and a vector of abundances by species, return the alpha term in the
#' Beverton-Holt equation. The order of species between the two should be the
#' same as no checks are done.
#' Typically the competition matrix is an euclidean trait distance matrix
#' between species. The closer the species are the higher the combination.
#' The term is computed as follow:
#'
#' \deqn{
#'     \alpha_i = \sum_{j = 1, j \neq i}^{S} N_{t, j, x} \times (1 -
#'                                                               \delta_{ij})
#' }{
#'     alpha_i = sum_{j = 1, j != i}^{S} Ntjx * (1 - delta_ij)
#' },
#' where alpha_i is the competition term of species i; Ntjx the abundance of
#' species j, at time t, in patch x and delta_ij the functional distance between
#' species i and species j.
#'
#' @param distance  dissimilarity matrix between species
#' @param Nts       vector of abundances of species at time t
#' @param A         scalar for the inter-specific competition
#' @param B         scalar for the intra-specific competition
#' @param di_thresh dissimilary threshold above which species are considered
#'                  maximally dissimilar
#' @export
alphaterm <- function(distance, Nts, A, B, di_thresh) {

    # From a certain distance, di_thresh, species are considered maximally
    # dissimilar
    # if di_thresh = max(distance) then only true maximally dissimilar are
    # considered maximally dissimilar
    distance[distance >= di_thresh] = max(distance)

    # We compute similarity matrix: when species are not distant they
    # are maximally similar (similarity = max(distance))
    similarity = max(distance) - distance
    # Similarity is for Inter-specific competition only
    diag(similarity) = 0

    A * (Nts %*% similarity) + B * (Nts)
}

#' Beverton-Holt function
#'
#' To simulate growth easily, use the Beverton-Holt equation. Which is:
#'
#' \deqn{
#'    N_{t+1, i, x} = \frac{R_{i, x} \times N_{t, i, x}}{1 + A \times \alpha}
#' }{
#'    N(t+1) = (R * N(t))/(1 + A * alpha)
#' }
#' where t is time, i is the species of interest and x the patch it occupies.
#'
#' @param R     a numeric vector of species growth rates
#' @param N     a numeric vector of species population sizes
#' @param alpha the competition coefficient see [alphaterm()] for its
#'              computation
#'
#' @export
bevHoltFct <- function(R, N, alpha){
    (R * N)/(1 + alpha)
}

#' Species growth rate for a given trait and environment
#'
#' Using traits that affect growth rate and specified environments, this
#' function returns a numeric value of expected growth rates given the traits
#' and the environmental value. The total growth rate is then the average of the
#' growth rates computed with each trait.
#'
#' For the moment the environmental filter follows a Gaussian distribution:
#'
#' \deqn{
#'     R_{i, x} = k \times \exp(- \frac{(trait_i - env_x)^2}{2\times {width}^2})
#' }{
#'     R_ix = k * exp((t_i - env_x)^2/(2 * width^2))
#' },
#' where t_i is trait of species i, env_x the environmental value in patch x, k
#' a scalar giving the maximal growth rate and width the environmental breadth
#' of species.
#'
#' @param trait_values  a numeric vector of species trait values
#' @param env_value     a single numeric value giving the environmental variable
#' @param trait_weights data.frame with at least three columns equal to `trait`
#'                      (giving the name of the concerned traits in `traits`
#'                      df), `growth_weight` the relative weight of the trait in
#'                      growth and `compet_weight` the relative weight of the
#'                      trait in competition, both hierarchical and based on
#'                      limiting similarity.
#' @param k      a scalar giving the maximum growth rate in optimal environment
#' @param width  a numeric for niche breadth, constant in gaussian function
#'
#' @export
env_curve <- function(trait_values, env_value, trait_weights, k = 2,
                      width = 0.5) {

    if (length(width) != 1 & length(width) != length(env_value)) {
        stop("There are either too many or not enough values for width")
    }

    if (identical(sum(trait_weights$growth_weight, na.rm = TRUE), 0)) {
        R <- k
    } else {
        fitness_traits <- trait_weights[trait_weights$growth_weight != 0 &
                                         !is.na(trait_weights$growth_weight),]
        given_values <- trait_values[fitness_traits$trait]

        # Weight traits to obtain a composite traits
        composite_trait <- weighted.mean(given_value,
                                         fitness_traits$growth_weight)

        # Each trait has similar impact
        R <- k * exp(-((composite_trait - env_value)^2)/(2*width^2))
    }
}

#' Check trait weights data.frame
#'
#' This is an internal help function to check the trait weights data.frame
#' (used in [multigen()]). The structure of the given trait weights is fixed:
#' one column should be named `trait`` with the names of the traits in it. The
#' two other columns `growth_weight` and `compet_weight` contains respectively
#' the relative weight of the trait in growth and competition.
#' The function also additionnally check that the traits in the `trait` column
#' are in the `traits` data.frame.
#'
#' @param trait_weights data.frame with at least three columns equal to `trait`
#'                      (giving the name of the concerned traits in `traits`
#'                      df), `growth_weight` the relative weight of the trait in
#'                      growth and `compet_weight` the relative weight of the
#'                      trait in competition.
#' @param traits        a species-traits data.frame with species as rownames and
#'                      traits as numeric columns with names matching
#'                      `trait_weights` column `trait`
#'
#' @return nothing if data.frame passes the checks, stops early otherwise.
#' @export
#' @examples
#' # Working trait weights data.frame
#' traits = data.frame(trait1 = 1, trait2 = 2, trait3 = 3)
#'
#' weight_1 = data.frame(
#'    trait = c("trait1", "trait2", "trait3"),
#'    growth_weight    = c(0.5, 0.5, 0),
#'    compet_weight    = c(0,   0.5, 0.5),
#'    hierarchy_weight = c(0,   0,   0))
#'
#' # Silent function
#' check_trait_weights(weight_1, traits)
#' \dontrun{
#' # Not valid trait weights data.frame
#' not_valid = data.frame(trait = c("trait1", "trait2", "trait3"),
#'     growth_weight = c(0.5, 0.8, 0),
#'     compet_weight = c(0, 0.5, 0.9))
#'
#' # Stop and error
#' check_trait_weights(not_valid, traits)
#' }
check_trait_weights = function(trait_weights, traits) {

    if (!is.data.frame(trait_weights) & !is.matrix(trait_weights)) {
        stop("trait_weights should be a data.frame or a matrix")
    }

    # Check that trait_weightts contains the needed columns
    needed_cols = c("trait", "growth_weight", "compet_weight",
                    "hierarchy_weight")

    if (!all(needed_cols %in% colnames(trait_weights))) {

        absent_columns = setdiff(needed_cols, colnames(trait_weights))

        stop("Column(s) ", paste(absent_columns, collapse = ", "),
             " is (are) not in trait_weights")

    }

    # Check that traits in trait_weights are in traits data.frame
    if (!all(trait_weights$trait %in% colnames(traits))) {

        absent_traits = setdiff(trait_weights$trait, colnames(traits))

        stop("Trait(s) ", paste(absent_traits, collapse = ", "),
             " (is/are) not in provided traits data.frame")
    }
}

#' Compute trait distance between species
#'
#' This function compute trait distance between species using a trait matrix and
#' a trait weights data.frame. For all the traits with competition weights not
#' equal to zero, it computes a weighted 'composite trait' that is then used to
#' compute euclidean trait distance between species. Trait distance is first
#' exponentiated then standardized between 0 and 1.
#'
#' @inheritParams check_trait_weights
#' @param exponent \[`numeric(1)`\] (default: 1)\cr{}
#'                 The exponent used before standardizing the distance
#'
#' @return an euclidean distance matrix (of type matrix)
#' @importFrom stats dist weighted.mean
#' @export
compute_compet_distance = function(trait_weights, traits, exponent = 1) {

    # If there is no defined competition trait, there is no competition
    if (sum(trait_weights$compet_weight, na.rm = TRUE) == 0) {

        disttraits <- matrix(1, nrow = nrow(traits), ncol = nrow(traits),
                             dimnames = list(rownames(traits),
                                             rownames(traits)))

    } else {

        # Extract weights of traits contributing to competition
        compet_weights <- trait_weights[!is.na(trait_weights$compet_weight) &
                                            trait_weights$compet_weight != 0,]

        compet_traits <- traits[, compet_weights$trait, drop = FALSE]


        # Transform columns prior computing distance
        scaled_compet_traits <- sweep(compet_traits, 2,
                                      compet_weights$compet_weight,
                                      function(x,y) x * sqrt(y))

        # Compute distance matrix
        disttraits <- as.matrix(dist(scaled_compet_traits))^exponent

        disttraits <- disttraits/max(disttraits)
    }

    return(disttraits)
}

#' Scale distance or matrix between 0 and 1
#'
#' dist - min(dist) / (max(dist) - min(dist))
#' @param dist_matrix a matrix or a dist object
scale_distance = function(dist_matrix) {
    denom_standard = 1

    if (diff(range(dist_matrix)) != 0) {
        denom_standard = diff(range(dist_matrix))
    }

    dist_matrix <- (dist_matrix - min(dist_matrix)) / denom_standard
}


#' Compute Hierarchical Competition coefficient at each time step
#'
#' Outputs a matrix of additional growth per patch per species given by
#' hierarchical competition. The values are considered
#' @param composition_given_time_step composition matrix at a given time step
#'                                    (a site-species matrix with sites in rows)
#' @param trait_values a trait matrix
#' @param trait_weights a scenario data.frame
#' @param H the hierarchical competition scalar
#' @inheritParams compute_compet_distance
#'
#' @export
compute_hierarchical_compet = function(composition_given_time_step,
                                       trait_values, trait_weights,
                                       H, exponent = 1) {

    # Add to R effect of hierarchical traits
    # (so far same traits as limiting similarity traits)
    hierarchical_trait <- trait_weights[
        trait_weights$hierarchy_weight != 0 &
            !is.na(trait_weights$hierarchy_weight),]

    if (nrow(hierarchical_trait) == 0) {
        Rh <- matrix(0,
                     nrow = nrow(composition_given_time_step[,,1]),
                     ncol = ncol(composition_given_time_step[,,1]),
                     dimnames = dimnames(composition_given_time_step[,,1]))
    } else {
        hierarchical_values <- trait_values[, hierarchical_trait$trait,
                                            drop = FALSE]

        # For each trait compute a hierarchical competition value
        traits_Rh = lapply(as.data.frame(hierarchical_values),function(x) {
            species = rownames(hierarchical_values)

            # Compute the difference of traits between all pairs of species
            single_trait_diff = outer(x, x, "-") / (diff(range(x)))

            # Smaller species have no effect in hierarchical competition
            single_trait_diff[single_trait_diff > 0] = 0

            # Because all differences are negative we need to take absolute
            # to exponentiate distances
            single_trait_diff = -(abs(single_trait_diff)^exponent)

            rownames(single_trait_diff) = species
            colnames(single_trait_diff) = species

            # Weight the difference by the abundance of species
            single_trait_Rh = composition_given_time_step %*%
                t(single_trait_diff)
            single_trait_Rh[composition_given_time_step == 0] = 0

            single_trait_Rh
        })

        # Weight the different hierarchical growth by the importance of
        # the trait for hierarchical competition
        Rh = Reduce(`+`,
                    Map(`*`, traits_Rh, hierarchical_trait$hierarchy_weight))

        Rh[composition_given_time_step == 0] = 0L
    }

    return(H * Rh)
}
