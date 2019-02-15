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
#' @param distance dissimilarity matrix between species
#' @param Nts      vector of abundances of species at time t
#' @param A        scalar for the inter-specific competition
#' @param B        scalar for the intra-specific competition
#' @export
alphaterm <- function(distance, Nts, A, B) {
    A * (Nts %*% (1 - distance)) + B * (Nts)
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
#' @param H      a numeric for hierarchical competition such as H/k <= 1
#' @param th_max a numeric for the hierarchical trait value maximizing
#'               hierarchical competiton
#' @param h_fun  a character name of a function that describes how hierarchical
#'               is combined to environmental-based growth (default: `"sum"`)
#'
#' @export
env_curve <- function(trait_values, env_value, trait_weights, k = 2,
                      width = 0.5, H = 0, th_max = 25, h_fun = "sum") {

    if (length(width) != 1 & length(width) != length(env_value)) {
        stop("There are either too many or not enough values for width")
    }

    if (identical(sum(trait_weights$growth_weight, na.rm = TRUE), 0)) {
        R <- k
    } else {
        fitness_traits <- trait_weights[trait_weights$growth_weight != 0 &
                                         !is.na(trait_weights$growth_weight),]

        given_values <- trait_values[fitness_traits$trait]

        # Each trait has similar impact
        R <- k * exp(-((given_values - env_value)^2)/(2*width^2))

        # Weigh each trait function of contribution to growth
        R <- weighted.mean(R, fitness_traits$growth_weight)

        # Add to R effect of hierarchical traits
        # (so far same traits as limiting similarity traits)
        hierarchical_trait <- trait_weights[trait_weights$hierarchy_weight != 0 &
                                         !is.na(trait_weights$hierarchy_weight),]

        hierarchical_values <- trait_values[hierarchical_trait$trait]

        if (any(hierarchical_values > th_max)) {
            stop("th_max must be superior to any hierarchical trait value in a",
                 " directional filter perspective")
        }

        R_h <- H * (hierarchical_values/(th_max))
        # Weigh each trait function of contribution to competition
        R_h <- weighted.mean(R_h, hierarchical_trait$hierarchy_weight)

        # Combining environmental-based growth and additional growth due to
        # hierarchical competition
        R <- get(h_fun)(R, R_h, na.rm = TRUE)
    }
}


#' Function to run the simulation
#'
#' Using specified parameters this function run the simulation
#'
#' @inheritParams compute_compet_distance
#' @param env         a vector of environmental values
#' @param time        an integer giving the total number of generations
#' @param species     an integer giving the total number of species to simulate
#' @param patches     an integer giving the total number of patches to simulate
#' @param composition the actual array containing species abundances per site
#'                    over time, giving the initial populations of each species
#' @param A           the scalar of inter-specific competition coefficient
#'                    (see [bevHoltFct()])
#' @param B           the scalar for intra-specific competition coefficient
#'                    by default B = A
#' @param d           a numeric value between 0 and 1 giving the percentage of
#'                    dispersal occurring across all patches
#' @param k           a scalar giving the maximum growth rate in optimal
#'                    environment
#' @param width       a numeric giving niche breadth of all species
#' @param H           a numeric for hierarchical competition such as H/k <= 1
#' @param th_max      a numeric for the hierarchical trait value maximizing
#'                    hierarchical competition
#' @param h_fun       a function that describes how hierarchical is combined to
#'                    environmental-based growth (default: `sum()`)
#'
#' @export
multigen <- function(traits, trait_weights, env, time, species, patches,
                     composition, A = A, B = B, d, k, width, H, th_max,
                     dist_power = dist_power, h_fun = "sum") {

    # Check k dimensions
    if ((length(k) != 1 & length(k) != species)) {
        stop("k should be either length one (one k for all species) or ",
             "k should be a matrix with one row per species")
    }

    # Check assumptions on trait_weights data.frame
    check_trait_weights(trait_weights, traits)

    # Calculate dist trait
    disttraits <- compute_compet_distance(trait_weights, traits, dist_power)

    traits_k_and_H <- cbind(traits, k, H)

    # Calculate fitness term (R = growth)
    env_param <- cbind(env, width, th_max)
    Rmatrix <- apply(traits_k_and_H, 1, function(x) { # Loop over the species
        apply(env_param, 1, function(y){ # Loop over env, k, width combinations
            env_curve(trait_values  = x[-c(length(x) - 1, length(x))],
                      env_value     = y[1],
                      trait_weights = trait_weights,
                      k             = x[length(x) - 1],
                      width         = y[2],
                      H             = x[length(x)],
                      th_max = y[3],
                      h_fun  = h_fun)
        })
    })

    # Compute only the environmental part of growth_rate
    Rmatrix_env <- apply(traits_k_and_H, 1, function(x) { # Loop over the species
        apply(env_param, 1, function(y){ # Loop over env, k, width combinations
            env_curve(trait_values  = x[-c(length(x) - 1, length(x))],
                      env_value     = y[1],
                      trait_weights = trait_weights,
                      k             = x[length(x) - 1],
                      width         = y[2],
                      H             = 0,
                      th_max        = y[3],
                      h_fun         = h_fun)
        })
    })

    # List of alphaterms
    alphalist = list()

    for (m in seq(1, time - 1)) {

        # Calculate niche term (alpha) including carrying capacity
        alpha <- alphaterm(disttraits, composition[,,m], A = A, B = B)

        alphalist[[m]] <- alpha

        composition[,, m + 1] <- bevHoltFct(Rmatrix, composition[,,m], alpha)

        # threshold number of individuals
        composition[,, m + 1] <- ifelse(composition[,, m + 1] < 2, 0,
                                        composition[,, m + 1])
        ## Dispersal
        # Probability of dispersal proportional to number of individuals
        # in patch given a certain probability 'd'
        migrate <- d * composition[,, m + 1]

        stay <- composition[,, m + 1] - migrate

        # All immigrants move evenly to all patches
        immigrants <- apply(migrate, 2, "sum")/patches
        # /!\ vector recycling when summing both matrices
        total <- t(t(stay) + immigrants)

        # Final composition
        composition[,, m + 1] <- total

        # threshold number of individuals
        composition[,, m + 1] <- ifelse(composition[,, m + 1] < 2, 0,
                                        composition[,, m + 1])
    }
    return(list(compo   = composition,
                alpha   = alphalist,
                rmatrix = Rmatrix,
                rmatenv = Rmatrix_env))
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
#' compute euclidean trait distance between species.
#'
#' @inheritParams check_trait_weights
#' @param dist_power power used when weighting distance exponentially
#'                   (`NULL` by default = euclidean distance)
#'
#' @return an euclidean distance matrix (of type matrix)
#' @importFrom stats dist weighted.mean
#' @export
compute_compet_distance = function(trait_weights, traits, dist_power = NULL) {

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
        disttraits <- as.matrix(dist(scaled_compet_traits))

        # Scale trait distance to balance growth
        disttraits <- scale_distance(disttraits)

        # If required weight distance exponentially
        if (!is.null(dist_power)) {
        disttraits <- exp(disttraits^dist_power)

        disttraits <- scale_distance(disttraits)
        }
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
