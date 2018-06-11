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
#' @param distance competiton matrix of species
#' @param Nts      vector of abundances of species at time t
#' @export
alphaterm <- function(distance, Nts) {
    Nts %*% (1 - distance)
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
#' @param A     a scalar for competition coefficient
#' @param alpha the competition coefficient see [alphaterm()] for its
#'              computation
#'
#' @export
bevHoltFct <- function(R, N, A, alpha){
    (R * N)/(1 + A * alpha)
}

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
#' @param trti   a numeric vector of species trait values
#' @param envx   a
#' @param types  a character vector of trait types indicating if traits should be
#'               used to compute growth rate (types `RA` or `R`)
#' @param k      a scalar giving the maximum growth rate in optimal environment
#' @param width  a numeric for niche breadth, constant in gaussian function
#' @param weight a numeric vector indicating the contribution of each growth
#'               trait to growth
#'
#' @export
env_curve <- function(trti, envx, types, k = 2, width = 0.5, weight = NULL) {

    if (length(width) != 1 & length(width) != length(envx)) {
        stop("There are either too many or not enough values for width")
    }

    if (!is.null(weight) & !is.numeric(weight)) {
        stop("Weight(s) should be numeric")
    }

    if (is.null(weight)) {
        weight = rep(1, sum(c("R", "RA") %in% types))
    }

    if (length(weight) != sum(c("R", "RA") %in% types)) {
        stop("Please specify a weight for each trait")
    }

    if (!(sum(c("R", "RA") %in% types))) {
        R <- k
    } else {
        fitness_traits <- trti[types == "R" | types == "RA"]

        # Each trait has similar impact
        R <- k * exp(-((fitness_traits - envx)^2)/(2*width^2))
    }

    # Weigh each trait function of contribution to growth
    Rfinal <- weighted.mean(R, weight)

    return(Rfinal)
}


#' Function to run the simulation
#'
#' Using specified parameters this function run the simulation
#'
#' @param traits      a species-traits data.frame with species as rownames and
#'                    traits as numeric columns with names matching `trait_type`
#' @param trait_type  a character vector indicating trait types (`R`=contributing
#'                    to growth only, `A` = contributing to competition only,
#'                    `RA` = contributing to both growth and competition or
#'                    `N` = trait not contributing to growth nor competition)
#' @param env         a vector of environmental values
#' @param time        an integer giving the total number of generations
#' @param species     an integer giving the total number of species to simulate
#' @param patches     an integer giving the total number of patches to simulate
#' @param composition the actual array containing species abundances per site
#'                    over time, giving the initial populations of each species
#' @param A           the scalar of competition coefficent (see [bevHoltFct()])
#' @param d           a numeric value between 0 and 1 giving the percentage of
#'                    dispersal occurring across all patches
#' @param k           a scalar giving the maximum growth rate in optimal
#'                    environment
#' @param width       a numeric giving niche breadth of all species
#'
#' @importFrom stats dist
#' @export
multigen <- function(traits, trait_type, env, time, species, patches,
                     composition, A, d, k, width) {

    # Calculate dist trait
    if (!(sum(c("A", "RA") %in% trait_type))) {

        # If there is no defined competition trait, there is no competition
        disttraits <- matrix(0, nrow = species, ncol = species)
    } else {

        disttraits <- as.matrix(dist(traits[, trait_type == "A" |
                                                trait_type == "RA"]))

        # Scale trait distance to balance growth
        disttraits <- (disttraits - min(disttraits)) / diff(range(disttraits))
    }

    # Calculate fitness term (R)
    env_param <- cbind(env, k, width)
    Rmatrix <- apply(traits, 1, function(x) { # Loop over the species traits
        apply(env_param, 1, function(y){ # Loop over env, k, width combinations
            env_curve(x, y[1], types = trait_type, k = y[2], width = y[3])
        })
    })

    # List of alphaterms
    alphalist = list()

    for (m in seq(1, time - 1)) {

        # Calculate niche term (alpha)
        alpha <- alphaterm(disttraits, composition[,,m])

        alphalist[[m]] = alpha

        composition[,, m + 1] <- bevHoltFct(Rmatrix, composition[,,m], A, alpha)
        # threshold number of individuals
        composition[,, m + 1] <- ifelse(composition[,, m + 1] < 2, 0,
                                        composition[,, m + 1])
        ## Dispersal
        # Probability of dispersal proportional to number of individuals in patch
        # given a certain probability
        migrate <- d * composition[,, m + 1]

        stay <- composition[,, m + 1] - migrate

        # All immigrants move evenly to all patches
        immigrants <- apply(migrate, 2, "sum")/patches
        total <- t(t(stay) + immigrants) # beware of vector recycling when summing both matrices
        composition[,,m + 1] <- total
    }
    return(list(compo = composition,
                alpha = alphalist))
}

