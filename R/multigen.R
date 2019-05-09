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
#' @param th_min      a numeric for the hierarchical trait value minimizing
#'                    hierarchical competition
#' @param h_fun       a function that describes how hierarchical is combined to
#'                    environmental-based growth (default: `sum()`)
#' @param di_thresh   dissimilary threshold above which species are considered
#'                    maximally dissimilar
#'
#' @export
multigen <- function(traits, trait_weights, env, time, species, patches,
                     composition, A = A, B = B, d, k, width, H, th_max, th_min,
                     h_fun = "sum", di_thresh = th_max - th_min) {

    # Check k dimensions
    if ((length(k) != 1 & length(k) != species)) {
        stop("k should be either length one (one k for all species) or ",
             "k should be a matrix with one row per species")
    }

    # Check assumptions on trait_weights data.frame
    check_trait_weights(trait_weights, traits)

    # Calculate dist trait
    disttraits <- compute_compet_distance(trait_weights, traits)

    traits_k_and_H <- cbind(traits, k, H)

    # Calculate fitness term (R = growth)
    env_param <- cbind(env, width, th_max, th_min)
    Rmatrix <- apply(traits_k_and_H, 1, function(x) { # Loop over the species
        apply(env_param, 1, function(y){ # Loop over env, k, width combinations
            env_curve(trait_values  = x[-c(length(x) - 1, length(x))],
                      env_value     = y[1],
                      trait_weights = trait_weights,
                      k             = x[length(x) - 1],
                      width         = y[2],
                      H             = x[length(x)],
                      th_max = y[3],
                      th_min = y[4],
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
                      th_min        = y[4],
                      h_fun         = h_fun)
        })
    })

    # List of alphaterms
    alphalist = list()

    for (m in seq(1, time - 1)) {

        # Calculate niche term (alpha) including carrying capacity
        alpha <- alphaterm(disttraits, composition[,,m], A = A, B = B,
                           di_thresh = di_thresh)

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
