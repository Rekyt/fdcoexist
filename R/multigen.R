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
#' @param h_fun       a function that describes how hierarchical is combined to
#'                    environmental-based growth (default: `sum()`)
#' @param di_thresh   dissimilary threshold above which species are considered
#'                    maximally dissimilar
#' @param lim_sim_exponent exponent to use for limiting similarity distances
#' @param hierar_exponent  exponent to use for hierarchical compet. distances
#'
#' @export
multigen <- function(traits, trait_weights, env, time, species, patches,
                     composition, A = A, B = B, d, k, width, H, h_fun = "+",
                     di_thresh = 24, lim_sim_exponent = 2,
                     hierar_exponent = 0.5, K) {

    # Check k dimensions
    if ((length(k) != 1 & length(k) != species)) {
        stop("k should be either length one (one k for all species) or ",
             "k should be a matrix with one row per species")
    }

    # Check assumptions on trait_weights data.frame
    check_trait_weights(trait_weights, traits)

    # Calculate dist trait
    disttraits <- compute_compet_distance(trait_weights, traits,
                                          lim_sim_exponent)

    traits_k <- cbind(traits, k)

    # Calculate fitness term (R = growth)
    env_param <- cbind(env, width)
    Rmatrix <- apply(traits_k, 1, function(x) { # Loop over the species
        apply(env_param, 1, function(y){ # Loop over env, k, width combinations
            env_curve(trait_values  = x[-c(length(x))],
                      env_value     = y[1],
                      trait_weights = trait_weights,
                      k             = x[length(x)],
                      width         = y[2])
        })
    })

    # List of alphaterms
    alphalist = vector("list", time - 1)

    # List of R_h terms (extra growth from hierarchical competiton)
    rh_list = vector("list", time - 1)

    for (m in seq(1, time - 1)) {
        # Calculate niche term (alpha) including carrying capacity
        alpha <- alphaterm(disttraits, composition[,,m], A = A, B = B,
                           di_thresh = di_thresh)

        alphalist[[m]] <- alpha

        # Compute hierarchical competition
        R_h <- compute_hierarchical_compet(
            composition_given_time_step = composition[,,m],
            trait_values                = traits,
            trait_weights               = trait_weights,
            H                           = H,
            exponent                    = hierar_exponent)

        rh_list[[m]] <- R_h

        rh_list[[m]][composition[,,m] == 0] <- NA

        # Total growth
        R_tot <- get(h_fun)(Rmatrix, R_h)

        # Update composition
        composition[,, m + 1] <- bevHoltFct(R_tot, composition[,,m], alpha, K)

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
                rmatrix = Rmatrix,
                rhlist  = rh_list,
                ls_exponent = lim_sim_exponent,
                hierar_exponent = hierar_exponent))
}
