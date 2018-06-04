## Define algorithms and fcts for communtiy assembly
#' Function definition for deterministic model run with global dispersal

#' Compute alpha term in Beverton-Holt function
#'
#' From a competition matrix (for the moment the distance between species
#' traits) and a vector of abundances by species, return the alpha term in the
#' Beverton-Holt equation. The order of species between the two should be the
#' same as no checks are done.
#'
#' @param distance competiton matrix of species
#' @param Nts      vector of abundances of species at time t
#' @export
alphaterm <- function(distance, Nts) {
	Nts %*% distance
}

#' Beverton-Holt function
#'
#' To simulate growth easily, use the Beverton-Holt equation (ref?)
#'
#' @param R     a vector of species growth rates
#' @param N     a vector of species population sizes
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

    traits_matrix <- mapply(seq, from = min_val, to = max_val, length.out = n_species)
    row.names(traits_matrix) <- paste0("species", seq_len(n_species))
    colnames(traits_matrix) <- paste0("trait", seq_along(min_val))

    return(traits_matrix)

}

#' Species growth rate for a given trait and environment
#'
#' Using traits that affect growth rate and specified environments, this
#' functionr returns a data.frame of expected growth rates given the traits and
#' the environmental value. Suppose a gaussian relationship between growth rate
#' and environmental value. The total growth rate is then the average of the
#' growth rates computed with each trait.
#'
#' @param trti  a data.frame of species' traits
#' @param envx  a data.frame
#' @param types a vector of trait types indicating if traits should be used to
#'              compute growth rate (types `RA` or `R`)
#' @param k     scalar for growth rate
#' @param c     constant in gaussian function (standard deviation)
#'
#' @export
envtrtcurve <- function(trti, envx, types, k = 2, c = 0.5) {

    if (!(sum(c("R", "RA") %in% types))) {
        R <- k
    } else {
        fitness_traits <- trti[types == "R" | types == "RA"]

        # Each trait has similar impact
        R <- k * exp(-((fitness_traits - envx)^2)/2*c^2)
    }

    Rfinal <- mean(R)  # Each trait contributes equally to fitness

    return(Rfinal)
}


#' Function to run the simulation
#'
#' Using specified parameters this function run the simulation
#'
#' @param traits a species-traits data.frame
#' @param trait_type a character vector indicating trait types (`R`, `A`, `RA`
#'                   or `N`)
#' @param env a vector of environmenal values
#' @param time a integer giving the number of generations
#' @param species ?
#' @param patches the number of patches?
#' @param composition ?
#' @param A the scalar of competition coefficent (see [bevHoltFct()])
#' @param d ?
#' @param k a scalar for computation fo growth rate
#' @param c a constant to compute growth rates
#'
#' @importFrom stats dist
#' @export
multigen <- function(traits, trait_type, env, time, species, patches,
                     composition, A, d, k, c) {

    # Calculate dist trait
    if (!(sum(c("A", "RA") %in% trait_type))) {

        # If there is no defined competition trait, there is no competition
        disttraits <- matrix(0, nrow = species, ncol = species)
    } else {

        disttraits <- as.matrix(dist(traits[, trait_type == "A" |
                                                trait_type == "RA"]))

        # Scale trait distance to balance growth
        disttraits <- (disttraits - min(disttraits)) / (max(disttraits) - min(disttraits))
    }

	# Calculate fitness term (R)
	Rmatrix <- apply(traits, 1, function(x, types) {
	    sapply(env, envtrtcurve, trti = x, types = types, k = k, c = c)
	},
	types = trait_type)


	# List of alphaterms
	alphalist = vector("numeric", length = time - 1)

	for (m in seq(1, time - 1)) {

	    # Calculate niche term (alpha)
	    alpha <- alphaterm(disttraits, composition[,,m]) * k

	    alphalist[m] = alpha

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
	    total <- stay + immigrants
	    composition[,,m + 1] <- total
	}
    return(list(compo = composition,
                alpha = alphalist))
}

