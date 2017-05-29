## Define algorithms and fcts for communtiy assembly
#' Function definition for deterministic model run with global dispersal

alphaterm <- function(distance, Nts){#needs pop sizes at time t
	Nts %*% distance
}

bevHoltFct <- function(R, N, A, alpha){
    ((R * N) * 1)/(1 + A * alpha)
}

#create trait obj
#  function to generate traits
# n_species = number of species
# min_val = vector minimum values with length number of traits
# max_val = vector of maximum values with length number of traits
generate_traits <- function(n_species, min_val, max_val) {
    
    traits_matrix <- mapply(runif, n = n_species, min = min_val, max = max_val)
    row.names(traits_matrix) <- paste0("species", seq_len(n_species))
    colnames(traits_matrix) <- paste0("trait", seq_along(min_val))

    return(traits_matrix)

}

# Calculate R
# calculates R given some env and trait combo
#   trti = trait matrix
#   envx = environmental value (= optimal value),
#   k    = scalar to get normal growth rate,
#   c    = constant in gaussian function (standard deviation).
envtrtcurve <- function(trti, envx, types, k = 2, c = 0.5) {

    fitness_traits <- trti[types == "R" | types == "RA"]

    # Each trait has similar impact
    R <- k * exp(-((fitness_traits - envx)^2)/2*c^2)

    Rfinal <- mean(R)  # Each trait contributes equally to fitness

    return(Rfinal)
}


alphaterm <- function(distance, Nts){ # needs pop sizes at time t
    Nts %*% distance
}


multigen <- function(traits, trait_type, env, time, species, patches,
                     composition, A, d) {

    # Calculate dist trait
    disttraits <- as.matrix(dist(traits[, trait_type == "A" |
                                            trait_type == "RA"]))

	# Calculate fitness term (R)
	Rmatrix <- apply(traits, 1, function(x, types) {
	    sapply(env, envtrtcurve, trti = x, types = types)
	},
	types = trait_type)

	for (m in seq(1, time - 1)) {

	    # Calculate niche term (alpha)
	    alpha <- alphaterm(disttraits, composition[,,m])
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
    return(composition)
}

