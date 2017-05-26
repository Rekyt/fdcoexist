## Define algorithms and fcts for communtiy assembly
#' Function definition for deterministic model run with global dispersal

alphaterm <- function(distance, Nts){#needs pop sizes at time t
	Nts %*% distance
}
#calculates R given some env and trait combo
envtrtcurve <- function(trti, envx, k=2, c=0.5){ k*exp(-1*((trti-envx)^2)/2*c^2)}


bevHoltFct <- function(R, N, A, alpha){ ((R * N) * 1)/(1 + A * alpha)
	}
	
multigen <- function(traits, env, time, species, patches, composition, A, d) {
	#Calculate dist trait 
	disttraits <- as.matrix(dist(traits))
	#Calculate fitness term (R)
	Rmatrix <- sapply(traits, envtrtcurve, env)
	
	for (m in seq(1, time - 1)) {

	#Calculate niche term (alpha)
	alpha <- alphaterm(disttraits, composition[,,m])
	composition[,,m+1] <- bevHoltFct(Rmatrix, composition[,,m], A, alpha)	
	composition[,,m+1] <- ifelse(composition[,, m+1] <2, 0, composition[,, m+1]) # threshold number of individuals
		# dispersal
		migrate <- d*composition[,,m+1]
		stay <- composition[,,m+1] - migrate
		immigrants <- apply(migrate, 2, "sum")/patches
		total <- stay + immigrants
		composition[,,m+1] <- total
			}
    return(composition)
}

