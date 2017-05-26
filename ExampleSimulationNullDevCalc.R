#######################
### Code for example metacommunity simulation and beta-null deviation calculations
### with "Differentiating between niche and neutral assembly in metacommunities using
### null models of beta-diversity"
### Prepared May 14, 2014
### Authors Caroline Tucker, Lauren Shoemaker, Brett Melbourne
#######################


## Load required source files and libraries
library(reldist)
library(vegan)
library(bipartite)
library(mvtnorm)
source("MetacommunityDynamicsFctsOikos.r")

## Set number of patches, species, time
patches  <- 10   # Number of patches
species  <- 25   # Number of species
time     <- 150	 # Length of model run (generations)
initpop  <- 25   # initial population size
n_traits <- 3    # number of traits

composition <- array(NA, dim=c(patches, species, time), dimnames=list(paste0("patches", 1:patches), paste0("species", 1:species), paste0("time", 1:time))) #where values = N
composition[,,"time1"] <- initpop #N

#create env obj
env <- 1:patches
names(env) <- dimnames(composition)[[1]]

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

traits <- generate_traits(species,
                          rep(min(env), n_traits),
                          rep(max(env), n_traits))

# Vector that indicates trait contribution:
#   "R"  = trait contributes to fitness,
#   "A"  = trait contributes to niche difference,
#   "N"  = trait contributes to none,
#   "RA" = trait contributes to BOTH fitness and niche difference
trait_type <- rep("RA", n_traits)
names(trait_type) <- colnames(traits)

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

#plot(envtrtcurve(1, env, 2, 0.5)~env, type="l")
#points(envtrtcurve(6, env, 1, 0.5)~env, type="l")

# Computes R value for each combination of traits (per species) and environnment
Rmatrix <- apply(traits, 1, function(x, types) {
        sapply(env, envtrtcurve, trti = x, types = types)
    },
    types = trait_type)


# Calculate dist trait
disttraits <- as.matrix(dist(traits[, trait_type == "A" | trait_type == "RA"]))

alphaterm <- function(distance, Nts){#needs pop sizes at time t
	Nts %*% distance
}

alpha <- alphaterm(disttraits, composition[,,1])

A = 0.001

bevHoltFct <- function(R, N, A, alpha){ ((R * N) * 1)/(1 + A * alpha)
	}

for(m in seq(1, time - 1)){
	composition[,, m + 1] <- bevHoltFct(Rmatrix, composition[,,m], A, alpha)
}

plot(1:time, composition[5,13,], type="l")
plot(1:time, composition[5,1,], type="l")

final <- ifelse(composition[,,time] <2, 0, stoch[,]) # threshold number of individuals
