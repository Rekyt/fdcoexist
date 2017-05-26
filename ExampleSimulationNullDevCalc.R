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
source("MetacommunityDynamicsFctsOikos.r")

## Set number of patches, species, time
patches <- 10  # Number of patches
species <- 25  # Number of species
time <- 150	   # Length of model run (generations)
initpop <- 25 #initial population size

composition <- array(NA, dim=c(patches, species, time), dimnames=list(paste0("patches", 1:patches), paste0("species", 1:species), paste0("time", 1:time))) #where values = N
composition[,,"time1"] <- initpop #N

#create env obj
env <- 1:patches
names(env) <- dimnames(composition)[[1]]

#create trait obj
traits <- seq(from=1, to=patches, length.out=species)
names(traits) <- dimnames(composition)[[2]]

#Calculate dist trait 
disttraits <- as.matrix(dist(traits))

#Calculate R 
#calculates R given some env and trait combo
envtrtcurve <- function(trti, envx, k=2, c=0.5){ k*exp(-1*((trti-envx)^2)/2*c^2)}

#plot(envtrtcurve(1, env, 2, 0.5)~env, type="l")
#points(envtrtcurve(6, env, 1, 0.5)~env, type="l")

Rmatrix <- sapply(traits, envtrtcurve, envx=env)

alphaterm <- function(distance, Nts){#needs pop sizes at time t
	Nts %*% distance
}

alpha <- alphaterm(disttraits, composition[,,1])

A = 0.001

bevHoltFct <- function(R, N, A, alpha){ ((R * N) * 1)/(1 + A * alpha)
	}

for(m in seq(1, time - 1)){		
	composition[,,m+1] <- bevHoltFct(Rmatrix, composition[,,m], A, alpha)	
}

plot(1:time, composition[5,13,], type="l")
plot(1:time, composition[5,1,], type="l")

final <- ifelse(composition[,,time] <2, 0, stoch[,]) # threshold number of individuals
