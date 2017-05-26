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
source("multigen_function.r")

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

A = 0.001	#Alpha scalar
d = 0.05 #Dispersal percentage

results <- multigen(traits=traits, env=env, time=time, species=species, patches=patches, composition=composition, A=A, d=d)

plot(1:time, composition[5,13,], type="l")
plot(1:time, composition[5,1,], type="l")

final <- ifelse(composition[,,time] <2, 0, stoch[,]) # threshold number of individuals
