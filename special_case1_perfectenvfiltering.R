#need special case - 1 trait, sorted to match species, for perfect environmental filtering, and 1 case where neutral, and 1 case where competition


#create trait obj
#  function to generate traits
# n_species = number of species
# min_val = vector minimum values with length number of traits
# max_val = vector of maximum values with length number of traits
generate_traitsEF <- function(n_species, min_val, max_val) {
    
    traits_matrix <- mapply(seq, from = min_val, to = max_val, length.out = n_species)
    row.names(traits_matrix) <- paste0("species", seq_len(n_species))
    colnames(traits_matrix) <- paste0("trait", seq_along(min_val))

    return(traits_matrix)

}

#################


patches  <- 25   # Number of patches
species  <- 25   # Number of species
time     <- 150	 # Length of model run (generations)
initpop  <- 150   # initial population size
n_traits <- 1    # number of traits

# Competition + Dispersion coefficients
A =  0.0005	# Alpha scalar #from oikos paper
d = 0.01 # Dispersal percentage

# Generate environment
env <- 1:patches
names(env) <- dimnames(composition)[[1]]
composition <- array(NA, dim = c(patches, species, time),
                     dimnames = list(paste0("patches", 1:patches),
                                     paste0("species", 1:species),
                                     paste0("time", 1:time))) #where values = N

composition[,,"time1"] <- initpop # N
# Generate traits
traits <- generate_traitsEF(species,
                          rep(min(env), n_traits),
                          rep(max(env), n_traits))

# Vector that indicates trait contribution:
#   "R"  = trait contributes to fitness,
#   "A"  = trait contributes to niche difference,
#   "N"  = trait contributes to none,
#   "RA" = trait contributes to BOTH fitness and niche difference
trait_type <- rep("RA", n_traits)#c("RA", "RA") ## # #
names(trait_type) <- colnames(traits)


# Actual simulation ------------------------------------------------------------

results <- multigen(traits = traits, trait_type = trait_type, env = env,
                    time = time, species = species, patches = patches,
                    composition = composition, A = A, d = d, k=1.15, c=0.25)#1.45 from oikos paper

# threshold number of individuals
final <- ifelse(results[,,time] < 2, 0, results[,, time])


plot(1:time, results[1, 1, 1:time], type = "n", ylim=c(0, 2000))
apply(results, 2, function(x){points(1:time, x[1, 1:time], type = "l", ylim=c(0,max(final)), col=rainbow(25)[sample(1:25)])})
plot(1:time, results[2, 1, 1:time], type = "n", ylim=c(0, 2000))
apply(results, 2, function(x){points(1:time, x[9, 1:time], type = "l", ylim=c(0,max(final)), col=rainbow(25)[sample(1:25)])})


# Compute FD on last community -------------------------------------------------

# Vector to consider available traits
#   TRUE  = trait is available to compute FD,
#   FALSE = traits is unavailable.
available_traits = rep(TRUE, n_traits)

# Get names of the species that are present at least in a single patch
present_species = colnames(final)[colSums(final) != 0]

if(any(rowSums(final[, present_species])==0)){
	finalFD <- final[-which(rowSums(final[, present_species])==0),]
}else{
	finalFD <- final}
	

# Compute FD metrics
final_FD <- dbFD(traits[present_species, available_traits],
                 finalFD[, present_species],
                 calc.FRic = TRUE, stand.FRic = TRUE,
                 scale.RaoQ = TRUE,
                 calc.FDiv = TRUE)


