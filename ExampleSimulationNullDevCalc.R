# Packages ---------------------------------------------------------------------

library(FD)

# Source other functions -------------------------------------------------------

source("R/multigen_function.r")


# Initial values ---------------------------------------------------------------
## Set number of patches, species, time
patches  <- 10   # Number of patches
species  <- 25   # Number of species
time     <- 150	 # Length of model run (generations)
initpop  <- 25   # initial population size
n_traits <- 3    # number of traits

composition <- array(NA, dim = c(patches, species, time),
                     dimnames = list(paste0("patches", 1:patches),
                                     paste0("species", 1:species),
                                     paste0("time", 1:time))) #where values = N

composition[,,"time1"] <- initpop # N

# Generate environment
env <- 1:patches
names(env) <- dimnames(composition)[[1]]

# Competition + Dispersion coefficients
A = 0.001	# Alpha scalar
d = 0.05 # Dispersal percentage

# Generate tarits
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



# Actual simulation ------------------------------------------------------------

results <- multigen(traits = traits, trait_type = trait_type, env = env,
                    time = time, species = species, patches = patches,
                    composition = composition, A = A, d = d)

# threshold number of individuals
final <- ifelse(results[,,time] < 2, 0, results[,, time])


# Plots ------------------------------------------------------------------------

plot(1:time, composition[5,13,], type = "l")
plot(1:time, composition[5, 1,], type = "l")


# Compute FD on last community -------------------------------------------------

# Vector to consider available traits
#   TRUE  = trait is available to compute FD,
#   FALSE = traits is unavailable.
available_traits = rep(TRUE, n_traits)

# Get names of the species that are present at least in a single patch
present_species = colnames(final)[colSums(final) != 0]

# Compute FD metrics
final_FD <- dbFD(traits[present_species, available_traits],
                 final[, present_species],
                 calc.FRic = TRUE, stand.FRic = TRUE,
                 scale.RaoQ = TRUE,
                 calc.FDiv = TRUE)
