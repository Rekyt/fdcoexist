# Packages ---------------------------------------------------------------------

library(FD)

# Source other functions -------------------------------------------------------

source("R/multigen_function.r")


# Initial values ---------------------------------------------------------------
## Set number of patches, species, time
patches  <- 10   # Number of patches
species  <- 25   # Number of species
time     <- 150	 # Length of model run (generations)
initpop  <- 20   # initial population size
n_traits <- 1    # number of traits

composition <- array(NA, dim = c(patches, species, time),
                     dimnames = list(paste0("patches", 1:patches),
                                     paste0("species", 1:species),
                                     paste0("time", 1:time))) #where values = N

composition[,,"time1"] <- initpop # N

# Generate environment
env <- 1:patches
names(env) <- dimnames(composition)[[1]]

# Competition + Dispersion coefficients
A = 0.0005	# Alpha scalar
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
                    composition = composition, A = A, d = d, k=1.15, c=0.5)

# threshold number of individuals
final <- ifelse(results[,,time] < 2, 0, results[,, time])


# Plots ------------------------------------------------------------------------

plot(1:time, results[5, 1, 1:time], type = "l", ylim=c(0,150))
points(1:time, results[5, 4, 1:time], type = "l", col="red")
points(1:time, results[5, 8, 1:time], type = "l", col="blue")
points(1:time, results[5, 12, 1:time], type = "l", col="green")
points(1:time, results[5, 20, 1:time], type = "l", col="purple")
points(1:time, results[5, 25, 1:time], type = "l", col="grey")
#
plot(1:time, results[9, 1, 1:time], type = "l", ylim=c(0,150))
points(1:time, results[9, 4, 1:time], type = "l", col="red")
points(1:time, results[9, 8, 1:time], type = "l", col="blue")
points(1:time, results[9, 12, 1:time], type = "l", col="green")
points(1:time, results[9, 20, 1:time], type = "l", col="purple")
points(1:time, results[9, 25, 1:time], type = "l", col="grey")


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
