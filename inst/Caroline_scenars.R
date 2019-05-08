# The aim of this script is to be used on the cluster to make simulations
library("tidyverse")
library("furrr")
library("future.apply")
suppressMessages({
    devtools::load_all()
})

# Parameters -------------------------------------------------------------------
list_A = c(0, 10^-6)
list_k = c(1.2)
list_B = c(0, 10^-4)
list_H = c(0, 1)
n_seed = 5
n_patches = 25
n_species = 100
n_gen = 50


# Generate all scenarios
single_trait = data.frame(trait = paste0("trait", 1:2),
                          growth_weight    = c(0, 1),
                          compet_weight    = c(0, 1),
                          hierarchy_weight = c(0, 1))

two_traits_R = data.frame(trait = paste0("trait", 1:2),
                          growth_weight    = c(0.5, 0.5),
                          compet_weight    = c(0,   0),
                          hierarchy_weight = c(0,   0))

two_traits_each = data.frame(trait = paste0("trait", 1:2),
                             growth_weight    = c(1, 0),
                             compet_weight    = c(0, 1),
                             hierarchy_weight = c(0, 1))

two_traits_all  = data.frame(trait = paste0("trait", 1:2),
                             growth_weight    = c(0.5, 0.5),
                             compet_weight    = c(0.5, 0.5),
                             hierarchy_weight = c(0.5, 0.5))
scenar_list = list(R100A100H100 = single_trait,
                   R50A0H0      = two_traits_R,
                   R0A100H100   = two_traits_each,
                   R50A50H50    = two_traits_all)

# Generate trait matrices for each seed
fixed_coef = 0.3
trait_seeds = lapply(seq(n_seed), function(seed) {
    set.seed(seed)
    list(
        uncor  = generate_cor_traits(n_patches, n_species, 1, cor_coef = 0),
        poscor = generate_cor_traits(n_patches, n_species, 1,
                                     cor_coef = fixed_coef),
        negcor = generate_cor_traits(n_patches, n_species, 1,
                                     cor_coef = -fixed_coef)
    )
})

# multidimensional matrix
composition <- array(NA, dim = c(n_patches, n_species, n_gen),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen))))

composition[,,1] = 50

# Actual simulations -----------------------------------------------------------

n_slots = 64

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

plan(multicore, workers = 8)

tictoc::tic()
var_param = future_map(param_sets, function(x) {
    suppressMessages({
        devtools::load_all()
    })

    meta_simul(seed_number = x$run_n,
               given_k = x$k,
               given_A = x$A,
               given_B = x$B,
               given_scenars = scenar_list,
               given_H = x$H,
               given_traits = trait_seeds[[x$run_n]],
               given_h_fun = "sum",
               given_di_thresh = 24,
               given_env = 1:25,
               given_composition = composition,
               given_d = 0.05,
               given_env_width = 2)
}, .progress = TRUE)
tictoc::toc()

# Extract Performances & CWM ---------------------------------------------------

cat("\nExtracting Performance from each simul\n")

tictoc::tic()
caro_perfs = future_map_dfr(unlist(var_param, recursive = FALSE),
                          ~extract_performances_from_simul(.x, trait_seeds,
                                                           TRUE),
                          .progress = TRUE)
tictoc::toc()

full_trait_df = map_dfr(trait_seeds, ~.x %>%
            map_dfr(function(y) {
                y %>%
                    as.data.frame() %>%
                    rownames_to_column("species")
            }, .id = "trait_cor"), .id = "seed")

var_perfs = caro_perfs %>%
    mutate(seed = as.character(seed)) %>%
    inner_join(full_trait_df, by = c("species", "seed", "trait_cor"))

var_perfs %>%
    filter(H == 0, A == 0, seed == 1, trait_cor == "uncor", R_scenar == 100, A_scenar == 100, H_scenar == 100)
# Save files -------------------------------------------------------------------

saveRDS(var_param, file = paste0("inst/job_data/",
                                 "caroline_scenars.Rds"),
        compress = TRUE)
saveRDS(trait_seeds, file = paste0("inst/job_data/caroline_traits.Rds"))
saveRDS(caro_perfs, file = "inst/job_data/caroline_perfs.Rds", compress = TRUE)


