# The aim of this script is to be used on the cluster to make simulations
library("purrr")
library("future.apply")
suppressMessages({
    pkgload::load_all()
})


# Parameters -------------------------------------------------------------------
main_folder = "inst/job_data/perf_225384e/"

# A and H values to get a reduction of 20%, 40%, 60%, and 80% of growth compared
# to when B = 0 and they are respectively equal to 0
A_for_k_1.15 = c(0, 2.19e-4, 2.257e-4,  2.44e-4, 2.58e-4)
A_for_k_1.3  = c(0,   1e-7,     5e-7,   1.9e-6, 6.92e-6)
A_for_k_1.45 = c(0, 2.32e-9,  2.36e-8, 2.155e-7, 1.75e-6)
list_k = c(1.15, 1.3, 1.45)
list_B = c(0, 1.585e-4, 3.17e-4)
H_for_k_1.15 = c(0, 6.215e-4, 6.48e-4, 6.68e-4, 6.754e-4)
H_for_k_1.3  = c(0,     1e-6,    4e-6, 1.35e-5, 4.05e-5)
H_for_k_1.45 = c(0,  2.31e-8, 2.31e-7, 1.1e-6,  3e-6)
n_seed = 15
n_patches = 25
n_species = 100
n_gen = 50
n_traits = 2
init_pop = 50

set.seed(20190619)
seed_list = sample(1e6, size = n_seed)

# Functions --------------------------------------------------------------------
generate_param_comb = function(
   given_k, given_A, given_H, given_B = list_B, given_seed = seed_list) {
    list(
        run_n = given_seed,
        A     = given_A,
        k     = given_k,
        B     = given_B,
        H     = given_H) %>%
        purrr::cross()
}
# Sets of all parameters for each value
param_1.15 = generate_param_comb(
    given_A = A_for_k_1.15,
    given_k = 1.15,
    given_H = H_for_k_1.15)

param_1.3 = generate_param_comb(
    given_A = A_for_k_1.3,
    given_k = 1.3,
    given_H = H_for_k_1.3)

param_1.45 = generate_param_comb(
    given_A = A_for_k_1.45,
    given_k = 1.45,
    given_H = H_for_k_1.45)

param_sets = list(param_1.15, param_1.3, param_1.45) %>%
    purrr::flatten()

# Initial population matrix
composition = array(NA, dim = c(n_patches, n_species, n_gen),
                    dimnames = list(paste0("patches", seq(n_patches)),
                                    paste0("species", seq(n_species)),
                                    paste0("time", seq(n_gen))))
composition[,,1] = init_pop

# Traits & Contribution scenarios ----------------------------------------------
# Generate all trait contribution scenarios (trait contribution to growth,
# limiting similarity and hierarchical competition)
weights = round(seq(0, 100, length.out = 3), digits = 0)

scenar_list = cross(list(R = weights, A = weights, H = weights)) %>%
    map(~create_trait_weights(.x$R, .x$A, .x$H, 2)) %>%
    # Name scenarios
    set_names(nm = cross(list(R = weights, A = weights, H = weights)) %>%
                  map_chr(~paste0("R", .x$R, "A", .x$A, "H", .x$H)))

# Generate sets of traits fir each seed
used_trait_list = lapply(seed_list, function(given_seed) {
    set.seed(given_seed)
    uncor_traits = generate_cor_traits(n_patches, n_species, n_traits - 1,
                                       cor_coef = 0)
    set.seed(given_seed)
    poscor_traits = generate_cor_traits(n_patches, n_species,
                                        n_traits - 1,
                                        cor_coef = 0.3)
    set.seed(given_seed)
    negcor_traits = generate_cor_traits(n_patches, n_species,
                                        n_traits - 1,
                                        cor_coef = -0.3)

    # list(uncor  = uncor_traits,
    #      poscor = poscor_traits,
    #      negcor = negcor_traits)

    list(uncor  = uncor_traits)
})
names(used_trait_list) = seed_list

# Combine traits in a single data.frame
full_trait_df = map_dfr(used_trait_list,~.x %>%
                            map_dfr(function(x) {
                                x %>%
                                    as.data.frame() %>%
                                    tibble::rownames_to_column("species")
                            },
                            .id = "trait_cor"),
                        .id = "seed")

# Save Species Trait Information
saveRDS(full_trait_df, file = paste0(main_folder, "bigmem_trait_df.Rds"))

# Actual simulations -----------------------------------------------------------

# Setup parallelization scheme (all cores minus one)
plan(multisession, workers = future::availableCores() - 1)

# Actually run simulations (tictoc measures execution time)
tictoc::tic()
var_param = furrr::future_walk(seq_along(param_sets), function(given_n) {
    suppressMessages({
        pkgload::load_all()
    })

    cat("set: ", given_n, "\n")

    x = param_sets[[given_n]]

    # Proper simulation
    simul_list = meta_simul(seed_number = x$run_n,
               given_k = x$k,
               given_A = x$A,
               given_B = x$B,
               given_scenars = scenar_list,
               given_H = x$H,
               given_traits = used_trait_list[[as.character(x$run_n)]],
               given_h_fun = "+",
               given_di_thresh = 24,
               given_env = 1:25,
               given_composition = composition,
               given_d = 0.05,
               given_env_width = 2)

    # Save raw data file
    saveRDS(simul_list,
            paste0(main_folder, "simul_cat_", given_n, ".Rds"),
            compress = TRUE)

    # Extract species performance estimates
    simul_perf = map_dfr(simul_list, function(y) {
        extract_performances_from_simul(
            y, used_trait_list[[as.character(x$run_n)]], TRUE
    )
    })

    # Save species performance files
    saveRDS(simul_perf,
            paste0(main_folder, "perf_df_", given_n, ".Rds"),
            compress = TRUE)
}, .progress = TRUE)
tictoc::toc()
