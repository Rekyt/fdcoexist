# The aim of this script is to be used on the cluster to make simulations
library("purrr")
library("future.apply")
suppressMessages({
    devtools::load_all()
})

# Parameters -------------------------------------------------------------------
main_folder = "inst/job_data/perf_ef503c/"

list_A = c(0, 10^-(seq(1, 8, length.out = 6)))[c(1, 6)]
list_k = c(1.2, 1.3)
list_B = list_A
list_H = c(0, 10^-(seq(4, 5.5, length.out = 5)))
n_seed = 15
n_patches = 25
n_species = 100
n_gen = 50
n_traits = 2
init_pop = 50

set.seed(20190619)
seed_list = sample(1e6, size = n_seed)

# Sets of all parameters
param_sets = list(
    run_n = seed_list,
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

# Initial population matrix
composition = array(NA, dim = c(n_patches, n_species, n_gen),
                    dimnames = list(paste0("patches", seq(n_patches)),
                                    paste0("species", seq(n_species)),
                                    paste0("time", seq(n_gen))))
composition[,,1] = init_pop

# Traits & Contribution scenarios ----------------------------------------------
# Generate all trait contribution scenarios
weights = round(seq(0, 100, length.out = 3), digits = 0)

scenar_list = cross(list(R = weights, A = weights, H = weights)) %>%
    map(~create_trait_weights(.x$R, .x$A, .x$H, 2)) %>%
    # Name scenarios
    set_names(nm = cross(list(R = weights, A = weights, H = weights)) %>%
                  map_chr(~paste0("R", .x$R, "A", .x$A, "H", .x$H)))

# Generate sets of traits
used_trait_list = lapply(seed_list, function(given_seed) {
    set.seed(given_seed)
    given_traits = generate_cor_traits(n_patches, n_species, n_traits - 1,
                                       cor_coef = 0)
    set.seed(given_seed)
    cor_trait = given_traits = generate_cor_traits(n_patches, n_species,
                                                   n_traits - 1,
                                                   cor_coef = 0.3)
    set.seed(given_seed)
    negcor_trait = given_traits = generate_cor_traits(n_patches, n_species,
                                                      n_traits - 1,
                                                      cor_coef = -0.3)

    list(uncor  = given_traits,
         poscor = cor_trait,
         negcor = negcor_trait)
})
names(used_trait_list) = seed_list

full_trait_df = map_dfr(used_trait_list,~.x %>%
                            map_dfr(function(x) {
                                x %>%
                                    as.data.frame() %>%
                                    tibble::rownames_to_column("species")
                            },
                            .id = "trait_cor"),
                        .id = "seed")

saveRDS(full_trait_df, file = paste0(main_folder, "bigmem_trait_df.Rds"))

# Actual simulations -----------------------------------------------------------
plan(multiprocess, workers = future::availableCores())

tictoc::tic()
var_param = future_lapply(seq_along(param_sets), function(given_n) {
    suppressMessages({
        devtools::load_all()
    })

    cat("set: ", given_n, "\n")

    x = param_sets[[given_n]]

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

    simul_perf = map_dfr(simul_list, function(y) {
        extract_performances_from_simul(
            y, used_trait_list[[as.character(x$run_n)]], TRUE
    )
    })
    saveRDS(simul_perf,
            paste0(main_folder, "perf_df_", given_n, ".Rds"),
            compress = TRUE)
})
tictoc::toc()
