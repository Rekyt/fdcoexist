# The aim of this script is to be used on the cluster to make simulations
library("purrr")
library("future.apply")
suppressMessages({
    devtools::load_all()
})

# Parameters -------------------------------------------------------------------
list_A = c(0, 10^-(seq(1, 8, length.out = 6)))
list_k = seq(1, 1.5, length.out = 6)
list_B = list_A
list_H = seq(0, 0.2, length.out = 6)
n_seed = 30
n_patches = 25
n_species = 100
n_gen = 50
n_traits = 2
init_pop = 50

# Traits & Contribution scenarios ----------------------------------------------
# Generate all trait contribution scenarios
weights = round(seq(0, 100, length.out = 3), digits = 0)

scenar_list = cross(list(R = weights, A = weights, H = weights)) %>%
    map(~create_trait_weights(.x$R, .x$A, .x$H, 2)) %>%
    # Name scenarios
    set_names(nm = cross(list(R = weights, A = weights, H = weights)) %>%
                  map_chr(~paste0("R", .x$R, "A", .x$A, "H", .x$H)))

# Initial population matrix
composition = array(NA, dim = c(n_patches, n_species, n_gen),
                    dimnames = list(paste0("patches", seq(n_patches)),
                                    paste0("species", seq(n_species)),
                                    paste0("time", seq(n_gen))))
composition[,,1] = init_pop

# Generate sets of traits
used_trait_list = lapply(seq(n_seed), function(given_seed) {
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
names(used_trait_list) = seq(n_seed)

full_trait_df = map_dfr(used_trait_list,~.x %>%
                            map_dfr(function(x) {
                                x %>%
                                    as.data.frame() %>%
                                    tibble::rownames_to_column("species")
                            },
                            .id = "trait_cor"),
                        .id = "seed")

# Sets of all parameters
param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

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
               given_traits = used_trait_list[[x$run_n]],
               given_h_fun = "+",
               given_di_thresh = 24,
               given_env = 1:25,
               given_composition = composition,
               given_d = 0.05,
               given_env_width = 2)

    simul_perf = map_dfr(simul_list, function(y) {
        extract_performances_from_simul(y, used_trait_list[[x$run_n]], TRUE)
    })
    saveRDS(simul_perf, paste0("inst/job_data/perf_df_", given_n, ".Rds"), compress = TRUE)
})
tictoc::toc()

# Save files -------------------------------------------------------------------
# Save Trait data.frame
saveRDS(full_trait_df, file = "inst/job_data/bigmem_trait_df.Rds")
# Save performance extracted from simulations
# saveRDS(var_param, file = paste0("inst/job_data/perf_list_",
#                                 gsub("-", "_", Sys.Date()), ".Rds"),
#        compress = TRUE)
