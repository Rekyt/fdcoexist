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
list_H = seq(0, 1, length.out = 6)
n_seed = 30
n_patches = 25
n_species = 100
n_gen = 50


# Generate all scenarios
weights = round(seq(0, 100, length.out = 3), digits = 0)

scenar_df = expand.grid(R = weights, A = weights, H = weights) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(scenar_name   = paste0("R", R, "A", A, "H", H),
           trait_weights = purrr::pmap(list(R, A, H),
                                       ~create_trait_weights(..1, ..2, ..3)))
scenar_list = scenar_df$trait_weights
names(scenar_list) = scenar_df$scenar_name

# multidimensional matrix
composition <- array(NA, dim = c(n_patches, n_species, n_gen),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen))))

# Get trait list
used_trait_list = map(seq(n_seed), function(given_seed) {
    set.seed(given_seed)
    given_traits = generate_cor_traits(n_patches, n_species, n_traits - 1,
                                       cor_coef = 0)

    list(uncor = given_traits)
}) %>%
    set_names(nm = seq(n_seed))

full_trait_df = map_dfr(used_trait_list, ~.x$uncor %>%
                            as.data.frame() %>%
                            tibble::rownames_to_column("species"),
                        .id = "seed")

# Actual simulations -----------------------------------------------------------

n_slots = 64
job_task_id = 1

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    # Make all combinations but exclude cases where B > A
    cross(.filter = function(v, w, x, y, z) {y > w})

number_of_sets_per_task = 10010

# Return split sequence for a giving number of
f = function(a, b) {
    seq((a - 1) * b + 1, a * b, by = 1)
}

param_used = f(job_task_id, number_of_sets_per_task)

plan(multicore, workers = n_slots)

tictoc::tic()

var_param = future_lapply(param_sets[param_used], function(x) {
    suppressMessages({
        devtools::load_all()
    })

    meta_simul(seed_number = x$run_n,
               given_k = x$k,
               given_A = x$A,
               given_B = x$B,
               given_scenars = scenar_list,
               given_H = x$H,
               given_traits = NULL,
               given_h_fun = "sum",
               given_di_thresh = 24,
               given_env = 1:25,
               given_composition = composition,
               given_d = 0.05)
})
tictoc::toc()

# Extract Performances & CWM ---------------------------------------------------

cat("\nExtracting Performance from each simul\n")

tictoc::tic()
var_param_perfs = map_dfr(unlist(var_param, recursive = FALSE),
                          ~extract_performances_from_simul(.x, used_trait_list,
                                                           TRUE))
tictoc::toc()

# Save files -------------------------------------------------------------------

saveRDS(var_param, file = paste0("inst/job_data/var_param_bigmem_",
                                 min(param_used), "_", max(param_used),
                                 "_data.Rds"),
        compress = TRUE)
saveRDS(full_trait_df, file = "inst/job_data/bigmem_trait_df.Rds")
saveRDS(var_param_perfs, file = paste0("inst/job_data/bigmem_", min(param_used),
                                       "_", max(param_used),"perfs.Rds"))
