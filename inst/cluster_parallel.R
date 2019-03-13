# The aim of this script is to be used on the cluster to make simulations
library("purrr")
library("furrr")
suppressMessages({
    devtools::load_all()
})

# Parameters -------------------------------------------------------------------
list_A = c(0, 10^-(seq(1, 8, length.out = 6)))
list_k = seq(1, 1.5, length.out = 6)
list_B = list_A
list_H = seq(0, 1, length.out = 6)
n_seed = 50
n_patches = 25
n_species = 100
n_gen = 75


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

# Actual simulations -----------------------------------------------------------

cli_arguments = commandArgs(trailingOnly = TRUE)
job_task_id = cli_arguments[1]

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

number_of_jobs_per_task = 32

# Return split sequence for a giving number of
f = function(a, b) {
    seq((a - 1) * b + 1, a * b, by = 1)
}

param_used = f(job_task_id, number_of_jobs_per_task)


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
               given_traits = NULL,
               given_h_fun = "sum",
               given_di_thresh = 24,
               given_env = 1:25,
               given_composition = composition,
               given_d = 0.05)
})
tictoc::toc()

# Save files -------------------------------------------------------------------

saveRDS(var_param, file = paste0("~/projects/fdcoexist/inst/job_data/var_param_",
                                 job_task_id, ".Rds"),
        compress = TRUE)
