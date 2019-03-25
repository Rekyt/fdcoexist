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
n_traits = 2


# Generate all scenarios
weights = round(seq(0, 100, length.out = 3), digits = 0)

scenar_list = cross(list(R = weights, A = weights, H = weights)) %>%
    map(~create_trait_weights(.x$R, .x$A, .x$H, 2)) %>%
    # Name scenarios
    set_names(nm = cross(list(R = weights, A = weights, H = weights)) %>%
                  map_chr(~paste0("R", .x$R, "A", .x$A, "H", .x$H)))

# multidimensional matrix
composition = array(NA, dim = c(n_patches, n_species, n_gen),
                    dimnames = list(paste0("patches", seq(n_patches)),
                                    paste0("species", seq(n_species)),
                                    paste0("time", seq(n_gen))))

# Get trait list
used_trait_list = map(seq(n_seed), function(given_seed) {
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
}) %>%
    set_names(nm = seq(n_seed))

full_trait_df = map_dfr(used_trait_list,~.x %>%
                            map_dfr(function(x) {
                                x %>%
                                    as.data.frame() %>%
                                    tibble::rownames_to_column("species")
                            },
                            .id = "trait_cor"),
                        .id = "seed")

# Actual simulations -----------------------------------------------------------

n_slots = 64

cli_args = commandArgs(TRUE)
job_task_id = as.numeric(cli_args[1])

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    # Make all combinations but exclude cases where A > B
    cross(.filter = ~..2 > ..4)

number_of_sets_per_task = 1010

# Return split sequence for a giving number of
f = function(a, b) {
    seq((a - 1) * b + 1, a * b, by = 1)
}

param_used = f(job_task_id, number_of_sets_per_task)

if (max(param_used) > length(param_sets)) {
    param_used = seq(min(param_used), length(param_sets), by = 1)
}

plan(multicore)

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

# Save simuls
saveRDS(var_param, file = paste0("inst/job_data/var_param_bigmem_",
                                 min(param_used), "_", max(param_used),
                                 "_data.Rds"),
        compress = TRUE)

# Extract Performances & CWM ---------------------------------------------------

cat("\nExtracting Performance from each simul\n")

plan(multicore)

tictoc::tic()
var_param_perfs = future_lapply(unlist(var_param, recursive = FALSE),
                          function(x) {

                              suppressMessages({devtools::load_all})

                              extract_performances_from_simul(x, used_trait_list,
                                                           TRUE)
                          })
tictoc::toc()

# Save files -------------------------------------------------------------------

saveRDS(full_trait_df, file = "inst/job_data/bigmem_trait_df.Rds")
saveRDS(var_param_perfs, file = paste0("inst/job_data/bigmem_", min(param_used),
                                       "_", max(param_used),"perfs.Rds"),
        compress = TRUE)
