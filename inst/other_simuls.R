# Packages ---------------------------------------------------------------------
suppressPackageStartupMessages({
    suppressWarnings({
        library("tidyverse")
        library("furrr")

        devtools::load_all()
    })
})

list_A = c(0, 10^-(c(4, 5, 6)))
list_k = c(1.2, 1.5)
list_B = c(0, 10^-(c(2, 3, 4)))
list_H = c(0, 0.5)
n_seed = 10
n_patches = 25
n_species = 100
n_gen = 50

weights = c(0, 50, 100)
scenar_df = expand.grid(R = weights, A = weights, H = weights) %>%
    as_tibble() %>%
    mutate(scenar_name   = paste0("R", R, "A", A, "H", H),
           trait_weights = purrr::pmap(list(R, A, H),
                                       ~create_trait_weights(..1, ..2, ..3, 2)))
scenar_list = scenar_df$trait_weights
names(scenar_list) = scenar_df$scenar_name

# Generate trait matrices for each seed
fixed_coef = 0.3
trait_seeds = lapply(seq(n_seed), function(seed) {
    set.seed(seed)
    uncor  = generate_cor_traits(n_patches, n_species, 1, cor_coef = 0)
    set.seed(seed)
    poscor = generate_cor_traits(n_patches, n_species, 1,
                                 cor_coef = fixed_coef)
    set.seed(seed)
    negcor = generate_cor_traits(n_patches, n_species, 1,
                                 cor_coef = -fixed_coef)
    list(poscor = poscor,
         uncor  = uncor,
         negcor = negcor)
})

# multidimensional matrix
composition <- array(NA, dim = c(n_patches, n_species, n_gen),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen))))

composition[,,1] = 50

# Actual simulations -----------------------------------------------------------
param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

tictoc::tic()

plan(multiprocess, workers = 14)
var_param = future_map(param_sets[1:320], function(x) {
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
               given_d = 0.05)
}, .progress = TRUE)
tictoc::toc()

saveRDS(var_param, "job_data/other_simuls.Rds")
# All Traits -------------------------------------------------------------------

full_trait_df = map_dfr(trait_seeds, function(x) {
    map_dfr(x, ~.x %>%
                as.data.frame() %>%
                rownames_to_column("species"),
            .id = "trait_cor")
},.id = "seed") %>%
    mutate(seed = as.integer(seed))

# Extract performances indices -------------------------------------------------

tictoc::tic()
var_perfs = future_map_dfr(unlist(var_param, recursive = FALSE),
                           ~extract_performances_from_simul(.x, trait_seeds,
                                                            TRUE),
                           .progress = TRUE)
tictoc::toc()

saveRDS(trait_seeds, "job_data/other_simuls_traits.Rds")
saveRDS(var_perfs, "job_data/other_simuls_perfs.Rds")
