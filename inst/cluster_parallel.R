# The aim of this script is to be used on the cluster to make simulations
library("purrr")
library("snow")
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

tictoc::tic()

param_sets = list(
    run_n = seq(n_seed),
    A     = list_A,
    k     = list_k,
    B     = list_B,
    H     = list_H) %>%
    cross()

#lecture du fichier contenant la liste des machines
cli_arguments = commandArgs(TRUE)
pe_file = cli_arguments[1]

#construction de la liste des slots pour le "cluster"
node_list = read.table(pe_file, sep = " ", header = FALSE,
                       stringsAsFactors = FALSE)
node_names = node_list[,1]
slots_num  = node_list[,2]

workers = rep(node_names, slots_nm)

#We will run in parallel mode (socket) with
cl = makeSOCKcluster(workers)

var_param = parLapply(cl, param_sets, function(x) {
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
stopCluster(cl)
tictoc::toc()

# Trait data.frame -------------------------------------------------------------
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

# Save files -------------------------------------------------------------------

saveRDS(var_param, file = "~/projects/fdcoexist/var_param.Rds",
        compress = TRUE)
saveRDS(full_trait_df, file = "~/projects/fdcoexist/full_trait_df.Rds",
        compress = TRUE)
