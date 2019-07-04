# Simple scenarios to test the influence of hierarchical competition on the
# dynamics of simple communities
# Packages ---------------------------------------------------------------------
devtools::load_all()

# Simulation parameters --------------------------------------------------------
A = 1e-5
k = 1.25
B = 1e-2
H = 1e-2
n_patches = 5
n_species = 4
n_gen = 50
n_traits = 2


# Initial population matrix
composition = array(NA, dim = c(n_patches, n_species, n_gen),
                    dimnames = list(paste0("patches", seq(n_patches)),
                                    paste0("species", seq(n_species)),
                                    paste0("time", seq(n_gen))))
composition[1,,1] = c(50, 50,  0,  0)
composition[2,,1] = c( 0,  0, 50, 50)
composition[3,,1] = c(50,  0, 50, 50)
composition[4,,1] = c( 0, 50, 50, 50)
composition[5,,1] = c(50, 50, 50, 50)

# Traits & Contribution scenarios ----------------------------------------------
trait_weights = create_trait_weights(0, 0, 100, 2)


trait_mat = matrix(c(rep(1, 4), c(1, 4, 15, 16)), ncol = 2,
                   dimnames = list(paste0("species", 1:4),
                                   paste0("trait", 1:2)))

# Actual simulation ------------------------------------------------------------
simul = meta_simul(seed_number = 1,
           given_k = k,
           given_A = A,
           given_B = B,
           given_scenars = list(R0A0H100 = trait_weights),
           given_H = H,
           given_traits = list(uncor = trait_mat),
           given_h_fun = "+",
           given_di_thresh = 24,
           given_env = rep(1, 5),
           given_composition = composition,
           given_d = 0.05,
           given_env_width = 2)

# Compare plots ----------------------------------------------------------------

cowplot::plot_grid(plot_patch(simul[[1]]$compo[[1]], "patches2", 50) +
                       theme_bw(),
                   plot_patch(simul[[1]]$compo[[1]], "patches3", 50) +
                       theme_bw(),
                   plot_patch(simul[[1]]$compo[[1]], "patches4", 50) +
                       theme_bw(),
                   plot_patch(simul[[1]]$compo[[1]], "patches5", 50) +
                       theme_bw(),
                   ncol = 2, nrow = 2)
