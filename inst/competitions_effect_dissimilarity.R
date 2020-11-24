# Relationship between trait difference and effect of exponents
# in limiting similarity and hierarchical competition
library("dplyr")
library("ggplot2")

pkgload::load_all()

# Exponent of hierarchical compoent --------------------------------------------
hierar_df = data.frame(
    diff_trait = seq(-2, 2, length.out = 5000)
) %>%
    mutate(hierar_0.5  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^0.5),
           hierar_1.0  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^1),
           hierar_2.0  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^2),
           hierar_0.25 = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^0.25)) %>%
    tidyr::pivot_longer(hierar_0.5:hierar_0.25, names_to = "hierar_exp",
                        values_to = "hierar_effect") %>%
    tidyr::extract(hierar_exp, c("hierar_exp"), regex = ".*_(.*)$")

plot_hierar_effect = ggplot(
    hierar_df, aes(diff_trait, hierar_effect, color = hierar_exp)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_line(size = 1) +
    labs(x = "Difference in trait\n(focal species - other species)",
         y = "Addition to basal growth",
         color = "Hier. Comp.\nexponent",
         title = "Hierarchical Competition effect") +
    theme_bw() +
    theme(aspect.ratio = 1)

# Exponent of limiting similarity ----------------------------------------------

limsim_df = data.frame(
    diff_trait = seq(-2, 2, length.out = 5000)) %>%
    mutate(limsim_0.5  = abs(diff_trait)^0.5,
           limsim_1.0  = abs(diff_trait)^1,
           limsim_2.0  = abs(diff_trait)^2,
           limsim_0.25 = abs(diff_trait)^0.25) %>%
    tidyr::pivot_longer(limsim_0.5:limsim_0.25, names_to = "limsim_exp",
                        values_to = "limsim_effect") %>%
    tidyr::extract(limsim_exp, c("limsim_exp"), regex = ".*_(.*)$")

plot_limsim_effect = ggplot(
    limsim_df, aes(diff_trait, limsim_effect, color = limsim_exp)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_line(size = 1) +
    labs(x = "Difference in trait\n(focal species - other species)",
         y = "Divider to basal growth",
         color = "Lim. sim.\nexponent",
         title = "Limiting similarity effect") +
    theme_bw() +
    theme(aspect.ratio = 1)

cowplot::plot_grid(plot_hierar_effect, plot_limsim_effect, ncol = 1,
                   align = "hv")

# Limiting similarity behaviour -----------------------------------------------
# Expected behaviour: growth rate of a species should decrease strongly if a
# coexisting species has similar trait values, this influence should decrease
# as we move furhter away the species trait value

# 3 species with 1 trait contributing to both k and A
# species 1 and 2 very close, 1 and 3 functionally distant
# H = 0, low and strong A

set.seed(1) # set random seed
n_patches <- 25 # number of patches
n_species <- 3 # number of species
n_gen <- 50 # number of time steps
n_traits <- 2 # number of traits
init_pop <- 50 # number of individuals at t=0 for each species
d <- 0.05 # dispersal parameter
width <- 5 # standard deviation of the Gaussian environmental filtering
K <- 100 # carrying capacity (per species per patch)
# Initial population matrix, 50 individuals of each species
composition <- array(NA, dim = c(n_patches, n_species, n_gen*20),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen*20))))
composition[, , 1] <- init_pop
# trait value and contribution
trait_scenar <- data.frame(trait = "trait1", growth_weight = 1,
                           compet_weight = 1, hierarchy_weight = 0)
uncor_traits <- data.frame(trait1 = c(4, 4.1, 7))
rownames(uncor_traits) <- paste0("species", seq(1:3))

# parameter values
list_k <- 2
list_A <- c(1e-2, 1e-4, 1e-6, 0)
list_H <- 0

# multigen
simul_list <- list()
for(i in 1:length(list_A)){
    simul <- multigen(
        traits = uncor_traits, trait_weights = trait_scenar, env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = list_A[i], B = 1e-7, d = d, k = list_k, H = list_H,
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100)

    simul_list[[i]] <- simul
}

env_resp <- r_env(simul, n_patches = 25, sp = 3, plot = TRUE)
env_resp[[2]]

env_resp1 <- r_env(simul_list[[1]], n_patches = 25, sp = 3, plot = TRUE)

cowplot::plot_grid(env_resp1[[2]],
                   env_resp[[2]])
