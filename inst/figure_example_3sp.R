
library("dplyr")
library("ggplot2")
library("cowplot")
library("tidyr")

pkgload::load_all()

# 3 species with 1 trait contributing to both k and A
# species 1 and 2 very close, 1 and 3 functionally distant
# H = 0, low and strong A

## Model parameters -----------------------------------------------------------
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

## Only Lim.Sim. --------------------------------------------------------------
# parameter values
list_k <- 2
list_A <- c(5, 2, 1, 0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 0)
list_H <- 0

# multigen
simul_list <- list()
for(i in 1:length(list_A)){
    simul <- multigen(
        traits = uncor_traits, trait_weights = trait_scenar, env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = list_A[i], B = 1e-7, d = d, k = list_k, H = list_H,
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
        lim_sim_exponent = 1)

    simul_list[[i]] <- simul
}

# plot_patch and environmental response curves
env_resp <- r_env(simul, n_patches = 25, sp = 3, plot = TRUE)
env_resp_plot <- env_resp[[2]] +
    scale_color_brewer(palette = "Dark2") +
    guides(color = FALSE) +
    geom_vline(xintercept = 4, linetype = "dashed")

plotA_3sp <- plot_grid(
    plot_patch(simul_list[[5]]$compo, patch = 4, time =  50, equilibrium = FALSE) +
        scale_color_brewer(palette = "Dark2") +
        guides(color = FALSE) +
        geom_label(aes(label = ifelse(time == 50,
                                      gsub("species", "sp", as.character(sp)),
                                      NA))) +
        labs(title = paste0("Patch 4; A = ", list_A[[5]])) +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()),
    plot_patch(simul_list[[6]]$compo, patch = 4, time =  50, equilibrium = FALSE) +
        scale_color_brewer(palette = "Dark2") +
        guides(color = FALSE) +
        geom_label(aes(label = ifelse(time == 50,
                                      gsub("species", "sp", as.character(sp)),
                                      NA))) +
        labs(title = paste0("Patch 4; A = ", list_A[[6]])) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()),
    nrow = 1, rel_widths = c(1, 0.95))
plotA_3sp

# ## Only Hierarch. -----------------------------------------------------------
list_H <- c(1e-2, 1e-4)
list_A <- 0
trait_scenar <- data.frame(trait = "trait1", growth_weight = 1,
                           compet_weight = 0, hierarchy_weight = 1)


simul_list <- list()
for(i in 1:length(list_H)){
    simul <- multigen(
        traits = uncor_traits, trait_weights = trait_scenar, env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = list_A, B = 1e-7, d = d, k = list_k, H = list_H[i],
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
        lim_sim_exponent = 1)

    simul_list[[i]] <- simul
}

plotH_3sp <- plot_grid(
    plot_patch(simul_list[[1]]$compo, patch = 4, time =  50, equilibrium = FALSE) +
        scale_color_brewer(palette = "Dark2") +
        guides(color = FALSE) +
        geom_label(aes(label = ifelse(time == 50,
                                      gsub("species", "sp", as.character(sp)),
                                      NA))) +
        labs(title = paste0("Patch 4; H = ", list_H[[1]])),
    plot_patch(simul_list[[2]]$compo, patch = 4, time =  50, equilibrium = FALSE) +
        scale_color_brewer(palette = "Dark2") +
        guides(color = FALSE) +
        geom_label(aes(label = ifelse(time == 50,
                                      gsub("species", "sp", as.character(sp)),
                                      NA))) +
        labs(title = paste0("Patch 4; H = ", list_H[[2]])) +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank()),
    nrow = 1, rel_widths = c(1, 0.95))
plotH_3sp

## Plot -----------------------------------------------------------------------
plot_grid(env_resp_plot, plotA_3sp, plotH_3sp, nrow = 3,
          labels = c("a", "b", "c"))
