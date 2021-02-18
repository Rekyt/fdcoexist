
library("dplyr")
library("ggplot2")
library("cowplot")
library("ggrepel")
library("RColorBrewer")
library("tidyr")
library("gridExtra")

pkgload::load_all()

# A increasing, looking at trait variance

## Model parameters -----------------------------------------------------------
set.seed(1) # set random seed
n_patches <- 25 # number of patches
n_species <- 50 # number of species
n_gen <- 50 # number of time steps
n_traits <- 2 # number of traits
init_pop <- 50 # number of individuals at t=0 for each species
d <- 0.05 # dispersal parameter
width <- 5 # standard deviation of the Gaussian environmental filtering
K <- 100 # carrying capacity (per species per patch)
# Initial population matrix, 50 individuals of each species
composition <- array(NA, dim = c(n_patches, n_species, n_gen),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen))))
composition[, , 1] <- init_pop
# trait values
uncor_traits <- generate_cor_traits(n_patches, n_species, n_traits - 1,
                                    cor_coef = 0)

# parameter values: k not varying through the three examples
list_k <- 2

## Increasing A ---------------------------------------------------------------
scenarios <- list(
    scenar1 = list(list_A = 1e-8,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),

    scenar2 = list(list_A = 1e-7,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),

    scenar3 = list(list_A = 1e-6,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),
    scenar4 = list(list_A = 1e-5,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),
    scenar5 = list(list_A = 1e-4,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),
    scenar6 = list(list_A = 1e-3,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),
    scenar7 = list(list_A = 0.01,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0))
)

# multigen
simul_list <- list()
for(i in 1:length(scenarios)){
    simul <- multigen(
        traits = uncor_traits,
        trait_weights = scenarios[[i]]$trait_scenar,
        env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = scenarios[[i]]$list_A, H = scenarios[[i]]$list_H,
        B = 1e-7, d = d, k = list_k,
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
        lim_sim_exponent = 1)

    simul_list[[i]] <- simul
}

## Plots ----------------------------------------------------------------------
# Matching t50 with traits per scenario and per site
res <- c()
sp_tra <- data.frame(uncor_traits)
sp_tra$sp <- rownames(sp_tra)
A_scenar <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 0.01)

for(i in 1:length(simul_list)){
    # Abundances at t50
    tmp <- simul_list[[i]]$compo[, , 50]
    # Convert to data.frame
    tmp <- reshape2::melt(tmp)
    colnames(tmp) <- c("site", "sp", "ab")
    # Merging trait values
    tmp <- left_join(tmp, sp_tra, by = "sp")
    # Add A value
    tmp$A <- A_scenar[i]

    # Removing absent species
    tmp <- tmp[which(tmp$ab > 0), ]

    # Weighted variance by site
    tmp <- tmp %>%
        group_by(site) %>%
        mutate(var_w = wtd_var(trait1, mu = ab),
               var_w2 = wtd_var(trait2, mu = ab)) %>%
        ungroup() %>%
        as.data.frame()

    # Binding results
    res <- rbind(res, tmp)
}

# Plot
res$site_plot <- as.integer(gsub("patches", "", res$site))

cwv_A <-
    ggplot(res[!duplicated(res[, c("A", "site")]), ],
           aes(site_plot, var_w)) +
    geom_point(aes(fill = as.factor(A)), shape = 21, size = 2) +
    geom_line(aes(color = as.factor(A))) +
    scale_color_viridis_d("Limiting similarity (A)") +
    scale_fill_viridis_d("Limiting similarity (A)") +
    labs(x = "Site", y = "Weighted trait variance") +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA))

# Saving figure
ggsave2("inst/figures/paper_supp_fig4_A_cwv.png",
        cwv_A,
        width = 16.6, height = 16.6,
        units = "cm", dpi = 300)

