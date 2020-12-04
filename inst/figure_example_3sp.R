
library("dplyr")
library("ggplot2")
library("cowplot")
library("tidyr")
library("gridExtra")

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
# trait values
uncor_traits <- data.frame(trait1 = c(4, 4.1, 7))
rownames(uncor_traits) <- paste0("species", seq(1:3))

# parameter values: k not varying through the three examples
list_k <- 2

## Three scenarios ------------------------------------------------------------
# List of parameter values and trait contribution
three_scenarios <- list(
    scenar1 = list(list_A = 0,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 0, hierarchy_weight = 0)),
    scenar2 = list(list_A = 0.02,
                   list_H = 0,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 1, hierarchy_weight = 0)),
    scenar3 = list(list_A = 0,
                   list_H = 2e-2,
                   trait_scenar = data.frame(
                       trait = "trait1", growth_weight = 1,
                       compet_weight = 0, hierarchy_weight = 1))
)

# multigen
simul_list <- list()
for(i in 1:length(three_scenarios)){
    simul <- multigen(
        traits = uncor_traits,
        trait_weights = three_scenarios[[i]]$trait_scenar,
        env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = three_scenarios[[i]]$list_A, H = three_scenarios[[i]]$list_H,
        B = 1e-7, d = d, k = list_k,
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
        lim_sim_exponent = 1)

    simul_list[[i]] <- simul
}

## Plot -----------------------------------------------------------------------
# explicit legend (strong lim.sim (A = bla); strong Hierarch.)

# plot_patch and environmental response curves
env_resp <- r_env(simul, n_patches = 25, sp = 3, plot = TRUE)
env_resp_plot <- env_resp[[2]] +
    scale_color_brewer(palette = "Dark2") +
    guides(color = FALSE) +
    geom_vline(xintercept = 4, linetype = "dashed")

# Trait table
table_tra_theme <- ttheme_minimal(
    core = list(fg_params = list(cex = 1),
                bg_params = list(fill = c("grey90", "grey89", "grey70"),
                                 col = NA)), # blues9[1:3]
    colhead = list(fg_params = list(cex = 2)),
    rowhead = list(fg_params = list(cex = 2)))

sp_tra <- uncor_traits
colnames(sp_tra) <- "Trait"
rownames(sp_tra) <- paste("Sp.", seq(1:3))

tra_tbl <- tableGrob(sp_tra, theme = table_tra_theme)


# Species dynamics (no competition)
plotk_3sp <- plot_patch(simul_list[[1]]$compo, patch = 4, time =  20,
                        equilibrium = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 50)) +
    guides(color = FALSE) +
    geom_label(aes(label = ifelse(time == 20,
                                  gsub("species", "sp", as.character(sp)),
                                  NA))) +
    labs(title = paste0("Patch 4; No competition"))

# Species dynamics (lim.sim.)
plotA_3sp <- plot_patch(simul_list[[2]]$compo, patch = 4, time =  20,
                        equilibrium = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 50)) +
    guides(color = FALSE) +
    geom_label(aes(label = ifelse(time == 20,
                                  gsub("species", "sp", as.character(sp)),
                                  NA))) +
    labs(title = paste0("Patch 4; Limiting similarity")) +
    # \n(A=", three_scenarios[[2]]$list_A, ")")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

# Species dynamics (hierarchical competition)
plotH_3sp <- plot_patch(simul_list[[3]]$compo, patch = 4, time =  20,
                        equilibrium = FALSE) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0, 50)) +
    guides(color = FALSE) +
    geom_label(aes(label = ifelse(time == 20,
                                  gsub("species", "sp", as.character(sp)),
                                  NA))) +
    labs(title = paste0("Patch 4; Hierarchical competition")) +
    # (H=", three_scenarios[[3]]$list_H, ")")) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

# All plots together
paper_figure3spexample <-
    plot_grid(plot_grid(grid.arrange(tra_tbl, env_resp_plot, nrow = 1,
                                     as.table = TRUE, widths = c(0.5, 1))),
              plot_grid(plotk_3sp, plotA_3sp, plotH_3sp, nrow = 1,
                        rel_widths = c(1, 0.9, 0.9)),
              labels = c("a)", "b)"), nrow = 2)
paper_figure3spexample

# Saving figure
ggsave2("inst/figures/paper_figure_3spexample.png", paper_figure3spexample,
        width = 18, height = 10,
        units = "cm", dpi = 300, scale = 1.5)


