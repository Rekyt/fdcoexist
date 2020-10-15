# FDcoexist June 4 2020
# Authors: C. Tucker & M. Grenié
# Generate species level responses to site conditions. Only moderate value of A,
# H or both used.
# replicates of each required - currently ~100
# For plots 2 & 3 (update with st err)

# Load Packages ----------------------------------------------------------------
library(dplyr)
library(cowplot)
library(tidyr)
library(RColorBrewer)
library(sads)
devtools::load_all()  # Load all functions from package in local repo


# Functions --------------------------------------------------------------------
label_hierar_exp = function(string) {

  lapply(unlist(string), function(x) {
    paste0("hierarchical~distance^{", x, "}")
  })
}

# Initial Parameters -----------------------------------------------------------
set.seed(20201008) # set random seed
n_patches <- 25    # number of patches
n_species <- 50    # number of species
n_gen <- 100       # number of time steps
n_traits <- 2      # number of traits
init_pop <- 50     # number of individuals at t=0 for each species
d <- 0.05          # dispersal parameter
width <- 5         # standard deviation of the Gaussian environmental filtering
K <- 100           # carrying capacity (per species per patch)

# Initial population matrix, 50 individuals of each species
composition <- array(NA, dim = c(n_patches, n_species, n_gen),
                     dimnames = list(paste0("patches", seq(n_patches)),
                                     paste0("species", seq(n_species)),
                                     paste0("time", seq(n_gen))))

composition[, , 1] <- init_pop


# Trait scenario with one trait contributing to all processes
trait_scenar <- data.frame(trait            = "trait1",
                           growth_weight    = 1,
                           compet_weight    = 1,
                           hierarchy_weight = 1)

# Get random trait values per species within [0,25] range
uncor_traits <- generate_cor_traits(n_patches, n_species, n_traits - 1,
                                    cor_coef = 0)
uncor_traits <- uncor_traits[, 1, drop = FALSE]

# Number of combination of traits
trait_comb <- seq(100)


# Generate truly random trait distributions (trait 1 varies, instead of being
# similar for all sets)
for (i in trait_comb[-1]) {
    uncor_traits2 <- generate_cor_traits_rand(n_patches, n_species,
                                              n_traits - 1, cor_coef = 0)
    uncor_traits2 <- uncor_traits2[, 1, drop = FALSE]
    uncor_traits <- cbind(uncor_traits, uncor_traits2)
}

# Parameter values - reduced for speed
list_k <- c(2)
list_A <- c(0, 1e-3)
list_H <- c(0, 1e-3)
list_hierar_expo <- c(0.25, 0.5, 1, 2, 4)

# data.frame with all combinations of parameter values
comb <- data.frame(
  expand.grid(list_k, list_A, list_H, trait_comb, list_hierar_expo))

colnames(comb) <- c("k", "A", "H", "trait_comb", "hierar_exp")

comb <- distinct(comb) # remove duplicates


# Main Simulations -------------------------------------------------------------
simul <- vector("list", nrow(comb))    # list of simulations
mis_i <- vector("list", nrow(comb))    # list of data.frame with species
                                       # performance & patch of best performance
tra_env <- vector("list", nrow(comb))  # Contain CWM by environment

for (i in seq(nrow(comb))) {

    # Extract simulation trait data.frame
    traits <- data.frame(trait1 = uncor_traits[,comb[i, "trait_comb"]])
    simul_i <- multigen(
        traits = traits, trait_weights = trait_scenar, env = 1:n_patches,
        time = n_gen, species = n_species, patches = 25,
        composition = composition,
        A = comb[i, "A"], B = 1e-7, d = d, k = comb[i, "k"],
        H = comb[i, "H"],
        width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
        hierar_exponent = comb[i, "hierar_exp"])

    simul[[i]] <- simul_i

    ## Data frame of site vs. abundance at last time step
    tra_env_j <- reshape2::melt(simul_i$compo[, , n_gen])
    colnames(tra_env_j) <- c("site", "sp", "ab")

    ## Species trait data.frame
    sp_tra <- data.frame(sp = rownames(uncor_traits),
                         tra = uncor_traits[, comb[i, "trait_comb"]])
    tra_env_j <- suppressWarnings(left_join(tra_env_j, sp_tra, by = "sp"))

    ## Compute CWM
    tra_env_j <- tra_env_j %>%
        group_by(site) %>%
        summarise(cwm = weighted.mean(tra, w = ab, na.rm = TRUE),
                  .groups = "keep") %>%
        mutate(time = n_gen,
               k = comb[i,"k"], A = comb[i, "A"], H = comb[i, "H"],
               trait_comb = comb[i, "trait_comb"],
               param_comb = i, hierar_exp = comb[i, "hierar_exp"],
               env = as.integer(gsub("patches", "", x = site))) %>%
        as.data.frame()

    ## Compute growth rates
    max_time <- n_gen
    gw <- extract_growth_rates(simul = simul_i, chosen_time = n_gen)

    # Compute weighted growth rates
    gwm <- gw %>%
        inner_join(sp_tra, by = c(species = "sp")) %>%
        select(-final_abundance) %>%
        relocate(patch, species, tra) %>%
        group_by(patch) %>%
        # Scale growth rates relatively in each site
        mutate(across(avg_growth_rate:int_growth_rate,
                  ~scales::rescale(.x, c(0, 1)))) %>%
        # Weight trait values per site by growth rates
        summarise(across(avg_growth_rate:int_growth_rate,
                     ~wtd_mean(tra, w = .x, na.rm = TRUE)),
                  .groups = "keep") %>%
        as.data.frame()

    ## Bind patch level results
    tra_env[[i]] <- tra_env_j %>%
        inner_join(gwm, by = c(env = "patch"))

    ## Compute Species Performance and Patch of Best performance for all species
    # return information with growth from env. filtering only, growth from env.
    # filtering and hierarchical competition, real abundance, and abundance
    # based on environmental filtering only
    all_growth <- r_env_CT(simul_i, sp1 = n_species, n_patches = n_patches,
                           time = n_gen)

    ## Compute Species-level mismatches
    matches <- matrix(NA, ncol = 17, nrow = n_species)
    for (z in 1:n_species) {

        # Extract max performance and environmental of max performance for each
        # species
        allout <- extract_mismatchesCT(x = all_growth, z = z)

        sub <- na.omit(gw[which(gw$sp == paste("species", z, sep="")), ])
        if (nrow(sub) > 0) {
            patch_finalabund <- sub[which.max(sub$final_abundance), "patch"]
            patch_avg_growth_rate <- sub[which.max(sub$avg_growth_rate), "patch"]
            patch_max_growth_rate <- sub[which.max(sub$max_growth_rate), "patch"]
            patch_int_growth_rate <- sub[which.max(sub$int_growth_rate), "patch"]
            finalabund <- max(sub$final_abundance, na.rm=TRUE)
            avg_growth_rate <- max(sub$avg_growth_rate, na.rm=TRUE)
            max_growth_rate <- max(sub$max_growth_rate, na.rm=TRUE)
            int_growth_rate <- max(sub$int_growth_rate, na.rm=TRUE)

            matches[z,] <- c(z, allout,
                             patch_finalabund,      finalabund,
                             patch_avg_growth_rate, avg_growth_rate,
                             patch_max_growth_rate, max_growth_rate,
                             patch_int_growth_rate, int_growth_rate)
        } else {
            matches[z,] <- c(z, allout, NA, NA, NA, NA, NA, NA, NA, NA)
        }
    }

    colnames(matches) <- c("Species", "TrueRPatch", "TrueR", "ObsRPatch",
                           "ObsR", "TrueAbPatch", "TrueAb", "ObsAbPatch",
                           "ObsAb", "finalabundPatch", "finalabund",
                           "avgGRPatch", "avgGR", "maxGRPatch", "maxGR",
                           "intGRPatch", "intGR")

    matches <- cbind(matches, comb[rep(i, nrow(matches)), ], traits)

    mis_i[[i]] <- matches

    cat(100*i/nrow(comb), "\t")
}

# Compute Species Mismatches ---------------------------------------------------

# Calculate all relative species-level mismatches
# (either by performance or by patch)
allmis <- bind_rows(mis_i)


allmis <- allmis %>%
    mutate(
        ## Real final abudance vs. abundance with only environmental filtering
        MisAbP = finalabund - ObsAb,
        # Same but using patch of greatest abundance
        MisAbPatchP = finalabundPatch - ObsAbPatch,
        ## Real intrinsic growth rate (growth rate at first timesteps) vs.
        ## growth rate observed with environmental filtering only
        MisinstRP = intGR - ObsR,
        # Same but with patch of best growth
        MisinstRPatchP = intGRPatch - ObsRPatch,
        ## Average growth rate over all timesteps vs. growth rate observed with
        ## environmental filtering only
        MisavgGRP = avgGR - ObsR,
        # Same but with patch of best growth
        MisavgGRPatchP = avgGRPatch - ObsRPatch,
        ## Maximum growth rate over all timesteps vs. growth rate observed with
        ## environmental filtering only
        MismaxGRP = maxGR - ObsR,
        # Same but with patch of best growth
        MismaxGRPatchP = maxGRPatch - ObsRPatch,
        ## Concatenate parameter values
        comb = paste(k, A, H, hierar_exp, sep = "_")
    ) %>%
    # Make all mismatches relative to gradient length
    mutate_at(vars(starts_with("Mis")), ~./n_patches)

# Get average of mismatches per species
allmisA <- aggregate(allmis, by = list(allmis$Species, allmis$comb), "mean",
                     na.rm = TRUE)
allmisA$comb <- allmisA$Group.2
# Convert back comb column in parameter combinations without hierarchical exp.
allmisA <- allmisA %>%
    extract(comb, c("comb", "hierar_exp"), c("(.*)_(.+)$")) %>%
    mutate(hierar_exp = as.numeric(hierar_exp))

unique(allmisA$comb)
# subset into distinct scenarios
nocomp <- allmisA[which(allmisA$comb %in% c("2_0_0")), ]
Acomp <- allmisA[which(allmisA$comb %in% c("2_0.001_0")), ]
Hcomp <- allmisA[which(allmisA$comb %in% c("2_0_0.001")), ]
allcomp <- allmisA[which(allmisA$comb %in% c("2_0.001_0.001")), ]


# Fig. 2: species mismatches in function of competition type -------------------
all_minmax_mismatch = allmisA %>%
    select(Species, comb, hierar_exp, MisAbPatchP, MisinstRPatchP,
           MismaxGRPatchP, MisavgGRPatchP) %>%
    rowwise() %>%
    mutate(min_mismatch = min(MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
                              MisavgGRPatchP),
           max_mismatch = max(MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
                              MisavgGRPatchP)) %>%
    select(Species, min_mismatch, max_mismatch, comb, hierar_exp)

plot_species_mismatch = allmisA %>%
    filter(comb != "2_0.001_0.001", hierar_exp == 2) %>%
    select(Species, MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
           MisavgGRPatchP, comb, hierar_exp) %>%
    tidyr::gather("mismatch_name", "mismatch_value",
                  MisAbPatchP:MisavgGRPatchP) %>%
    inner_join(all_minmax_mismatch, by = c("Species", "comb", "hierar_exp")) %>%
    mutate(comb = factor(comb, levels = c("2_0_0", "2_0.001_0",
                                          "2_0_0.001"))) %>%
    ggplot(aes(mismatch_value, Species, color = mismatch_name,
               shape = mismatch_name)) +
    geom_segment(aes(x = min_mismatch, xend = max_mismatch,
                     yend = Species),
                 color = "grey", size = 2/3) +
    geom_vline(xintercept = 0, linetype = 1, size = 1/2, color = "black") +
    geom_point(size = 1.5) +
    facet_wrap(vars(comb), ncol = 3,
               labeller = labeller(
                   comb = c("2_0.001_0" = "+Limiting Similarity",
                            "2_0_0.001" = "+Hierarchical Competition",
                            "2_0_0" = "Environmental filtering\nonly"))) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_y_continuous(limits = c(1, 50), breaks = c(1, c(1,2,3,4,5)*10)) +
    scale_color_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate"
    )) +
    scale_shape_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate"
    )) +
    labs(x = "Relative Mismatch from True Optimal Patch (% of gradient)",
         shape = "Performance Measure",
         color = "Performance Measure") +
    theme_bw(12) +
    theme(aspect.ratio = 1,
          panel.grid.major.x = element_line(size = 1.3),
          panel.spacing.x = unit(4, "mm"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")

plot_species_mismatch


# Fig. 3: Deviation along the environment --------------------------------------

all_trait_env = bind_rows(tra_env) %>%
  tidyr::gather("cwm_type", "cwm_value", cwm,
                avg_growth_rate:int_growth_rate) %>%
  mutate(mismatch_value = cwm_value - env,
         comb = paste(k, A, H, sep = "_"))

plot_deviation_env = all_trait_env %>%
  filter(comb != "2_0.001_0.001", hierar_exp == 2) %>%
  mutate(comb = factor(comb, levels = c("2_0_0", "2_0.001_0",
                                        "2_0_0.001"))) %>%
  ggplot(aes(env, mismatch_value, color = cwm_type, shape = cwm_type)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_wrap(vars(comb), ncol = 3,
             labeller = labeller(
               comb = c("2_0.001_0" = "+Limiting Similarity",
                        "2_0_0.001" = "+Hierarchical Competition",
                        "2_0_0" = "Environmental filtering\nonly"))) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point", size = 1) +
  scale_color_discrete(labels = c(
    cwm = "CWM",
    avg_growth_rate = "Avg. Growth Rate Weighted-Mean",
    int_growth_rate = "Intrinsic Growth Rate Weighted-Mean",
    max_growth_rate = "Maximum Growth Rate Weighted-Mean"
  )) +
  scale_shape_discrete(labels = c(
    cwm = "CWM",
    avg_growth_rate = "Avg. Growth Rate Weighted-Mean",
    int_growth_rate = "Intrinsic Growth Rate Weighted-Mean",
    max_growth_rate = "Maximum Growth Rate Weighted-Mean"
  )) +
  labs(x = "Environment",
       y = "Deviation between Observed CWM vs. expected CWM",
       color = "CWM Type",
       shape = "CWM Type") +
  theme_bw(12) +
  theme(aspect.ratio = 1,
        panel.grid.major.x = element_line(size = 1.3),
        panel.spacing.x = unit(4, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_blank(),
        legend.position = "top")

plot_deviation_env

# Fig. 4: CWM & CWV in function of trait contribution to growth ----------------

# Generate a second trait value randomly
set.seed(20201015)
second_trait = generate_cor_traits_rand(n_patches, n_species,
                                        n_traits - 1, cor_coef = 0)[,2]

# Generate different scenarios where trait contributes differently to growth
var_growth_scenars = lapply(seq(0, 1, length.out = 6), function(x) {
  data.frame(trait            = c("trait1", "trait2"),
             growth_weight    = c(x, 1-x),
             compet_weight    = 0,
             hierarchy_weight = 0)
})

# Get a combinations of simulations parameters
var_growth_comb = expand.grid(trait_comb    = trait_comb,
                              growth_scenar = seq_along(var_growth_scenars))

# Initialize list of simulations
growth_simuls <- vector("list", nrow(var_growth_comb))

# Simulations of varying trait growth contribution
for (i in seq_len(nrow(var_growth_comb))) {

  # Get trait dataset
  traits <- data.frame(trait1 = uncor_traits[,var_growth_comb[i, "trait_comb"]],
                       trait2 = second_trait)

  # Get trait contribution scenario
  growth_scenar <- var_growth_scenars[[var_growth_comb[i, "growth_scenar"]]]

  # Actual Simulation
  growth_simul <- multigen(
    traits = traits, trait_weights = growth_scenar, env = 1:n_patches,
    time = n_gen, species = n_species, patches = 25,
    composition = composition,
    A = 0, B = 1e-7, d = d, k = 2,
    H = 0 ,
    width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
    hierar_exponent = 2)

  growth_simuls[[i]] <- growth_simul
}

var_growth_cwm <- purrr::map2_dfr(
  growth_simuls, seq_len(nrow(var_growth_comb)), function(x, y) {

    # Extract actual growth contribution
    growth_contribution <- var_growth_scenars[[var_growth_comb[y, "growth_scenar"]]]$growth_weight[[1]]

    # Extract final abundances
    abund_df <- x$compo[,,100] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("patch") %>%
      tidyr::gather("species", "abundance", -patch)

    # Get back trait data
    trait_df <- data.frame(trait1 = uncor_traits[,var_growth_comb[y, "trait_comb"]]) %>%
      tibble::rownames_to_column("species")

    # Compute Community-Weighted Mean and Community-Weighted Variance
    cwm_df <- abund_df %>%
      dplyr::inner_join(trait_df, by = "species") %>%
      group_by(patch) %>%
      summarise(cwm = wtd_mean(trait1, abundance),
                cwv = wtd_var(trait1, abundance),
                .groups = "drop") %>%
      mutate(growth_weight = growth_contribution,
             env = as.numeric(gsub("patches", "", patch)),
             param_comb = y)
})

plot_cwm_cwv_growth = var_growth_cwm %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(growth_weight))) +
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CWM", cwv = "CWV"))) +
  labs(x = "Environment",
       y = "CWM or CWV value") +
  scale_colour_manual(name = "Trait Contribution\nto growth",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_growth

# Supp. Fig. 1: Parameter space ------------------------------------------------

# Generate parameter space
param_space = expand.grid(
  k          = c(2),
  A          = c(0, 1e-7, 1e-5, 1e-3, 1e-1),
  H          = c(0, 1e-7, 1e-5, 1e-3, 1e-1),
  trait_comb = trait_comb[1:30]
)

param_space_simuls <- vector("list", nrow(var_growth_comb))

# Simulations of varying trait growth contribution
for (i in seq_len(nrow(param_space))) {

  # Get trait dataset
  traits <- data.frame(trait1 = uncor_traits[,param_space[i, "trait_comb"]])

  # Actual Simulation
  param_space_simul <- multigen(
    traits = traits, trait_weights = trait_scenar, env = 1:n_patches,
    time = n_gen, species = n_species, patches = 25,
    composition = composition,
    B = 1e-7, d = d,
    A = param_space[i, "A"],
    k = param_space[i, "k"],
    H = param_space[i, "H"],
    width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
    hierar_exponent = 2)

  param_space_simuls[[i]] <- param_space_simul
}

param_space_abund <- purrr::map2_dfr(
  param_space_simuls, seq_len(nrow(param_space)), function(x, y) {
    # Extract final abundances
    abund_df <- x$compo[,,100] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("patch") %>%
      tidyr::gather("species", "abundance", -patch) %>%
      mutate(k = param_space[y, "k"],
             A = param_space[y, "A"],
             H = param_space[y, "H"],
             trait_comb = param_space[y, "trait_comb"])

  })

param_space_avg = param_space_abund %>%
  group_by(k, A, H, trait_comb, patch) %>%
  summarise(n_sp = sum(abundance > 0),
            avg_ab = mean(abundance))

plot_param_space_rich = param_space_avg %>%
  ungroup() %>%
  mutate(A_chr = format(A, scientific = FALSE),
         H_chr = format(H, scientific = FALSE)) %>%
  ggplot(aes(A_chr, H_chr, z = n_sp)) +
  stat_summary_2d() +
  scale_fill_viridis_c(name = "Species\nRichness") +
  scale_x_discrete(labels = function(x) scales::scientific(as.numeric(x))) +
  scale_y_discrete(labels = function(x) scales::scientific(as.numeric(x))) +
  labs(x = "Limiting similarity intensity",
       y = "Hierarchical competition intensity") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        legend.position = "top")

plot_param_space_abund = param_space_avg %>%
  ungroup() %>%
  mutate(A_chr = format(A, scientific = FALSE),
         H_chr = format(H, scientific = FALSE)) %>%
  ggplot(aes(A_chr, H_chr, z = avg_ab)) +
  stat_summary_2d() +
  scale_fill_viridis_c(name = "Average\nAbundance") +
  scale_x_discrete(labels = function(x) scales::scientific(as.numeric(x))) +
  scale_y_discrete(labels = function(x) scales::scientific(as.numeric(x))) +
  labs(x = "Limiting similarity intensity",
       y = "Hierarchical competition intensity") +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        legend.position = "top")

plot_param_space = cowplot::plot_grid(plot_param_space_rich,
                                      plot_param_space_abund, labels = "AUTO")



# Supp. Fig. 2: effect of hierarchical exponent --------------------------------
plot_hierar_exp_species_mismatch = allmisA %>%
    filter(comb == "2_0_0.001") %>%
    select(Species, MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
           MisavgGRPatchP, comb, hierar_exp) %>%
    tidyr::gather("mismatch_name", "mismatch_value",
                  MisAbPatchP:MisavgGRPatchP) %>%
    inner_join(all_minmax_mismatch, by = c("Species", "comb", "hierar_exp")) %>%
    mutate(comb = factor(comb, levels = c("2_0_0", "2_0.001_0",
                                          "2_0_0.001"))) %>%
    ggplot(aes(mismatch_value, Species, color = mismatch_name,
               shape = mismatch_name)) +
    geom_segment(aes(x = min_mismatch, xend = max_mismatch,
                     yend = Species),
                 color = "grey", size = 2/3) +
    geom_vline(xintercept = 0, linetype = 1, size = 1/2, color = "black") +
    geom_point(size = 1.5) +
    facet_wrap(vars(hierar_exp), ncol = 3,
               labeller = as_labeller(label_hierar_exp, default = label_parsed)) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_y_continuous(limits = c(1, 50), breaks = c(1, c(1,2,3,4,5)*10)) +
    scale_color_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate"
    )) +
    scale_shape_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate"
    )) +
    labs(x = "Relative Mismatch from True Optimal Patch (% of gradient)",
         shape = "Performance Measure",
         color = "Performance Measure") +
    theme_bw(12) +
    theme(aspect.ratio = 1,
          panel.grid.major.x = element_line(size = 1.3),
          panel.spacing.x = unit(4, "mm"),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 8),
          legend.position = "top")

plot_hierar_exp_species_mismatch


# Supp. Fig. 3: contribution to competition effect on CWM and CWV --------------

var_comp_scenars = lapply(seq(0, 1, length.out = 6), function(x) {
  data.frame(trait            = c("trait1", "trait2"),
             growth_weight    = c(0, 1),
             compet_weight    = c(x, 1-x),
             hierarchy_weight = 0)
})

var_hierar_scenars = lapply(seq(0, 1, length.out = 6), function(x) {
  data.frame(trait            = c("trait1", "trait2"),
             growth_weight    = c(0, 1),
             compet_weight    = 0,
             hierarchy_weight = c(x, 1-x))
})

var_comp_scenars = c(var_comp_scenars, var_hierar_scenars)

# Get a combinations of simulations parameters
var_comp_comb = expand.grid(trait_comb    = trait_comb,
                            comp_scenar = seq_along(var_comp_scenars))

# Initialize list of simulations
comp_simuls <- vector("list", nrow(var_comp_comb))

# Simulations of varying trait growth contribution
for (i in seq_len(nrow(var_comp_comb))) {

  # Get trait dataset
  traits <- data.frame(trait1 = uncor_traits[,var_comp_comb[i, "trait_comb"]],
                       trait2 = second_trait)

  # Get trait contribution scenario
  comp_scenar <- var_comp_scenars[[var_comp_comb[i, "comp_scenar"]]]

  A = ifelse(comp_scenar$compet_weight[[1]] == 0, 0, 1e-3)
  H = ifelse(comp_scenar$hierarchy_weight[[1]] == 0, 0, 1e-3)

  # Actual Simulation
  comp_simul <- multigen(
    traits = traits, trait_weights = comp_scenar, env = 1:n_patches,
    time = n_gen, species = n_species, patches = 25,
    composition = composition,
    A = A, B = 1e-7, d = d, k = 2,
    H = H ,
    width = rep(width, n_patches), h_fun = "+", di_thresh = 24, K = 100,
    hierar_exponent = 2)

  comp_simuls[[i]] <- comp_simul
}

var_comp_cwm <- purrr::map2_dfr(
  comp_simuls, seq_len(nrow(var_comp_comb)), function(x, y) {

    # Extract trait contribution to competitions
    comp_contribution <- var_comp_scenars[[var_comp_comb[y, "comp_scenar"]]]$compet_weight[[1]]
    hier_contribution <- var_comp_scenars[[var_comp_comb[y, "comp_scenar"]]]$hierarchy_weight[[1]]

    # Extract final abundances
    abund_df <- x$compo[,,100] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("patch") %>%
      tidyr::gather("species", "abundance", -patch)

    # Get back trait data
    trait_df <- data.frame(trait1 = uncor_traits[,var_comp_comb[y, "trait_comb"]]) %>%
      tibble::rownames_to_column("species")

    # Compute Community-Weighted Mean and Community-Weighted Variance
    cwm_df <- abund_df %>%
      dplyr::inner_join(trait_df, by = "species") %>%
      group_by(patch) %>%
      summarise(cwm = wtd_mean(trait1, abundance),
                cwv = wtd_var(trait1, abundance),
                .groups = "drop") %>%
      mutate(comp_weight = comp_contribution,
             hier_weight = hier_contribution,
             env = as.numeric(gsub("patches", "", patch)),
             param_comb = y)
  })

plot_cwm_cwv_comp = var_comp_cwm %>%
  filter(hier_weight == 0) %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(comp_weight))) +
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CWM", cwv = "CWV"))) +
  labs(x = "Environment",
       y = "CWM or CWV value") +
  scale_colour_manual(name = "Trait Contribution\nto limiting similarity",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_hier = var_comp_cwm %>%
  filter(comp_weight == 0) %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(hier_weight))) +
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CWM", cwv = "CWV"))) +
  labs(x = "Environment",
       y = "CWM or CWV value") +
  scale_colour_manual(name = "Trait Contribution\nto hierarchical comp.",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_comp_contrib = cowplot::plot_grid(
  plot_cwm_cwv_comp,
  plot_cwm_cwv_hier,
  labels = "AUTO",
  ncol = 1, nrow = 2)
