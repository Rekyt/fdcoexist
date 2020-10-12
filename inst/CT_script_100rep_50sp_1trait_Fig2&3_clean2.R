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


# Computation & Outputs --------------------------------------------------------
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
    for (z in 1:n_species){

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

# Compute Mismatches for each Species ------------------------------------------

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


# Plot Figure 2: species mismatches in function of competition type ------------
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
label_hierar_exp = function(string) {

    lapply(unlist(string), function(x) {
        paste0("hierarchical~distance^{", x, "}")
    })
}

# Plot Figure 3: Deviation along the environment -------------------------------

# Plot Figure 4: CWM & CWV in function of trait contribution to growth ---------

# Supplementary Figure 1: effect of hierarchical exponent ----------------------
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
          strip.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")

plot_hierar_exp_species_mismatch


# Supplementary Figure 2: Parameter space --------------------------------------


# Correlation of Mismatches ----------------------------------------------------
###tables of correlations
#allcomp
cor(allcomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
               "MisavgGRPatchP")],
    method="spearman")

# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000     0.4164496     0.6729627     0.8061727
# MisPatchinstR  0.4164496     1.0000000     0.3293066     0.5395412
# MismaxGRPatch  0.6729627     0.3293066     1.0000000     0.8102633
# MisavgGRPatch  0.8061727     0.5395412     0.8102633     1.0000000

#H comp
cor(Hcomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
             "MisavgGRPatchP")],
    method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000    0.40114217    0.10229505    0.49460843
# MisPatchinstR  0.4011422    1.00000000    0.04985318    0.82697422
# MismaxGRPatch  0.1022950    0.04985318    1.00000000   -0.09408434
# MisavgGRPatch  0.4946084    0.82697422   -0.09408434    1.00000000

cor(Acomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
             "MisavgGRPatchP")], method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb     1.0000000     0.7729351     0.6976272     0.9281220
# MisPatchinstR  0.7729351     1.0000000     0.4811728     0.7313187
# MismaxGRPatch  0.6976272     0.4811728     1.0000000     0.7848516
# MisavgGRPatch  0.9281220     0.7313187     0.7848516     1.0000000

#no comp
cor(nocomp[,c("MisAbPatchP", "MisinstRPatchP", "MismaxGRPatchP",
              "MisavgGRPatchP")], method="spearman")
# MisPatchAb MisPatchinstR MismaxGRPatch MisavgGRPatch
# MisPatchAb             1            NA            NA            NA
# MisPatchinstR         NA             1            NA            NA
# MismaxGRPatch         NA            NA             1            NA
# MisavgGRPatch         NA            NA            NA             1


# Mismatches by patch (CWM/GWM) ------------------------------------------------

cwmdat <- bind_rows(tra_env)
cwmdat$comb <- paste(cwmdat$k, cwmdat$A, cwmdat$H, sep = "_")
cwmdat$patch <- as.numeric(gsub("[[:alpha:]]","", cwmdat$site))

site_avg_perf = cwmdat %>%
    select(comb, site, env, cwm, avg_growth_rate, max_growth_rate, int_growth_rate) %>%
    # Average CWM/GWM across all trait replicates
    group_by(comb, site, env) %>%
    summarise(across(cwm:int_growth_rate, mean))

plot_site_mismatch = site_avg_perf %>%
    # Remove limiting similarity only scenarios
    filter(comb != "2_0.001_0.001") %>%
    ungroup() %>%
    # Tidy data to get properly ordered for ggplot2
    tidyr::gather("perf_name", "perf_value",
                  cwm:int_growth_rate) %>%
    # Reorder factor levels
    mutate(comb = factor(comb, levels = c("2_0_0", "2_0.001_0",
                                          "2_0_0.001")),
           perf_name = factor(
               perf_name,
               levels = c("cwm", "avg_growth_rate", "int_growth_rate",
                          "max_growth_rate"))
    ) %>%
    ggplot(aes(env, (perf_value - env)/25, color = perf_name,
               shape = perf_name)) +
    geom_hline(yintercept = 0, linetype = 1, size = 1/2, color = "black") +
    geom_line() +
    geom_point(size = 1.5) +
    facet_wrap(vars(comb), ncol = 3,
               labeller = labeller(
                   comb = c("2_0.001_0" = "+Limiting Similarity",
                            "2_0_0.001" = "+Hierarchical Competition",
                            "2_0_0" = "Environmental filtering\nonly"))) +
    scale_x_continuous(breaks = c(1, (1:5)*5)) +
    scale_y_continuous(labels = scales::label_percent()) +
    scale_color_discrete(labels = c(
        cwm = "Abundance",
        avg_growth_rate = "Average Growth Rate",
        int_growth_rate = "Intrinsic Growth Rate",
        max_growth_rate = "Maximum Growth Rate"
    )) +
    scale_shape_discrete(labels = c(
        cwm = "Abundance",
        avg_growth_rate = "Average Growth Rate",
        int_growth_rate = "Intrinsic Growth Rate",
        max_growth_rate = "Maximum Growth Rate"
    )) +
    labs(x = "Site",
         y = "Relative Site-level Mismatch (% of gradient)",
         shape = "Performance Measure",
         color = "Performance Measure") +
    theme_bw(12) +
    theme(aspect.ratio = 1,
          panel.grid.major.y = element_line(size = 1.3),
          panel.spacing.x = unit(4, "mm"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")

plot_site_mismatch

##These need to have error bars added!
save.image(file="/Users/tuckerca/Dropbox/Documents/Caroline/fdcoexist/2020/Data/100rep_50sp_1trait_scenarios.RData")


# Compute distinctiveness ------------------------------------------------------

all_di = uncor_traits %>%
    apply(2, function(x) {
      traits = x

      sc_traits = (traits - min(traits))/(max(traits) - min(traits))

      names(sc_traits) = 1:50

      trait_dist = dist(sc_traits)

      funrar::distinctiveness_global(trait_dist, "reg_di")
    }) %>%
    purrr::set_names(nm = seq(ncol(uncor_traits))) %>%
  bind_rows(.id = "trait_comb")


mis_di = bind_rows(mis_i_pred) %>%
    filter(hierar_exp == 2.00) %>%  # Filter original hierarchical exponent value
    inner_join(all_di %>%
                   mutate(trait_comb = as.integer(trait_comb),
                          Species = as.integer(species)), by = c("trait_comb", "Species")) %>%
    group_by(finalabundPatch) %>%
    summarise(k = mean(k, A = mean(A), H = mean(H)))

patch_di = simul_i$compo[,,"time100"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("patch") %>%
  tidyr::gather("species", "abund", species1:species50) %>%
  inner_join(all_di %>%
               filter(trait_comb == 1) %>%
               mutate(species = paste0("species", species)),
             by = "species") %>%
  filter(abund > 0) %>%
  group_by(patch) %>%
  summarise(avg_di = mean(reg_di),
            cwm_di = weighted.mean(reg_di, abund))

avg_di_df = purrr::map2(simul, tra_env, function(x, y) {

  x$compo[,,"time100"] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("patch") %>%
    tidyr::gather("species", "abund", species1:species50)  %>%
    inner_join(all_di %>%
                 filter(trait_comb == unique(y$trait_comb)) %>%
                 mutate(species = paste0("species", species)),
               by = "species") %>%
    filter(abund > 0) %>%
    group_by(patch) %>%
    summarise(avg_di = mean(reg_di),
              cwm_di = weighted.mean(reg_di, abund),
              .groups = "drop") %>%
    mutate(k = unique(y$k),
           A = unique(y$A),
           H = unique(y$H))
}) %>%
  bind_rows()

