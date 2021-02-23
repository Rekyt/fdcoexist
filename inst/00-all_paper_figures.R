# FDcoexist June 4 2020
# Authors: C. Tucker & M. Grenié
# Generate species level responses to site conditions. Only moderate value of A,
# H or both used.
# replicates of each required - currently ~100

# Load Packages ----------------------------------------------------------------
library(dplyr)
library(cowplot)
library(tidyr)
library(RColorBrewer)
devtools::load_all()  # Load all functions from package in local repo


# Functions --------------------------------------------------------------------
label_hierar_exp = function(string) {

  lapply(unlist(string), function(x) {
    paste0("hierarchical~distance^{", x, "}")
  })
}

# Initial Parameters -----------------------------------------------------------
set.seed(20201008) # To be fixed to replicate our results
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
list_A <- c(0, 1e-4)
list_H <- c(0, 5e-3)
list_hierar_expo <- c(0.5, 1, 2)



# data.frame with all combinations of parameter values
comb <- expand.grid(k = list_k, A = list_A, H = list_H, trait_comb = trait_comb,
                    hierar_exp = list_hierar_expo) %>%
  data.frame()

comb <- distinct(comb) # remove duplicates

# Short handles for different scenarios
scenarios = comb %>%
  mutate(comb = paste(k, A, H, sep = "_")) %>%
  distinct(comb) %>%
  pull()

names(scenarios) = c("none", "limsim", "hier", "both")

legend_comb =  c("+Limiting Similarity", "+Hierarchical Competition")
names(legend_comb) = c(scenarios[["limsim"]],
                       scenarios[["hier"]])

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
        hierar_exponent = comb[i, "hierar_exp"], lim_sim_exponent = 1)

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
        ## Real final abundance vs. abundance with only environmental filtering
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
    mutate_at(vars(starts_with("Mis")), list(relative = ~./n_patches))

# Get average of mismatches per species
allmisA <- aggregate(allmis, by = list(allmis$Species, allmis$comb), "mean",
                     na.rm = TRUE) %>%
  select(-Species, -trait_comb, -comb) %>%
  rename(Species = Group.1,
         comb = Group.2) %>%
  select(k, A, H, hierar_exp, comb, Species, everything())
# Convert back comb column in parameter combinations without hierarchical exp.
allmisA <- allmisA %>%
    extract(comb, c("comb", "hierar_exp"), c("(.*)_(.+)$")) %>%
    mutate(hierar_exp = as.numeric(hierar_exp))

unique(allmisA$comb)
# subset into distinct scenarios
nocomp <- allmisA[which(allmisA$comb %in% scenarios[["none"]]), ]
Acomp <- allmisA[which(allmisA$comb %in% scenarios[["limsim"]]), ]
Hcomp <- allmisA[which(allmisA$comb %in% scenarios[["hier"]]), ]
allcomp <- allmisA[which(allmisA$comb %in% scenarios[["both"]]), ]


saveRDS(allmisA, "inst/all_mismatches.Rds")

# Fig. 1: 3 species example of the effect of processes -------------------------
# Please refer to script inst/01-figure1_example_3sp.R

# Fig. 2: species mismatches in function of competition type -------------------

# Get for each species in each parameter combination the minimum and maximum
# mistmatch values observed across all Patch mismatches
all_minmax_mismatch = allmisA %>%
    select(Species, comb, hierar_exp, MisAbPatchP, MisinstRPatchP,
           MismaxGRPatchP, MisavgGRPatchP) %>%
    rowwise() %>%
    mutate(min_mismatch = min(MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
                              MisavgGRPatchP),
           max_mismatch = max(MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
                              MisavgGRPatchP)) %>%
    select(Species, min_mismatch, max_mismatch, comb, hierar_exp)

# Combine and keep only relevant information for figure 2
mismatch_extract = allmisA %>%
  select(k, A, H, Species, comb, hierar_exp, MisAbPatchP, MisinstRPatchP,
         MismaxGRPatchP, MisavgGRPatchP) %>%
  tidyr::pivot_longer(
    c(MisAbPatchP, MisinstRPatchP, MismaxGRPatchP, MisavgGRPatchP),
    names_to = "mismatch_name", values_to = "mismatch_value") %>%
  inner_join(all_minmax_mismatch, by = c("Species", "comb", "hierar_exp")) %>%
  filter(comb != scenarios[["none"]], comb != scenarios[["both"]],
         hierar_exp == 0.5)

# Get the correlation between mismatches
mismatch_cor = mismatch_extract %>%
  select(comb, Species, mismatch_name, mismatch_value) %>%
  tidyr::pivot_wider(names_from = mismatch_name,
                     values_from = mismatch_value) %>%
  tidyr::nest(mismatches = c(Species, MisAbPatchP, MisinstRPatchP,
                             MisavgGRPatchP, MismaxGRPatchP)) %>%
  mutate(cor_mat = purrr::map(
    mismatches, ~cor(.x[, -1], method = "spearman", use = "complete.obs") %>%
                                round(2)))

## Convert correlation table to image
# Limiting similarity scenario
limsim_table = mismatch_cor$cor_mat[[1]][, 1:3]
limsim_table[upper.tri(limsim_table)] = ""
limsim_table = as.data.frame(limsim_table, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("mismatch_name") %>%
  mutate(mismatch_name = case_when(
    mismatch_name == "MisinstRPatchP" ~ "&#9632;",
    mismatch_name == "MisavgGRPatchP" ~ "&#9650;",
    mismatch_name == "MismaxGRPatchP" ~ "+",
    TRUE ~ "",
  )) %>%
  mutate(
    mismatch_name = kableExtra::cell_spec(
      mismatch_name, color = c("white", "#00BFC4", "#7CAE00", "#C77CFF"),
      bold = TRUE, escape = FALSE, align = "c")
  )
limsim_table[1, 2:4] = c(
  '<span style=" font-weight: bold;    color: #F8766D !important;" >&#9679;</span>',
  '<span style=" font-weight: bold;    color: #00BFC4 !important;" >&#9632;</span>',
  '<span style=" font-weight: bold;    color: #7CAE00 !important;" >&#9650;</span>')
limsim_table[2, 3] = ""
limsim_table[3, 4] = ""

limsim_table %>%
  kableExtra::kable(escape = FALSE, col.names = NULL) %>%
  kableExtra::kable_minimal(full_width = FALSE, font_size = 14) %>%
  kableExtra::as_image(file = "inst/figures/paper_figure2_limsimtable.png")
limsim_png = png::readPNG("inst/figures/paper_figure2_limsimtable.png",
                          native = TRUE)

# Hierarchical Competition Scenario
hiercomp_table = mismatch_cor$cor_mat[[2]][, 1:3]
hiercomp_table[upper.tri(hiercomp_table)] = ""
hiercomp_table = as.data.frame(hiercomp_table, stringsAsFactors = FALSE) %>%
  tibble::rownames_to_column("mismatch_name") %>%
  mutate(mismatch_name = case_when(
    mismatch_name == "MisinstRPatchP" ~ "&#9632;",
    mismatch_name == "MisavgGRPatchP" ~ "&#9650;",
    mismatch_name == "MismaxGRPatchP" ~ "+",
    TRUE ~ "",
  )) %>%
  mutate(
    mismatch_name = kableExtra::cell_spec(
      mismatch_name, color = c("white", "#00BFC4", "#7CAE00", "#C77CFF"),
      bold = TRUE,
      escape = FALSE)
  )
hiercomp_table[1, 2:4] = c(
  '<span style=" font-weight: bold;    color: #F8766D !important;" >&#9679;</span>',
  '<span style=" font-weight: bold;    color: #00BFC4 !important;" >&#9632;</span>',
  '<span style=" font-weight: bold;    color: #7CAE00 !important;" >&#9650;</span>')
hiercomp_table[2, 3] = ""
hiercomp_table[3, 4] = ""

hiercomp_table %>%
  kableExtra::kable(escape = FALSE, col.names = NULL) %>%
  kableExtra::row_spec(1, background = "white") %>%
  kableExtra::row_spec(2, background = "white") %>%
  kableExtra::row_spec(3, background = "white") %>%
  kableExtra::row_spec(4, background = "white") %>%
  kableExtra::kable_minimal(full_width = FALSE, font_size = 14) %>%
  kableExtra::as_image(file = "inst/figures/paper_figure2_hiercomptable.png")
hiercomp_png = png::readPNG("inst/figures/paper_figure2_hiercomptable.png",
                            native = TRUE)

# Actual Main plot
plot_species_mismatch = mismatch_extract %>%
  mutate(comb = factor(comb, levels = c(scenarios[["limsim"]],
                                        scenarios[["hier"]]))) %>%
  ggplot(aes(mismatch_value, Species, color = mismatch_name,
             shape = mismatch_name)) +
  geom_segment(aes(x = min_mismatch, xend = max_mismatch,
                   yend = Species),
               color = "grey", size = 1/3) +
  geom_vline(xintercept = 0, linetype = 1, size = 1/2, color = "black") +
  geom_point(size = 1.5) +
  facet_wrap(vars(comb), ncol = 2, labeller = labeller(comb = legend_comb)) +
  scale_y_continuous(limits = c(1, 50), breaks = c(1, c(1,2,3,4,5)*10)) +
  scale_color_discrete(labels = c(
    MisAbPatchP = "Abundance",
    MisavgGRPatchP = "Average Growth Rate",
    MisinstRPatchP = "Intrinsic Growth Rate",
    MismaxGRPatchP = "Maximum Growth Rate"
  ), guide = guide_legend(nrow = 2)) +
  scale_shape_discrete(labels = c(
    MisAbPatchP = "Abundance",
    MisavgGRPatchP = "Average Growth Rate",
    MisinstRPatchP = "Intrinsic Growth Rate",
    MismaxGRPatchP = "Maximum Growth Rate"
  ), guide = guide_legend(nrow = 2)) +
  labs(x     = "Patch Mismatch",
       shape = "Performance\nMeasure",
       color = "Performance\nMeasure") +
  theme_bw(10) +
  theme(aspect.ratio = 1,
        panel.grid.major.x = element_line(size = 1.3),
        panel.spacing.x = unit(4, "mm"),
        strip.background = element_blank(),
        panel.border = element_blank(),
        legend.position = "top")

paper_figure2 = plot_species_mismatch +
  patchwork::inset_element(limsim_png,   0.15, 0.5, 0.25,  0.7,
                           align_to = "panel") +
  patchwork::inset_element(hiercomp_png, 0.75, 0.5, 0.85, 0.7,
                           align_to = "panel") +
  theme_void()

paper_figure2

saveRDS(paper_figure2, "inst/figures/paper_figure2.Rds")

ggsave2("inst/figures/paper_figure2.pdf", paper_figure2,
        width = 16.6, height = 8.5,
        units = "cm", dpi = 300, scale = 1.5)

ggsave2("inst/figures/paper_figure2.svg", paper_figure2,
        width = 16.6, height = 8.5,
        units = "cm", dpi = 300, scale = 1.5)

ggsave2("inst/figures/paper_figure2.png", paper_figure2,
        width = 16.6, height = 8.5,
        units = "cm", dpi = 300, scale = 1.5)

# Fig. 3: Deviation along the environment --------------------------------------

all_trait_env = bind_rows(tra_env) %>%
  tidyr::gather("cwm_type", "cwm_value", cwm,
                avg_growth_rate:int_growth_rate) %>%
  mutate(mismatch_value = cwm_value - env,
         comb = paste(k, A, H, sep = "_"))

saveRDS(all_trait_env, "inst/all_trait_env.Rds")

legend_comb[3] = "Environmental filtering\nonly"
names(legend_comb)[3] = scenarios[["none"]]

plot_deviation_env = all_trait_env %>%
  filter(comb != scenarios[["both"]], hierar_exp == 2) %>%
  mutate(comb = factor(comb,
                       levels = c(scenarios[["none"]], scenarios[["limsim"]],
                                  scenarios[["hier"]]))) %>%
  ggplot(aes(env, mismatch_value, color = cwm_type, shape = cwm_type)) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_wrap(vars(comb), ncol = 3, labeller = labeller(comb = legend_comb)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point", size = 1) +
  scale_color_discrete(labels = c(
    cwm = "CWM",
    avg_growth_rate = "Avg. Growth Rate\nWeighted-Mean",
    int_growth_rate = "Intrinsic Growth Rate\nWeighted-Mean",
    max_growth_rate = "Maximum Growth Rate\nWeighted-Mean"
  )) +
  scale_shape_discrete(labels = c(
    cwm = "CWM",
    avg_growth_rate = "Avg. Growth Rate\nWeighted-Mean",
    int_growth_rate = "Intrinsic Growth Rate\nWeighted-Mean",
    max_growth_rate = "Maximum Growth Rate\nWeighted-Mean"
  )) +
  labs(x = "Environment",
       y = "Deviation (Observed CWM vs. Expected CWM)",
       color = "CWM\nType",
       shape = "CWM\nType") +
  theme_bw(10) +
  theme(aspect.ratio = 1,
        panel.grid.major.x = element_line(size = 1.3),
        panel.spacing.x = unit(4, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.border = element_blank(),
        legend.position = "top")

plot_deviation_env

saveRDS(plot_deviation_env, "inst/figures/paper_figure3.Rds")

ggsave2("inst/figures/paper_figure3.pdf", plot_deviation_env,
        width = 16.6, height = 8.8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_figure3.png", plot_deviation_env,
        width = 16.6, height = 8.8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_figure3.svg", plot_deviation_env,
        width = 16.6, height = 8.8,
        units = "cm", dpi = 300)

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

saveRDS(var_growth_cwm, "inst/var_growth_cwm.Rds")

plot_cwm_cwv_growth = var_growth_cwm %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(growth_weight))) +
  # 1:1 line only on CWM facet
  geom_abline(data = data.frame(cwm_name = "cwm", slope = 1, intercept = 0),
              aes(slope = slope, intercept = 0), linetype = 2) +
  # Rest of the geoms
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CW Mean", cwv = "CW Variance"))) +
  labs(x = "Environment",
       y = "Community-Weighted Value") +
  scale_colour_manual(name = "Trait Contribution\nto growth",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw(10) +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_growth

saveRDS(plot_cwm_cwv_growth, "inst/figures/paper_figure4.Rds")

ggsave2("inst/figures/paper_figure4.pdf", plot_cwm_cwv_growth,
        width = 12.7, height = 8.8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_figure4.png", plot_cwm_cwv_growth,
        width = 12.7, height = 8.8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_figure4.svg", plot_cwm_cwv_growth,
        width = 12.7, height = 8.8,
        units = "cm", dpi = 300)

# Table 2: model the influence of each parameter -------------------------------

cwm_env_model = all_trait_env %>%
  filter(hierar_exp == 0.5) %>%
  select(-site, -time, -trait_comb, -param_comb, -hierar_exp, -comb) %>%
  mutate(A = factor(A),
         H = factor(H)) %>%
  tidyr::nest(cwm_env = c(k, A, H, cwm_value, mismatch_value, env)) %>%
  mutate(
    cwm_mod = purrr::map(
      cwm_env, function(x) lm(cwm_value ~ env + A*H, data = x)),
    mis_mod = purrr::map(
      cwm_env, function(x) lm(mismatch_value ~ env + A*H, data = x)),
    cwm_coef = purrr::map(cwm_mod, broom::tidy),
    mis_coef = purrr::map(cwm_mod, broom::tidy))

table2_figure = sjPlot::plot_models(cwm_env_model$cwm_mod,
                                    m.labels = cwm_env_model$cwm_type) +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = "top")

# Supp. Fig. 1: Parameter space ------------------------------------------------

# Generate parameter space
param_space = expand.grid(
  k          = c(2),
  A          = c(0, 1e-7, 1e-5, 1e-4, 3e-4, 5e-4, 7e-4, 1e-3),
  H          = c(0, 1e-4, 1e-3, 1e-2, 5e-2, 1e-1, 5e-1),
  trait_comb = trait_comb[1:30]
)

param_space_simuls <- vector("list", nrow(param_space))

# Simulations of varying trait growth contribution
for (i in seq_len(nrow(param_space))) {

  # Get trait dataset
  traits <- data.frame(trait1 = uncor_traits[, param_space[i, "trait_comb"]])

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
    lim_sim_exponent = 1, hierar_exponent = 2)

  param_space_simuls[[i]] <- param_space_simul
}

param_space_abund <- purrr::map2_dfr(
  param_space_simuls, seq_len(nrow(param_space)), function(x, y) {
    # Extract final abundances
    x$compo[,,100] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("patch") %>%
      tidyr::gather("species", "abundance", -patch) %>%
      mutate(k = param_space[y, "k"],
             A = param_space[y, "A"],
             H = param_space[y, "H"],
             trait_comb = param_space[y, "trait_comb"])

  })

saveRDS(param_space_abund, "inst/param_space_abund.Rds")

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
  theme_bw(10) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45),
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
  theme_bw(10) +
  theme(aspect.ratio = 1,
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45),
        legend.position = "top")

plot_param_space = cowplot::plot_grid(plot_param_space_rich,
                                      plot_param_space_abund, labels = "AUTO")


saveRDS(plot_param_space, "inst/figures/paper_supp_fig1_param_space.Rds")

ggsave2("inst/figures/paper_supp_fig1_param_space.pdf", plot_param_space,
        width = 16.6, height = 8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig1_param_space.png", plot_param_space,
        width = 16.6, height = 8,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig1_param_space.svg", plot_param_space,
        width = 16.6, height = 8,
        units = "cm", dpi = 300)

# Supp. Fig. 2: effect of hierarchical exponent --------------------------------
plot_hierar_exp_species_mismatch = allmisA %>%
    filter(comb == scenarios[["hier"]]) %>%
    select(Species, MisAbPatchP, MisinstRPatchP, MismaxGRPatchP,
           MisavgGRPatchP, comb, hierar_exp) %>%
    tidyr::gather("mismatch_name", "mismatch_value",
                  MisAbPatchP:MisavgGRPatchP) %>%
    inner_join(all_minmax_mismatch, by = c("Species", "comb", "hierar_exp")) %>%
    mutate(comb = factor(comb,
                         levels = c(scenarios[["none"]], scenarios[["limsim"]],
                                    scenarios[["hier"]]))) %>%
    ggplot(aes(mismatch_value, Species, color = mismatch_name,
               shape = mismatch_name)) +
    geom_segment(aes(x = min_mismatch, xend = max_mismatch,
                     yend = Species),
                 color = "grey", size = 2/3) +
    geom_vline(xintercept = 0, linetype = 1, size = 1/2, color = "black") +
    geom_point(size = 1.5) +
    facet_wrap(vars(hierar_exp), ncol = 3,
               labeller = as_labeller(label_hierar_exp, default = label_parsed)) +
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
    labs(x = "Patch Mismatch",
         shape = "Performance\nMeasure",
         color = "Performance\nMeasure") +
    theme_bw(10) +
    theme(aspect.ratio = 1,
          panel.grid.major.x = element_line(size = 1.3),
          panel.spacing.x = unit(4, "mm"),
          panel.border = element_blank(),
          strip.background = element_blank(),
          axis.text = element_text(size = 8),
          legend.position = "top")

plot_hierar_exp_species_mismatch

saveRDS(plot_hierar_exp_species_mismatch,
        "inst/figures/paper_supp_fig2_hierar_exponent.Rds")

ggsave2("inst/figures/paper_supp_fig2_hierar_exponent.pdf",
        plot_hierar_exp_species_mismatch,
        width = 16.6, height = 12, scale = 1.5,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig2_hierar_exponent.png",
        plot_hierar_exp_species_mismatch,
        width = 16.6, height = 12, scale = 1.5,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig2_hierar_exponent.svg",
        plot_hierar_exp_species_mismatch, scale = 1.5,
        width = 16.6, height = 12,
        units = "cm", dpi = 300)

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

# Simulations of varying trait contribution to competition
for (i in seq_len(nrow(var_comp_comb))) {

  # Get trait dataset
  traits <- data.frame(trait1 = uncor_traits[,var_comp_comb[i, "trait_comb"]],
                       trait2 = second_trait)

  # Get trait contribution scenario
  comp_scenar <- var_comp_scenars[[var_comp_comb[i, "comp_scenar"]]]

  A = ifelse(comp_scenar$compet_weight[[1]] == 0, 0, list_A[[2]])
  H = ifelse(comp_scenar$hierarchy_weight[[1]] == 0, 0, list_H[[2]])

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

saveRDS(var_comp_cwm, "inst/var_comp_cwm.Rds")

plot_cwm_cwv_comp = var_comp_cwm %>%
  filter(hier_weight == 0) %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(comp_weight))) +
  # 1:1 line only on CWM facet
  geom_abline(data = data.frame(cwm_name = "cwm", slope = 1, intercept = 0),
              aes(slope = slope, intercept = 0), linetype = 2) +
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CWM", cwv = "CWV"))) +
  labs(x = "Environment",
       y = "CWM or CWV value") +
  scale_colour_manual(name = "Trait Contribution\nto limiting similarity",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw(10) +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_hier = var_comp_cwm %>%
  filter(comp_weight == 0) %>%
  tidyr::gather("cwm_name", "cwm_value", cwm, cwv) %>%
  ggplot(aes(env, cwm_value, color = as.factor(hier_weight))) +
  # 1:1 line only on CWM facet
  geom_abline(data = data.frame(cwm_name = "cwm", slope = 1, intercept = 0),
              aes(slope = slope, intercept = 0), linetype = 2) +
  stat_summary(geom = "smooth") +
  facet_wrap(vars(cwm_name), scales = "free_y",
             labeller = labeller(cwm_name = c(cwm = "CWM", cwv = "CWV"))) +
  labs(x = "Environment",
       y = "CWM or CWV value") +
  scale_colour_manual(name = "Trait Contribution\nto hierarchical comp.",
                      values = RColorBrewer::brewer.pal(7, "YlOrRd")[2:7],
                      labels = function(x) scales::percent(as.numeric(x))) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw(10) +
  theme(aspect.ratio = 1,
        legend.position = "top",
        strip.background = element_blank(),
        panel.grid = element_blank())

plot_cwm_cwv_comp_contrib = cowplot::plot_grid(
  plot_cwm_cwv_comp,
  plot_cwm_cwv_hier,
  labels = "AUTO",
  ncol = 1, nrow = 2)

saveRDS(plot_cwm_cwv_comp_contrib,
        "inst/figures/paper_supp_fig3_comp_contrib.Rds")

ggsave2("inst/figures/paper_supp_fig3_comp_contrib.pdf",
        plot_cwm_cwv_comp_contrib,
        width = 16.6, height = 16.6,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig3_comp_contrib.png",
        plot_cwm_cwv_comp_contrib,
        width = 16.6, height = 16.6,
        units = "cm", dpi = 300)

ggsave2("inst/figures/paper_supp_fig3_comp_contrib.svg",
        plot_cwm_cwv_comp_contrib,
        width = 16.6, height = 16.6,
        units = "cm", dpi = 300)


# Supp. Fig 4: Effect of limiting similarity strength of trait variance -------
# Please refer to script inst/02-figure_s4_limsim.R

