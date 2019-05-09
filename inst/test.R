var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 100, A_scenar == 0, H_scenar == 0,
           k == "fixed = 1.3", H == 0.4, A == 1e-5, N150 > 0,
           species %in% paste0("species", c(1, 25, 50, 100))) %>%
    ggplot(aes(env_growth_rate, max_growth_rate)) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = 2) +
        facet_wrap(~species)

var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99))) %>%
    gather("gr_type", "gr_value", env_growth_rate, max_growth_rate) %>%
    ggplot(aes(patch, gr_value, color = as.factor(seed), group = seed)) +
    geom_line(size = 1) +
    facet_grid(vars(gr_type), vars(species), scales = "free_y",
               labeller = labeller(gr_type =
                                       c(env_growth_rate = "Theoretical Growth Rate",
                                         max_growth_rate = "Realized Growth Rate"))) +
    labs(x = "Environment",
         y = "Growth Rate") +
    theme(legend.position = "none",
          aspect.ratio = 1)

var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99))) %>%
    ggplot(aes(env_growth_rate, max_growth_rate, color = as.factor(seed), group = seed)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_grid(NULL, vars(species), scales = "free_y") +
    labs(x = "Theoretical Growth Rate",
         y = "Realized Growth Rate") +
    theme(legend.position = "none",
          aspect.ratio = 1)

var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99)), N150 > 0) %>%
    ggplot(aes(distance_to_optimum, N150, color = as.factor(seed))) +
    geom_point() +
    facet_grid(NULL, vars(species), scales = "free_y") +
    labs(x = "Distance to Optimum",
         y = "Abundance") +
    theme(legend.position = "none",
          aspect.ratio = 1)


null_simuls = meta_simul_A_B_k_and_H(seed_number = 1,
                                     given_k = 1.3,
                                     given_A =  0,
                                     given_B = 1e-4,
                                     given_H =  0,
                                     given_scenars = scenar_list,
                                     given_h_fun = "sum",
                                     given_di_thresh = 25)
null_simul_perfs = null_simuls %>%
    map_dfr(
        ~extract_performances_from_simul(.x, used_trait_list[1], TRUE))

null_simul_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0,
           species %in% paste0("species", c(1, 25, 50, 99)), N150 > 0) %>%
    ggplot(aes(distance_to_optimum, N150)) +
    geom_point() +
    facet_grid(NULL, vars(species), scales = "free_y") +
    labs(x = "Distance to Optimum",
         y = "Abundance") +
    theme(legend.position = "none",
          aspect.ratio = 1) +
    scale_y_log10()


real_dist_abund = var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99)), N150 > 0) %>%
    ggplot(aes(distance_to_optimum, N150)) +
    geom_point(aes(color = as.factor(seed))) +
    geom_line(data = null_simul_perfs %>%
                  filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0,
                         species %in% paste0("species", c(1, 25, 50)), N150 > 0) %>%
                  select(species, N150, distance_to_optimum),
              color = "black", size = 1) +
    geom_point(data = null_simul_perfs %>%
                   filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0,
                          species %in% paste0("species", c(1, 25, 50)), N150 > 0) %>%
                   select(species, N150, distance_to_optimum),
               color = "black", size = 1) +
    facet_grid(NULL, vars(species)) +
    labs(x = "Distance to Optimum",
         y = "Abundance at t = 75",
         color = "seed") +
    theme(aspect.ratio = 1) +
    scale_y_log10()

null_simul_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0,
           species %in% paste0("species", c(1, 25, 50)), N150 > 0) %>%
    select(species, N150, distance_to_optimum) -> null_simul_dist

var_scenars_perfs %>%
    filter(di_thresh == 25, R_scenar == 0, A_scenar == 0, H_scenar == 0, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99)), N150 > 0) %>%
    select(seed, species, N150, distance_to_optimum) -> real_simul_dist

real_simul_dist %>%
    inner_join(null_simul_dist %>%
                   rename(null_abund = N150),
               by = c("species", "distance_to_optimum")) %>%
    mutate(th_obs_diff = null_abund - N150) %>%
    ggplot(aes(distance_to_optimum, th_obs_diff, color = as.factor(seed))) +
    geom_point() +
    facet_grid(NULL, vars(species))

var_scenars_perfs %>%
    filter(di_thresh == 25, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           species %in% paste0("species", c(1, 25, 50, 99)),
           H_scenar == 0, A_scenar == 0) %>%
    ggplot(aes(patch, env_growth_rate, color = as.factor(seed))) +
    geom_line(size = 1) +
    facet_grid(vars(R_scenar), vars(species), labeller = label_both) +
    labs(x = "Environment",
         y = "Theoretical Growth Rate (environmental)",
         color = "Seed")

var_scenars_perfs %>%
    filter(di_thresh == 25, H == 0.4, A == 1e-5, k == "fixed = 1.3",
           H_scenar == 0, A_scenar == 0, R_scenar == 0) %>%
    ggplot(aes(patch, max_growth_rate, group = interaction(species, seed), color = as.factor(seed))) +
    geom_line() +
    facet_wrap(~seed) +
    labs(x = "Environment",
         y = "Realized Growth Rate",
         color = "Seed")

# Answering Caroline's email ---------------------------------------------------

# Abundance vs. Environment per species in polyculture
var_scenars_perf %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5,
           species %in% paste0("species", c(1, 25, 50, 75, 99)),
           R_scenar == 0, A_scenar == 0, H_scenar == 0) %>%
    ggplot(aes(patch, N150, color = as.factor(seed))) +
    geom_line(size = 1) +
    facet_wrap(~species) +
    scale_y_log10()


# Estimating optimal environment using abundance
patch_optim = var_scenars_perf %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5,
           N150 > 0) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, species) %>%
    filter(N150 == max(N150)) %>%
    summarise(patch_optim = patch)

# Estimating species optimal environment using growth rate
patch_optim_growth = var_scenars_perf %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5,
           N150 > 0) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, species) %>%
    filter(max_growth_rate == max(max_growth_rate)) %>%
    summarise(patch_optim = patch)

# Plot of species optimal environment vs. trait values
# + CWM
patch_optim %>%
    ungroup() %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, species) %>%
    mutate(growth_trait = weighted.mean(c(trait1, trait2),
                                        c(100 - R_scenar, R_scenar))) %>%
    gather("trait_name", "trait_value", trait2, growth_trait) %>%
    ggplot(aes(patch_optim, trait_value, group = interaction(seed, trait_name,
                                                             H_scenar),
               linetype = trait_name, color = as.factor(H_scenar))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    stat_smooth(geom = "line", method = "lm", se = FALSE, alpha = 1/10,
                color = "black") +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    geom_smooth(data = var_scenars_cwm %>%
                    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3",
                           A == 1e-5),
                size = 1, method = "lm", se = FALSE,
                aes(x = patch, y = trait2_cwm, group = NULL, linetype = "trait2_cwm")) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.3,
                          aes(label = ..eq.label.., group = NULL),
                          label.x.npc = "right", label.y.npc = "bottom") +
    facet_grid(vars(R_scenar), vars(A_scenar), labeller = label_both) +
    labs(x = "Environment with max. abundance",
         y = "Trait Value")

# Same but for optimal environment given by growth rate
patch_optim_growth %>%
    ungroup() %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, species) %>%
    mutate(growth_trait = weighted.mean(c(trait1, trait2),
                                        c(100 - R_scenar, R_scenar))) %>%
    gather("trait_name", "trait_value", trait2, growth_trait) %>%
    ggplot(aes(patch_optim, trait_value, group = interaction(seed, trait_name,
                                                             H_scenar),
               linetype = trait_name, color = as.factor(H_scenar))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    # stat_smooth(geom = "line", method = "lm", se = FALSE, alpha = 1/10,
    #             color = "black") +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.3,
                          aes(label = ..eq.label.., group = NULL),
                          label.x.npc = "right", label.y.npc = "bottom") +
    facet_grid(vars(R_scenar), vars(A_scenar), labeller = labeller(
        R_scenar = function(x) paste0("Trait 2 Contrib.\n",
                                      "to Growth : ", x, "%"),
        A_scenar = function(x) paste0("Trait 2 Contrib.\nto ",
                                      "Compet. : ", x, "%"),
        .default = label_both)) +
    scale_linetype_discrete(name = "Trait Category",
                            labels = c(growth_trait = "Composite Trait",
                                       trait2 = "Trait 2 Polyculture")) +
    labs(x = "Environment with max. growth rate",
         y = "Trait Value",
         color = "Contrib. to\nHierarchical Compet.") +
    theme(aspect.ratio = 1)

# Extract values of environmental growth
patch_optim_env = var_scenars_perf %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, species) %>%
    filter(env_growth_rate == max(env_growth_rate)) %>%
    summarise(patch_optim = patch) %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    group_by(species, add = TRUE) %>%
    mutate(growth_trait = weighted.mean(c(trait1, trait2),
                                        c(100 - R_scenar, R_scenar))) %>%
    gather("trait_name", "trait_value", trait2, growth_trait)

patch_optim_env %>%
    ggplot(aes(patch_optim, trait_value,
               group = interaction(seed, trait_name, H_scenar),
               linetype = trait_name, color = as.factor(H_scenar))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(aes(shape = trait_name)) +
    stat_smooth(geom = "line", method = "lm", se = FALSE, alpha = 1/10,
                color = "black") +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.3,
                          aes(label = ..eq.label.., group = NULL),
                          label.x.npc = "right", label.y.npc = "bottom") +
    facet_grid(vars(R_scenar), vars(A_scenar), labeller = label_both) +
    labs(x = "Environment with max. growth rate",
         y = "Trait Value")

# Instead of selecting the environment where each species exhibit its maximum
# growth we can select the species that has the highest growth rate for each
# environment
species_optim = var_scenars_perf %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5, N150 > 0) %>%
    group_by(R_scenar, A_scenar, H_scenar, seed, patch) %>%
    filter(env_growth_rate == max(env_growth_rate)) %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    group_by(patch, species, add = TRUE) %>%
    mutate(growth_trait = weighted.mean(c(trait1, trait2),
                                        c(100 - R_scenar, R_scenar))) %>%
    gather("trait_name", "trait_value", trait2, growth_trait)

species_optim %>%
    ggplot(aes(patch, trait_value,
               group = interaction(seed, trait_name, H_scenar),
               linetype = trait_name, color = as.factor(H_scenar))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    stat_smooth(geom = "line", method = "lm", se = FALSE, alpha = 1/10,
                color = "black") +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.3,
                          aes(label = ..eq.label.., group = NULL),
                          label.x.npc = "right", label.y.npc = "bottom") +
    facet_grid(vars(R_scenar), vars(A_scenar), labeller = label_both) +
    labs(x = "Environment with max. growth rate",
         y = "Trait Value")

# Now focusing on the change in parameter set ----------------------------------

param_patch_optim = var_scenars_perf %>%
    filter(di_thresh == 25, R_scenar == 50, A_scenar == 50, H_scenar == 50,
           N150 > 0) %>%
    group_by(A, H, k, seed, species) %>%
    filter(N150 == max(N150)) %>%
    summarise(R_scenar = unique(R_scenar),
              patch_optim = patch) %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    group_by(A, H, k, seed, species) %>%
    mutate(growth_trait = weighted.mean(c(trait1, trait2),
                                        c(100 - R_scenar, R_scenar))) %>%
    gather("trait_name", "trait_value", trait2, growth_trait)

param_patch_optim_env =  var_scenars_perf %>%
    filter(di_thresh == 25, R_scenar == 50, A_scenar == 50, H_scenar == 50,
           N150 > 0) %>%
    group_by(A, H, k, seed, species) %>%
    filter(env_growth_rate == max(env_growth_rate)) %>%
    summarise(R_scenar = unique(R_scenar), patch_optim = patch) %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    rename(trait2_monocult = trait2) %>%
    gather("trait_name", "trait_value", trait2_monocult)

fig_param_polycultures = param_patch_optim %>%
    bind_rows(param_patch_optim_env) %>%  # Join Trait 2 monocultre
    bind_rows(var_scenars_cwm %>%
                  filter(di_thresh == 25, R_scenar == 50, A_scenar == 50,
                         H_scenar == 50) %>%
                  select(A, H, k, seed, R_scenar, patch, trait2_cwm) %>%
                  rename(patch_optim = patch) %>%
                  gather("trait_name", "trait_value", trait2_cwm) %>%
                  mutate(trait1 = NA_real_,
                         trait3 = NA_real_,
                         trait4 = NA_real_,
                         species = NA_character_)) %>%
    ungroup() %>%
    # Reorder trait categories for easier reading
    mutate(trait_name = fct_relevel(trait_name, "growth_trait",
                                    "trait2_monocult", "trait2",
                                    "trait2_cwm")) %>%
    # Actual plot
    ggplot(aes(patch_optim, trait_value, group = interaction(seed, trait_name),
               linetype = trait_name, color = trait_name)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    # Trait 2 polyculture & Composite Trait (individual seed)
    stat_smooth(geom = "line", method = "lm", se = FALSE, alpha = 1/10,
                color = "black") +
    # Trait 2 polyculture & Composite Trait (average across seeds)
    geom_smooth(size = 1.25, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE, size = 2.8,
                          aes(label = ..eq.label.., group = NULL),
                          label.x.npc = "right", label.y.npc = "bottom") +
    facet_grid(vars(k), vars(A, H), labeller = label_both) +
    scale_linetype_manual(name = "Trait Category",
                          values = c(growth_trait = 1,
                                     trait2_monocult = 3,
                                     trait2 = 5,
                                     trait2_cwm = 2),
                          labels = c(growth_trait = "Composite Trait",
                                     trait2_monocult = "Trait 2 Monoculture",
                                     trait2 = "Trait 2 Polyculture",
                                     trait2_cwm = "CWM of Trait 2")) +
    scale_color_discrete(name = "Trait Category",
                         labels = c(growth_trait = "Composite Trait",
                                    trait2_monocult = "Trait 2 Monoculture",
                                    trait2 = "Trait 2 Polyculture",
                                    trait2_cwm = "CWM of Trait 2")) +
    labs(x = "Environment",
         y = "Trait Value") +
    theme(aspect.ratio = 1)

# Compare individual growth ----------------------------------------------------
var_scenar_perfs %>%
    filter(di_thresh == 25, H == 0.4, k == "fixed = 1.3", A == 1e-5, seed == 1)  %>%
    filter(A_scenar == 0, H_scenar == 0) %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)),
               by = c("seed", "species")) %>%
    mutate(growth_trait = pmap_dbl(list(trait1, trait2, R_scenar),
                                   ~weighted.mean(c(..1, ..2),
                                                  c(100 - ..3, ..3))),
           nona_growth = ifelse(is.na(max_growth_rate),
                                0, max_growth_rate)) %>%
    group_by(R_scenar, patch) %>%
    summarise_at(vars(trait1, trait2, growth_trait),
                 list(env_growth = ~weighted.mean(., env_growth_rate),
                      max_growth = ~weighted.mean(., nona_growth))) %>%
    gather("growth_type", "growth_value", ends_with("growth")) %>%
    ggplot(aes(patch, growth_value, color = as.factor(R_scenar))) +
    geom_line(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    facet_wrap(~growth_type, scales = "free")

# Make figures as needed by CT -------------------------------------------------
pkgload::load_all()
library("tidyverse")

# caro_scenars = readRDS("inst/job_data/caroline_scenars.Rds")
caro_traits = readRDS("inst/job_data/caroline_traits.Rds")
caro_perfs = readRDS("inst/job_data/caroline_perfs.Rds")


full_trait_df = caro_traits %>%
    map_dfr(~map_dfr(.x, function(x) {
        x %>%
            as.data.frame() %>%
            rownames_to_column("species")
    }, .id = "trait_cor"), .id = "seed")

caro_full = caro_perfs %>%
    inner_join(full_trait_df %>%
                   mutate(seed = as.integer(seed)), by = c("species", "seed",
                                                           "trait_cor"))

caro_sum = caro_full %>%
    group_by(k, A, B, H, trait_cor, R_scenar, A_scenar, H_scenar, seed,
             patch) %>%
    summarise_at(
        vars(matches("trait[12]")),
        .funs = list(cwm       = ~weighted.mean(., N150,
                                                na.rm = TRUE),
                     w_max     = ~weighted.mean(., max_growth_rate,
                                                na.rm = TRUE),
                     w_max_cap = ~weighted.mean(., max_growth_rate_per_capita,
                                                na.rm = TRUE)))
# median scenar & no trait correlation
caro_med = caro_sum %>%
    filter(R_scenar == 50, A_scenar == 505, H_scenar == 0,
           trait_cor == "uncor") %>%
    ungroup()


# gather("trait_type", "trait_value", matches("trait[12].*")) %>%
#     separate("trait_type", c("trait_name", "measure_type"), sep = "_", extra = "merge") %>%

caro_indiv = caro_med %>%
    filter(B == 0, A == 0, H == 0) %>%
    select(k, seed, patch,
           trait1_w_max, trait2_w_max,
           trait1_w_max_cap, trait2_w_max_cap) %>%
    rename_at(vars(contains("w_max")), .funs = list(~gsub("w_max", "indiv", .)))

caro_mono = caro_med %>%
    filter(B != 0, A == 0, H == 0)

caro_poly = caro_med %>%
    filter(B != 0, A != 0) %>%
    inner_join(caro_indiv, by = c("k", "seed", "patch"))

caro_perf_env_plot = caro_poly %>%
    gather("trait_type", "measure_value", matches("trait[12].*")) %>%
    separate("trait_type", c("trait_name", "measure_type"),
             sep = "_", extra = "merge") %>%
    filter(trait_name == "trait2", B == 1e-6, A %in% c(1e-6, 1e-4),
           H %in% c(0, 1), !grepl("cap", measure_type, fixed = TRUE)) %>%
    ggplot(aes(patch, measure_value, color = measure_type,
               group = interaction(seed, measure_type))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(alpha = 1/5) +
    geom_smooth(method = "lm", se = FALSE, aes(group = measure_type)) +
    facet_grid(vars(k), vars(A, H), labeller = labeller(
        A = function(x) paste0("Limiting Sim. = ", x),
        H = function(x) paste0("Hierarch. Compet. = ", x),
        k = function(x) paste0("Max Growth = ", x)
    )) +
    scale_color_discrete(labels = c(cwm = "Trait 2 CWM",
                                    indiv = "Individual",
                                    w_max = "Polyculture Growth Rate")) +
    labs(x = "Environment",
         y = "Performance Measure",
         caption = paste0("B = 1e-6; n = 10; Median scenario ",
                          "(each trait contributing to all processes)"),
         title = "Performance Estimates vs. Environment") +
    theme(aspect.ratio = 1,
          strip.background = element_blank())

# Species richnes figures ------------------------------------------------------

sp_rich = caro_var_perfs %>%
    group_by(k, A, B, H, scenario, trait_cor, seed,
             patch) %>%
    summarise(species_rich = sum(N150 > 0, na.rm = TRUE))

sp_rich %>%
    filter(trait_cor == "uncor") %>%
    ggplot(aes(as.factor(k), as.factor(A))) +
    stat_summary_2d(aes(z = species_rich), fun = "mean", geom = "raster") +
    facet_grid(vars(H, B), vars(scenario), labeller = label_both) +
    scale_fill_viridis_c() +
    labs(fill = "Species Richness") +
    theme(legend.position = "top",
          legend.text = element_text(angle = 45),
          aspect.ratio = 1)

# Assemble all performances of simulations -------------------------------------

all_perfs = list.files("inst/job_data", "other_simuls_.*_perfs.Rds",
                       full.names = TRUE) %>%
    .[c(1, 11, 2:10, 12:13)] %>%
    purrr::walk(function(x) {
        given_data = readRDS(x)

        fst_name = gsub(".Rds", ".fst", x, fixed = TRUE)

        fst::write_fst(given_data, fst_name, compress = 100)
    })


full_trait_df = readRDS("inst/job_data/other_simuls_traits.Rds") %>%
    purrr::map_dfr(function(x) {
        purrr::map_dfr(x, ~.x %>%
                           as.data.frame() %>%
                           tibble::rownames_to_column("species"),
                       .id = "trait_cor")
    }, .id = "seed") %>%
    dplyr::mutate(seed = as.integer(seed))

all_perfs = list.files("inst/job_data", "other_simuls_env_2.*_perfs.Rds",
                       full.names = TRUE) %>%
    .[c(1, 11, 2:10, 12:13)] %>%
    purrr::map(function(x) {
        given_dt = readRDS(x) %>%
            dplyr::filter(N150 > 0)
    }) %>%
    dplyr::bind_rows() %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(full_trait_df, by = c("seed", "trait_cor", "species"))

all_cwm = all_perfs %>%
    group_by(k, A, B, H, R_scenar, A_scenar, H_scenar, trait_cor, seed, patch) %>%
    summarise_at(vars(dplyr::matches("trait[12]")),
                 .funs = list(cwm = ~weighted.mean(., N150, na.rm = TRUE)))

all_cwv = all_perfs %>%
    group_by(k, A, B, H, R_scenar, A_scenar, H_scenar, trait_cor, seed, patch) %>%
    summarise_at(vars(dplyr::matches("trait[12]")),
                 .funs = list(cwv = ~Hmisc::wtd.var(., N150, na.rm = TRUE)))

all_rich = all_perfs %>%
    filter(N150 > 0) %>%
    group_by(k, A, B, H, R_scenar, A_scenar, H_scenar, trait_cor, seed,
             patch) %>%
    summarise(species_rich = n()) %>%
    summarise(species_rich = mean(species_rich))

all_cwm %>%
    filter(k == 1.2, A == 1e-4, B == 1e-2, H == 0.5, trait_cor == "uncor") %>%
    ggplot(aes(patch, trait2_cwm, color = as.factor(H_scenar))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(aes(group = interaction(H_scenar, seed)), alpha = 1/5) +
    geom_smooth(aes(group = H_scenar), method = "lm", se = FALSE) +
    facet_grid(vars(R_scenar), vars(A_scenar),
               labeller = labeller(
                   R_scenar = function(x) paste0("Contrib. to\nGrowth: ", x, "%"),
                   A_scenar = function(x) paste0("Contrib. to\nLimit. Sim.: ", x, "%"))) +
    ylim(1, 25) +
    scale_color_discrete(labels = function(x) paste0(x, "%")) +
    theme_bw(14) +
    theme(aspect.ratio = 1, legend.position = "top") +
    labs(x = "Environment",
         y = "Trait CWM",
         color = "Contrib. to\nHierach. Compet.",
         caption = "10 seeds; k = 1.2; A = 1e-4; B = 1e-2; H = 0.5; uncorrelated traits")

all_cwm %>%
    filter(R_scenar == 50, A_scenar == 50, H_scenar == 50, B != 1e-4, A != 1e-5,
           trait_cor == "uncor") %>%
    ggplot(aes(patch, trait2_cwm, color = as.factor(H))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(aes(group = interaction(H, seed)), alpha = 1) +
    geom_smooth(aes(group = H), method = "lm", se = FALSE) +
    facet_grid(vars(k, B), vars(A), labeller = label_both) +
    ylim(1, 25) +
    theme_bw(14) +
    theme(aspect.ratio = 1, legend.position = "top") +
    labs(x = "Environment",
         y = "Trait CWM",
         color = "Hierach. Compet.\nIntensity",
         caption = "1 seed; R50A50H50; uncorrelated traits")

### Coexistence plot
param_labs = c(k = "Maximum Growth Rate (k)",
               A = "Limit. Sim. Intensity (A)",
               B = "Intra-sp. Compet. Intensity (B)",
               H = "Hierarchical Compet. Intensity (H)")

plot_rich_simul = function(rich_df, x_var, y_var, given_labs = param_labs) {
    rich_df %>%
        filter(trait_cor == "uncor") %>%
        ggplot(aes_string(x = paste0("as.factor(", x_var, ")"),
                          y = paste0("as.factor(", y_var, ")"),
                          z = "species_rich")) +
        stat_summary_2d(fun = mean) +
        scale_fill_viridis_c(limits = c(0, 70)) +
        labs(x = given_labs[x_var],
             y = given_labs[y_var],
             fill = "SR") +
        theme(aspect.ratio = 1)
}

fig_k_A = all_rich %>%
    plot_rich_simul("k", "A")

fig_k_H = all_rich %>%
    plot_rich_simul("k", "H")

fig_A_B = all_rich %>%
    plot_rich_simul("A", "B")

fig_A_H = all_rich %>%
    plot_rich_simul("A", "H")

fig_k_B = all_rich %>%
    plot_rich_simul("k", "B")

fig_B_H = all_rich %>%
    plot_rich_simul("B", "H")

library("cowplot")

plot_grid(fig_k_A, fig_k_H, fig_k_B, fig_B_H, fig_A_H, fig_A_B,
          ncol = 3, nrow = 2)

### Community level trait distribution

three_cwm = all_perfs %>%
    filter(k == 1.5, A == 1e-6, B == 1e-2, H == 0.5, trait_cor == "uncor",
           patch %in% c(1, 13, 25), seed == 1) %>%
    filter(N150 > 0) %>%
    mutate(scenario = paste0("R", R_scenar, "A", A_scenar, "H", H_scenar)) %>%
    filter(scenario %in% c("R100A0H0", "R0A100H0", "R0A0H100")) %>%
    mutate(scenario = fct_relevel(scenario, c("R100A0H0", "R0A100H0",
                                              "R0A0H100"))) %>%
    group_by(patch, scenario) %>%
    summarise(trait2_cwm = weighted.mean(trait2, N150))

all_perfs %>%
    filter(k == 1.5, A == 1e-6, B == 1e-2, H == 0.5, trait_cor == "uncor",
           patch %in% c(1, 13, 25), seed == 1) %>%
    filter(N150 > 0) %>%
    mutate(scenario = paste0("R", R_scenar, "A", A_scenar, "H", H_scenar)) %>%
    filter(scenario %in% c("R100A0H0", "R0A100H0", "R0A0H100")) %>%
    mutate(scenario = fct_relevel(scenario, c("R100A0H0", "R0A100H0",
                                              "R0A0H100"))) %>%
    ggplot(aes(trait2)) +
    geom_histogram(color = "white") +
    geom_vline(data = three_cwm, aes(xintercept = trait2_cwm, color = "CWM"),
               linetype = 2, size = 1) +
    facet_grid(vars(scenario), vars(patch),
               labeller = labeller(patch = label_both)) +
    labs(x = "Trait value",
         y = "Number of Species",
         colour = "",
         caption = "1 seed; k = 1.5; A = 1e-6; B = 1e-2; H = 0.5; uncorrelated traits") +
    scale_color_manual(values = c(CWM = "darkblue")) +
    theme(legend.position = "top",
          aspect.ratio = 1)

### CWV graph

plot_cwm_param = function(cwv_df, chosen_param, given_labs = param_labs) {
    cwv_df %>%
        ungroup() %>%
        filter(B == 1e-3, trait_cor == "uncor",
               seed == 1) %>%
        mutate(scenario = paste0("R", R_scenar, "A", A_scenar, "H", H_scenar)) %>%
        filter(scenario %in% c("R100A0H0", "R0A100H0", "R0A0H100")) %>%
        mutate(scenario = fct_relevel(scenario, c("R100A0H0", "R0A100H0",
                                                  "R0A0H100"))) %>%
        ggplot(aes_string(paste0("as.factor(", chosen_param, ")"),
                          "trait2_cwv")) +
        #geom_abline(slope = 1, intercept = 0, linetype = 2) +
        geom_boxplot() +
        facet_grid(cols = vars(scenario)) +
        labs(x = given_labs[chosen_param],
             y = "Trait CWV")

}

fig_cwv = c("A", "H", "k") %>%
    map(~plot_cwm_param(all_cwv, .x))

plot_grid(plotlist = fig_cwv, ncol = 1, nrow = 3, align = "hv")


### CWM with trait correlations
all_cwm %>%
    ungroup() %>%
    mutate(scenario = paste0("R", R_scenar, "A", A_scenar, "H", H_scenar)) %>%
    filter(k == 1.5, A == 1e-6, B == 1e-2, H == 0.5, R_scenar == 0,
           A_scenar == 0, H_scenar == 0) %>%
    mutate(trait_cor = fct_relevel(trait_cor, "negcor", "uncor", "poscor")) %>%
    ggplot(aes(patch, trait2_cwm, color = trait_cor)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_smooth(method = "lm", se = FALSE, size = 1) +
    scale_color_discrete(labels = c(uncor = "r = 0",
                                    poscor = "r = 0.3",
                                    negcor = "r = -0.3")) +
    labs(x = "Environment",
         y = "Trait CWM",
         color = "Trait\nCorrelation") +
    theme(aspect.ratio = 1,
          legend.position = "top")


### Abundance – Distance to Optimum figure
theme_set(theme_bw(12))

th_optim = all_perfs %>%
    filter(trait_cor == "uncor", B == 1e-2) %>%
    distinct(k, B) %>%
    mutate(th_abund = purrr::map2(
        k, B, ~data.frame(distance_to_optimum = seq(0, 25, length.out = 100)) %>%
            mutate(N150 = (1/.y) * (.x * exp(-(distance_to_optimum^2) /
                                               (2*2^2)) - 1)) %>%
            mutate(N150_corr = ifelse(N150 < 0, 0, N150))))

th_optim %>%
    unnest(th_abund) %>%
    gather("abund_type", "abund", N150, N150_corr) %>%
    ggplot(aes(distance_to_optimum, abund, color = as.factor(k))) +
    geom_line(aes(group = interaction(abund_type, k)), size = 1) +
    facet_grid(vars(abund_type))

all_perfs %>%
    filter(trait_cor == "uncor", B == 1e-2) %>%
    ggplot(aes(distance_to_optimum, N150)) +
    geom_point(shape = ".") +
    geom_line(data = th_optim %>%
                  unnest(th_abund), aes(y = N150_corr),
              color = "darkblue", size = 1) +
    facet_grid(vars(k, H), vars(A), labeller = label_both)

all_perfs %>%
    select(seed, trait_cor, h_fun, di_thresh, k, A, B, H, patch, species, N150, distance_to_optimum) %>%
    filter(B == 1e-2,trait_cor == "uncor") %>%
    mutate(th_n150 = (1/B) * (k * exp(-(distance_to_optimum^2) / (2 * 2^2)) - 1),
           th_n150_corr = ifelse(th_n150 < 0, 0, th_n150)) %>%
    mutate(sq_deviation = (N150 - th_n150_corr)^2) %>%
    group_by(k, A, B, H, trait_cor, seed) %>%
    summarise(mse = mean(sq_deviation), rmse = sqrt(mse)) %>%
    ungroup() %>%
    ggplot(aes(as.factor(k), rmse)) +
    ggbeeswarm::geom_quasirandom(width = 0.1) +
    labs(x = "Maximum Basal Growth Rate (k)",
         y = "RMSE of abundance vs. dopt") +
    theme(legend.position = "top",
          aspect.ratio = 1)


# Actual figures ---------------------------------------------------------------

single_simul %>%
    filter(species %in% paste0("species", c(1, 40, 50, 99))) %>%
    ggplot(aes(patch, env_growth_rate, color = species)) +
    geom_line(size = 1) +
    geom_vline(data = single_simul %>%
                   filter(species %in% paste0("species", c(1, 40, 50, 99))) %>%
                   group_by(species) %>%
                   filter(env_growth_rate == max(env_growth_rate)), aes(xintercept = patch), linetype = 2)

# Obtaining individual T by E
single_simul = var_perfs %>%
    filter(trait_cor == "uncor", R_scenar == 100, A_scenar == 100,
           H_scenar == 100, seed == 1, k == 1.2, H == 0, B == 1e-4, A == 0)

base_simul = var_perfs %>%
    filter(trait_cor == "uncor", R_scenar == 100, A_scenar == 100,
           H_scenar == 100, seed == 1, k == 1.2, H == 0, B == 0, A == 0)

plot_individual = single_simul %>%
    group_by(species) %>%
    filter(env_growth_rate == max(env_growth_rate)) %>%
    ungroup() %>%
    ggplot(aes(patch, trait2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    labs(x = "Environment",
         y = "Focal Trait")

plot_monoculture = single_simul %>%
    group_by(species) %>%
    filter(N150 > 0, N150 == max(N150)) %>%
    ungroup() %>%
    ggplot(aes(patch, trait2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    labs(x = "Environment",
         y = "Focal Trait")

plot_max_growth = single_simul %>%
    group_by(species) %>%
    filter(max_growth_rate == max(max_growth_rate)) %>%
    ungroup() %>%
    ggplot(aes(patch, trait2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    labs(x = "Environment",
         y = "Focal Trait")

single_simul %>%
    bind_rows(base_simul) %>%
    filter(species %in% paste0("species", c(1, 40, 50, 99))) %>%
    ggplot(aes(patch, N150, color = species)) +
    geom_line() +
    facet_wrap(~B, scales = "free_y")

single_simul %>%
    bind_rows(base_simul) %>%
    group_by(patch, B) %>%
    summarise(cwm = weighted.mean(trait2, N150)) %>%
    ggplot(aes(patch, cwm, color = as.factor(B), group = B)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line() +
    lims(y = c(1, 25))

var_perfs %>%
    filter(trait_cor == "uncor", seed == 1, R_scenar == 100, A_scenar == 100,
           H_scenar == 100, k == 1.2, B == 1e-4) %>%
    group_by(A, H, patch) %>%
    summarise(cwm = weighted.mean(trait2, N150)) %>%
    ungroup() %>%
    ggplot(aes(patch, cwm, color = interaction(A, H))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(size = 0.9) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    scale_color_discrete(
        name = "CWM Type",
        labels = c("0.0" = "Monoculture",
                   "1e-06.0" = "Community\nOnly Limit. Sim.",
                   "0.1" = "Community\nOnly Hierarch. Compet.",
                   "1e-06.1" = "Community\nLimit. Sim. + Hierarch. Compet."))
var_perfs %>%
    filter(trait_cor == "uncor", R_scenar == 100, A_scenar == 100,
           H_scenar == 100, k == 1.2, B == 1e-4) %>%
    group_by(seed, A, H, patch) %>%
    summarise(cwm = weighted.mean(trait2, N150)) %>%
    ungroup() %>%
    ggplot(aes(patch, cwm, color = interaction(A, H))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(size = 0.9) +
    geom_smooth(data = var_perfs %>%
                    filter(trait_cor == "uncor", R_scenar == 100, A_scenar == 100,
                           H_scenar == 100, k == 1.2, B == 1e-4, A != 0 | H != 0) %>%
                    group_by(seed, A, H, species) %>%
                    filter(N150 == max(N150)) %>%
                    select(patch, trait2),
                aes(patch, trait2), color = "darkred", size = 0.9, method = "lm", se = FALSE) +
    ggpmisc::stat_poly_eq(formula = y ~ x, parse = TRUE) +
    scale_color_discrete(
        name = "CWM Type",
        labels = c("0.0" = "Monoculture",
                   "1e-06.0" = "Community\nOnly Limit. Sim.",
                   "0.1" = "Community\nOnly Hierarch. Compet.",
                   "1e-06.1" = "Community\nLimit. Sim. + Hierarch. Compet.")) +
    facet_wrap(~seed)
