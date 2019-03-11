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
