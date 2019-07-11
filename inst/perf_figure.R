# Script make performance comparison figure
# Packages ---------------------------------------------------------------------
library("dplyr")

# Load Data --------------------------------------------------------------------
# Select data for k = 1.3 A = 2.5e-7, B = 6.3e-6 and variable H (not 0 nor 1e-8)
all_perf_df = paste0("inst/job_data/perf_2fd398/perf_data/perf_df_",
                     all_param_df %>%
                         filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B > 6.3e-6,
                                B < 6.4e-6) %>%
                         pull(file_number), ".Rds") %>%
    purrr::map_dfr(readRDS)

all_trait = readRDS("inst/job_data/perf_2fd398/bigmem_trait_df.Rds")


# Subset data ------------------------------------------------------------------

all_perf_df = all_perf_df %>%
    inner_join(all_trait %>%
                   mutate(seed = as.integer(seed)),
               by = c("species", "seed", "trait_cor"))

# Single Parameter Set CWMs
single_cwm = all_cwm %>%
    filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B > 6.3e-6, B < 6.4e-6, H > 6.3e-6,
           H < 6.4e-6, trait_cor == "uncor", patch >= 5, patch <= 20)

seed_df = all_perf_df %>%
    distinct(seed) %>%
    slice(1:10)


# Figure 2: Abundance environment curve ----------------------------------------
# Effect of simulation parameters on single species abundance along the
# environment.
# (Abundance vs. Patch in function of values of parameters)
# Four panels figure: color code for parameter values k, A, B, H
subset_perf = all_perf_df %>%
    filter(seed == 227470, trait_cor == "uncor", species == "species55",
           R_scenar == 50, A_scenar == 50, H_scenar == 50)

all_perf_df %>%
    filter(seed == 227470, trait_cor == "uncor", species == "species48",
           R_scenar == 0, A_scenar == 50, H_scenar == 0) %>%
    ggplot(aes(patch, N150, color = as.factor(H))) +
    geom_line(size = 1, alpha = 2/3) +
    facet_wrap(~species, scales = "free_y") +
    scale_color_viridis_d() +
    labs(x = "Environment",
         y = "Final Abundance",
         color = "Hierarch. Compet.\nIntensity") +
    theme(aspect.ratio = 1,
          legend.position = "top")


# Abundance of species in a patch in function of H
all_perf_df %>%
    filter(seed == 227470, trait_cor == "uncor",
           R_scenar == 0, A_scenar == 50, H_scenar == 0, patch == 12,
           N150 > 0)  %>%
    select(H, species, N150) %>%
    arrange(H, desc(N150)) %>%
    ggplot(aes(H, N150, color = species)) +
    geom_line(size = 1) +
    scale_x_log10() +
    scale_y_log10()

# Figure 3: CWM-Environment figure ---------------------------------------------
# Effect of Trait contribution scenario on CWM_Trait2-Environment Relationship
# (CWM vs. Patch, 9 panels, R_scenar on rows, A_scenar on columns,
# H_scenar color)
single_cwm %>%
    filter(patch >= 5, patch <= 20) %>%
    ggplot(aes(patch, trait2_cwm, color = as.factor(H_scenar),
               # specify interaction otherwise it messes with groups
               group = interaction(as.factor(H_scenar), seed))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(size = 1, alpha = 1/20) +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(aes(label = ..eq.label.., group = NULL),
                          formula = y ~ x, parse = TRUE, size = 3,
                          label.y.npc = "bottom", label.x.npc = "right",
                          coef.digits = 2) +
    facet_grid(vars(R_scenar), vars(A_scenar),
               labeller = labeller(
                   R_scenar = function(x) paste0("Trait 2 Contrib.\n",
                                                 "to Growth : ", x, "%"),
                   A_scenar = function(x) paste0("Trait 2 Contrib.\nto ",
                                                 "Compet. : ", x, "%"),
                   .default = label_both)) +
    scale_color_viridis_d() +
    labs(x = "Environmental Variable",
         y = "CWM of Trait 2",
         color = "Trait 2 contrib. to Hierarch. Compet.",
         subtitle = expression("Trait 2 CWM <-> Environment (t"["opt"] *
                                   ") relationship")) +
    theme(aspect.ratio = 1,
          legend.position = "top")

# Test the influence of parameters
all_cwm %>%
    filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B > 6.3e-6, B < 6.4e-6, H > 6.3e-6,
           trait_cor == "uncor", patch >= 5, patch <= 20) %>%
    ggplot(aes(patch, trait2_cwm, color = as.factor(H_scenar),
               # specify interaction otherwise it messes with groups
               group = interaction(as.factor(H_scenar), seed))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_line(size = 1, alpha = 1/20) +
    geom_smooth(size = 1, method = "lm", se = FALSE, aes(group = NULL)) +
    ggpmisc::stat_poly_eq(aes(label = ..eq.label.., group = NULL),
                          formula = y ~ x, parse = TRUE, size = 3,
                          label.y.npc = "bottom", label.x.npc = "right",
                          coef.digits = 2) +
    facet_grid(vars(R_scenar), vars(H, A_scenar),
               labeller = labeller(
                   R_scenar = function(x) paste0("Trait 2 Contrib.\n",
                                                 "to Growth : ", x, "%"),
                   A_scenar = function(x) paste0("Trait 2 Contrib.\nto ",
                                                 "Compet. : ", x, "%"),
                   .default = label_both)) +
    scale_color_viridis_d() +
    labs(x = "Environmental Variable",
         y = "CWM of Trait 2",
         color = "Trait 2 contrib. to Hierarch. Compet.",
         subtitle = expression("Trait 2 CWM <-> Environment (t"["opt"] *
                                   ") relationship")) +
    theme(aspect.ratio = 1,
          legend.position = "top")

# Compute Different Performances indices ---------------------------------------

perf_subset = all_perf_df %>%
    filter(patch >= 5, patch <= 20, trait_cor == "uncor", R_scenar == 50,
           A_scenar == 50, H_scenar == 50) %>%
    semi_join(seed_df, by = "seed") %>%
    group_by(k, A, B, H, R_scenar, A_scenar, H_scenar, seed, patch)

tictoc::tic()
perf_ind = perf_subset %>%
    summarise(monocult_trait2 = weighted.mean(trait2, env_growth_rate,
                                              na.rm = TRUE),
              polycult_trait2 = weighted.mean(
                  trait2,
                  ifelse(is.na(max_growth_rate), 0, max_growth_rate),
                  na.rm = TRUE),
              cwm_trait2      = weighted.mean(trait2, N150, na.rm = TRUE))
tictoc::toc()

# Other definition of monoculture performance index
perf_best_growth = perf_subset %>%
    filter(env_growth_rate == max(env_growth_rate)) %>%
    rename(monobest_trait2 = trait2) %>%
    select(monobest_trait2)

perf_poly_growth = perf_subset %>%
    filter(max_growth_rate == max(max_growth_rate, na.rm = TRUE)) %>%
    rename(polybest_trait2 = trait2) %>%
    select(polybest_trait2)

perf_ind_all = perf_ind %>%
    full_join(perf_best_growth) %>%
    full_join(perf_poly_growth)
# Figure 4: Difference between performance indices -----------------------------
# Difference in performance indices affected by the intensity of ecological
# processes
# Difference between CWM and monoculture (trait of species with highest purely
# environmental growth rate per site)
# Difference between CWM and polyculture (trait of species with highest realized
# growth rate per site)
# Figure of estimates of mixed model explaining these differences against
# simulation parameters (k, A, B, H) (+random effect on seed)

perf_diff = perf_ind_all %>%
    group_by(patch, add = TRUE) %>%
    summarise(diff_cwm_mono  = cwm_trait2 - monobest_trait2,
              diff_cwm_poly  = cwm_trait2 - polybest_trait2)


# Difference between CWM to Monoculture
mod_diff_cwm_mono = perf_diff %>%
    ungroup() %>%
    mutate_at(vars(A, B, H),
              .funs = list(log = ~as.numeric(scale(log(. + 1e-16))))) %>%
    lme4::lmer(diff_cwm_mono ~ k + A + B + H_log + (1|seed), data = .)

sjPlot::plot_model(mod_diff_cwm_mono, show.values = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#001c00", size = 1,
               alpha = 1/2) +
    labs(title = "Parameters Effect on Difference between CWM - Monoculture")

# Difference between CWM to Polyculture
mod_diff_cwm_poly = perf_diff %>%
    ungroup() %>%
    mutate_at(vars(A, B, H),
              .funs = list(log = ~as.numeric(scale(log(. + 1e-16))))) %>%
    lme4::lmer(diff_cwm_poly ~ k + A + B + H_log + (1|seed), data = .)

sjPlot::plot_model(mod_diff_cwm_poly, show.values = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#001c00", size = 1,
               alpha = 1/2) +
    labs(title = "Parameters Effect on Difference between CWM - Polyculture")

# Plot Figures -----------------------------------------------------------------
perf_ind_all %>%
    tidyr::gather("perf_name", "perf_value", ends_with("trait2")) %>%
    ggplot(aes(patch, perf_value, color = perf_name)) +
    facet_wrap(~perf_name, scales = "free_y") +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 1/3)

perf_ind_all %>%
    filter(polycult_trait2 > -40, polycult_trait2 < 40) %>%
    ggplot(aes(polycult_trait2, polybest_trait2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point()

perf_ind_all %>%
    ggplot(aes(monocult_trait2, monobest_trait2)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point()

single_cwm %>%
    ggplot(aes(patch, trait1_cwm, color = as.factor(seed))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    stat_smooth(method = "lm", geom = "line", size = 1)

perf_diff %>%
    tidyr::gather("diff_name", "diff_value", starts_with("diff")) %>%
    filter(diff_name == "diff_cwm_mono") %>%
    ggplot(aes(interaction(R_scenar, H_scenar), diff_value)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_point(alpha = 1/5, position = position_jitter(0.2)) +
    facet_wrap(~diff_name, scales = "free_y") +
    ylim(-10, 10)


