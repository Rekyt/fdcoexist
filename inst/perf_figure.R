# Script make performance comparison figure
# Packages ---------------------------------------------------------------------
library("dplyr")
library("furrr")
library("cowplot")
pkgload::load_all()

plan(multiprocess, workers = 15)

# Functions ---------------------------------------------------------------------
theme_set(theme_bw(14))
# Function so that facet labels are properly formatted
scientific_notation = function(x, label = "H") {
    paste0(label, ": ", scales::scientific_format()(as.numeric(x)))
}

plot_env_abund = function(perf_df, variable = "H", legend_label) {
    subset_perf %>%
        ggplot(aes_string("patch", "N150",
                          color = paste0("as.factor(", variable, ")"))) +
        geom_line(size = 1, alpha = 1/3) +
        #geom_line(size = 1.3, alpha = 1/3, aes(group = NULL)) +
        facet_wrap(~seed, scales = "free_y") +
        scale_color_viridis_d(labels = function(x) x %>%
                                  as.numeric() %>%
                                  scales::scientific()) +
        labs(x = "Environment",
             y = "Final Abundance",
             color = legend_label) +
        theme(aspect.ratio = 1,
              legend.position = "top")
}

# Load Data --------------------------------------------------------------------
# Select data for k = 1.3 A = 2.5e-7, B = 6.3e-6 and variable H (not 0 nor 1e-8)

main_folder = "inst/job_data/perf_15162ff/"

all_trait = readRDS(paste0(main_folder, "/bigmem_trait_df.Rds"))

tictoc::tic()
list.files(main_folder, "perf_df_*", full.names = TRUE) %>%
    future_walk(function(x) {

        pkgload::load_all()

        x %>%
            readRDS() %>%
            compute_weighted_performance(all_trait) %>%
            saveRDS(gsub("perf_df_", "cwm_df_", x), compress = TRUE)
    },
    .progress = TRUE)
tictoc::toc()

all_sp_perf_df = list.files(main_folder, "perf_df_*", full.names = TRUE) %>%
    purrr::map_dfr(readRDS)

all_perf_df = list.files(main_folder, "cwm_df_*", full.names = TRUE) %>%
    purrr::map_dfr(readRDS)

tidy_perf = all_perf_df %>%
    filter(trait_cor == "uncor") %>%
    tidyr::gather("comperf_name", "comperf_value", matches("trait[12]")) %>%
    tidyr::separate("comperf_name", c("comperf_name", "trait"), sep = "_") %>%
    filter(trait == "trait2")

# Subset data ------------------------------------------------------------------

all_perf_df = all_perf_df %>%
    inner_join(all_trait %>%
                   mutate(seed = as.integer(seed)),
               by = c("species", "seed", "trait_cor"))

# Single Parameter Set CWMs
single_cwm = tidy_perf %>%
    filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B > 2.5e-7, B < 2.6e-7, H == 1e-4,
           trait_cor == "uncor", patch >= 5, patch <= 20,
           comperf_name == "cwm")

seed_df = all_perf_df %>%
    distinct(seed) %>%
    slice(1:10)

# Compare Avg. Growth Rate wit various param -----------------------------------
# Compare values of average growth rate per patch when param. equals 0 or ≠ than
# 0
target_param = rlang::sym("H")

perf_growth = all_sp_perf_df %>%
    filter(R_scenar == 0, A_scenar == 0, H_scenar == 0) %>%
    group_by(B, seed, trait_cor, patch, species) %>%
    select(groups(), max_growth_rate, !!target_param) %>%
    ungroup() %>%
    filter(trait_cor == "uncor") %>%
    mutate(!!target_param := as.numeric(!!target_param)) %>%
    arrange(B, seed, patch, species, !!target_param, max_growth_rate) %>%
    group_by(B, seed, trait_cor, patch, species) %>%
    mutate(percent_growth = max_growth_rate/first(max_growth_rate)) %>%
    group_by(B, trait_cor, !!target_param) %>%
    summarise(rel_growth = mean(percent_growth, na.rm = TRUE))


# Figure 2: Abundance environment curve ----------------------------------------
# Effect of simulation parameters on single species abundance along the
# environment.
# (Abundance vs. Patch in function of values of parameters for a fixed trait
# contribution scenario)
# Four panels figure: color code for parameter values k, A, B, H
subset_perf = all_sp_perf_df %>%
    filter(trait_cor == "uncor", species == "species67",
           R_scenar == 0, A_scenar == 0, H_scenar == 0)

fig_env_abund_k = plot_env_abund(subset_perf, "k", "Max. Growth Rate")
fig_env_abund_A = plot_env_abund(subset_perf, "A", "Limiting. Sim.\nIntensity")
fig_env_abund_B = plot_env_abund(subset_perf, "B", "Intrasp. Compet.\nIntensity")
fig_env_abund_H = plot_env_abund(subset_perf, "H", "Hierach. Compet.\nIntensity")

fig_env_abund = plot_grid(fig_env_abund_k, fig_env_abund_A,
                          fig_env_abund_B, fig_env_abund_H,
                          ncol = 2, labels = c("k", "A", "B", "H"))

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
# H_scenar color, for a fixed simulation parameter set)
single_cwm %>%
    filter(patch >= 5, patch <= 20, trait == "trait2") %>%
    ggplot(aes(patch, comperf_value, color = as.factor(H_scenar),
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
all_perf_df %>%
    filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B >  2.5e-7, B <  2.6e-7,
           H > 3.16e-5, H < 3.17e-5, trait_cor == "uncor", patch >= 5,
           patch <= 20) %>%
    ggplot(aes(patch, cwm_trait2, color = as.factor(H_scenar),
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

# Figure 4: Effect of Param. Values on Estimators of TxE relationship -----------
# Difference in performance indices affected by the intensity of ecological
# processes

tidy_perf %>%
    filter(R_scenar == 100, A_scenar == 100, H_scenar == 100, k == 1.3,
           !(comperf_name %in% c("monocult", "polycult"))) %>%
    ggplot(aes(patch, comperf_value, color = comperf_name)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 1/50) +
    stat_smooth(se = FALSE, geom = "line", size = 1,
                alpha = 1/10, aes(group = interaction(seed, comperf_name))) +
    stat_smooth(se = FALSE, geom = "line", size = 1.3,
                alpha = 1/2) +
    facet_grid(vars(A, B), vars(H),
               labeller = labeller(H = function(x) scientific_notation(x, "H"),
                                   A = function(x) scientific_notation(x, "A"),
                                   B = function(x) scientific_notation(x, "B"))) +
    guides()
    labs(x = "Environment",
         y = "Performance Value") +
    theme_bw() +
    theme(aspect.ratio = 1)
# Figure 5: Higher-order moments CWV & CWS -------------------------------------

# Figure 6: Single Patch Trait-Perforamance relationship -----------------------
fig_trait_perf_A = all_sp_perf_df %>%
    filter(R_scenar == 100, A_scenar == 100, H_scenar == 100, k == 1.3,
           B > 2e-7, H == 1e-4, seed == 227470, patch == 5) %>%
    inner_join(all_trait %>%
                   mutate(seed = as.numeric(seed))) %>%
    mutate(N150 = log10(N150 + 1)) %>%
    tidyr::gather("perf_index", "perf_value", N150, max_growth_rate,
                  env_growth_rate) %>%
    ggplot(aes(trait2, perf_value, color = as.factor(A))) +
    geom_point() +
    geom_line() +
    facet_wrap(~perf_index, scales = "free_y") +
    scale_color_viridis_d() +
    labs(x = "Trait Value",
         y = "Performance Value")

fig_trait_perf_B = all_sp_perf_df %>%
    filter(R_scenar == 100, A_scenar == 100, H_scenar == 100, k == 1.3,
           A > 2e-7, H == 1e-4, seed == 227470, patch == 5) %>%
    inner_join(all_trait %>%
                   mutate(seed = as.numeric(seed))) %>%
    mutate(N150 = log10(N150 + 1)) %>%
    tidyr::gather("perf_index", "perf_value", N150, max_growth_rate,
                  env_growth_rate) %>%
    ggplot(aes(trait2, perf_value, color = as.factor(B))) +
    geom_point() +
    geom_line() +
    facet_wrap(~perf_index, scales = "free_y") +
    scale_color_viridis_d() +
    labs(x = "Trait Value",
         y = "Performance Value")

fig_trait_perf_H = all_sp_perf_df %>%
    filter(R_scenar == 100, A_scenar == 100, H_scenar == 100, k == 1.3,
           A > 2e-7, B > 2e-7, seed == 227470, patch == 5) %>%
    inner_join(all_trait %>%
                   mutate(seed = as.numeric(seed))) %>%
    mutate(N150 = log10(N150 + 1)) %>%
    tidyr::gather("perf_index", "perf_value", N150, max_growth_rate,
                  env_growth_rate) %>%
    ggplot(aes(trait2, perf_value, color = as.factor(H))) +
    geom_point() +
    geom_line() +
    facet_wrap(~perf_index, scales = "free_y") +
    scale_color_viridis_d() +
    labs(x = "Trait Value",
         y = "Performance Value")

plot_grid(fig_trait_perf_A, fig_trait_perf_B, fig_trait_perf_H,
          nrow = 3)

# Other Figures ----------------------------------------------------------------


