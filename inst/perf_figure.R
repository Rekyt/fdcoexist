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

main_folder = "inst/job_data/perf_c4a9018/"

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

all_perf_df = list.files(main_folder, "^cwm_df_*", full.names = TRUE) %>%
    purrr::map_dfr(readRDS)




tidy_perf = all_perf_df %>%
    filter(trait_cor == "uncor") %>%
    tidyr::gather("comperf_name", "comperf_value", matches("trait[12]")) %>%
    tidyr::separate("comperf_name", c("trait", "comperf_name"), sep = "_",
                    extra = "merge") %>%
    filter(trait == "trait2")

saveRDS(tidy_perf, paste0(main_folder, "tidy_perf_c4a9018.Rds"),
        compress = TRUE)
# Extract performances in first generations ------------------------------------
trait_list = split(all_trait, list(all_trait$seed)) %>%
    lapply(function(z) {
        z %>%
            split(list(z$trait_cor)) %>%
            lapply(function(x) {
                x %>%
                    select(-seed, -trait_cor) %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames("species") %>%
                    as.matrix()
            })
})

tictoc::tic()
list.files(main_folder, "simul_cat_*", full.names = TRUE) %>%
    future_walk(function(x) {

        pkgload::load_all()

        perf_df = x %>%
            readRDS() %>%
            purrr::map_dfr(function(y) {

                extract_performances_from_simul(
                    y, trait_list[[as.character(y$seed)]],  chosen_time = 4
                )})

        saveRDS(perf_df, gsub("simul_cat_", "t10_perf_df_", x),
                compress = TRUE)

        perf_df %>%
            compute_weighted_performance(all_trait) %>%
            saveRDS(gsub("simul_cat_", "t10_cwm_df_", x), compress = TRUE)
    }, .progress = TRUE)
tictoc::toc()

all_perf_t4_df = list.files(main_folder, "^t10_cwm_df_*", full.names = TRUE) %>%
    purrr::map_dfr(readRDS)

tidy_perf_t4 = all_perf_t4_df %>%
    filter(trait_cor == "uncor") %>%
    tidyr::gather("comperf_name", "comperf_value", matches("trait[12]")) %>%
    tidyr::separate("comperf_name", c("trait", "comperf_name"), sep = "_",
                    extra = "merge") %>%
    filter(trait == "trait2")

saveRDS(tidy_perf_t4, paste0(main_folder, "tidy_perf_t4_c4a9018.Rds"),
        compress = TRUE)

# Combine Estimates at all times -----------------------------------------------

tidy_perf_all = bind_rows(tidy_perf, tidy_perf_t4)

# Compare Avg. Growth Rate wit various param -----------------------------------
# Compare values of average growth rate per patch when param. equals 0 or ≠ than
# 0
target_param = rlang::sym("A")

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

# Reproducing preliminary figures ----------------------------------------------

perf_estimate = c(
    cwm                 = "CWM",
    cwv                 = "CWV",
    cws                 = "CWS",
    cwk                 = "CWK",
    pure_env            = "Weighted Th. GR",
    weighted_avg_growth = "Weighted Avg. GR",
    weighted_int_growth = "Weighted. Int. GR (exp. model)",
    weighted_max_growth = "Weighted Max.GR"
)

A_for_k_1.15 = c(0, 2.19e-4, 2.257e-4,  2.44e-4, 2.58e-4)
A_for_k_1.3  = c(0,   1e-7,     5e-7,   1.9e-6, 6.92e-6)
A_for_k_1.45 = c(0, 2.32e-9,  2.36e-8, 2.155e-7, 1.75e-6)
list_k = c(1.15, 1.3, 1.45)
list_B = c(0, 1.585e-4, 3.17e-4)
H_for_k_1.15 = c(0, 6.215e-4, 6.48e-4, 6.68e-4, 6.754e-4)
H_for_k_1.3  = c(0,     1e-6,    4e-6, 1.35e-5, 4.05e-5)
H_for_k_1.45 = c(0,  2.31e-8, 2.31e-7, 1.1e-6,  3e-6)

a_values = data.frame(
    k = rep(c(1.15, 1.3, 1.45), each = 5),
    A = c(A_for_k_1.15, A_for_k_1.3, A_for_k_1.45),
    A_red = rep(c(0, 0.2, 0.4, 0.6, 0.8), 3)
)

h_values = data.frame(
    k = rep(c(1.15, 1.3, 1.45), each =  5),
    H = c(H_for_k_1.15, H_for_k_1.3, H_for_k_1.45),
    H_red = rep(c(0, 0.2, 0.4, 0.6, 0.8), 3)
)

base_scenario = tidy_perf_all %>%
    filter(R_scenar == 100, A_scenar == 100, H_scenar == 100,
           trait_cor == "uncor") %>%
    mutate(perf_minus_expect = comperf_value - patch,
           pred_trunc_gaussian = purrr::map_dbl(
               patch,
               ~truncated_gaussian(.x, 2, 1, 25)),
           perf_minus_trunc_gaussian = comperf_value - pred_trunc_gaussian)

# Base figure: No competition --------------------------------------------------
# Figure in the absence of any competition
fig_no_competition = base_scenario %>%
    filter(A == 0, B == 0, H == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(time), labeller = label_both) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "No competition",
         caption = "A = 0; B = 0; H = 0; uncorrelated traits; N = 15; CI±95%")

# Figure Intrasp. competition ---------------------------------------------------
# Figure showing the influence of intra-specific competition only
fig_intra_competition = base_scenario %>%
    filter(A == 0, H == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    ggplot(aes(patch, perf_minus_expect, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(B, time), labeller = label_both) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Intra-specific competition only",
         caption = "A = 0; H = 0; uncorrelated traits; N = 15 CI ±95%")


# Figure Limiting Similarity ---------------------------------------------
# Figure showing the influence of limiting similarity only
fig_lim_sim = base_scenario %>%
    filter(B == 0, H == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(a_values, by = c("k", "A")) %>%
    ggplot(aes(patch, perf_minus_expect, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(A_red, time), labeller = labeller(
        k     = label_both,
        A_red = function(x) paste0("Growth Red.: ",
                                   scales::percent(as.numeric(x), accuracy = 1))
    )) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Limiting similarity only",
         caption = "B = 0; H = 0; uncorrelated traits; N = 15; CI ±95%")

# Figure Hierarchical competition ----------------------------------------------
# Figure showing the effect of hierarchical competition
fig_hierarch_comp = base_scenario %>%
    filter(B == 0, A == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(h_values, by = c("k", "H")) %>%
    ggplot(aes(patch, perf_minus_expect, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(H_red, time), labeller = labeller(
        k     = label_both,
        H_red = function(x) paste0("Growth Red.: ",
                                   scales::percent(as.numeric(x), accuracy = 1))
    )) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Hierarchical Competition only",
         caption = "B = 0; A = 0; uncorrelated traits; N = 15; CI ±95%")

# Figure All competitions ------------------------------------------------------
fig_both_comp = base_scenario %>%
    filter(B == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(h_values, by = c("k", "H")) %>%
    inner_join(a_values, by = c("k", "A")) %>%
    ggplot(aes(patch, perf_minus_expect, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(A_red, time), vars(k, H_red), labeller = labeller(
        k     = label_both,
        H_red = function(x) paste0("Growth Red. (H): ",
                                   scales::percent(as.numeric(x), acc = 1)),
        A_red = function(x) paste0("Growth Red. (A): ",
                                   scales::percent(as.numeric(x), acc = 1))
    ), scales = "free_y") +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1,
          legend.position = "top") +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Hierarchical Competition only",
         caption = "B = 0; uncorrelated traits; N = 15; CI ±95%")

# Figures Multi-traits ---------------------------------------------------------

multi_base = tidy_perf_all %>%
    filter(R_scenar == 50, A_scenar == 100, H_scenar == 100,
           trait_cor == "uncor") %>%
    mutate(perf_minus_expect = comperf_value - patch,
           pred_trunc_gaussian = purrr::map_dbl(
               patch,
               ~truncated_gaussian(.x, 2, 1, 25)),
           perf_minus_trunc_gaussian = comperf_value - pred_trunc_gaussian)

fig_multi_no_competition = multi_base %>%
    filter(A == 0, B == 0, H == 0) %>%
    filter(!grepl("cw[vsk]|best", comperf_name)) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(time), labeller = label_both) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Multi-traits (R50%) – No competition",
         caption = "A = 0; B = 0; H = 0; uncorrelated traits; N = 15; CI ±95%")

fig_multi_intra_competition = multi_base %>%
    filter(A == 0, H == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(B, time), labeller = label_both) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Multi-traits (R50%) – Intra-specific competition only",
         caption = "A = 0; H = 0; uncorrelated traits; N = 15 CI ±95%")

fig_multi_lim_sim = multi_base %>%
    filter(B == 0, H == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(a_values, by = c("k", "A")) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(A_red, time), labeller = labeller(
        k     = label_both,
        A_red = function(x) paste0("Growth Red.: ",
                                   scales::percent(as.numeric(x), accuracy = 1))
    )) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Multi-traits (R50%) – Limiting similarity only",
         caption = "B = 0; H = 0; uncorrelated traits; N = 15; CI ±95%")

fig_multi_hierarch_comp = multi_base %>%
    filter(B == 0, A == 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(h_values, by = c("k", "H")) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(k), vars(H_red, time), labeller = labeller(
        k     = label_both,
        H_red = function(x) paste0("Growth Red.: ",
                                   scales::percent(as.numeric(x), accuracy = 1))
    )) +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1) +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Multi-traits (R50%) – Hierarchical Competition only",
         caption = "B = 0; A = 0; uncorrelated traits; N = 15; CI ±95%")

fig_multi_both_comp = multi_base %>%
    filter(B == 0, A != 0, H != 0, !grepl("cw[vsk]|best", comperf_name)) %>%
    inner_join(h_values, by = c("k", "H")) %>%
    inner_join(a_values, by = c("k", "A")) %>%
    ggplot(aes(patch, perf_minus_trunc_gaussian, color = comperf_name)) +
    geom_hline(yintercept = 0, linetype = 2) +
    stat_summary(fun.y = mean, geom = "line") +
    stat_summary(fun.data = mean_cl_boot, alpha = 1/5) +
    facet_grid(vars(A_red, time), vars(k, H_red), labeller = labeller(
        k     = label_both,
        H_red = function(x) paste0("Growth Red. (H): ",
                                   scales::percent(as.numeric(x), acc = 1)),
        A_red = function(x) paste0("Growth Red. (A): ",
                                   scales::percent(as.numeric(x), acc = 1))
    ), scales = "free_y") +
    scale_color_discrete(labels = perf_estimate) +
    theme(aspect.ratio = 1,
          legend.position = "top") +
    labs(x = "Environment",
         y  = "Deviation from expectation",
         color = "Performance Estimates",
         title = "Multi-traits (R50%) – Hierarchical Competition only",
         caption = "B = 0; uncorrelated traits; N = 15; CI ±95%")
