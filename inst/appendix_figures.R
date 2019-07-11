# Script to assemble all performances from all fdcoexist simulations and output
# final figures
# Author: Matthias Grenié <matthias dot grenie at gmail dot com>
# Packages ---------------------------------------------------------------------
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("cowplot")

# Functions --------------------------------------------------------------------
plot_param_space = function(provided_df, x_var, y_var,
                            z_var = "estimate") {

    x_name = paste0("as.factor(", x_var, ")")
    y_name = paste0("as.factor(", y_var, ")")

    ggplot(provided_df, aes_string(x_name, y_name, z = z_var)) +
        stat_summary_2d(fun = "mean", geom = "tile",
                        fun.args = list(na.rm = FALSE),
                        drop = TRUE) +
        scale_fill_viridis_c(na.value = "gray65") +
        scale_x_discrete(labels = function(x) {
            x %>%
                as.character() %>%
                as.numeric() %>%
                format(digits = 2)
        }) +
        scale_y_discrete(labels = function(x) {
            x %>%
                as.character() %>%
                as.numeric() %>%
            format(digits = 2)
        }) +
        theme_bw() +
        theme(aspect.ratio = 1,
              legend.position = "top",
              axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_sr_param_space = function(provided_df, x_var, y_var) {

    plot_param_space(provided_df, x_var, y_var, z_var = "species_rich") +
        scale_fill_viridis_c(limits = c(0, 28))
}

# Parameter space --------------------------------------------------------------
weights = c(0, 50, 100)
list_A = c(0, 10^-(seq(1, 8, length.out = 6)))
list_k = seq(1, 1.5, length.out = 6)
list_B = list_A
list_H = list_A

set.seed(20190619)
seed_list = sample(1e6, size = 30)

# Data Frame of combination of all parameters
all_param_df = expand.grid(seed = seed_list,
                               A = list_A, k = list_k,
                               B = list_B, H = list_H) %>%
    mutate(file_number = row_number())

# Load data --------------------------------------------------------------------

all_cwm = readRDS("inst/job_data/perf_2fd398/all_cwm.Rds")

all_slopes = readRDS("inst/job_data/perf_2fd398/all_slopes.Rds")

# Extract scenario -------------------------------------------------------------

# Get all combination of scenarios
all_comb = all_cwm %>%
    select(k, A, B, H, h_fun, di_thresh, R_scenar, A_scenar, H_scenar,
           trait_cor, seed) %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50) %>%
    distinct() %>%
    full_join(crossing(trait_cor = c("negcor", "uncor", "poscor"),
                       cwm_name  = c("trait1_cwm", "trait2_cwm")),
              by = "trait_cor")

# Select median scenario and add missing combination
median_scenario = all_slopes %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50) %>%
    # Add missing combinations
    full_join(all_comb %>%
                  filter(trait_cor == "uncor"),
              by = colnames(all_comb)) %>%
    filter(cwm_name == "trait2_cwm")

# Compute average slope per parameter combination
median_scenario_mean_slope = median_scenario %>%
    group_by(k, A, B, H) %>%
    summarise(mean_slope = mean(estimate, na.rm = TRUE))


# Supplementary Figure 1: Slope parameter space --------------------------------

fig_k_A = plot_param_space(median_scenario, "k", "A")

fig_k_B = plot_param_space(median_scenario, "k", "B")

fig_k_H = plot_param_space(median_scenario, "k", "H")

fig_A_B = plot_param_space(median_scenario, "A", "B")

fig_A_H = plot_param_space(median_scenario, "A", "H")

fig_B_H = plot_param_space(median_scenario, "B", "H")

estimate_legend = cowplot::get_legend(fig_k_A)

param_space = cowplot::plot_grid(
    estimate_legend,
    list(fig_k_A, fig_A_B, fig_k_B, fig_A_H, fig_k_H, fig_B_H) %>%
        lapply(function(x) x +
                   theme(legend.position = "none")) %>%
        cowplot::plot_grid(plotlist = ., nrow = 3, align = "hv"),
    rel_heights = c(0.05, 0.95), nrow = 2, ncol = 1)

ggsave("fig_param_space.png", plot = param_space, width = 14, height = 21,
       units = "cm")


# Supplementary Figure 2: SR parameter space -----------------------------------

all_sr = all_cwm %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50, patch == 13)

fig_sr_k_A = plot_sr_param_space(all_sr, "k", "A")

fig_sr_k_B = plot_sr_param_space(all_sr, "k", "B")

fig_sr_k_H = plot_sr_param_space(all_sr, "k", "H")

fig_sr_A_B = plot_sr_param_space(all_sr, "A", "B")

fig_sr_A_H = plot_sr_param_space(all_sr, "A", "H")

fig_sr_B_H = plot_sr_param_space(all_sr, "B", "H")

sr_legend = cowplot::get_legend(fig_sr_k_A)

param_sr_space = cowplot::plot_grid(
    sr_legend,
    list(fig_sr_k_A, fig_sr_A_B, fig_sr_k_B, fig_sr_A_H, fig_sr_k_H,
         fig_sr_B_H) %>%
        lapply(function(x) x +
                   theme(legend.position = "none")) %>%
        cowplot::plot_grid(plotlist = ., nrow = 3, align = "hv"),
    rel_heights = c(0.05, 0.95), nrow = 2, ncol = 1)

ggsave("fig_param_space_sr.png", plot = param_sr_space, width = 14, height = 21,
       units = "cm")
