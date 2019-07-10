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
        stat_summary_2d(fun = "mean", geom = "tile", na.rm = TRUE,
                        drop = TRUE) +
        scale_fill_viridis_c() +
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
# Load data --------------------------------------------------------------------

all_cwm = readRDS("inst/job_data/perf_2fd398/all_cwm.Rds")

all_slopes = readRDS("inst/job_data/perf_2fd398/all_slopes.Rds")

subset_cwm = all_cwm %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50, patch >= 5, patch <= 20)
# Isolate problematic releves --------------------------------------------------
# Some combinations of parameters leads to more than 25 values (across the range
# of patches) It may be because of anomaly during Assemblage of CWMs
all_cwm %>%
    count() %>%
    filter(n > 25)

# Parameter space --------------------------------------------------------------
weights = c(0, 50, 100)
list_A = c(0, 10^-(seq(1, 8, length.out = 6)))
list_k = seq(1, 1.5, length.out = 6)
list_B = list_A
list_H = seq(0, 0.2, length.out = 6)
all_param_df = crossing(h_fun = "+", di_thresh = 24, k = list_k, A = list_A,
                        B = list_B, H = list_H,
                        trait_cor = c("negcor", "uncor", "poscor"),
                        R_scenar = weights, A_scenar = weights,
                        H_scenar = weights, seed = 1:30)

# Compute slopes ---------------------------------------------------------------
median_scenario = all_slopes %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50)

subset_data = subset_cwm %>%
    tidyr::nest(patch, trait1_cwm, trait2_cwm, trait1_cwv, trait2_cwv,
                species_rich) %>%
    mutate(all_na = purrr::map_lgl(data, ~all(is.na(.$trait1_cwm))))

subset_slopes = subset_data %>%
    filter(!all_na) %>%
    mutate(lm_mod = purrr::map(data, ~lm(trait1_cwm ~ patch,
                                         data = .x, na.action = na.exclude)),
           lm_sum = purrr::map(lm_mod, broom::tidy)) %>%
    unnest(lm_sum) %>%
    filter(term == "patch") %>%
    full_join(subset_data)

# Plot parameter space ---------------------------------------------------------

fig_k_A = plot_param_space(subset_slopes, "k", "A")

fig_k_B = plot_param_space(subset_slopes, "k", "B")

fig_k_H = plot_param_space(subset_slopes, "k", "H")

fig_A_B = plot_param_space(subset_slopes, "A", "B")

fig_A_H = plot_param_space(subset_slopes, "A", "H")

fig_B_H = plot_param_space(subset_slopes, "B", "H")

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

# Figure CWM-Environment -------------------------------------------------------

fig_all_cwm = subset_cwm %>%
    ungroup() %>%
    mutate(A_H = paste("A = ", format(A, digits = 2, scientific = TRUE),
                       "; H = ", format(H, digits = 2, scientific = TRUE),
                       sep = "")) %>%
    ggplot(aes(patch, trait1_cwm, color = A_H,
               group = interaction(A_H, seed))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, size = 0.8) +
    geom_point(alpha = 1/5) +
    stat_smooth(geom = "line", size = 1, alpha = 1/3) +
    facet_grid(vars(k), vars(B), labeller = labeller(B = function(x) {
        x %>%
            as.numeric() %>%
            format(digits = 2, scientific = TRUE) %>%
            paste0("B = ", .)
    }, k = label_both)) +
    scale_color_discrete() +
    theme_bw() +
    theme(aspect.ratio = 1)

ggsave("fig_all_cwm.png", fig_all_cwm, width = 29.7, height = 21,
       units = "cm")

# Make figures using species richness ------------------------------------------

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
