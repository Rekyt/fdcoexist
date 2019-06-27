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
plot_param_space = function(provided_df, x_var, y_var) {

    x_name = paste0("as.factor(", x_var, ")")
    y_name = paste0("as.factor(", y_var, ")")

    ggplot(provided_df, aes_string(x_name, y_name, z = "estimate")) +
        stat_summary_2d(fun = mean, geom = "tile", na.rm = TRUE,
                        drop = FALSE) +
        scale_fill_viridis_c(limits = c(0, 1)) +
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

# Load data --------------------------------------------------------------------

all_cwm = readRDS("inst/job_data/older_hierarch/all_cwm.Rds")

all_slopes = readRDS("inst/job_data/older_hierarch/all_slopes.Rds")


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
    filter(cwm_name == "trait1_cwm") %>%
    ungroup() %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50)

# Plot parameter space ---------------------------------------------------------

fig_k_A = plot_param_space(median_scenario, "k", "A")

fig_k_B = plot_param_space(median_scenario, "k", "B")

fig_k_H = plot_param_space(median_scenario, "k", "H")

fig_A_B = plot_param_space(median_scenario, "A", "B")

fig_A_H = plot_param_space(median_scenario, "A", "H")

fig_B_H = plot_param_space(median_scenario, "B", "H")

estimate_legend = cowplot::get_legend(fig_k_A)

param_space = cowplot::plot_grid(
    estimate_legend,
    list(fig_k_A, fig_A_B, fig_k_B, fig_A_H, fig_k_B, fig_B_H) %>%
        lapply(function(x) x +
                   theme(legend.position = "none")) %>%
        cowplot::plot_grid(plotlist = ., nrow = 3, align = "hv"),
    rel_heights = c(0.05, 0.95), nrow = 2, ncol = 1)

ggsave("fig_param_space.png", plot = param_space, width = 14, height = 21,
       units = "cm")

# Figure CWM-Environment -------------------------------------------------------

subset_cwm = all_cwm %>%
    filter(trait_cor == "uncor", R_scenar == 0, A_scenar == 0,
           H_scenar == 0)

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
