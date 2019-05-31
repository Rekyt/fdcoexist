list.files("inst/job_data/perfs_52920/cwm/", full.names = TRUE) %>%
    purrr::map_dfr(readRDS) %>%
    nest() %>%
    mutate(mod_cwm1 = map(data, ~lm(trait1_cwm ~ patch, data = .x)),
           mod_cwm2 = map(data, ~lm(trait2_cwm ~ patch, data = .x)),
           tidy_cwm1 = map(mod_cwm1, broom::tidy),
           tidy_cwm2 = map(mod_cwm2, broom::tidy)) -> ji

weights = c(0, 50, 100)

all_param_df = crossing(h_fun = "+", di_thresh = 24, k = list_k, A = list_A,
                        B = list_B, H = list_H,
                        trait_cor = c("negcor", "uncor", "poscor"),
                        R_scenar = weights, A_scenar = weights,
                        H_scenar = weights, seed = 1:30)

all_cwm = ji %>%
    unnest(tidy_cwm1) %>%
    filter(term == "patch") %>%
    full_join(all_param_df)

median_scenario = all_cwm %>%
    filter(trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50)

plot_param_space = function(provided_df, x_var, y_var) {

    x_name = paste0("as.factor(", x_var, ")")
    y_name = paste0("as.factor(", y_var, ")")

    ggplot(provided_df, aes_string(x_name, y_name, z = "estimate")) +
        stat_summary_2d(fun = mean, geom = "tile", na.rm = TRUE,
                        drop = FALSE) +
        scale_fill_viridis_c(limits = c(0, 1)) +
        theme_bw() +
        theme(aspect.ratio = 1,
              legend.position = "top")
}

fig_k_A = plot_param_space(median_scenario, "k", "A")

fig_k_B = plot_param_space(median_scenario, "k", "B")

fig_k_H = plot_param_space(median_scenario, "k", "H")

fig_A_B = plot_param_space(median_scenario, "A", "B")

fig_A_H = plot_param_space(median_scenario, "A", "H")

fig_B_H = plot_param_space(median_scenario, "B", "H")

estimate_legend = cowplot::get_legend(fig_k_A)

cowplot::plot_grid(
    estimate_legend,
    list(fig_k_A, fig_A_B, fig_k_B, fig_A_H, fig_k_B, fig_B_H) %>%
        lapply(function(x) x +
                   theme(legend.position = "none")) %>%
        cowplot::plot_grid(plotlist = ., nrow = 3),
    rel_heights = c(0.05, 0.95), nrow = 2, ncol = 1)
