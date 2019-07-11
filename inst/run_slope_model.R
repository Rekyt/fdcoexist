library("dplyr")

all_slopes = readRDS("inst/job_data/perf_2fd398/all_slopes.Rds")

filtered_slopes = all_slopes %>%
    mutate(seed = as.factor(seed),
           trait_cor = as.factor(trait_cor)) %>%
    mutate_at(vars(k, H, R_scenar, A_scenar, H_scenar),
              .funs = list(scale = ~as.numeric(scale(.)))) %>%
    mutate_at(vars(A, B, H),
              .funs = list(log = ~as.numeric(scale(log(. + 1e-16)))))
tictoc::tic()
full_model = lme4::lmer(estimate ~ k_scale + A_log + B_log + H_scale +
                            R_scenar_scale + A_scenar_scale + H_scenar_scale +
                            trait_cor + (1 | seed),
                        data = filtered_slopes)
tictoc::toc()
saveRDS(full_model, "inst/job_data/perf_2fd398/mixed_model.Rds")

fig_mod_estimates = sjPlot::plot_model(full_model,
                   terms = c("k_scale", "A_log", "B_log", "H_scale",
                             "R_scenar_scale", "A_scenar_scale",
                             "H_scenar_scale", "trait_cor"),
                   show.values = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, colour = "#001c00", size = 1,
               alpha = 1/2) +
    scale_x_discrete(
        labels = c(k_scale        = "Max. environmental\n growth (k)",
                   A_log          = "Interspecific\ncompet. (A)",
                   B_log          = "Intraspepcific\ncompet. (B)",
                   H_scale        = "Hierarchical\ncompet. (H)",
                   R_scenar_scale = "Trait contrib. to\ngrowth\n(R_scenar)",
                   A_scenar_scale = "Trait contrib. to\nintersp. compet.\n(A_scenar)",
                   H_scenar_scale = "Trait contrib. to\nhierarch. compet.\n(H_scenar)"),
        limits = c("trait_cor",
                   "H_scenar_scale",
                   "A_scenar_scale",
                   "R_scenar_scale",
                   "H_scale",
                   "B_log",
                   "A_log",
                   "k_scale")) +
    labs(title = "Simulation Parameters Effect on CWM <-> Environment Slope") +
    theme_bw(12) +
    theme(aspect.ratio = 1,
          axis.text = element_text(colour = "#001c00"))

fig_mod_estimates
