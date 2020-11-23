# Relationship between trait difference and effect of exponents
# in limiting similarity and hierarchical competition
library("dplyr")
library("ggplot2")

# Exponent of hierarchical compoent --------------------------------------------
hierar_df = data.frame(
    diff_trait = seq(-2, 2, length.out = 5000)
) %>%
    mutate(hierar_0.5  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^0.5),
           hierar_1.0  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^1),
           hierar_2.0  = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^2),
           hierar_0.25 = ifelse(diff_trait >= 0, 0, -abs(diff_trait)^0.25)) %>%
    tidyr::pivot_longer(hierar_0.5:hierar_0.25, names_to = "hierar_exp",
                        values_to = "hierar_effect") %>%
    tidyr::extract(hierar_exp, c("hierar_exp"), regex = ".*_(.*)$")

plot_hierar_effect = ggplot(
    hierar_df, aes(diff_trait, hierar_effect, color = hierar_exp)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_line(size = 1) +
    labs(x = "Difference in trait\n(focal species - other species)",
         y = "Addition to basal growth",
         color = "Hier. Comp.\nexponent",
         title = "Hierarchical Competition effect") +
    theme_bw() +
    theme(aspect.ratio = 1)

# Exponent of limiting similarity ----------------------------------------------

limsim_df = data.frame(
    diff_trait = seq(-2, 2, length.out = 5000)
) %>%
    mutate(limsim_0.5  = abs(diff_trait)^0.5,
           limsim_1.0  = abs(diff_trait)^1,
           limsim_2.0  = abs(diff_trait)^2,
           limsim_0.25 = abs(diff_trait)^0.25) %>%
    tidyr::pivot_longer(limsim_0.5:limsim_0.25, names_to = "limsim_exp",
                        values_to = "limsim_effect") %>%
    tidyr::extract(limsim_exp, c("limsim_exp"), regex = ".*_(.*)$")

plot_limsim_effect = ggplot(
    limsim_df, aes(diff_trait, limsim_effect, color = limsim_exp)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_line(size = 1) +
    labs(x = "Difference in trait\n(focal species - other species)",
         y = "Divider to basal growth",
         color = "Lim. sim.\nexponent",
         title = "Limiting similarity effect") +
    theme_bw() +
    theme(aspect.ratio = 1)

cowplot::plot_grid(plot_hierar_effect, plot_limsim_effect, ncol = 1, align = "hv")
