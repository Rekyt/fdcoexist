# FDcoexist October 16 2020
# Authors: P. Denelle
# Generate correlation tables between performance metrics
# they will be reproduced on Inskcape and directly annotated on Figure 2
# temporary

# Packages ----
library(ggplot2)
library(tidyr)

## Loading figure 2 ----
fig2 <- readRDS("inst/figures/paper_figure2.Rds")

fig2_p <- fig2$data
# Removing environmental filtering only
fig2_p <- fig2_p[which(fig2_p$comb != "2_0_0"), ]

# Correlation tables ----
# Limiting similarity
lim_sim <- spread(fig2_p[which(fig2_p$comb == "2_0.001_0"), ],
                  mismatch_name, mismatch_value)

mat_cor <- cor(lim_sim[, c("MisAbPatchP", "MisavgGRPatchP", "MisinstRPatchP",
                           "MismaxGRPatchP")], method = "spearman")
mat_cor <- round(mat_cor, 2)
diag(mat_cor) <- ""
mat_cor[upper.tri(mat_cor)] <- ""
rownames(mat_cor) <- colnames(mat_cor) <- c("ab", "avg", "int", "max")

# Hierarchical competition
hie <- spread(fig2_p[which(fig2_p$comb == "2_0_0.001"), ],
              mismatch_name, mismatch_value)

mat_cor <- cor(hie[, c("MisAbPatchP", "MisavgGRPatchP", "MisinstRPatchP",
                       "MismaxGRPatchP")], method = "spearman")
mat_cor <- round(mat_cor, 2)
diag(mat_cor) <- ""
mat_cor[upper.tri(mat_cor)] <- ""
rownames(mat_cor) <- colnames(mat_cor) <- c("ab", "avg", "int", "max")

# Figure 2 ----
ggplot(fig2_p, aes(mismatch_value, Species, color = mismatch_name,
                   shape = mismatch_name)) +
    geom_segment(aes(x = min_mismatch, xend = max_mismatch,
                     yend = Species),
                 color = "grey", size = 2/3) +
    geom_vline(xintercept = 0, linetype = 1, size = 1/2, color = "black") +
    geom_point(size = 1.5) +
    facet_wrap(vars(comb), ncol = 3,
               labeller = labeller(
                   comb = c("2_0.001_0" = "+Limiting Similarity",
                            "2_0_0.001" = "+Hierarchical Competition",
                            "2_0_0" = "Environmental filtering\nonly"))) +
    scale_x_continuous(labels = scales::label_percent()) +
    scale_y_continuous(limits = c(1, 50), breaks = c(1, c(1,2,3,4,5)*10)) +
    scale_color_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate")) +
    scale_shape_discrete(labels = c(
        MisAbPatchP = "Abundance",
        MisavgGRPatchP = "Average Growth Rate",
        MisinstRPatchP = "Intrinsic Growth Rate",
        MismaxGRPatchP = "Maximum Growth Rate")) +
    labs(x = "Relative Mismatch from True Optimal Patch (% of gradient)",
         shape = "Performance\nMeasure",
         color = "Performance\nMeasure") +
    theme_bw(10) +
    theme(aspect.ratio = 1,
          panel.grid.major.x = element_line(size = 1.3),
          panel.spacing.x = unit(4, "mm"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "top")
