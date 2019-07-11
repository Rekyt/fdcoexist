# Script make performance comparison figure
# Packages ---------------------------------------------------------------------
library("dplyr")

# Load Data --------------------------------------------------------------------
all_perf = list.files("inst/job_data/perf_2fd398/perf_data",  "perf_df_411*",
                      full.names = TRUE)

all_perf_df = purrr::map_dfr(all_perf, readRDS)

all_trait = readRDS("inst/job_data/perf_2fd398/bigmem_trait_df.Rds")


# Subset data ------------------------------------------------------------------

all_perf_df = all_perf_df %>%
    inner_join(all_trait %>%
                   mutate(seed = as.integer(seed)),
               by = c("species", "seed", "trait_cor"))

# Single Parameter Set CWMs
single_cwm = all_cwm %>%
    filter(k == 1.3, A > 2.5e-7, A < 2.6e-7, B > 6.3e-6, B < 6.4e-6, H > 6.3e-6,
           H < 6.4e-6, trait_cor == "uncor", R_scenar == 50, A_scenar == 50,
           H_scenar == 50, patch >= 5, patch <= 20)

# Compute Different Performances indices ---------------------------------------

perf_ind = all_perf_df %>%
    filter(patch >= 5, patch <= 20,
           R_scenar == 100, A_scenar == 100, H_scenar == 100) %>%
    group_by(seed, patch) %>%
    summarise(monocult_trait2 = weighted.mean(trait2, env_growth_rate,
                                              na.rm = TRUE),
              polycult_trait2 = weighted.mean(
                  trait2,
                  ifelse(is.na(max_growth_rate), 0, max_growth_rate),
                  na.rm = TRUE),
              cwm_trait2      = weighted.mean(trait2, N150, na.rm = TRUE))

# Plot Figures -----------------------------------------------------------------
perf_ind %>%
    tidyr::gather("perf_name", "perf_value", ends_with("trait2")) %>%
    ggplot(aes(patch, perf_value, color = perf_name)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point(alpha = 1/3)

single_cwm %>%
    ggplot(aes(patch, trait1_cwm, color = as.factor(seed))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    geom_point() +
    stat_smooth(method = "lm", geom = "line", size = 1)
