#' Community Performance Estimators
#'
#' Compute for each community all traits:
#' * trait of best monoculture
#' * trait of best polyculture
#' * trait CWM
#' * weighted trait by monoculture growth
#' * weighted trait by polyculture growth
#' @param perf_df  species performance data frame output by
#'                 [extract_performances_from_simul()]
#' @param trait_df trait data.frame containing columns `trait_cor`, `seed` and
#'                 `species`
#'
#' @import dplyr
#' @export
compute_weighted_performance = function(perf_df, trait_df) {

    # Add trait information to performances
    perf_df = perf_df %>%
        inner_join(trait_df %>%
                              mutate(seed = as.integer(seed)),
                   by = c("trait_cor", "seed", "species")) %>%
        group_by(k, A, B, H, R_scenar, A_scenar, H_scenar, trait_cor, seed,
                 patch)

    species_rich = perf_df %>%
        summarise(species_rich = sum(N150 > 0, na.rm = TRUE))

    # Compute monoculture estimator using best environmental growth per patch
    best_mono_growth = perf_df %>%
        filter(env_growth_rate == max(env_growth_rate)) %>%
        rename_at(vars(matches("trait[0-9]+")), ~ paste0("monobest_", .)) %>%
        select_at(vars(group_cols(), starts_with("monobest_"))) %>%
        ungroup()

    # Compute polycuture estimator using best realized growth per patch
    best_poly_growth = perf_df %>%
        filter(max_growth_rate == max(max_growth_rate, na.rm = TRUE)) %>%
        rename_at(vars(matches("trait[0-9]+")), ~ paste0("polybest_", .)) %>%
        select_at(vars(group_cols(), starts_with("polybest_"))) %>%
        ungroup()

    # Combine all estimators
    perf_df %>%
        # Weighted Estimators
        summarise_at(
            vars(matches("trait[0-9]+")),
            list(monocult = ~weighted.mean(., env_growth_rate, na.rm = TRUE),
                 polycult = ~weighted.mean(
                     ., ifelse(is.na(max_growth_rate), 0, max_growth_rate),
                     na.rm = TRUE),
                 cwm          = ~weighted.mean(., N150, na.rm = TRUE),
                 cwv          = ~wtd_var(., N150, na.rm = TRUE),
                 cws          = ~wtd_skewness(., N150, na.rm = TRUE),
                 cwk          = ~wtd_kurtosis(., N150, na.rm = TRUE))) %>%
        # Reverse suffix to prefix
        rename_at(vars(ends_with("monocult")),
                  ~paste0("monocult_",
                          gsub("_monocult", "", .,fixed = TRUE))) %>%
        rename_at(vars(ends_with("polycult")),
                  ~paste0("polycult_",
                          gsub("_polycult", "", .,fixed = TRUE))) %>%
        rename_at(vars(ends_with("cwm")),
                  ~paste0("cwm_",
                          gsub("_cwm", "", .,fixed = TRUE))) %>%
        # Add estimators from best growth
        inner_join(best_mono_growth, by = c(group_vars(.), "patch")) %>%
        inner_join(best_poly_growth, by = c(group_vars(.), "patch")) %>%
        inner_join(species_rich,     by = c(group_vars(.), "patch")) %>%
        ungroup()
}
