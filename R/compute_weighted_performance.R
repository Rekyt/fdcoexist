#' Community Performance Estimators
#'
#' Compute for each community all traits:
#' * weighted mean of trait by species purely environmental growth rate
#' `pur_env`,
#' * weighted mean of trait by species realized growth rate `weighted_growth`,
#' * weighted mean of trait by species abundance `cwm`,
#' * trait of species with highest growth rate per community `best_growth`,
#' * trait of species with highest abundance per community `best_abund`,
#' * other community weighted moments such as variance, skewness and kurtosis.
#'
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

    # Species Richness per community
    species_rich = perf_df %>%
        summarise(species_rich = sum(N150 > 0, na.rm = TRUE))

    # Pure Environmental Filtering ---------------------------------------------
    envbest_growth = perf_df %>%
        filter_top_perf_per_trait(env_growth_rate, envbest_growth)

    # Estimate performance using best perf. species ----------------------------
    best_growth = perf_df %>%
        filter_top_perf_per_trait(max_growth_rate, best_growth)

    best_abund = perf_df %>%
        filter_top_perf_per_trait(N150, best_abund)

    # Combine all estimators ---------------------------------------------------
    full_perf = perf_df %>%
        # Weighted Estimators
        summarise_at(
            vars(matches("trait[0-9]+")),
            list(
                # Consider Only Environmental Fitting
                pure_env        = ~wtd_mean(., env_growth_rate, na.rm = TRUE),
                # Weight traits using growth rate (analogous of CWM but with GR)
                # need to rescale growth rate between 0 and 1 to account for
                # negative growth rate in some patches
                weighted_growth = ~wtd_mean(., max_growth_rate %>%
                                            scales::rescale(c(0,1)),
                                            na.rm = TRUE),
                # Community Weighted Moments
                cwm             = ~wtd_mean(.,     N150, na.rm = TRUE),
                cwv             = ~wtd_var(.,      N150, na.rm = TRUE),
                cws             = ~wtd_skewness(., N150, na.rm = TRUE),
                cwk             = ~wtd_kurtosis(., N150, na.rm = TRUE)))
    list(full_perf,
         envbest_growth,
         best_growth,
         best_abund,
         species_rich) %>%
        {Reduce(function(x, y) full_join(x, y, by = c(group_vars(x), "patch")),
                .)} %>%
        group_by(patch, add = TRUE) %>%
        ungroup()
}

filter_top_perf_per_trait = function(df, wt, given_name) {

    df %>%
        tidyr::gather("trait_name", "trait_value", matches("trait[0-9]+")) %>%
        group_by(trait_name, add = TRUE) %>%
        # Select first values with maximum weight
        top_n(1, !!enquo(wt)) %>%
        select(group_vars(.), trait_name, trait_value) %>%
        # In the highly improbable case of ex-aequo per community compute the
        # average trait
        summarise(trait_value = mean(trait_value, na.rm = TRUE)) %>%
        rename(!!enquo(given_name) := trait_value) %>%
        # Renaming the column properly
        tidyr::gather("perf_type", "perf_value", !!enquo(given_name)) %>%
        mutate(perf_type = paste0(trait_name, "_", perf_type)) %>%
        select(-trait_name) %>%
        tidyr::spread(perf_type, perf_value)
}
