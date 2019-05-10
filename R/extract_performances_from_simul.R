#' Extract Performance Estimators From Simulations
#'
#' @param simul \[`list(1+)`\]\cr{}
#'              single simulation from [`multigen()`]
#' @param trait_list \[`list(1+)`\]\cr{}
#'                   a nested list of trait ordered by seeds used as well as
#'                   trait correlation name
#' @param realized_growth_rate \[`boolean(1)`\]\cr{}
#'                             should realized growth rate be computed
#'
#' @export
extract_performances_from_simul = function(simul, trait_list,
                                           realized_growth_rate = TRUE) {

    # Extract Compositional Data
    given_compo = simul$compo[[1]]

    max_time = dim(given_compo)[3]

    trait_df = trait_list[[simul$traits]]

    n_other = ncol(trait_df) - 2

    # Get the contribution of focal trait to each process
    contrib = extract_trait_contrib_from_scenar(simul$scenario)

    given_k = simul$k

    if (!length(given_k) == 1) {
        given_k = "gaussian"
    }

    # Compute Growth Rate ------------------------------------------------------
    env_growth_rate = simul$rmatrix %>%
        tibble::as_tibble() %>%
        tibble::rownames_to_column("patch") %>%
        dplyr::mutate(patch = as.numeric(patch)) %>%
        tidyr::gather("species", "env_growth_rate", -patch)

    # Get performance for all species and all sites
    lapply(seq_len(nrow(given_compo)), function(site_index) {

        site_abund = given_compo[site_index,,]

        # Manual derivative to get growth rate
        max_growth = apply(site_abund, 1, function(given_abund) {

            # Get first moment where species goes extinct
            time_before_extinct = which(given_abund == 0)[1] - 1

            # When species doesn't go extinct consider maximum time
            if (is.na(time_before_extinct)) {
                time_before_extinct = length(given_abund)
            }

            growth_rate = NA_real_

            if (time_before_extinct > 10) {
                given_time = 1:time_before_extinct

                growth_num = given_abund[seq(2, time_before_extinct)] -
                    given_abund[seq(time_before_extinct - 1)]

                growth_denom = given_abund[seq(time_before_extinct - 1)]
                growth_rate = max(growth_num/growth_denom, na.rm = TRUE)
            }

            ifelse(is.null(growth_rate) | is.infinite(growth_rate),
                   NA_real_, growth_rate)
        }) %>%
            tibble::enframe("species", "max_growth_rate")

        # Get optimality: distance to optimum environment given by weighted
        # average with contribution of first trait & second trait
        optim_dist = apply(trait_df, 1, function(given_traits) {
            opt_dist = weighted.mean(abs(site_index - given_traits),
                                     c(c(100 - contrib$R, contrib$R,
                                         rep(0, n_other))))

            return(opt_dist)
        }) %>%
            tibble::enframe("species", "distance_to_optimum")

        # Abundance at final time step
        sp_abund = site_abund[, max_time] %>%
            tibble::enframe("species", "N150") %>%
            dplyr::mutate(patch = site_index) %>%
            dplyr::select(patch, dplyr::everything())


        # Combine all data
        sp_abund %>%
            dplyr::inner_join(optim_dist, by = "species") %>%
            dplyr::inner_join(max_growth, by = "species")
    }) %>%
        dplyr::bind_rows() %>%
        dplyr::inner_join(env_growth_rate, by = c("patch", "species")) %>%
        dplyr::mutate(seed      = simul$seed,
                      trait_cor = simul$traits,
                      h_fun     = simul$h_fun,
                      di_thresh = simul$di_thresh,
                      k         = given_k,
                      A         = simul$A,
                      B         = simul$B,
                      H         = simul$H,
                      R_scenar  = contrib$R,
                      A_scenar  = contrib$A,
                      H_scenar  = contrib$H) %>%
        dplyr::mutate_at(dplyr::vars(dplyr::contains("growth_rate")),
                         list(per_capita = ~ . / N150))
}
