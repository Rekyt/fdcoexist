#' Extract Performance Estimators From Simulations
#'
#' @param simul \[`list(1+)`\]\cr{}
#'              single simulation from [`multigen()`]
#' @param trait_list \[`list(1+)`\]\cr{}
#'                   a nested list of trait ordered by seeds used as well as
#'                   trait correlation name
#' @param chosen_time \[`integer(1)`\]\cr{}
#'                    generation which should be used to get final abundance
#'
#'
#' @importFrom stats coef lm
#' @export
extract_performances_from_simul = function(simul, trait_list,
                                           chosen_time = NULL) {

    # Extract information from simulation data ---------------------------------
    # Extract Compositional Data
    given_compo = simul$compo[[1]]

    if (is.null(chosen_time)) {
        chosen_time = dim(given_compo)[3]
    }
    max_time = chosen_time  # maximum time to consider

    trait_df = trait_list[[simul$traits]]

    n_other = ncol(trait_df) - 2

    # Get the contribution of focal trait to each process
    contrib = extract_trait_contrib_from_scenar(simul$scenario)

    given_k = simul$k

    # Environmental Growth Rate ------------------------------------------------
    env_growth_rate = simul$rmatrix %>%
        tibble::as_tibble() %>%
        tibble::rownames_to_column("patch") %>%
        dplyr::mutate(patch = as.numeric(patch)) %>%
        tidyr::gather("species", "env_growth_rate", -patch)

    # Realized Growth Rate -----------------------------------------------------
    # Get performance for all species and all sites
    sp_perf = lapply(seq_len(nrow(given_compo)), function(site_index) {

        site_abund = given_compo[site_index,, 1:max_time]

        # Compute all types of growth rate
        growth_rates = apply(site_abund, 1, function(given_abund) {

            # Get first moment where species goes extinct
            time_before_extinct = which(given_abund == 0)[1] - 1

            # When species doesn't go extinct consider maximum time
            if (is.na(time_before_extinct)) {
                time_before_extinct = length(given_abund)
            }

            avg_growth_rate = NA_real_
            max_growth_rate = NA_real_
            int_growth_rate = NA_real_

            if (time_before_extinct >= ifelse(max_time > 10, 10, max_time)) {

                # Get abundance at t + 1
                growth_num = given_abund[seq(2, time_before_extinct)]

                # Get abundance at t
                growth_denom = given_abund[seq(time_before_extinct - 1)]

                # Take out cases where difference in abundances is smaller than
                # 0
                consider_cases = which(abs(log(growth_num/growth_denom)) >
                                           .Machine$double.eps)

                # Compute average growth rate
                avg_growth_rate = mean(growth_num[consider_cases] /
                                       growth_denom[consider_cases],
                                   na.rm = TRUE)

                # Naive maximum growth rate (maximum growth rate over all gen.)
                max_growth_rate = max(growth_num[consider_cases] /
                                          growth_denom[consider_cases],
                                      na.rm = TRUE)

                # Estimate of intrinsic growth rate using model
                abund_model = lm(
                    log_abund ~ time,
                    data = data.frame(
                        log_abund = log(given_abund[1:time_before_extinct]),
                        time = 1:time_before_extinct))

                int_growth_rate = exp(coef(abund_model)[2])

            }

            # Final data.frame
            data.frame(
                avg_growth_rate = na_if_null_or_inf(avg_growth_rate),
                max_growth_rate = na_if_null_or_inf(max_growth_rate),
                int_growth_rate = na_if_null_or_inf(int_growth_rate)
            )
        }) %>%
            bind_rows(.id = "species")

        # Abundance at final time step
        sp_abund = site_abund[, max_time] %>%
            tibble::enframe("species", "final_abundance") %>%
            dplyr::mutate(patch = site_index) %>%
            dplyr::select(patch, dplyr::everything())


        # Combine all data
        sp_abund %>%
            dplyr::inner_join(growth_rates, by = "species")
    })

    sp_perf %>%
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
                      H_scenar  = contrib$H,
                      time = max_time)
}
