#' Simulations with various parameters
#'
#' @param seed_number       the actual seed used to generate traits
#' @param given_k           the maximum achievable growth
#' @param given_A           the intensity of limiting similarity
#' @param given_B           the intensity of intra-specific competition
#' @param given_H           the intensity of hierarchical competition
#' @param given_scenars     if provided the data.frame of scenarios
#' @param given_traits      if provided the traits data.frame to use otherwise
#'        generate uncorrelated traits
#' @param given_h_fun       function to use to combine hierarchical competition
#' @param given_di_thresh   the threshold distance at which species are
#'        considered maximally dissimilar and do not compete anymore in
#'        limiting similarity
#' @param given_env         the environmental width of the different patches
#' @param given_composition the empty matrix with good dimensions to be used
#'        to follow population dynamics
#' @param given_d           the index of dispersal
#'
#' @export
meta_simul = function(seed_number, given_k = k, given_A = A,
                      given_B = B,
                      given_H = H, given_scenars = NULL,
                      given_traits = NULL,
                      given_h_fun  = "sum",
                      given_di_thresh = 24,
                      given_env,
                      given_composition,
                      given_d) {

    our_scenars = list(R50A50H50 = create_trait_weights(50, 50, 50))

    if (!is.null(given_scenars)) {
        our_scenars = given_scenars
    }

    compo_dim = dim(given_composition)
    given_gen = compo_dim[3]
    given_species = compo_dim[2]
    given_patches = compo_dim[1]


    set.seed(seed_number)
    uncor_traits = generate_cor_traits(given_patches, given_species, 3,
                                       cor_coef = 0)
    all_cor = list(uncor = uncor_traits)

    if (!is.null(given_traits)) {
        all_cor = given_traits
    }

    guessed_th_max = lapply(all_cor, max)
    guessed_th_max = max(unlist(guessed_th_max))
    guessed_th_min = lapply(all_cor, min)
    guessed_th_min = max(unlist(guessed_th_min))

        all_env = list(constant = rep(given_patches, 5))

    all_compet = list(compet = list(A        = given_A,
                                    B        = given_B,
                                    H        = given_H,
                                    B_over_A = given_B/given_A))


    all_facets = expand.grid(compet_status = names(all_compet),
                             env_width     = names(all_env),
                             cor_level     = names(all_cor),
                             scenario      = names(our_scenars))



    apply(all_facets, 1, function(given_row) {
        simul = multigen(
            traits        = all_cor[[given_row[["cor_level"]]]],
            trait_weights = our_scenars[[given_row[["scenario"]]]],
            env           = given_env,
            time          = given_gen,
            species       = given_species,
            patches       = given_patches,
            composition   = given_composition,
            A             = all_compet[[given_row[["compet_status"]]]][["A"]],
            B             = all_compet[[given_row[["compet_status"]]]][["B"]],
            d             = given_d,
            k             = given_k,
            H             = all_compet[[given_row[["compet_status"]]]][["H"]],
            th_max        = guessed_th_max,
            th_min        = guessed_th_min,
            width         = all_env[[given_row[["env_width"]]]],
            h_fun         = given_h_fun,
            di_thresh     = given_di_thresh)

        return(
            list(
                k             = given_k,
                A             = all_compet[[given_row[["compet_status"]]]][["A"]],
                B             = all_compet[[given_row[["compet_status"]]]][["B"]],
                B_over_A      = all_compet[[given_row[["compet_status"]]]][["B_over_A"]],
                H             = all_compet[[given_row[["compet_status"]]]][["H"]],
                seed          = seed_number,
                compet_status = given_row[["compet_status"]],
                env           = given_row[["env_width"]],
                traits        = given_row[["cor_level"]],
                scenario      = given_row[["scenario"]],
                compo         = list(simul$compo),
                rmatrix       = simul$rmatrix,
                rmatenv       = simul$rmatenv,
                alphaterm     = simul$alpha,
                h_fun         = given_h_fun,
                di_thresh     = given_di_thresh))
    })
}

#' @export
extract_performances_from_simul = function(simul, trait_list,
                                           realized_growth_rate = TRUE) {

    given_compo = simul$compo[[1]]

    max_time = dim(given_compo)[3]

    trait_df = trait_list[[simul$seed]][[simul$traits]]

    n_other = ncol(trait_df) - 2

    contrib = extract_trait_contrib_from_scenar(simul$scenario)

    given_k = simul$k

    if (!length(given_k) == 1) {
        given_k = "gaussian"
    }

    th_growth_rate = simul$rmatrix %>%
        tibble::as_tibble() %>%
        tibble::rownames_to_column("patch") %>%
        dplyr::mutate(patch = as.numeric(patch)) %>%
        tidyr::gather("species", "th_growth_rate", -patch)

    env_growth_rate = simul$rmatenv %>%
        tibble::as_tibble() %>%
        tibble::rownames_to_column("patch") %>%
        dplyr::mutate(patch = as.numeric(patch)) %>%
        tidyr::gather("species", "env_growth_rate", -patch)

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
            enframe("species", "max_growth_rate")

        # Get optimality
        optim_dist = apply(trait_df, 1, function(given_traits) {
            opt_dist = weighted.mean(abs(site_index - given_traits),
                                     c(c(100 - contrib$R, contrib$R,
                                         rep(0, n_other))))

            return(opt_dist)
        }) %>%
        tibble::enframe("species", "distance_to_optimum")

        # Abundance
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
        dplyr::inner_join(th_growth_rate, by = c("patch", "species")) %>%
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
        dplyr::mutate_at(dplyr::vars(dplyr::contains("growth_rate")), list(per_capita = ~ . / N150))
}

