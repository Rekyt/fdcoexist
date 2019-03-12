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

    set.seed(seed_number)
    uncor_traits = generate_cor_traits(n_patches, n_species, n_traits - 1,
                                       cor_coef = 0)
    all_cor = list(uncor = uncor_traits)

    if (!is.null(given_traits)) {
        all_cor = given_traits
    }

    guessed_th_max = lapply(all_cor, max)
    guessed_th_max = max(unlist(guessed_th_max))
    guessed_th_min = lapply(all_cor, min)
    guessed_th_min = max(unlist(guessed_th_min))

    compo_dim = dim(given_composition)
    given_gen = compo_dim[3]
    given_species = compo_dim[2]
    given_patches = compo_dim[1]

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
