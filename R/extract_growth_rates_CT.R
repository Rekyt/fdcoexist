
extract_growth_rates = function(simul, chosen_time = NULL){

    if(!is.null(chosen_time)){
        if(chosen_time < 2){
            stop("chosen_time must be an integer superior to 1.")
        }
        simul$compo <- simul$compo[, , 1:chosen_time]
    }

    gw <- lapply(seq_len(nrow(simul$compo)), function(site_index) {
        site_abund = simul$compo[site_index,, 1:chosen_time]

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

            if (time_before_extinct >= ifelse(chosen_time > 10, 10,
                                              chosen_time)) {

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
        })  %>%
            bind_rows(.id = "species")

        # Abundance at final time step
        sp_abund = site_abund[, chosen_time] %>%
            tibble::enframe("species", "final_abundance") %>%
            dplyr::mutate(patch = site_index) %>%
            dplyr::select(patch, dplyr::everything())


        # Combine all data
        sp_abund %>%
            dplyr::inner_join(growth_rates, by = "species")
    })

    # Bind list into one data frame
    gw <- as.data.frame(bind_rows(gw))

    return(gw)
}

# Graphical function to plot environmental response curves of species

#' Plot environmental response curves of species in all the patches
#'
#' Automatically assigns one color per species.
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param sp integer for the number of species to consider
#' @param plot a boolean determining whether to plot or not the response curves
#'
#' @import ggplot2
#' @export
#'
#'
r_env <- function(simul, n_patches, sp, time, plot = TRUE){
    sp <- 1:sp

    env_growth <- data.frame(simul[["rmatrix"]][, sp])
    env_growth$env <- 1:n_patches
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        rename(max_r_env = env) %>%
        as.data.frame()

    env_growth <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")],
                            by = "sp")

    if(plot == TRUE){
        env_plot <- ggplot(env_growth, aes(env, r_env)) +
            geom_line(aes(color = as.factor(sp))) +
            geom_label(aes(label = ifelse(env == max_r_env, sp, NA)),
                       size = 3) +
            labs(x = "Patch", y = expression("r"["env"]),
                 title = "Environmental response curve of species") +
            theme_classic() +
            guides(color = FALSE) +
            theme(panel.border = element_rect(fill = NA))

        return(list(env_growth, env_plot))
    }else {
        return(env_growth)
    }
}

r_all <- function(simul, n_patches, sp, time=50, plot = TRUE){
    sp <- 1:sp
	
    env_growth <- data.frame(simul[["r_tot"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        rename(max_r_env = env) %>%
        as.data.frame()

    env_growth <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")
        return(env_growth)
    }
    
r_compo <- function(simul, n_patches, sp, time=50, plot = TRUE){
    sp <- 1:sp

    env_growth <- data.frame(simul[["compo"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        rename(max_r_env = env) %>%
        as.data.frame()

    env_growth <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")
        return(env_growth)
    }    
    
r_envab <- function(simul, n_patches, sp, time=50, plot = TRUE){
    sp <- 1:sp

    env_growth <- data.frame(simul[["env_ab"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- gather(env_growth, "sp", "r_env", contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        group_by(sp) %>%
        top_n(1, r_env) %>%
        rename(max_r_env = env) %>%
        as.data.frame()

    env_growth <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")
        return(env_growth)
    }      