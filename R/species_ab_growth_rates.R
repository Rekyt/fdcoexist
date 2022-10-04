# Graphical function to plot environmental response curves of species

#' Plot environmental response curves of species in all the patches
#'
#' Automatically assigns one color per species.
#'
#' @param simul results array of a given simulation (a [multigen()] output)
#' @param n_patches an integer indicating the patch number
#' @param sp integer for the number of species to consider
#' @param time integer determining the time step of interest
#'
#' @import ggplot2
#' @export
#'
#'
sp_ab_gr <- function(simul, n_patches, sp, time){
    sp <- 1:sp # vector of species of interest

    # Growth rate determined solely by the environmental filtering
    env_growth <- data.frame(simul[["rmatrix"]][, sp])
    # Growth rate determined by the numerator of the eq. (env filt + hierarch)
    env_hierarch_growth <- data.frame(simul[["r_tot"]][, sp, time])
    # Observed abundances at the time step of interest
    obs_ab <- data.frame(simul[["compo"]][, sp, time])
    # Abundances determined solely by the environmental filtering
    env_ab <- data.frame(simul[["env_ab"]][, sp, time])

    # Long-format dataframes
    env_growth$env <- 1:n_patches
    env_growth <- tidyr::gather(env_growth, "sp", "r_env",
                                tidyr::contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    env_hierarch_growth$env <- 1:n_patches
    env_hierarch_growth <- tidyr::gather(
        env_hierarch_growth, "sp", "r_env_hierarch", tidyr::contains("species")
    )
    env_hierarch_growth$sp <- gsub("species", "sp", env_hierarch_growth$sp)

    obs_ab$env <- 1:n_patches
    obs_ab <- tidyr::gather(obs_ab, "sp", "obs_ab", tidyr::contains("species"))
    obs_ab$sp <- gsub("species", "sp", obs_ab$sp)

    env_ab$env <- 1:n_patches
    env_ab <- tidyr::gather(env_ab, "sp", "env_ab", tidyr::contains("species"))
    env_ab$sp <- gsub("species", "sp", env_ab$sp)

    # Patch where each species has its maximal growth rate/obs ab/...
    max_r_env <- env_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env) %>%
        dplyr::rename(max_r_env = env) %>%
        as.data.frame()

    max_r_env_hierarch <- env_hierarch_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env_hierarch) %>%
        dplyr::rename(max_r_env_hierarch = env) %>%
        as.data.frame()

    max_obs_ab <- obs_ab %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, obs_ab) %>%
        dplyr::rename(max_obs_ab = env) %>%
        as.data.frame()

    max_env_ab <- env_ab %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, env_ab) %>%
        dplyr::rename(max_env_ab = env) %>%
        as.data.frame()

    # Adding the maximums to the corresponding dataframes
    env_growth <- dplyr::left_join(
        env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp"
    )

    env_hierarch_growth <- dplyr::left_join(
        env_hierarch_growth,
        max_r_env_hierarch[, c("sp", "max_r_env_hierarch")], by = "sp"
    )

    obs_ab <- dplyr::left_join(
        obs_ab, max_obs_ab[, c("sp", "max_obs_ab")], by = "sp"
    )

    env_ab <- dplyr::left_join(
        env_ab, max_env_ab[, c("sp", "max_env_ab")], by = "sp"
    )

    # Merging all the dataframes together
    all_dat <- dplyr::left_join(
        env_growth, env_hierarch_growth, by = c("env", "sp")
    ) %>%
        dplyr::left_join(all_dat, obs_ab, by = c("env", "sp")) %>%
        dplyr::left_join(all_dat, env_ab, by = c("env", "sp"))

    return(all_dat)

}
