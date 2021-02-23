#' Extract Species by Patch Growth Rates and Optimal Patches
#'
#' @param simul     a simulation output from `multigen()`
#' @param n_patches the number of patches in the simulation
#' @param sp1       the number of species
#' @param time      the time slice at which e
#'
#' @return
#' Outputs 4 data.frame of species by patch mismatches each with four columns:
#'
#' * `env_all` -> Environmental Growth only
#' * `r_all`   -> Environmental + Hierarchical competition Growth only
#' * `r_comp`  -> Abundance
#' * `r_envab` -> Abundance based on environmental filtering only
#'
#' The 4 data.frame are:
#'
#' * `env`       -> Patch
#' * `sp`        -> Species
#' * `r_env`     -> Observed statistic (environmental growth, total growth,
#'                   abundance, environmental filtering abundance)
#' * `max_r_env` -> Patch of maximum observed statistic
#
r_env_CT <- function(simul, n_patches, sp1, time) {
    sp <- 1:sp1
    env_growth <- data.frame(simul[["rmatrix"]][, sp])
    env_growth$env <- 1:n_patches
    env_growth <- tidyr::gather(env_growth, "sp", "r_env",
                                tidyr::contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate based on environment
    # only
    max_r_env <- env_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env) %>%
        dplyr::rename(max_r_env = env) %>%
        as.data.frame()

    renv <- dplyr::left_join(env_growth, max_r_env[, c("sp", "max_r_env")],
                             by = "sp")

    # Growth of species based on total growth termsÒ
    # (environment + hierarchical competition)
    sp <- 1:sp1
    env_growth <- data.frame(simul[["r_tot"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- tidyr::gather(env_growth, "sp", "r_env",
                                tidyr::contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal total
    # (env + hierarch. competition) growth rate
    max_r_env <- env_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env) %>%
        dplyr::rename(max_r_env = env) %>%
        as.data.frame()

    rtot <- dplyr::left_join(env_growth, max_r_env[, c("sp", "max_r_env")],
                             by = "sp")

    #r_comp
    sp <- 1:sp1
    env_growth <- data.frame(simul[["compo"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- tidyr::gather(env_growth, "sp", "r_env",
                                tidyr::contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env) %>%
        dplyr::rename(max_r_env = env) %>%
        as.data.frame()

    rcomp <- dplyr::left_join(env_growth, max_r_env[, c("sp", "max_r_env")],
                              by = "sp")


    #r_envab
    sp <- 1:sp1
    env_growth <- data.frame(simul[["env_ab"]][, sp, (time-1)])
    rownames(env_growth) <- rownames(data.frame(simul[["rmatrix"]][, sp]))
    env_growth$env <- 1:n_patches
    env_growth <- tidyr::gather(env_growth, "sp", "r_env",
                                tidyr::contains("species"))
    env_growth$sp <- gsub("species", "sp", env_growth$sp)

    # Patch where each species has its maximal growth rate
    max_r_env <- env_growth %>%
        dplyr::group_by(sp) %>%
        dplyr::top_n(1, r_env) %>%
        dplyr::rename(max_r_env = env) %>%
        as.data.frame()

    envab <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")

    return(list(env_all = renv, r_all = rtot, r_comp = rcomp, r_envab = envab))

}
