#Combination of all my 'special' functions

## Equilibrium function

# Function to compute time before equilibrium
equilibrium <- function(simul, n_gen, patch){
    nb_ind_df <- data.frame(time_step = 1,
                            nb_ind = round(sum(simul$compo[patch, , 1]), 0),
                            diff = NA)
    for(i in 2:n_gen) {
        nb_ind_df <- rbind(
            nb_ind_df,
            data.frame(time_step = i,
                       nb_ind = round(sum(simul$compo[patch, , i]), 0),
                       diff = sum(abs(round(simul$compo[patch, , i], 0) -
                                          round(simul$compo[patch, , (i-1)], 0)))))
    }

    first_0 <- nb_ind_df[which(nb_ind_df$diff == 0), "time_step"][1]

    # Control for potential later outbreak of a species
    if (!is.na(first_0)) {
        # If more than 10 individuals difference during one time step,
        # search for another equilibrium after that 10 diff
        if (sum(abs(nb_ind_df[first_0:nrow(nb_ind_df), "diff"]) > 10) > 0) {
            last_10 <- nb_ind_df[which(abs(nb_ind_df$diff) > 10),
                                 "time_step"]
            last_10 <- last_10[length(last_10)]

            nb_ind_df <- nb_ind_df[last_10:nrow(nb_ind_df), ]

            first_0 <- nb_ind_df[which(nb_ind_df$diff == 0), "time_step"][1]
        }

    }

    require(vegan)

    first_0_NA <- ifelse(is.na(first_0), n_gen, first_0)

    return(data.frame(
        rich = length(simul$compo[patch, , first_0_NA][simul$compo[patch, , first_0_NA] > 0]),
        even = diversity(round(simul$compo[patch,, first_0_NA], 0))/
            log(specnumber(round(simul$compo[patch,, first_0_NA], 0))),
        nb_ind = round(sum(simul$compo[patch, , first_0_NA])),
        eq = first_0))
}

###
# Generate random traits but how different is it from our method?
# Compared to generate_cor_traits() introduce a little of variability in first
# trait as instead of being directly determined by the species number it adds
# little white noise to it and scale it to a minimum of 0 if negative and
# maximum of 25 if maximum value is over 25.
generate_cor_traits_rand <- function (number_patches, number_species,
                                      number_other = 9, cor_coef = 0.7,
                                      min_value = 1)
 {
     x1 = seq(min_value, number_patches, length.out = number_species)
     x1 = x1 + rnorm(length(x1), mean = 0, sd = 1)
     x1 = sort(ifelse(x1 < 0, 0, x1))
     x1 = ifelse(x1 < number_patches, x1, number_patches)
     cov_mat = matrix(c(1, cor_coef, 0, sqrt(1 - cor_coef^2)), nrow = 2)
     n_other <- number_other
     res <- lapply(seq(n_other), function(x) {
         x2 <- stats::runif(number_species, 1, number_patches)
         cov_mat %*% matrix(c(x1, x2), nrow = 2, byrow = TRUE)
     })
     res <- do.call(rbind, lapply(res, function(x) x[2, ]))
     res <- rbind(x1, res)
     res <- t(res)
     trait_mat <- cbind(res[, 1], sapply(2:ncol(res), function(x) {
         val <- res[, x]
         (val - min(val))/(max(val) - min(val)) * (diff(range(x1))) +
             min(x1)
     }))
     colnames(trait_mat) <- paste0("trait", seq(number_other +
         1))
     rownames(trait_mat) <- paste0("species", seq(number_species))
     return(trait_mat)
 }



###extract growth rates
r_env <- function(simul, n_patches, sp1, time) {
    sp <- 1:sp1
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

renv <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")

#r_all
    sp <- 1:sp1
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

rtot <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")

#r_comp
    sp <- 1:sp1
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

    rcomp <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")


#r_envab
    sp <- 1:sp1
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

    envab <- left_join(env_growth, max_r_env[, c("sp", "max_r_env")], by = "sp")

return(list(env_all= renv, r_all= rtot, r_comp= rcomp, r_envab= envab))

}

#
extract_mismatchesCT <- function(x=all_growth, z=z){
	dat <- x[["env_all"]]
	a <- dat[which(dat$sp==paste("sp", z, sep="")), ]
	env <- which.max(a$r_env)
	r <- max(a$r_env, na.rm=TRUE)

	dat <- x[["r_all"]]
	a <- dat[which(dat$sp==paste("sp", z, sep="")), ]
	env1 <- which.max(a$r_env)
	r1 <- max(a$r_env, na.rm=TRUE)

	dat <- x[["r_comp"]]
	a <- dat[which(dat$sp==paste("sp", z, sep="")), ]
	env2 <- which.max(a$r_env)
	r2 <- max(a$r_env, na.rm=TRUE)

	dat <- x[["r_envab"]]
	a <- dat[which(dat$sp==paste("sp", z, sep="")), ]
	env3 <- which.max(a$r_env)
	r3 <- max(a$r_env, na.rm=TRUE)

	return(c(env, r, env1, r1, env2, r2, env3, r3))
}



