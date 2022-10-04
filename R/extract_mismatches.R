#' Extract species mismatches
#'
#' @param x Simulation object from [multigen()]
#' @param z Number of the species
#'
#' @export
extract_mismatches <- function(x, z){
	dat <- x[["env_all"]]
	a <- dat[which(dat$sp == paste("sp", z, sep = "")), ]
	env <- which.max(a$r_env)
	r <- max(a$r_env, na.rm = TRUE)

	dat <- x[["r_all"]]
	a <- dat[which(dat$sp == paste("sp", z, sep = "")), ]
	env1 <- which.max(a$r_env)
	r1 <- max(a$r_env, na.rm = TRUE)

	dat <- x[["r_comp"]]
	a <- dat[which(dat$sp == paste("sp", z, sep = "")), ]
	env2 <- which.max(a$r_env)
	r2 <- max(a$r_env, na.rm = TRUE)

	dat <- x[["r_envab"]]
	a <- dat[which(dat$sp == paste("sp", z, sep = "")), ]
	env3 <- which.max(a$r_env)
	r3 <- max(a$r_env, na.rm = TRUE)

	return(c(env, r, env1, r1, env2, r2, env3, r3))
}
