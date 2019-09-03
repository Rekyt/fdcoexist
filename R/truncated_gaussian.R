#' Truncated Gaussian Function
#'
#' Returns the expected CWM of a community based on a truncated Gaussian as
#' proposed in Denelle et al. 2019 <doi:10.1111/oik.05851> (see Appendix 1,
#' OIK-05851). Under the hypothesis of a Gaussian environmental filter.
#' For each argument of the function the user can provide one or multiple values
#' as long as the vectors have the same lengths
#' @param t_opt     \[`numeric(1+)`\]\cr{}
#'                  community optimal trait
#' @param sigma_opt \[`numeric(1+)`\]\cr{}
#'                  community optimal standard deviation
#' @param a         \[`numeric(1+)`\]\cr{}
#'                  minimum value of trait in trait range
#' @param b         \[`numeric(1+)`\]\cr{}
#'                  maximum value of trait in trait range
#'
#' @return the predicted CWM(s)
#' @export
truncated_gaussian = function(t_opt, sigma_opt, a, b) {

    arg_lengths = lapply(list(t_opt, sigma_opt, a, b), length)

    if (!all(sapply(arg_lengths, function(x) x == length(b)))) {
        stop("Not all arguments are of the same size", call. = FALSE)
    }

    alpha = (a - t_opt)/sigma_opt
    beta  = (b - t_opt)/sigma_opt

    z = pnorm(beta) - pnorm(alpha)

    t_opt + ((dnorm(alpha) - dnorm(beta))/z) * sigma_opt
}
