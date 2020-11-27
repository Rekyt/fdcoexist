# Script with functions to compute weighted moments of trait distribution

#' Weighted mean that allows NA in values and weights
#'
#' @param x     values whose weighted mean is to be computed
#' @param w     weights of the same length as `x` to be used against x
#' @param na.rm a logical values indicating whether `NA` values in both `x` and
#'              `w` should be stripped before computation
#'
#' @export
wtd_mean = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    weighted.mean(x, w)
}

#' Weighted Variance
#'
#' @inheritParams wtd_mean
#' @export
wtd_var = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.var(x, w)
}

#' Weighted Skewness with na.rm
#'
#' Compute weighted Skewness using [`Weighted.Desc.Stat::w.skewness`] but add an
#' option to remove `NA` values
#'
#' @inheritParams wtd_mean
#'
#' @export
wtd_skewness = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.skewness(x, w)
}

#' Weighted Kurtosis with na.rm
#'
#' Compute weighted kurtosis using [`Weighted.Desc.Stat::w.kurtosis`] but add an
#' option to remove `NA` values
#'
#' @inheritParams wtd_mean
#'
#' @export
wtd_kurtosis = function(x, w, na.rm = TRUE) {

    if (na.rm) {
        s <- !is.na(x + w)
        x <- x[s]
        w <- w[s]
    }

    Weighted.Desc.Stat::w.kurtosis(x, w)
}
