na_if_null_or_inf = function(vect) {
    ifelse(is.null(vect) | is.infinite(vect), NA_real_, vect)
}
