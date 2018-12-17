context("test support functions")

# Preamble ---------------------------------------------------------------------

sp_traits = matrix(c(1, 2, 1, 1, 2, 0), ncol = 2)
dimnames(sp_traits) = list(species = paste0("sp", 1:3),
                           trait   = paste0("trait", 1:2))
trait_weights = data.frame(trait = paste0("trait", 1:2),
                           growth_weight    = c(1, 0),
                           compet_weight    = c(0, 1),
                           hierarchy_weight = c(0, 1),
                           stringsAsFactors = FALSE)


N0 = matrix(rep(c(50, 10), 3, each = TRUE), ncol = 3,
            dimnames = list(site    = paste0("patches", 1:2),
                            species = paste0("sp", 1:3)))

given_compo = simplify2array(lapply(seq(1, 0.25, length.out = 5),
                                    function(x) x * N0))
# Tests ------------------------------------------------------------------------

test_that("landscape_df() works", {

    expect_silent(landscape_df(list(compo = given_compo), 5))
})
