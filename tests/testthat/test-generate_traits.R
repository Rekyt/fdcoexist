context("test functions in generate_traits.R")

test_that("Trait generation works", {

    n_sp = sample(25, size = 1)
    n_patch = sample(2:25, size = 1)
    n_traits = sample(25, size = 1)

    expect_silent(generate_cor_traits(n_patch, n_sp, n_traits))


    # Test values
    given_traits = generate_cor_traits(n_patch, n_sp, n_traits)

    expect_is(given_traits, "matrix")

    expect_equal(dimnames(given_traits),
                 list(paste0("species", seq(n_sp)),
                      paste0("trait", seq(n_traits + 1))))

    expect_equivalent(given_traits[, 1], seq(1, n_patch, length.out = n_sp))

    expect_equal(range(given_traits[, 2]), c(1, n_patch))

    expect_equal(dim(given_traits), c(n_sp, n_traits + 1))


    # Test correlation

    low_cor = generate_cor_traits(n_patch, 1000, 1, cor_coef = 0.3)
    mid_cor = generate_cor_traits(n_patch, 1000, 1, cor_coef = 0.7)
    high_cor = generate_cor_traits(n_patch, 1000, 1, cor_coef = 0.9)

    expect_equivalent(cor(low_cor),
                      matrix(c(1, 0.3, 0.3, 1), nrow = 2, ncol = 2,
                             dimnames = list(paste0("trait", 1:2),
                                             paste0("trait", 1:2))),
                      tolerance = 5e-2)

    expect_equivalent(cor(mid_cor),
                      matrix(c(1, 0.7, 0.7, 1), nrow = 2, ncol = 2,
                             dimnames = list(paste0("trait", 1:2),
                                             paste0("trait", 1:2))),
                      tolerance = 5e-2)

    expect_equivalent(cor(high_cor),
                      matrix(c(1, 0.9, 0.9, 1), nrow = 2, ncol = 2,
                             dimnames = list(paste0("trait", 1:2),
                                             paste0("trait", 1:2))),
                      tolerance = 5e-2)
})

test_that("Trait Weights generation works", {

    expect_error(create_trait_weights("a", "nb", "k"),
                 "All provided weights should be integer", fixed = TRUE)

    expect_error(create_trait_weights(-5, 6, 76),
                 "All weights should be in [0; 100]", fixed = TRUE)

    expect_equal(create_trait_weights(5, 0, 0),
                 data.frame(
                     trait            = paste0("trait", 1:3),
                     growth_weight    = c(0.95, 0.05, 0),
                     compet_weight    = c(   0,    0, 1),
                     hierarchy_weight = c(   0,    0, 0)
                 ))

    expect_equal(create_trait_weights(5, 0, 50),
                 data.frame(
                     trait            = paste0("trait", 1:3),
                     growth_weight    = c(0.95, 0.05, 0),
                     compet_weight    = c(   0,    0, 1),
                     hierarchy_weight = c(   0,  0.5, 0)
                 ))
})
