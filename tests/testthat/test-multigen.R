context("Test multigen functions")

# Preamble ---------------------------------------------------------------------

sp_traits = matrix(c(1, 2, 1, 1, 2, 0), ncol = 2)
dimnames(sp_traits) = list(species = paste0("sp", 1:3),
                           trait   = paste0("trait", 1:2))
trait_weights = data.frame(trait = paste0("trait", 1:2),
                           growth_weight = c(1, 0),
                           compet_weight = c(0, 1))

# Actual tests -----------------------------------------------------------------

test_that("testing check_trait_weights()", {

    ## Input errors
    # Incorrect trait weights input
    expect_error(check_trait_weights(c(1, 0, 1, 0), sp_traits),
                 "trait_weights should be a data.frame or a matrix")

    # Correct Input
    expect_silent(check_trait_weights(trait_weights, sp_traits))

    ## Missing columns
    # Forgotten single colunms
    expect_error(check_trait_weights(trait_weights[, 1:2], sp_traits),
                 "Column(s) compet_weight is (are) not in trait_weights",
                 fixed = TRUE)
    expect_error(check_trait_weights(trait_weights[, c(1, 3)], sp_traits),
                 "Column(s) growth_weight is (are) not in trait_weights",
                 fixed = TRUE)
    expect_error(check_trait_weights(trait_weights[, 2:3], sp_traits),
                 "Column(s) trait is (are) not in trait_weights",
                 fixed = TRUE)
    # Forgot several columns
    expect_error(check_trait_weights(trait_weights[, 1, drop = FALSE],
                                     sp_traits),
                 paste0("Column(s) growth_weight, compet_weight is (are) not ",
                        "in trait_weights"), fixed = TRUE)

    ## Missing traits in weights df that are in traits df
    expect_error(check_trait_weights(trait_weights, sp_traits[, 1, drop = FALSE]),
                 "Trait(s) trait2 (is/are) not in provided traits data.frame",
                 fixed = TRUE)
})


test_that("compute_compet_distance() works", {
    # Consider distances using a single trait
    expect_equal(compute_compet_distance(trait_weights, sp_traits),
                 matrix(c(0, 1/2, 1/2, 1/2, 0, 1, 1/2, 1, 0), ncol = 3,
                        dimnames = list(paste0("sp", 1:3), paste0("sp", 1:3))))

    # Multitrait distance
    multi_trait_dist = as.matrix(dist(sp_traits))
    multi_trait_dist = multi_trait_dist/(diff(range(multi_trait_dist)))

    expect_equal(compute_compet_distance(data.frame(trait = paste0("trait", 1:2),
                                                    growth_weight = c(1, 0),
                                                    compet_weight = c(0.5, 0.5)),
                                         sp_traits),
                 multi_trait_dist)
    expect_equal(compute_compet_distance(data.frame(trait = paste0("trait", 1:2),
                                                    growth_weight = c(1, 0),
                                                    compet_weight = c(0.3, 0.3)),
                                         sp_traits),
                 multi_trait_dist)
    expect_equal(compute_compet_distance(data.frame(trait = paste0("trait", 1:2),
                                                    growth_weight = c(1, 0),
                                                    compet_weight = c(1, 1)),
                                         sp_traits),
                 multi_trait_dist)
})


test_that("env_curve() works as expected", {
    # Single trait growth
    expect_equal(env_curve(sp_traits[1,, drop = FALSE], 1, trait_weights,
                           1.15, 5),
                 1.15)
    expect_equal(env_curve(sp_traits[1,, drop = FALSE], 1, trait_weights,
                           1.25, 5),
                 1.25)
    expect_equal(env_curve(sp_traits[1,, drop = FALSE], 1, trait_weights,
                           1.25, 10),
                 1.25)
    expect_equal(env_curve(sp_traits[2,, drop = FALSE], 1, trait_weights,
                           1.25, 10),
                 1.25 * exp(-1/200))

    # Single trait growth without any traits defined
    expect_equal(env_curve(sp_traits[2,, drop = FALSE], 1,
                           data.frame(trait = paste0("trait", 1:2),
                                      growth_weight = c(0, 0),
                                      compet_weight = c(1, 1)),
                           1.25, 10),
                 1.25)

    # Multiple trait growth
    expect_equal(env_curve(sp_traits[3,, drop = FALSE], 1,
                           data.frame(trait = paste0("trait", 1:2),
                                      growth_weight = c(0.5, 0.5),
                                      compet_weight = c(1, 1)),
                           1.25, 10),
                 mean(c(1.25 * exp(-1/200), 1.25)))
})
