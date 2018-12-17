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


    # When no trait contribute to competition there is no competition

    expect_equal(compute_compet_distance(data.frame(trait = paste0("trait", 1:2),
                                                    growth_weight = c(1, 0),
                                                    compet_weight = c(0, 0)),
                                         sp_traits),
                 matrix(0, ncol = 3, nrow = 3,
                        dimnames = list(rownames(sp_traits),
                                        rownames(sp_traits))))

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

test_that("Computing trait distance with exponents", {
    expect_equal(compute_compet_distance(trait_weights, sp_traits,
                                         dist_power = 0),
                 matrix(rep(0, 9), ncol = 3,
                        dimnames = list(rownames(sp_traits),
                                        rownames(sp_traits))))

    exp05 = (exp(0.5) - 1) / (exp(1) - 1)

    expect_equal(compute_compet_distance(trait_weights, sp_traits,
                                         dist_power = 1),
                 matrix(c(0, exp05, exp05, exp05, 0, 1, exp05, 1, 0), ncol = 3,
                        dimnames = list(rownames(sp_traits),
                                        rownames(sp_traits))))

    exp025 = (exp(0.25) - 1) / (exp(1) - 1)

    expect_equal(compute_compet_distance(trait_weights, sp_traits,
                                         dist_power = 2),
                 matrix(c(0, exp025, exp025, exp025, 0, 1, exp025, 1, 0),
                        ncol = 3, dimnames = list(rownames(sp_traits),
                                                  rownames(sp_traits))))
})

test_that("bevHoltFct() returns good result", {

    R = runif(1, 1, 2)

    N = rbinom(1, 80, 1/2)

    expect_equal(bevHoltFct(0, 0, 0), 0)
    expect_equal(bevHoltFct(R, N, 0), R * N)
    expect_equal(bevHoltFct(R, 0, 0), 0)
    expect_equal(bevHoltFct(0, N, 0), 0)
    expect_equal(bevHoltFct(R, N, R * N - 1), 1)
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


    # Change environmental filter width
    expect_equal(env_curve(sp_traits[3,, drop = FALSE], 1,
                           data.frame(trait = paste0("trait", 1:2),
                                      growth_weight = c(0.5, 0.5),
                                      compet_weight = c(1, 1)),
                           1.25, 1),
                 mean(c(1.25 * exp(-1/2), 1.25)))

    expect_error(env_curve(sp_traits[2,, drop = FALSE], 1,
                           data.frame(trait = paste0("trait", 1:2),
                                      growth_weight = c(0, 0),
                                      compet_weight = c(1, 1)),
                           1.25, c(10, 5)),
                 "There are either too many or not enough values for width",
                 fixed = TRUE)

    expect_error(env_curve(sp_traits[2,, drop = FALSE], c(1, 2, 3),
                           data.frame(trait = paste0("trait", 1:2),
                                      growth_weight = c(0, 0),
                                      compet_weight = c(1, 1)),
                           1.25, c(10, 5)),
                 "There are either too many or not enough values for width",
                 fixed = TRUE)

    # Add Hierarchical Competition
})

test_that("alphaterm() works as expected", {

    # Initial population
    N0 = matrix(c(10, 5, 10, 5), dimnames = list(c("p1", "p2"), c("s1", "s2")),
                ncol = 2)

    # Species distance
    given_dist = matrix(c(0, 1, 1, 0), ncol = 2,
                        dimnames = list(c("s1", "s2"), c("s1", "s2")))

    expect_equal(alphaterm(given_dist, N0, 1, 0), N0)
    expect_equal(alphaterm(given_dist, N0, 1, 1), N0 * 2)
    expect_equal(alphaterm(given_dist, N0, 1, 1/2), N0 * 3/2)
})
