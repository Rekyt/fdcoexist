context("test support functions")

# Preamble ---------------------------------------------------------------------

test_that("na_if_null_or_inf works", {

    expect_equal(na_if_null_or_inf(Inf), NA_real_)
    expect_equal(na_if_null_or_inf(-Inf), NA_real_)
    expect_equal(na_if_null_or_inf(NA_real_), NA_real_)
    expect_equal(na_if_null_or_inf(0), 0)

    expect_equal(na_if_null_or_inf(c(1, NA, Inf, -Inf, 2)),
                 c(1, NA_real_, NA_real_, NA_real_, 2))
})

