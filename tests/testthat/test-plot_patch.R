context("test plot_patch()")

test_that("plot_patch() works", {
    N0 = matrix(c(50, 50, 10, 10), ncol = 2,
                dimnames = list(site    = paste0("p", 1:2),
                                species = paste0("s", 1:2)))

    given_compo = simplify2array(lapply(seq(1, 0.25, length.out = 5),
                                        function(x) x * N0))

    expect_silent(plot_patch(given_compo, "p1", 5))

    p1_plot = plot_patch(given_compo, "p1", 5)

    expect_is(p1_plot, "ggplot")
    expect_equal(p1_plot$labels$title, "Patch p1")
    expect_equal(p1_plot$labels$x, "Time")
    expect_equal(p1_plot$labels$y, "Number of individuals")
    expect_is(p1_plot$layers[[1]]$geom, "GeomLine")
})
