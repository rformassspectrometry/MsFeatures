test_that("corRows works", {
    x <- rbind(
        c(12, 34, 231, 234, 9, 5, 7),
        c(900, 900, 800, 10, 12, 9, 4),
        c(25, 70, 400, 409, 15, 8, 4),
        c(12, 13, 14, 15, 16, 17, 18),
        c(14, 36, 240, 239, 12, 7, 8),
        c(100, 103, 80, 2, 3, 1, 1)
    )

    res <- corRows(x)
    expect_true(nrow(res) == nrow(x))

    ## Support additional params
    res <- corRows(x, other_opt = 4)
    expect_true(nrow(res) == nrow(x))
})
