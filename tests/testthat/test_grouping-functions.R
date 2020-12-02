test_that("groupClosest works", {
    x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7, 9, 9.5, 15)
    res <- groupClosest(x)
    expect_equal(res, c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 6, 7))

    res <- groupClosest(x, maxDiff = 0.3)
    expect_equal(res, c(1, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 11, 12))

    idx <- sample(seq_along(res))
    res_2 <- groupClosest(x[idx], maxDiff = 0.3)
    expect_equal(res[idx], res_2)

    a <- c(4.9, 5.2, 5.4)
    res <- groupClosest(a, maxDiff = 0.3)
    expect_equal(res, c(1, 2, 2))
})

test_that("groupSimilarityMatrix works", {
    x <- rbind(
        c(1, 0.9, 0.6, 0.8, 0.5),
        c(0.9, 1, 0.7, 0.92, 0.8),
        c(0.6, 0.7, 1, 0.91, 0.7),
        c(0.8, 0.92, 0.91, 1, 0.9),
        c(0.5, 0.8, 0.7, 0.9, 1)
    )
    expect_error(groupSimilarityMatrix(x[1:4, ]), "symmetric matrix")

    res <- groupSimilarityMatrix(x, threshold = 0.9)
    expect_equal(res, c(2, 1, 3, 1, 4))

    res <- groupSimilarityMatrix(x, threshold = 0)
    expect_true(all(res == 1))

    ## Add also a correlation between 3 and 2
    x[2, 3] <- 0.9
    x[3, 2] <- 0.9
    res <- groupSimilarityMatrix(x, threshold = 0.9)
    expect_equal(res, c(2, 1, 1, 1, 3))

    ## Add a higher correlation between 4 and 5
    x[4, 5] <- 0.99
    x[5, 4] <- 0.99
    res <- groupSimilarityMatrix(x, threshold = 0.9)
    expect_equal(res, c(2, 2, 3, 1, 1))

    ## Increase correlation between 2 and 3
    x[2, 3] <- 0.92
    x[3, 2] <- 0.92
    res <- groupSimilarityMatrix(x, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1))

    ## 3 and 5 above threshold
    x[3, 5] <- 0.9
    x[5, 3] <- 0.9
    res <- groupSimilarityMatrix(x, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1))

    ## Real data
    load(system.file("extdata/cors.RData", package = "MsFeatures"))
    res <- groupSimilarityMatrix(cors, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1, 4, 1))

    res <- groupSimilarityMatrix(cors, threshold = 0.8)
    expect_equal(res, c(2, 1, 1, 1, 1, 3, 1))

    res <- groupSimilarityMatrix(cors, threshold = 0.7)
    expect_equal(res, c(2, 1, 1, 1, 1, 1, 1))
})
