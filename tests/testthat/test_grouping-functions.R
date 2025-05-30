test_that("groupConsecutive works", {
    x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7, 9,
           9.5, 15)
    res <- groupConsecutive(x)
    expect_equal(res, c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 6, 6, 7))

    res <- groupConsecutive(x, maxDiff = 0.3)
    expect_equal(res, c(1, 2, 2, 3, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 11, 12))

    idx <- sample(seq_along(res))
    res_2 <- groupConsecutive(x[idx], maxDiff = 0.3)
    expect_equal(res[idx], res_2)

    a <- c(4.9, 5.2, 5.4)
    res <- groupConsecutive(a, maxDiff = 0.3)
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

test_that("groupClosest works", {
    data(se)
    rts <- rowData(se)$rtmed
    ## Check that all differences within each group are below threshold!
    res <- groupClosest(rts, 10, FUN = groupSimilarityMatrixTree)
    resl <- split(rts, res)
    comp_fun <- function(z) {
        if (length(z) == 1)
            TRUE
        else all(diff(z) <= 10)
    }
    expect_true(all(vapply(resl, comp_fun, logical(1))))

    res_2 <- groupClosest(rts, 10, FUN = groupSimilarityMatrix)
    resl <- split(rts, res_2)
    comp_fun <- function(z) {
        if (length(z) == 1)
            TRUE
        else all(diff(z) <= 10)
    }
    expect_true(all(vapply(resl, comp_fun, logical(1))))

    res <- groupClosest(rts, 60, FUN = groupSimilarityMatrixTree)
    resl <- split(rts, res)
    comp_fun <- function(z) {
        if (length(z) == 1)
            TRUE
        else all(diff(z) <= 60)
    }
    expect_true(all(vapply(resl, comp_fun, logical(1))))

    res_2 <- groupClosest(rts, 60, FUN = groupSimilarityMatrix)
    resl <- split(rts, res_2)
    comp_fun <- function(z) {
        if (length(z) == 1)
            TRUE
        else all(diff(z) <= 60)
    }
    expect_true(all(vapply(resl, comp_fun, logical(1))))
})


test_that("groupSimilarityMatrixTree works", {
  data(se)
  rts <- rowData(se)$rtmed

  # Make distances
  dists <- dist(rts, method = "manhattan")

  # Check if all groups are clustered if maxdiff is equal to the maximum distance
  res <- groupSimilarityMatrixTree(dists, maxDiff = max(dists))
  expect_true(length(unique(res)) == 1)

  # Some features in the rts have the exact same rt. Check if these are clustered.
  res <- groupSimilarityMatrixTree(dists, maxDiff = min(dists))
  zero_distance_groups <- names(table(res)[table(res)>1])
  expect_true(all(vapply(zero_distance_groups, function(x) {
    all(as.matrix(dists)[res == x, res == x] == 0)},
    logical(1))))


  # Check if groups are clustered correctly
  rts <- c(1, 1.5, 2,   # group 1
           3, 3.1, 3.2, # group 2
           6, 7, 8      # group 3
  )
  dists <- dist(rts, method = "manhattan")
  res <- groupSimilarityMatrixTree(dists, maxDiff = 2)
  expect_equal(as.factor(as.character(res)), as.factor(c(1, 1, 1, 4, 4, 4, 7, 7, 7)))

  # Ensure that it is a one dimensional vector (and that euclidean and manhattan gives the same result)
  res1 <- groupSimilarityMatrixTree(dist(rts, method = "manhattan"), maxDiff = 2)
  res2 <- groupSimilarityMatrixTree(dist(rts, method = "euclidean"), maxDiff = 2)
  expect_equal(res1, res2)
})