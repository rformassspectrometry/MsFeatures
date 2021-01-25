test_that("AbundanceSimilarityParam works", {
    res <- AbundanceSimilarityParam(threshold = 1, other_stuff = "b",
                                    method = "c")
    expect_equal(res@threshold, 1)
    expect_equal(res@dots, list(other_stuff = "b", method = "c"))
    expect_error(AbundanceSimilarityParam(threshold = 1:2), "length 1")
})

test_that("groupFeatures,AbundanceSimilarityParam,matrix works", {
    x <- rbind(
        c(12, 34, 231, 234, 9, 5, 7),
        c(900, 900, 800, 10, 12, 9, 4),
        c(25, 70, 400, 409, 15, 8, 4),
        c(12, 13, 14, 15, 16, 17, 18),
        c(14, 36, 240, 239, 12, 7, 8),
        c(100, 103, 80, 2, 3, 1, 1)
    )

    res <- groupFeatures(x, AbundanceSimilarityParam())
    expect_equal(res, c(1, 2, 1, 3, 1, 2))
    res <- groupFeatures(x, AbundanceSimilarityParam(threshold = 0.95))
    expect_equal(res, c(1, 2, 1, 3, 1, 2))

    res <- groupFeatures(x, AbundanceSimilarityParam(method = "kendal",
                                                     threshold = 0.94))
    expect_equal(res, c(2, 1, 3, 4, 5, 1))

    ## Doesn't fail with additional parameter...
    res <- groupFeatures(x, AbundanceSimilarityParam(method = "kendal",
                                                     threshold = 0.94,
                                                     other_opts = 4))
    expect_equal(res, c(2, 1, 3, 4, 5, 1))

    expect_error(groupFeatures(x, AbundanceSimilarityParam(simFun = cor)),
                 "symmetric numeric matrix")

    expect_error(groupFeatures(x, AbundanceSimilarityParam(subset = 29)),
                 "out of bounds")

    res <- groupFeatures(x, AbundanceSimilarityParam(subset = c(1, 2, 6, 7)))
    expect_equal(res, c(2, 1, 2, 3, 2, 1))

    x <- matrix("a", nrow = 3, ncol = 2)
    expect_error(groupFeatures(x, AbundanceSimilarityParam(simFun = cor)),
                 "numeric matrix")

})

test_that("groupFeatures,AbundanceSimilarityParam,SummarizedExperiment works", {
    x <- rbind(
        c(12, 34, 231, 234, 9, 5, 7),
        c(900, 900, 800, 10, 12, 9, 4),
        c(25, 70, 400, 409, 15, 8, 4),
        c(12, 13, 14, 15, 16, 17, 18),
        c(14, 36, 240, 239, 12, 7, 8),
        c(100, 103, 80, 2, 3, 1, 1)
    )

    tmp <- SummarizedExperiment::SummarizedExperiment(assays = x)
    res <- groupFeatures(tmp, AbundanceSimilarityParam())
    expect_true(is(res, "SummarizedExperiment"))
    expect_equal(featureGroups(res), c("FG.001", "FG.002", "FG.001", "FG.003",
                                       "FG.001", "FG.002"))

    res <- groupFeatures(tmp, AbundanceSimilarityParam(method = "kendal",
                                                       threshold = 0.94))
    expect_equal(featureGroups(res), c("FG.002", "FG.001", "FG.003", "FG.004",
                                       "FG.005", "FG.001"))

    expect_error(groupFeatures(tmp, AbundanceSimilarityParam(), i = c(1, 2)),
                 "more than")

    featureGroups(tmp) <- c("1", "1", "2", "2", "1", "1")
    res <- groupFeatures(tmp, AbundanceSimilarityParam())
    expect_equal(featureGroups(res), c("1.001", "1.002", "2.001", "2.002",
                                       "1.001", "1.002"))

    featureGroups(tmp) <- c(NA, NA, 1, 1, 1, NA)
    res <- groupFeatures(tmp, AbundanceSimilarityParam())
    expect_equal(featureGroups(res), c(NA, NA, "1.001", "1.002", "1.001", NA))
})