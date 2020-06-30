test_that("SimillarRtimeParam works", {
    res <- SimilarRtimeParam(4)
    expect_true(res@diffRt == 4)

    expect_error(SimilarRtimeParam(1:2), "positive numeric")
    expect_error(SimilarRtimeParam(-1), "positive numeric")

    expect_error(SimilarRtimeParam(method = "some"), "should be one")
    expect_error(new("SimilarRtimeParam", method = c("groupClosest", "greedy")),
                 "should be either")
})

test_that(".group_similar_rtime works", {
    rts <- c(4, 3, 9, 12, 10, 11, 3, 20, 21, 25)
    res <- .group_similar_rtime(rts, tol = 1)
    expect_equal(res, paste0("FG.0", c(1, 1, 2, 3, 2, 3, 1, 4, 4, 5)))
    res <- .group_similar_rtime(rts, tol = 1, method = "greedy")
    expect_equal(res, paste0("FG.0", c(1, 1, 2, 2, 2, 2, 1, 3, 3, 4)))

    f <- c(NA, 1, 1, 1, NA, 2, 2, 2, 2, NA)
    res <- .group_similar_rtime(rts, f = f, tol = 1)
    expect_true(all(is.na(res[c(1, 5, 10)])))
    expect_equal(res, c(NA, "1.1", "1.2", "1.3", NA, "2.2", "2.1",
                        "2.3", "2.3", NA))
})

test_that("groupFeatures,numeric,SimilarRtimeParam works", {
    rts <- c(4, 3, 45, 6, 4, 23, 24, 25)
    res <- groupFeatures(rts, SimilarRtimeParam(3))
    expect_equal(res, c("FG.1", "FG.1", "FG.3", "FG.1", "FG.1", "FG.2", "FG.2", "FG.2"))
})

test_that("groupFeatures,XCMSnExp,SimilarRtimeParam works", {
    prm <- SimilarRtimeParam(3)
    xod <- xcms::dropFeatureDefinitions(xdata)

    expect_error(groupFeatures(xod, prm), "No feature definitions")
    res <- groupFeatures(xdata, prm)
    expect_true(any(colnames(featureDefinitions(res)) == "feature_group"))

    res <- groupFeatures(res, prm)
    expect_equal(featureDefinitions(res)$feature_group[1], "FG.025.1")

    tmp <- xdata
    featureDefinitions(tmp)$ms_level[c(1:3, 5)] <- 2
    res_2 <- groupFeatures(tmp, prm)
    fgs_2 <- featureDefinitions(res_2)$feature_group
    expect_true(all(is.na(fgs_2[c(1:3, 5)])))

    expect_error(groupFeatures(res, prm, msLevel = 1:2), "Currently only")
})

test_that("groupFeatures,SummarizedExperiment,SimilarRtimeParam works", {
    se <- xcms::quantify(xdata)
    res <- groupFeatures(se, param = SimilarRtimeParam(3), column = "rtmed")
    expect_true(any(colnames(SummarizedExperiment::rowData(res)) ==
                    "feature_group"))
    res_2 <- se
    SummarizedExperiment::rowData(res_2)$feature_group <- "FG"
    res_2 <- groupFeatures(res_2, param = SimilarRtimeParam(3), column = "rtmed")
    expect_equal(SummarizedExperiment::rowData(res),
                 SummarizedExperiment::rowData(res_2))
    res_2 <- se
    SummarizedExperiment::rowData(res_2)$feature_group <- "FG"
    SummarizedExperiment::rowData(res_2)$feature_group[c(1, 4, 5)] <- NA
    res <- groupFeatures(res_2, param = SimilarRtimeParam(3), column = "rtmed")
    expect_true(all(is.na(SummarizedExperiment::rowData(res)$feature_group[c(1, 4, 5)])))
})
