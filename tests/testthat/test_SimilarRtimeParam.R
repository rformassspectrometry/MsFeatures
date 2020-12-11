test_that("SimillarRtimeParam works", {
    res <- SimilarRtimeParam(4)
    expect_true(res@diffRt == 4)
    expect_equal(res@groupFun, groupSorted)

    expect_error(SimilarRtimeParam(1:2), "positive numeric")
    expect_error(SimilarRtimeParam(-1), "positive numeric")
    expect_error(SimilarRtimeParam(groupFun = 4), "function")

    prm <- SimilarRtimeParam(3)
})

test_that("groupFeatures,SimilarRtimeParam,numeric works", {
    prm <- SimilarRtimeParam(1, groupFun = MsCoreUtils::group)

    x <- c(2, 3, 4, 5, 10, 11, 12, 14, 15)
    res <- groupFeatures(x, prm)
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1, 1, 1, 1, 2, 2, 2, 3, 3)))

    prm <- SimilarRtimeParam(2, groupFun = MsCoreUtils::group)
    res <- groupFeatures(x, prm)
    expect_equal(res, factor(c(1, 1, 1, 1, 2, 2, 2, 2, 2)))

    prm <- SimilarRtimeParam(2, groupFun = groupSorted)
    res <- groupFeatures(x, prm)
    expect_equal(res, factor(c(1, 1, 2, 2, 3, 3, 3, 4, 4)))
})

test_that(".format_id works", {
    vals <- 1:4
    res <- .format_id(vals)
    expect_equal(res, c("001", "002", "003", "004"))
    res <- .format_id(1:1000)
    expect_equal(res[1], "0001")
})

test_that("groupFeatures,SummarizedExperiment,SimilarRtimeParam works", {
    data(se)
    prm <- SimilarRtimeParam(10, groupFun = groupSorted)
    rts <- rowData(se)$rtmed
    res <- groupFeatures(se, prm, rtime = "rtmed")
    expect_true(!any(is.na(featureGroups(res))))

    tmp <- split(rts, featureGroups(res))
    test_fun <- function(z) {
        if (length(z))
            all(diff(sort(z)) < 10)
        else TRUE
    }
    expect_true(all(vapply(tmp, test_fun, logical(1))))

    ## Pre-defined features.
    fgs <- c(rep("1", 5), rep("2", 20), rep("4", 10), rep("2", 20),
             rep("1", 20), rep("3", 30), rep("4", 120))
    featureGroups(se) <- fgs
    res <- groupFeatures(se, prm, rtime = "rtmed")
    tmp <- strsplit(featureGroups(res), split = ".", fixed = TRUE)
    expect_equal(fgs, vapply(tmp, function(z) z[1], character(1)))
    tmp <- split(rts, featureGroups(res))
    expect_true(all(vapply(tmp, test_fun, logical(1))))

    ## How does this work with NAs?
    fgs[c(5, 14, 67)] <- NA
    featureGroups(se) <- fgs
    res <- groupFeatures(se, param = prm, rtime = "rtmed")
    expect_true(all(is.na(featureGroups(res)[c(5, 14, 67)])))
    tmp <- split(rts, featureGroups(res))
    expect_true(all(vapply(tmp, test_fun, logical(1))))

    expect_error(groupFeatures(se, prm), "numeric values")
    expect_error(groupFeatures(se, prm, rtime = "peakidx"), "numeric values")

    ## Same with groupClosest
    rowData(se)$feature_group <- NULL
    prm <- SimilarRtimeParam(10, groupFun = groupClosest)
    rts <- rowData(se)$rtmed
    res <- groupFeatures(se, prm, rtime = "rtmed")
    expect_true(!any(is.na(featureGroups(res))))

    tmp <- split(rts, featureGroups(res))
    test_fun <- function(z) {
        if (length(z))
            all(diff(z) < 10)
        else TRUE
    }
    expect_true(all(vapply(tmp, test_fun, logical(1))))
})
