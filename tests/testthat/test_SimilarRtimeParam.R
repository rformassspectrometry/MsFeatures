test_that("SimillarRtimeParam works", {
    res <- SimilarRtimeParam(4)
    expect_true(res@diffRt == 4)
    expect_equal(res@groupFun, groupClosest)

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

    prm <- SimilarRtimeParam(2, groupFun = groupClosest)
    res <- groupFeatures(x, prm)
    expect_equal(res, factor(c(1, 1, 2, 2, 3, 3, 3, 4, 4)))
})
