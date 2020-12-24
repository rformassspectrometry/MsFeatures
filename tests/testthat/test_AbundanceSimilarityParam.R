test_that("AbundanceSimilarityParam works", {
    res <- AbundanceSimilarityParam(threshold = 1, other_stuff = "b",
                                    method = "c")
    expect_equal(res@threshold, 1)
    expect_equal(res@dots, list(other_stuff = "b", method = "c"))
    expect_error(AbundanceSimilarityParam(threshold = 1:2), "length 1")
})
