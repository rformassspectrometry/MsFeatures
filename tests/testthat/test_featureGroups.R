test_that("featureGroups,SummarizedExperiment works", {
    data(se)
    res <- featureGroups(se)
    expect_true(is.character(res))
    expect_true(all(is.na(res)))

    tmp <- se
    rowData(tmp)$feature_group <- seq_len(nrow(rowData(tmp)))
    expect_true(is.character(featureGroups(tmp)))
    expect_equal(featureGroups(tmp), as.character(seq_len(nrow(rowData(tmp)))))
})

test_that("featureGroups<-,SummarizedExperiment works", {
    data(se)
    res <- featureGroups(se)
    expect_true(is.character(res))
    expect_true(all(is.na(res)))

    featureGroups(se) <- 1
    expect_true(is.character(rowData(se)$feature_group))

    expect_error(featureGroups(se) <- 1:2, "length")
})
