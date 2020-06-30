test_that(".group_logic_matrix works", {
    xmat <- rbind(c(TRUE, FALSE, FALSE, FALSE),
                  c(FALSE, TRUE, FALSE, FALSE),
                  c(FALSE, FALSE, TRUE, FALSE),
                  c(FALSE, FALSE, FALSE, TRUE))
    expect_error(.group_logic_matrix(xmat[1:3, ]))
    res <- .group_logic_matrix(xmat)
    expect_true(length(res) == nrow(xmat))
    expect_equal(res, list(1, 2, 3, 4))

    xmat <- rbind(c(TRUE, FALSE, FALSE, FALSE, TRUE),
                  c(FALSE, TRUE, FALSE, FALSE, FALSE),
                  c(FALSE, FALSE, TRUE, FALSE, TRUE),
                  c(FALSE, FALSE, FALSE, TRUE, FALSE),
                  c(TRUE, FALSE, TRUE, FALSE, TRUE))
    res <- .group_logic_matrix(xmat)
    expect_equal(res, list(c(1, 3, 5), 2, 4))

    xcor <- matrix(FALSE, ncol = 13, nrow = 13)
    for (i in 1:13)
        xcor[i, i] <- TRUE
    xcor[8, 6] <- TRUE
    xcor[8, 7] <- TRUE
    xcor[9, 7] <- TRUE
    xcor[11, 7] <- TRUE
    xcor[6, 8] <- TRUE
    xcor[7, 8] <- TRUE
    xcor[10, 8] <- TRUE
    xcor[13, 8] <- TRUE
    xcor[7, 9] <- TRUE
    xcor[8, 10] <- TRUE
    xcor[7, 11] <- TRUE
    xcor[12, 11] <- TRUE
    xcor[11, 12] <- TRUE
    xcor[8, 13] <- TRUE
    res <- .group_logic_matrix(xcor)
    expect_equal(res, list(1, 2, 3, 4, 5, c(6:13)))
    
    xcor <- matrix(FALSE, ncol = 10, nrow = 10)
    for (i in seq_len(ncol(xcor))) {
        xcor[i, i] <- TRUE
    }
    xcor[1, 4] <- TRUE
    xcor[4, 1] <- TRUE
    xcor[2, 8] <- TRUE
    xcor[8, 2] <- TRUE
    xcor[3, 9] <- TRUE
    xcor[9, 3] <- TRUE
    xcor[8, 9] <- TRUE
    xcor[9, 8] <- TRUE
    res <- .group_logic_matrix(xcor)
    expect_equal(res, list(c(1, 4), c(2, 3, 8, 9), 5, 6, 7, 10))
})

test_that(".index_list_to_factor works", {
    x <- list(c(1, 5, 2), c(3, 4), c(6), 7)
    res <- .index_list_to_factor(x)
    expect_equal(res, factor(c(1, 1, 2, 2, 1, 3, 4)))
})

test_that("groupByCorrelation works", {
    x <- rbind(c(1, 2, 3, 4),
               c(2, 4, 6, 8),
               c(0, 2, 1, 2),
               c(1, 3, 4, 5))
    res <- groupByCorrelation(x)
    expect_true(is.factor(res))
    expect_equal(res, factor(c(1, 1, 2, 1)))

    res_2 <- groupByCorrelation(x, greedy = TRUE)
    expect_equal(res, res_2)
    
    expect_error(groupByCorrelation(x, threshold = c(0.4, 0.3)), "length 1")

    x <- rbind(x,
               c(2, 4, 6, 9),
               c(1, 4, 1, 4),
               c(1, 2, 3, 4))
    f <- c(1, 2, 2, 1, 1, 2, 2)
    res <- groupByCorrelation(x, f = f)
    expect_equal(res, factor(c("1.1", "2.1", "2.2", "1.1", "1.1", "2.2", "2.1")))

    f <- c(1, 2, NA, NA, 1, 2, 2)
    res <- groupByCorrelation(x, f = f)
    expect_equal(res, factor(c("1.1", "2.1", NA, NA, "1.1", "2.2", "2.1")))
    
    expect_error(groupByCorrelation(x, f = 3), "its length has to ")
})

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

test_that(".group_correlation_matrix works", {
    x <- rbind(
        c(1, 0.9, 0.6, 0.8, 0.5),
        c(0.9, 1, 0.7, 0.92, 0.8),
        c(0.6, 0.7, 1, 0.91, 0.7),
        c(0.8, 0.92, 0.91, 1, 0.9),
        c(0.5, 0.8, 0.7, 0.9, 1)
    )
    expect_error(.group_correlation_matrix(x[1:4, ]), "symmetric matrix")
    
    res <- .group_correlation_matrix(x, threshold = 0.9)
    expect_equal(res, c(2, 1, 3, 1, 4))   

    res <- .group_correlation_matrix(x, threshold = 0)
    expect_true(all(res == 1))
    
    ## Add also a correlation between 3 and 2
    x[2, 3] <- 0.9
    x[3, 2] <- 0.9
    res <- .group_correlation_matrix(x, threshold = 0.9)
    expect_equal(res, c(2, 1, 1, 1, 3))
    
    ## Add a higher correlation between 4 and 5
    x[4, 5] <- 0.99
    x[5, 4] <- 0.99
    res <- .group_correlation_matrix(x, threshold = 0.9)
    expect_equal(res, c(2, 2, 3, 1, 1))
    
    ## Increase correlation between 2 and 3
    x[2, 3] <- 0.92
    x[3, 2] <- 0.92
    res <- .group_correlation_matrix(x, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1))

    ## 3 and 5 above threshold
    x[3, 5] <- 0.9
    x[5, 3] <- 0.9
    res <- .group_correlation_matrix(x, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1))

    ## Real data
    load(system.file("extdata/cors.RData", package = "MsFeatures"))
    res <- .group_correlation_matrix(cors, threshold = 0.9)
    expect_equal(res, c(3, 2, 2, 1, 1, 4, 1))

    res <- .group_correlation_matrix(cors, threshold = 0.8)
    expect_equal(res, c(2, 1, 1, 1, 1, 3, 1))

    res <- .group_correlation_matrix(cors, threshold = 0.7)
    expect_equal(res, c(2, 1, 1, 1, 1, 1, 1))
})

test_that(".format_groups works", {
    res <- .format_groups(1:3)
    expect_equal(res, c("1", "2", "3"))
    res <- .format_groups(1:10)
    expect_equal(res, c("01", "02", "03", "04", "05", "06", "07", "08",
                        "09", "10"))
})
