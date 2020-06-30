#' @description
#'
#' Utility function to group rows/columns in an n x n logical matrix that are
#' `TRUE`. Note that this groups also rows/columns together that are not
#' direclty *linked* with a `TRUE`, but that are linked via other rows/columns
#' they have in common (i.e. if between row 2 and 4 is a `TRUE`, but also
#' between 3 and 4, all 3 of them are joined together, even if they are not
#' directly linked).
#'
#' @param x `logical` `matrix` with same number of rows and columns. See below
#'      for examples.
#'
#' @return `list` of `integer` indexes of rows (columns) that are grouped.
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Find groups of rows with a correlation higher than 0.9
#' x <- rbind(
#'     c(1, 3, 2, 5),
#'     c(2, 6, 4, 7),
#'     c(1, 1, 3, 1),
#'     c(1, 3, 3, 6),
#'     c(0, 4, 3, 1))
#' xcor <- cor(t(x))
#'
#' .group_logic_matrix(xcor > 0.9)
#'
#' xcor <- matrix(FALSE, ncol = 13, nrow = 13)
#' for (i in 1:13)
#'     xcor[i, i] <- TRUE
#' xcor[8, 6] <- TRUE
#' xcor[8, 7] <- TRUE
#' xcor[9, 7] <- TRUE
#' xcor[11, 7] <- TRUE
#' xcor[6, 8] <- TRUE
#' xcor[7, 8] <- TRUE
#' xcor[10, 8] <- TRUE
#' xcor[13, 8] <- TRUE
#' xcor[7, 9] <- TRUE
#' xcor[8, 10] <- TRUE
#' xcor[7, 11] <- TRUE
#' xcor[12, 11] <- TRUE
#' xcor[11, 12] <- TRUE
#' xcor[8, 13] <- TRUE
#' .group_logic_matrix(xcor)
#'
#' xcor <- matrix(FALSE, ncol = 10, nrow = 10)
#' for (i in seq_len(ncol(xcor))) {
#'     xcor[i, i] <- TRUE
#' }
#' xcor[1, 4] <- TRUE
#' xcor[4, 1] <- TRUE
#' xcor[2, 8] <- TRUE
#' xcor[8, 2] <- TRUE
#' xcor[3, 9] <- TRUE
#' xcor[9, 3] <- TRUE
#' xcor[8, 9] <- TRUE
#' xcor[9, 8] <- TRUE
#' .group_logic_matrix(xcor)
.group_logic_matrix <- function(x) {
    x <- unname(x)
    nr <- nrow(x)
    if (nr != ncol(x))
        stop("'x' is supposed to be a symmetric matrix")
    groups <- vector("list", nr)
    idx <- which(x, arr.ind = TRUE)
    is_in_group <- rep(FALSE, nr)
    for (i in seq_len(nr)) {
        if (is_in_group[i])
            next
        elms <- idx[idx[, 1] == i, 2]
        if (length(elms) == 1) {
            groups[[i]] <- i
            is_in_group[i] <- TRUE
        } else {
            while(!all(is_in_group[elms])) {
                is_in_group[elms] <- TRUE
                ## Get all rows containing these elements
                elms <- unique(idx[idx[, 1] %in% elms, 2])
            }
            groups[[i]] <- elms
        }
    }
    groups[lengths(groups) > 0]
}

#' @note this expects a `list` of `integer` were each index is only present
#'     **once** and in addition no indices should be missing!
#'
#' @noRd
.index_list_to_factor <- function(x) {
    len <- sum(lengths(x))
    res <- integer(len)
    for (i in seq_along(x))
        res[x[[i]]] <- i
    as.factor(res)
}

#' @description
#'
#' This function creates clusters of highest correlations between elements in
#' which **all** elements have a correlation `>= threshold` with each other.
#' This is in contrast to
#' `.group_logic_matrix(x >= threshold)` that creates clusters of elements that
#' have a correlation `>= threshold` to any other (i.e. at least one other)
#' element of a cluster. Thus, this function creates smaller clusters of
#' highly correlated elements. Note however that with this approach single
#' elements in one cluster could also have a correlation `>= threshold` to
#' another element in another cluster. The average correlation to its own
#' cluster will however be larger to that of the other.
#'
#' The algorithm is defined as follows:
#' - all pairs of values in `x` which are `>= threshold` are identified and
#'   sorted decreasingly.
#' - starting with the pair with the highest correlation groups are defined:
#' - if none of the two is in a group, both are put into the same new group.
#' - if one of the two is already in a group, the other is put into the same
#'   group if **all** correlations of it to that group are `>= threshold`
#'   (and are not `NA`).
#' - if both are already in the same group nothing is done.
#' - if both are in different groups: an element is put into the group of the
#'   other if a) all correlations of it to members of the other's group 
#'   are not `NA` and `>= threshold` **and** b) the average correlation to the
#'   other group is larger than the average correlation to its own group.
#'
#' This ensures that groups are defined in which all elements have a correlation
#' `>= threshold` with each other and the correlation between members of the
#' same group is maximized.
#' 
#' @param x symmetrix `numeric` `matrix`.
#'
#' @param threshold `numeric(1)` above which rows in `x` should be grouped.
#'
#' @param full `logical(1)` whether the full matrix should be considered, or
#'     just the upper triangular matrix (including the diagonal).
#'
#' @noRd
#' 
#' @examples
#' 
#' x <- rbind(
#'     c(1, 0.9, 0.6, 0.8, 0.5),
#'     c(0.9, 1, 0.7, 0.92, 0.8),
#'     c(0.6, 0.7, 1, 0.91, 0.7),
#'     c(0.8, 0.92, 0.91, 1, 0.9),
#'     c(0.5, 0.8, 0.7, 0.9, 1)
#'     )
#'
#' .group_correlation_matrix(x, threshold = 0.9)
#'
#' .group_correlation_matrix(x, threshold = 0.1)
#' 
#' ## In contrast to the "greedy" version that groups all elements together
#' ## that have a correlation >= 0.9 with any other element in a cluster.
#' CompMetaboTools:::.group_logic_matrix(x >= 0.9)
#'
#' ## Add also a correlation between 3 and 2
#' x[2, 3] <- 0.9
#' x[3, 2] <- 0.9
#' x
#' .group_correlation_matrix(x, threshold = 0.9)
#'
#' ## Add a higher correlation between 4 and 5
#' x[4, 5] <- 0.99
#' x[5, 4] <- 0.99
#' x
#' .group_correlation_matrix(x, threshold = 0.9)
#'
#' ## Increase correlation between 2 and 3
#' x[2, 3] <- 0.92
#' x[3, 2] <- 0.92
#' x
#' .group_correlation_matrix(x, threshold = 0.9) ## Don't break previous cluster!
.group_correlation_matrix <- function(x, threshold = 0.9, full = TRUE) {
    if (!full)
        x[lower.tri(x)] <- NA
    nr <- nrow(x)
    if (nr != ncol(x))
        stop("'x' should be a symmetric matrix")
    res <- rep(NA_integer_, nr)
    x[cbind(1:nr, 1:nr)] <- NA
    idx_pairs <- which(x >= threshold, arr.ind = TRUE)
    idx_pairs <- idx_pairs[order(x[idx_pairs], decreasing = TRUE), , drop = FALSE]
    grp_id <- 1
    the_other <- c(2, 1)
    for (i in seq_len(nrow(idx_pairs))) {
        got_grp <- res[idx_pairs[i, ]]
        nas <- is.na(got_grp)
        if (any(nas)) {
            ## at least one of them is not in a group
            if (sum(nas) == 2) {
                ## none of the two is in a group yet. Check:
                ## correlation above threshold for both with any other
                ## complete group?
                grps <- unique(res[!is.na(res)])
                mean_cor_to_grp <- integer(length(grps))
                names(mean_cor_to_grp) <- grps
                idx <- idx_pairs[i, ]
                for (grp in grps) {
                    idx_grp <- which(res == grp)
                    cor_to_grp <- x[idx, idx_grp]
                    if (full)
                        cor_to_grp <- c(cor_to_grp, x[idx_grp, idx])
                    if (!(any(is.na(cor_to_grp)) || any(cor_to_grp < threshold)))
                        mean_cor_to_grp[grp] <- mean(cor_to_grp)
                }
                mean_cor_to_grp <- mean_cor_to_grp[mean_cor_to_grp > 0]
                if (length(mean_cor_to_grp)) {
                    ## yes: put them into the group with which both have the
                    ## highest correlation.
                    res[idx] <- as.integer(
                        names(sort(mean_cor_to_grp, decreasing = TRUE)))
                } else {
                    ## no: add them as new group
                    res[idx] <- grp_id
                    grp_id <- grp_id + 1
                }
            } else {
                ## One is not in a group. Put that into the group of the
                ## other if a) no cor is NA and b) cor to all of the group
                ## are >= threshold.
                idx <- idx_pairs[i, nas]
                idx_grp <- which(res == got_grp[!nas])
                cor_to_grp <- x[idx, idx_grp]
                if (full)
                    cor_to_grp <- c(cor_to_grp, x[idx_grp, idx])
                if (!(any(is.na(cor_to_grp)) || any(cor_to_grp < threshold)))
                    res[idx] <- got_grp[!nas]
            }
        } else {
            ## both are in a group
            if (length(unique(got_grp)) > 1) {
                grp_1 <- which(res == got_grp[1])
                grp_2 <- which(res == got_grp[2])
                cor_1_1 <- x[idx_pairs[i, 1], grp_1]
                cor_1_2 <- x[idx_pairs[i, 1], grp_2]
                cor_2_1 <- x[idx_pairs[i, 2], grp_1]
                cor_2_2 <- x[idx_pairs[i, 2], grp_2]
                if (full) {
                    cor_1_1 <- c(cor_1_1, x[grp_1, idx_pairs[i, 1]])
                    cor_1_2 <- c(cor_1_2, x[grp_2, idx_pairs[i, 1]])
                    cor_2_1 <- c(cor_2_1, x[grp_1, idx_pairs[i, 2]])
                    cor_2_2 <- c(cor_2_2, x[grp_2, idx_pairs[i, 2]])
                }
                mcor_1_1 <- mean(cor_1_1, na.rm = TRUE)
                mcor_1_2 <- mean(cor_1_2)
                mcor_2_1 <- mean(cor_2_1)
                mcor_2_2 <- mean(cor_2_2, na.rm =TRUE)
                ## Put the elements into the group of the other, if its
                ## correlation to its group is larger than to its own group
                if (is.finite(mcor_1_2) && is.finite(mcor_1_1) &&
                    !any(cor_1_2 < threshold) && mcor_1_2 > mcor_1_1)
                    res[idx_pairs[i, 1]] <- got_grp[2]
                if (is.finite(mcor_2_1) && is.finite(mcor_2_2) &&
                    !any(cor_2_1 < threshold) && mcor_2_1 > mcor_2_2)
                    res[idx_pairs[i, 2]] <- got_grp[1]
                ## if (!(is.na(mcor_1_2) || any(cor_1_2 < threshold))
                ##     && is.finite(mcor_1_1) && mcor_1_2 > mcor_1_1)
                ##     res[idx_pairs[i, 1]] <- got_grp[2]
                ## if (!(is.na(mcor_2_1) || any(cor_2_1 < threshold))
                ##     && is.finite(mcor_2_2) && mcor_2_1 > mcor_2_2)
                ##     res[idx_pairs[i, 2]] <- got_grp[1]
            } # else nothing to do - they are already in the same group
        }
    }
    nas <- is.na(res)
    if (any(nas))
        res[nas] <- seq(grp_id, length.out = sum(nas))
    res
}

#' @title Group rows in a matrix based on their correlation
#'
#' @description
#'
#' The `groupByCorrelation` allows to group rows in a numeric matrix based on
#' their correlation with each other.
#'
#' Two types of groupings are available:
#' 
#' - `greedy = FALSE` (the default): the algorithm creates small groups of
#'   highly correlated members, all of which have a correlation with each other
#'   that are `>= threshold`. Note that with this algorithm, rows in `x` could
#'   still have a correlation `>= threshold` with one or more elements of a
#'   group they are not part of. See notes below for more information.
#' - `greedy = TRUE`: the algorithm creates large groups containing rows that
#'   have a correlation `>= threshold` with at least one element of that group.
#'   For example, if row 1 and 3 have a correlation above the threshold and
#'   rows 3 and 5 too (but correlation between 1 and 5 is below the threshold)
#'   all 3 are grouped into the same group (i.e. rows 1, 3 **and** 5).
#'
#' Note that with parameter `f` it is also possible to pre-define groups of
#' rows that should be further sub-grouped based on correlation with each other.
#' In other words, if `f` is provided, correlations are calculated only between
#' rows with the same value in `f` and hence these pre-defined groups of rows
#' are further sub-grouped based on pairwise correlation. The returned `factor`
#' is then `f` with the additional subgroup appended (and separated with a
#' `"."`). See examples below.
#'
#' @note
#'
#' Implementation note of the grouping algorithm:
#'
#' - all correlations between rows in `x` which are `>= threshold` are
#'   identified and sorted decreasingly.
#' - starting with the pair with the highest correlation groups are defined:
#' - if none of the two is in a group, both are put into the same new group.
#' - if one of the two is already in a group, the other is put into the same
#'   group if **all** correlations of it to that group are `>= threshold`
#'   (and are not `NA`).
#' - if both are already in the same group nothing is done.
#' - if both are in different groups: an element is put into the group of the
#'   other if a) all correlations of it to members of the other's group 
#'   are not `NA` and `>= threshold` **and** b) the average correlation to the
#'   other group is larger than the average correlation to its own group.
#'
#' This ensures that groups are defined in which all elements have a correlation
#' `>= threshold` with each other and the correlation between members of the
#' same group is maximized.
#' 
#' @param x `numeric` `matrix` where rows should be grouped based on
#'     correlation of their values across columns being larger than `threshold`.
#'
#' @param method `character(1)` with the method to be used for correlation. See
#'     [corr()] for options.
#'
#' @param use `character(1)` defining which values should be used for the
#'     correlation. See [corr()] for details.
#'
#' @param threshold `numeric(1)` defining the cut of value above which
#'     rows are considered to be correlated and hence grouped.
#'
#' @param f optional vector of length equal to `nrow(x)` pre-defining groups
#'     of rows in `x` that should be further sub-grouped. See description for
#'     details.
#'
#' @param greedy `logical(1)` whether a version of the grouping algorithm should
#'     be used that leads to larger, more loosely correlated, groups. The
#'     default is `greedy = FALSE`. See description for more information.
#' 
#' @return `factor` with same length than `nrow(x)` with the group each row
#'     is assigned to.
#'
#' @author Johannes Rainer
#'
#' @importFrom stats cor
#' 
#' @family grouping operations
#'
#' @export
#'
#' @examples
#'
#' x <- rbind(
#'     c(1, 3, 2, 5),
#'     c(2, 6, 4, 7),
#'     c(1, 1, 3, 1),
#'     c(1, 3, 3, 6),
#'     c(0, 4, 3, 1),
#'     c(1, 4, 2, 6),
#'     c(2, 8, 2, 12))
#'
#' ## define which rows have a high correlation with each other
#' groupByCorrelation(x)
#'
#' ## assuming we have some prior grouping of rows, further sub-group them
#' ## based on pairwise correlation.
#' f <- c(1, 2, 2, 1, 1, 2, 2)
#' groupByCorrelation(x, f = f)
groupByCorrelation <- function(x, method = "pearson",
                               use = "pairwise.complete.obs",
                               threshold = 0.9, f = NULL, greedy = FALSE) {
    if (length(threshold) > 1)
        stop("'threshold' has to be of length 1")
    if (!is.null(f)) {
        if (length(f) != nrow(x))
            stop("If 'f' is provided its length has to be equal to 'nrow(x)'")
        if (!is.factor(f))
            f <- factor(f, levels = unique(f))
        fnew <- rep(NA_character_, length(f))
        for (fg in levels(f)) {
            idx <- which(f == fg)
            idxl <- length(idx)
            if (idxl > 1) {
                cors <- cor(t(x[idx, ]), method = method, use = use)
                if (greedy) {
                    ## Ensure diagonal matrix is TRUE so that even if some
                    ## features have a correlation value of NA they are not
                    ## dropped
                    cors <- cors >= threshold
                    cors[cbind(1:idxl, 1:idxl)] <- TRUE
                    fids <- .index_list_to_factor(.group_logic_matrix(cors))
                } else
                    fids <- .group_correlation_matrix(
                        cors, threshold = threshold)
                fnew[idx] <- paste0(fg, ".", .format_groups(fids))
            } else
                fnew[idx] <- paste0(fg, ".1")
        }
        as.factor(fnew)
    } else {
        cors <- cor(t(x), method = method, use = use)
        if (greedy)
            .index_list_to_factor(.group_logic_matrix(cors >= threshold))
        else
            as.factor(.group_correlation_matrix(cors, threshold = threshold))
    }
}

#' @title Grouping of values into sets with smallest differences
#'
#' @description
#'
#' `groupClosest` groups values in `x` for which the difference is smaller than
#' `maxDiff`. As a result, the mean value between the groups will always be
#' larger than `maxDiff`. Values which would be assigned to more than one
#' group, are assigned to the one with the smallest difference to the group
#' mean value. A use case for this would be to group MS features based on their
#' retention times into groups of ~ co-eluting features.
#'
#' In detail, from the sorted `x`, the function starts from the smallest value
#' defining the first group as the one containing all
#' values in `x` with a difference to this first value which is `<= maxDiff`.
#' The next group is the defined based on the next larger value not being part
#' of the first group and includes all values with a difference `<= maxDiff` to
#' this value. For values fulfilling this criteria but being already part of
#' a previous group, the differences to the mean value of the current group
#' and to the mean of previous groups are compared and values are assigned to
#' the group to which they have the smallest difference.
#'
#' Example: values `1.1, 1.9, 2.2` should be grouped with a `maxDiff = 1`. The
#' first group is defined to include all values for which the difference to the
#' first value (`1.1`) is smaller `maxDiff`. Thus, the first group is defined
#' to contain values `1.1 and 1.9`. Then the next group is defined based on the
#' next larger value not part of any group, `2.2`. This group contains values
#' `1.9` and `2.2` with the value `1.9` being already assigned to the first
#' group. The difference between this value `1.9` and the mean of the
#' current group (`mean(c(1.9, 2.2)`) is then compared to the difference of
#' `1.9` to the mean value of the group `1.9` is already part of
#' (which is `mean(c(1.1, 1.9))`). Since the difference to the second group is
#' smaller, `1.9` is removed from the first group and assigned to the second
#' one.
#'
#' @note
#'
#' This grouping approach, in contrast to [xcms::groupOverlaps()], creates
#' smaller groups and values might not be included into the same group even
#' if their difference is smaller than `maxDiff` (see examples below).
#' 
#' @param x `numeric` of values that should be grouped.
#'
#' @param maxDiff `numeric(1)` defining the threshold for difference between
#'     values in `x` to be grouped into the same group.
#'
#' @return `integer` with the group assignment (values grouped together have
#'     the same value).
#'
#' @author Johannes Rainer
#' 
#' @family grouping operations
#'
#' @export
#' 
#' @examples
#'
#' ## The example described above
#' x <- c(1.1, 1.9, 2.2)
#' groupClosest(x)
#' 
#' x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7, 9, 9.5, 15)
#'
#' groupClosest(x)
#' ## value 5.2 was initially grouped with 4.3 (because their difference is
#' ## smaller 1, but then re-grouped together with 5.4 because the difference
#' ## between 5.4 (the next value outside the group of 4.3) and 5.2 is smaller
#' ## than its difference to the mean value of the group for value 4.3
#'
#' ## Example for a case in which values are NOT grouped into the same group
#' ## even if the difference between them is <= maxDiff
#' a <- c(4.9, 5.2, 5.4)
#' groupClosest(a, maxDiff = 0.3)
groupClosest <- function(x, maxDiff = 1) {
    if (is.unsorted(x)) {
        idx <- order(x)
        x <- x[idx]
    } else idx <- integer()
    x_len <- length(x)
    x_groups <- rep(NA_integer_, x_len)
    i <- 1
    group_id <- 1
    while(any(is.na(x_groups))) {
        grp <- which(abs(x - x[i]) <= maxDiff)
        ## Check if they are already part of a previous group
        not_in_prev_grp <- is.na(x_groups[grp])
        in_prev_grp <- grp[!not_in_prev_grp]
        ## grp <- grp[not_in_prev_grp]
        if (length(in_prev_grp)) {
            ## compare difference to current x[i] to mean of previous group(s)
            ## they are part of and assign them to the group with the closest
            ## difference
            ## i_diff <- abs(x[in_prev_grp] - x[i])
            ## Compare to the average of the current group.
            i_diff <- abs(x[in_prev_grp] - mean(x[grp]))
            prev_grp <- x_groups[in_prev_grp]
            to_rem <- rep(FALSE, length(in_prev_grp))
            for (j in unique(prev_grp)) {
                j_diff <- abs(x[in_prev_grp] - mean(x[which(x_groups == j)]))
                to_rem <- to_rem | j_diff < i_diff
            }
            grp <- c(in_prev_grp[!to_rem], grp[not_in_prev_grp])
        }
        x_groups[grp] <- group_id
        group_id <- group_id + 1
        i <- which.max(is.na(x_groups))
    }
    x_groups[idx] <- x_groups
    x_groups
}

.format_groups <- function(x) {
    digits <- ceiling(log10(length(x) + 1L))
    sprintf(paste0("%0", digits, "d"), as.integer(x))
}
