#' @title Grouping of sorted values into sets with smallest differences
#'
#' @description
#'
#' `groupConsecutive` groups **sorted** values in `x` for which the difference
#' is smaller than `maxDiff`. As a result, the mean difference between the
#' groups will always be larger than `maxDiff`, but difference between
#' individual values within the same group (e.g. between the first and last)
#' can be larger `maxDiff`.
#'
#' In detail, from the sorted `x`, the function starts from the smallest value
#' defining the first group as the one containing all values in `x` with a
#' difference to this first value which is `<= maxDiff`.
#' The next group is the defined based on the next larger value that is not part
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
#' The difference between consecutive (ordered) values within a defined group
#' is always `<= maxDiff`, but the difference between e.g. the first and the
#' last of the (ordered) values can be larger than `maxDiff`. See
#' [groupClosest()] for a more stringent grouping function.
#'
#' @param x `numeric` of values that should be grouped.
#'
#' @param maxDiff `numeric(1)` defining the threshold for difference between
#'     values in `x` to be grouped into the same group.
#'
#' @return `integer` with the group assignment (values grouped together have
#'     the same return value).
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
#' groupConsecutive(x)
#'
#' x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7,
#'     9, 9.5, 15)
#'
#' groupConsecutive(x)
#' ## value 5.2 was initially grouped with 4.3 (because their difference is
#' ## smaller 1, but then re-grouped together with 5.4 because the difference
#' ## between 5.4 (the next value outside the group of 4.3) and 5.2 is smaller
#' ## than its difference to the mean value of the group for value 4.3
#'
#' ## Example for a case in which values are NOT grouped into the same group
#' ## even if the difference between them is <= maxDiff
#' a <- c(4.9, 5.2, 5.4)
#' groupConsecutive(a, maxDiff = 0.3)
groupConsecutive <- function(x, maxDiff = 1) {
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

#' @title Group rows of a diagonal matrix using a threshold
#'
#' @description
#'
#' This function groups elements (rows or columns) of a diagonal matrix, such as
#' a pairwise correlation matrix or similarity matrix, with a value `>=
#' threshold`. This creates clusters of elements in which **all** elements have
#' a value `>= threshold` with **any** other element in that cluster. On a
#' correlation matrix (such as created with `cor`) it will generate small
#' clusters of highly correlated elements. Note however that single elements in
#' one cluster could also have a correlation `>= threshold` to another element
#' in another cluster. The average similarity to its own cluster will however
#' be higher to that of the other.
#'
#' @details
#'
#' The algorithm is defined as follows:
#' - all pairs of values in `x` which are `>= threshold` are identified and
#'   sorted decreasingly.
#' - starting with the pair with the highest correlation, groups are defined:
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
#' @param ... ignored.
#'
#' @return `integer` same length than `nrow(x)`, grouped elements (rows) defined
#'     by the same value.
#'
#' @author Johannes Rainer
#'
#' @family grouping operations
#'
#' @export groupSimilarityMatrix
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
#' groupSimilarityMatrix(x, threshold = 0.9)
#'
#' groupSimilarityMatrix(x, threshold = 0.1)
#'
#' ## Add also a correlation between 3 and 2
#' x[2, 3] <- 0.9
#' x[3, 2] <- 0.9
#' x
#' groupSimilarityMatrix(x, threshold = 0.9)
#'
#' ## Add a higher correlation between 4 and 5
#' x[4, 5] <- 0.99
#' x[5, 4] <- 0.99
#' x
#' groupSimilarityMatrix(x, threshold = 0.9)
#'
#' ## Increase correlation between 2 and 3
#' x[2, 3] <- 0.92
#' x[3, 2] <- 0.92
#' x
#' groupSimilarityMatrix(x, threshold = 0.9) ## Don't break previous cluster!
groupSimilarityMatrix <- function(x, threshold = 0.9, full = TRUE, ...) {
    nr <- nrow(x)
    if (nr != ncol(x))
        stop("'x' should be a symmetric matrix")
    if (!full)
        x[lower.tri(x)] <- NA
    res <- rep(NA_integer_, nr)
    sl <- seq_len(nr)
    x[cbind(sl, sl)] <- NA
    idx_pairs <- which(x >= threshold, arr.ind = TRUE)
    idx_pairs <- idx_pairs[order(x[idx_pairs], decreasing = TRUE), ,
                           drop = FALSE]
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
                    if (!(any(is.na(cor_to_grp)) ||
                          any(cor_to_grp < threshold)))
                        mean_cor_to_grp[grp] <- mean(cor_to_grp)
                }
                mean_cor_to_grp <- mean_cor_to_grp[mean_cor_to_grp > 0]
                if (length(mean_cor_to_grp)) {
                    ## yes: put them into the group with which both have the
                    ## highest correlation.
                    res[idx] <- as.integer(
                        names(sort(mean_cor_to_grp, decreasing = TRUE)))[1L]
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
                    !any(cor_1_2 < threshold) && mcor_1_2 >= mcor_1_1)
                    res[idx_pairs[i, 1]] <- got_grp[2]
                else if (is.finite(mcor_2_1) && is.finite(mcor_2_2) &&
                         !any(cor_2_1 < threshold) && mcor_2_1 >= mcor_2_2)
                    res[idx_pairs[i, 2]] <- got_grp[1]
            } # else nothing to do - they are already in the same group
        }
    }
    nas <- is.na(res)
    if (any(nas))
        res[nas] <- seq(grp_id, length.out = sum(nas))
    res
}


#' @title Finds all leaves of a cluster in a tree data structure.
#'
#' @description
#' Auxillary function for `groupSimilarityMatrixTree()` that searches through a
#' complete clustering tree for branches with internal distances less
#' than maxDiff.
#'
#' @details
#'
#' The algorithm recursively seeks the tree for the members of a cluster.
#'
#' @param hc is a hclust object.
#'
#' @param cluster_index is The index of the cluster in the hclust object.
#'
#' @return `integers` representing all leafs (retention time indexes) of a given
#'  node (and sub-nodes) in the tree
#'
#' @author Johan Lassen
#'
#' @family grouping operations
#' @noRd
get_cluster_members <- function(hc, cluster_index) {
  #' Recursively gets the members of a cluster in an hclust object.
  if (cluster_index < 0) {
    return(-cluster_index)
  } else {
    merge_row <- hc$merge[cluster_index, ]
    left_members <- get_cluster_members(hc, merge_row[1])
    right_members <- get_cluster_members(hc, merge_row[2])
    return(c(left_members, right_members))
  }
}

#' @title Identifies groups with intragroup pairwise distances less than maxDiff
#'
#' @description
#' Searches through a complete clustering tree for branches with internal
#' distances less than maxDiff.
#'
#' @details
#'
#' The search is performed top-down so that every node is evaluated for
#'  qualifying as a group.
#' If the node does not qualify the algorithm iterates to the next node in the
#' tree. If the node qualifies, all members are added as a group to the *groups
#' vector* and all following child nodes are skipped.
#'
#' @param dists is a matrix object with pairwise retention time distances.
#'
#' @param maxDiff is the maximum pairwise distance between members of a group.
#'
#' @return `integers` representing groups of peaks based on retention times.
#' Should be of the same length as the feature definitions.
#'
#' @author Johan Lassen
#' @importFrom stats hclust
#' @family grouping operations
groupSimilarityMatrixTree <- function(dists, maxDiff) {
  hc <- hclust(as.dist(dists), method = "complete")
  # Convert dist object to vector
  dist_matrix <- as.matrix(dists)
  n_clusters <- nrow(hc$merge)
  already_grouped <- c()
  groups <- seq_len(nrow(dist_matrix))
  for(i in rev(seq_len(n_clusters))){
    # Extract cluster members (top-down approach)
    members <- get_cluster_members(hc, i)
    if (max(dist_matrix[members, members]) > maxDiff) next
    if (any(members %in% already_grouped)) next
    # Save grouped members
    groups[members] <- min(groups[members])
    already_grouped <- c(already_grouped, members)
    if (length(already_grouped) == nrow(dist_matrix)) break
  }
  return(groups)
}

#' @title Group values with differences below threshold
#'
#' @description
#'
#' Group values with a difference between them being smaller than a user
#' defined threshold. This function uses the [groupSimilarityMatrix()] function
#' to create groups with smallest differences between its members. Differences
#' between **all** members of one group are below the user defined threshold
#' `maxDiff`. This is a more stringent grouping than what [groupConsecutive()]
#' performs leading thus to smaller groups (with smaller differences between
#' its members).
#'
#' @param x `numeric` of values that should be grouped.
#'
#' @param maxDiff `numeric(1)` defining the threshold for difference between
#'     values in `x` to be grouped into the same group.
#'
#' @param FUN supported similarity calculation function. Can be either
#'     [groupSimilarityMatrixTree()] (the default) or [groupSimilarityMatrix()].
#'
#' @return `integer` with the group assignment (values grouped together have
#'     the same return value).
#'
#' @author Johannes Rainer
#' @author Johan Lassen
#'
#' @family grouping operations
#'
#' @export
#'
#' @importFrom stats dist
#' @importFrom stats as.dist
#'
#' @examples
#'
#' x <- c(1.1, 1.9, 2.2)
#' groupClosest(x)
#' ## Although the difference between the 1st and 2nd element would be smaller
#' ## than the threshold, they are not grouped because the difference between
#' ## the 2nd and 3rd element is even smaller. The first element is also not
#' ## put into the same group, because it has a difference > diffRt to the 3rd
#' ## element.
#'
#' x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7,
#'     9, 9.5, 15)
#'
#' groupClosest(x)
groupClosest <- function(x, maxDiff = 1, FUN = groupSimilarityMatrixTree) {
  dists <- as.matrix(dist(x, method = "manhattan"))
  FUN(dists, maxDiff)
}
