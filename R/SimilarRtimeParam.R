#' @title Group features based on similar retention times
#'
#' @name groupFeatures-similar-rtime
#'
#' @description
#'
#' Group features based on similar retention time. This method is supposed to be
#' used as an initial *crude* grouping of co-eluting features based on their
#' retention. All features with a difference in their retention time
#' `<=` parameter `diffRt` of the `SimilarRtimeParam` object
#' are grouped together. If this method is called on an `XCMSnExp` object or an
#' `SummarizedExperiment` and a column `"feature_group"` is present in
#' its [xcms::featureDefinitions()] or [SummarizedExperiment::rowData] data
#' frame, respectively, these feature groups are further sub-grouped by this
#' method.
#'
#' Two different grouping methods are available:
#'
#' - `method = "greedy"`: this approach consecutively groups elements together
#'   if their difference in retention time is smaller than `diffRt`. If two
#'   features are grouped into one group, also all other features with a
#'   retention time within the defined window to any of the two features are
#'   also included into the feature group. This grouping is recursively
#'   expanded which can lead, depending on `diffRt` to very large feature
#'   groups spanning a large retention time window.
#'
#' - `method = "groupClosest"`: this approach uses the [groupClosest()] function
#'   that groups values together if their difference is smaller than `diffRt`.
#'   If the difference of a feature to more than one group is smaller `diffRt`
#'   it is assigned to the group to which its retention time is closest (most
#'   similar) to the mean retention time of that group. This leads to smaller
#'   group sizes. See [groupClosest()] for details and examples.
#'
#' @param column for `groupFeatures,SummarizedExperiment`: `character(1)`
#'     defining the column in the object's `rowData` containing the retention
#'     times on which the features (rows) should be grouped. Defaults to
#'     `column = "rtime"`.
#'
#' @param diffRt `numeric(1)` defining the retention time window within which
#'     features should be grouped. All features with a rtime difference
#'     smaller or equal than `diffRt` are grouped.
#'
#' @param method `character(1)` defining which grouping approach should be
#'     taken. Allowed values are `method = "groupClosest"` (the default) and
#'     `method = "greedy"`. See description for details.
#'
#' @param msLevel for `groupFeatures,XCMSnExp`: `integer(1)` defining the MS
#'     level on which the features should be grouped.
#'
#' @param object can be a
#'     - [XCMSnExp()] object containing correspondence results (i.e. feature
#'       definitions).
#'     - `numeric` with the retention times.
#'     - [SummarizedExperiment()] object were each row is considered to
#'       represent one *feature*. One column in its `rowData` is supposed to
#'       contain the retention times.
#'
#' @param param `SimilarRtimeParam` object with the settings for the method.
#'
#' @return
#'
#' depending on the input `object`:
#'
#' - `XCMSnExp`: the input object with feature groups added in column
#'   `"feature_group"` of its `featureDefinitions` data frame.
#' - `SummarizedExperiment`: the input object with feature groups added in
#'   column `"feature_group"` in the object's `rowData`.
#' - `numeric`: an `integer` with the group assignment of each value.
#'
#' @family feature grouping methods
#'
#' @rdname groupFeatures-similar-rtime
#'
#' @importClassesFrom xcms Param
#'
#' @exportClass SimilarRtimeParam
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ### --------------------------------------------- ###
#' ## Grouping a simple numeric vector
#'
#' x <- c(3, 5, 3, 4, 12, 4, 13, 11)
#' groupFeatures(x, SimilarRtimeParam(1))
#'
#'
#' ### --------------------------------------------- ###
#' ## Feature grouping on a result object from xcms
#'
#' library(xcms)
#' ## Loading the example data
#' data("xdata")
#' ## Update the path to the files for the local system
#' xcms::dirname(xdata) <- paste0(system.file("cdf/", package = "faahKO"),
#'                                rep(c("KO", "WT"), each = 4))
#'
#' ## The feature definitions
#' featureDefinitions(xdata)
#'
#' res <- groupFeatures(xdata, SimilarRtimeParam(25))
#'
#' ## The feature groups
#' table(featureDefinitions(res)$feature_group)
NULL

setClass("SimilarRtimeParam",
         slots = c(diffRt = "numeric",
                   method = "character"),
         contains = "Param",
         prototype = prototype(
             diffRt = 1,
             method = "groupClosest"
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@diffRt) != 1 || object@diffRt < 0)
                 msg <- c("'diffRt' has to be a positive numeric of length 1")
             if (length(object@method) != 1 || !object@method %in%
                 c("greedy", "groupClosest"))
                 msg <- c(
                     msg,
                     "'method' should be either \"groupClosest\" or \"greedy\"")
             msg
         })

#' @rdname groupFeatures-similar-rtime
#'
#' @export
SimilarRtimeParam <- function(diffRt = 1, method = c("groupClosest", "greedy")) {
    method <- match.arg(method)
    new("SimilarRtimeParam", diffRt = diffRt, method = method)
}

#' @description
#'
#' Helper function that performs the grouping based on the provided retention
#' times.
#'
#' @param x `numeric` with the retention times.
#'
#' @param f `factor` or `character` optionally providing pre-defined feature
#'     groups.
#'
#' @return `character` with the feature groups.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @importFrom MsCoreUtils group
.group_similar_rtime <- function(x, f = rep("FG", length(x)),
                                 method = "groupClosest",
                                 tol = 2) {
    f_new <- as.character(f)
    if (is.factor(f))
        f <- droplevels(f)
    else
        f <- factor(f, levels = unique(f))
    for (fg in levels(f)) {
        idx <- which(f == fg)
        idxl <- length(idx)
        if (idxl > 1) {
            if (method == "greedy")
                fids <- group(x[idx], tolerance = tol)
            if (method == "groupClosest")
                fids <- groupClosest( x[idx], maxDiff = tol)
            f_new[idx] <- paste0(fg, ".", .format_groups(fids))
        } else
            f_new[idx] <- paste0(fg, ".1")
    }
    f_new
}

#' @rdname groupFeatures-similar-rtime
setMethod(
    "groupFeatures",
    signature(object = "numeric", param = "SimilarRtimeParam"),
    function(object, param) {
        .group_similar_rtime(object, method = param@method,
                             tol = param@diffRt)
    })

#' @rdname groupFeatures-similar-rtime
#'
#' @importMethodsFrom xcms hasFeatures featureDefinitions featureDefinitions<-
#'
#' @importMethodsFrom MSnbase fileNames
#'
#' @importClassesFrom xcms XCMSnExp XProcessHistory
setMethod(
    "groupFeatures",
    signature(object = "XCMSnExp", param = "SimilarRtimeParam"),
    function(object, param, msLevel = 1L) {
        if (!hasFeatures(object))
            stop("No feature definitions present. Please run ",
                 "first 'groupChromPeaks'")
        if (length(msLevel) > 1)
            stop("Currently only grouping of features from a single MS level",
                 " is supported.")
        if (any(colnames(featureDefinitions(object)) == "ms_level"))
            is_msLevel <- featureDefinitions(object)$ms_level == msLevel
        else is_msLevel <- rep(TRUE, nrow(featureDefinitions(object)))
        if (any(colnames(featureDefinitions(object)) == "feature_group"))
            f <- featureDefinitions(object)$feature_group
        else
            f <- rep("FG", nrow(featureDefinitions(object)))
        f[!is_msLevel] <- NA
        f_new <- .group_similar_rtime(featureDefinitions(object)$rtmed, f = f,
                                      method = param@method,
                                      tol = param@diffRt)
        featureDefinitions(object)$feature_group <- f_new
        xph <- new("XProcessHistory", param = param, date = date(),
                   type = xcms:::.PROCSTEP.FEATURE.GROUPING,
                   fileIndex = 1:length(fileNames(object)),
                   msLevel = as.integer(msLevel))
        object@.processHistory[[(length(object@.processHistory) + 1)]] <- xph
        object
    })

#' @rdname groupFeatures-similar-rtime
#'
#' @importMethodsFrom SummarizedExperiment rowData rowData<-
setMethod(
    "groupFeatures",
    signature(object = "SummarizedExperiment", param = "SimilarRtimeParam"),
    function(object, param, column = "rtime") {
        if (any(colnames(rowData(object)) == "feature_group"))
            f <- rowData(object)$feature_group
        else
            f <- rep("FG", nrow(rowData(object)))
        f_new <- .group_similar_rtime(rowData(object)[, column], f = f,
                                      method = param@method,
                                      tol = param@diffRt)
        rowData(object)$feature_group <- f_new
        object
    })
