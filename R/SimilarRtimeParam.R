#' @include grouping-functions.R

#' @title Group features based on approximate retention times
#'
#' @name groupFeatures-similar-rtime
#'
#' @description
#'
#' Group features based on similar retention time. This method is supposed to be
#' used as an initial *crude* grouping of LC-MS features based on the median
#' retention time of all their chromatographic peaks. All features with a
#' difference in their retention time which is `<=` parameter `diffRt` of the
#' parameter object are grouped together. If a column `"feature_group"` is
#' found in
#' [SummarizedExperiment::colData()] this is further sub-grouped by this method.
#'
#' Parameter `groupFun` allows to specify the function that should be used for
#' the actual grouping. Two possible choices are:
#'
#' - `groupFun = groupClosest` (the default): the [groupClosest()] function
#'   groups values together if their difference is smaller than `diffRt`.
#'   If the difference of a feature to more than one group is smaller `diffRt`
#'   it is assigned to the group to which its retention time is closest (most
#'   similar) to the mean retention time of that group. This leads to smaller
#'   group sizes. See [groupClosest()] for details and examples.
#'
#' - `groupFun = MsCoreUtils::group`: this function consecutively groups
#'   elements together if their difference in retention time is smaller than
#'   `diffRt`. If two features are grouped into one group, also all other
#'   features with a retention time within the defined window to any of the two
#'   features are also included into the feature group. This grouping is
#'   recursively expanded which can lead, depending on `diffRt`, to very large
#'   feature groups spanning a large retention time window. See
#'   [MsCoreUtils::group()] for details.
#'
#' Other grouping functions might be added in future. Alternatively it is also
#' possible to provide a custom function for the grouping operation.
#'
#' @param diffRt `numeric(1)` defining the retention time window within which
#'     features should be grouped. All features with a rtime difference
#'     smaller or equal than `diffRt` are grouped.
#'
#' @param groupFun `function` that can be used to group values. Defaults to
#'     `groupFun = groupClosest`. See description for details and alternatives.
#'
#' @param object input object that provides the retention times that should be
#'     grouped. The `MsFeatures` package defines a method for `object` being
#'     either a `numeric` or a `SummarizedExperiment`.
#'
#' @param param `SimilarRtimeParam` object with the settings for the method.
#'
#' @param ... additional parameters passed to the `groupFun` function.
#'
#' @return
#'
#' Depending on parameter `object`:
#'
#' - for `object` being a `numeric`: returns a `factor` defining the feature
#'   groups.
#'
#' @family feature grouping methods
#'
#' @seealso [groupFeatures()] for the general concept of feature grouping.
#'
#' @rdname groupFeatures-similar-rtime
#'
#' @importClassesFrom ProtGenerics Param
#'
#' @exportClass SimilarRtimeParam
#'
#' @importFrom MsCoreUtils group
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Simple grouping of a numeric vector.
#' ##
#' ## Define a numeric vector that could represent retention times of features
#' x <- c(2, 3, 4, 5, 10, 11, 12, 14, 15)
#'
#' ## Group the values using a `group` function. This will create larger
#' ## groups.
#' groupFeatures(x, param = SimilarRtimeParam(2, MsCoreUtils::group))
#'
#' ## Group values using the default `groupClosest` function. This creates
#' ## smaller groups in which all elements have a difference smaller than the
#' ## defined `diffRt` with each other.
#' groupFeatures(x, param = SimilarRtimeParam(2, groupClosest))
NULL

setClass("SimilarRtimeParam",
         slots = c(diffRt = "numeric",
                   groupFun = "function"),
         contains = "Param",
         prototype = prototype(
             diffRt = 1,
             groupFun = groupClosest
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@diffRt) != 1 || object@diffRt < 0)
                 msg <- c("'diffRt' has to be a positive numeric of length 1")
             if (!is.function(object@groupFun))
                 msg <- c(msg, "'groupFun' should be a function")
             msg
         })

#' @rdname groupFeatures-similar-rtime
#'
#' @importFrom methods new
#'
#' @export
SimilarRtimeParam <- function(diffRt = 1, groupFun = groupClosest) {
    new("SimilarRtimeParam", diffRt = diffRt, groupFun = groupFun)
}

#' @rdname groupFeatures-similar-rtime
#'
#' @importFrom MsCoreUtils group
setMethod(
    "groupFeatures",
    signature(object = "numeric", param = "SimilarRtimeParam"),
    function(object, param, ...) {
        factor(param@groupFun(object, param@diffRt, ...))
    })

#' @rdname groupFeatures-similar-rtime
setMethod(
    "groupFeatures",
    signature(object = "SummarizedExperiment", param = "SimilarRtimeParam"),
    function(object, param, ...) {

    })

#' Simple helper function to add a (minimum) number of leading 0s to a numeric
#' value to create *better looking* identifiers.
#'
#' @noRd
.format_id <- function(x, length = 3L) {
    digits <- max(ceiling(log10(length(x) + 1L)), length)
    sprintf(paste0("%0", digits, "d"), as.integer(x))
}
