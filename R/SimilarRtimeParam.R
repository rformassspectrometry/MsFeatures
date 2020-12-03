## #' @title Group features based on approximate retention times
## #'
## #' @name groupFeatures-approximate-rtime
## #'
## #' @description
## #'
## #' Group features based on similar retention time. This method is supposed to be
## #' used as an initial *crude* grouping of features based on the median retention
## #' time of all their chromatographic peaks. All features with a difference in
## #' their retention time which is `<=` parameter `diffRt` of the parameter object
## #' are grouped together. If a column `"feature_group"` is found in
## #' [xcms::featureDefinitions()] this is further sub-grouped by this method.
## #'
## #' Two different grouping methods are available:
## #'
## #' - `method = "greedy"`: this approach consecutively groups elements together
## #'   if their difference in retention time is smaller than `diffRt`. If two
## #'   features are grouped into one group, also all other features with a
## #'   retention time within the defined window to any of the two features are
## #'   also included into the feature group. This grouping is recursively
## #'   expanded which can lead, depending on `diffRt` to very large feature
## #'   groups spanning a large retention time window.
## #'
## #' - `method = "groupClosest"`: this approach uses the [groupClosest()] function
## #'   that groups values together if their difference is smaller than `diffRt`.
## #'   If the difference of a feature to more than one group is smaller `diffRt`
## #'   it is assigned to the group to which its retention time is closest (most
## #'   similar) to the mean retention time of that group. This leads to smaller
## #'   group sizes. See [groupClosest()] for details and examples.
## #'
## #' @param diffRt `numeric(1)` defining the retention time window within which
## #'     features should be grouped. All features with a rtime difference
## #'     smaller or equal than `diffRt` are grouped.
## #'
## #' @param method `character(1)` defining which grouping approach should be
## #'     taken. Allowed values are `method = "groupClosest"` (the default) and
## #'     `method = "greedy"`. See description for details.
## #'
## #' @param msLevel `integer(1)` defining the MS level on which the features
## #'     should be grouped.
## #'
## #' @param object [XCMSnExp()] object containing also correspondence results.
## #'
## #' @param param `SimilarRtimeParam` object with the settings for the method.
## #'
## #' @return input `XCMSnExp` with feature groups added (i.e. in column
## #'     `"feature_group"` of its `featureDefinitions` data frame.
## #'
## #' @family feature grouping methods
## #'
## #' @seealso feature-grouping for a general overview.
## #'
## #' @rdname groupFeatures-approximate-rtime
## #'
## #' @importClassesFrom xcms Param
## #'
## #' @importFrom MsFeatures groupClosest
## #'
## #' @exportClass SimilarRtimeParam
## #'
## #' @author Johannes Rainer
## #'
## #' @examples
## #'
## #' ## Performing a quick preprocessing of a test data set.
## #' library(faahKO)
## #' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
## #'         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
## #'         system.file('cdf/WT/wt19.CDF', package = "faahKO"))
## #'
## #' od <- readMSData(fls, mode = "onDisk")
## #' xod <- findChromPeaks(
## #'     od, param = CentWaveParam(noise = 10000, snthresh = 40,
## #'                               prefilter = c(3, 10000)))
## #' pdp <- PeakDensityParam(sampleGroups = c(1, 1, 2))
## #' xodg <- groupChromPeaks(xod, param = pdp)
## #'
## #' ## Group features based on similar retention time (i.e. difference <= 2 seconds)
## #' xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 2))
## #'
## #' ## Feature grouping get added to the featureDefinitions in column "feature_group"
## #' head(featureDefinitions(xodg_grp)$feature_group)
## #'
## #' table(featureDefinitions(xodg_grp)$feature_group)
## #' length(unique(featureDefinitions(xodg_grp)$feature_group))
## #'
## #' ## Using the "greedy" method to create larger groups
## #' xodg_grp <- groupFeatures(xodg,
## #'     param = SimilarRtimeParam(diffRt = 2, method = "greedy"))
## #'
## #' length(unique(featureDefinitions(xodg_grp)$feature_group))
## NULL

## setClass("SimilarRtimeParam",
##          slots = c(diffRt = "numeric",
##                    method = "character"),
##          contains = "Param",
##          prototype = prototype(
##              diffRt = 1,
##              method = "groupClosest"
##          ),
##          validity = function(object) {
##              msg <- NULL
##              if (length(object@diffRt) != 1 || object@diffRt < 0)
##                  msg <- c("'diffRt' has to be a positive numeric of length 1")
##              if (length(object@method) != 1 || !object@method %in%
##                  c("greedy", "groupClosest"))
##                  msg <- c(
##                      msg,
##                      "'method' should be either \"groupClosest\" or \"greedy\"")
##              msg
##          })
