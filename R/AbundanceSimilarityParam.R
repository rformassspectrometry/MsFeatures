#' @include AllGenerics.R
#' @include grouping-functions.R

#' @title Group features based on abundance similarities across samples
#'
#' @name groupFeatures-similar-abundance
#'
#' @description
#'
#' Group features based on similar abundances (i.e. *feature values*) across
#' samples. Parameter `subset` allows to define a sub set of samples on which
#' the similarity calculation should be performed. It might for example be
#' better to exclude QC samples from the analysis because feature values  are
#' supposed to be constant in these samples.
#'
#' The function first calculates a nxn similarity matrix with n being the
#' number of features and subsequently groups features for which the similarity
#' is higher than the user provided threshold. Parameter `simFun` allows
#' to specify the function to calculate the pairwise similarities on the feature
#' values (eventually transformed by the function specified with parameter
#' `transform`). `simFun` defaults to a function that uses `cor` to calculate
#' similarities between rows in `object` but any function that calculates
#' similarities between rows and that returns a (symmetric) numeric similarity
#' matrix can be used.
#'
#' Parameter `groupFun` allows to specify the function to group the features
#' based on the similarity function. It defaults to `groupSimilarityMatrix`. See
#' [groupSimilarityMatrix()] for details.
#'
#' Additional settings for the `groupFun` and `simFun` functions can be passed
#' to the **parameter object**  with the `...` in the `AbundanceSimilarityParam`
#' constructor function. Other additional parameters specific for the type
#' of `object` can be passed *via* `...` in the `groupFeatures` call.
#'
#' @param groupFun `function` to group features based on the calculated
#'     similarity matrix. Defaults to `groupFun = groupSimilarityMatrix`. See
#'     [groupSimilarityMatrix()] for details.
#'
#' @param object object containing the feature abundances on which features
#'     should be grouped.
#'
#' @param param `AbundanceSimilarityParam` defining the settings for the
#'     grouping based on feature values.
#'
#' @param simFun `function` to be used to calculate the (pairwise)
#'     similarities. Defaults to `simFun = function(x, use =
#'     "pairwise.complete.obs", ...) cor(t(x), use = use, ...)` which
#'     calculates similarities between **rows** in `object` using the `cor`
#'     function. See description for more details.
#'
#' @param subset `integer` or `logical` defining a subset of samples (at least
#'     2) on which the similarity calculation should be performed. By default
#'     the calculation is performed on all samples.
#'
#' @param threshold `numeric(1)` defining the (similarity) threshold to be used
#'     for the feature grouping. This parameter is passed to the `groupFun`
#'     function.
#'
#' @param transform `function` to be used to transform feature abundances prior
#'     to the similarity calculation. Defaults to `transform = identity`.
#'     Alternatively, values could e.g. transformed into log2 scale with
#'     `transform = log2`.
#'
#' @param ... for `AbundanceSimilarityParam`: optional parameters to be passed
#'     along to `simFun` and `groupFun`. For `groupFeatures`: optional
#'     parameters for the extraction/definition of the feature values from
#'     `object`.
#'
#' @family feature grouping methods
#'
#' @seealso [groupFeatures()] for the general concept of feature grouping.
#'
#' @seealso [featureGroups()] for the function to extract defined feature
#'     groups from a `SummarizedExperiment`.
#'
#' @rdname groupFeatures-similar-abundance
#'
#' @exportClass AbundanceSimilarityParam
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' ## Define a simple numeric matrix on which we want to group the rows
#' x <- rbind(
#'     c(12, 34, 231, 234, 9, 5, 7),
#'     c(900, 900, 800, 10, 12, 9, 4),
#'     c(25, 70, 400, 409, 15, 8, 4),
#'     c(12, 13, 14, 15, 16, 17, 18),
#'     c(14, 36, 240, 239, 12, 7, 8),
#'     c(100, 103, 80, 2, 3, 1, 1)
#'     )
#'
#' ## Group rows based on similarity calculated with Pearson's correlation
#' ## on the actual data values (without transforming them).
#' res <- groupFeatures(x, AbundanceSimilarityParam())
#' res
#'
#' ## Use Spearman's rho to correlate rows of the log2 transformed x matrix
#' res <- groupFeatures(x, AbundanceSimilarityParam(method = "spearman",
#'     transform = log2))
#' res

setClass("AbundanceSimilarityParam",
         slots = c(threshold = "numeric",
                   simFun = "function",
                   groupFun = "function",
                   subset = "integer",
                   transform = "function",
                   dots = "list"),
         contains = "Param",
         prototype = prototype(
             threshold = 0.9,
             simFun = function(x, use = "pairwise.complete.obs", ...)
                 cor(t(x), use = use, ...),
             groupFun = groupSimilarityMatrix,
             subset = integer(),
             transform = identity,
             dots = list()
         ),
         validity = function(object) {
             msg <- NULL
             if (length(object@threshold) != 1)
                 msg <- "'threshold' has to be of length 1"
             msg
         })

#' @rdname groupFeatures-similar-abundance
#'
#' @importFrom stats cor
#'
#' @export
AbundanceSimilarityParam <-
    function(threshold = 0.9,
             simFun = function(x, use = "pairwise.complete.obs", ...)
                 cor(t(x), use = use, ...),
             groupFun = groupSimilarityMatrix,
             subset = integer(), transform = identity,
             ...) {
        if (is.logical(subset))
            subset <- which(subset)
        if (is.numeric(subset))
            subset <- as.integer(subset)
        if (!is.integer(subset))
            stop("'subset' has to be either a logical or an integer vector")
        new("AbundanceSimilarityParam", threshold = threshold, simFun = simFun,
            groupFun = groupFun, subset = subset, transform = transform,
            dots = list(...))
    }

#' @rdname groupFeatures-similar-abundance
#'
#' @export
setMethod(
    "groupFeatures",
    signature(object = "matrix", param = "AbundanceSimilarityParam"),
        function(object, param, ...) {
            if (!is.numeric(object))
                stop("'object' needs to be a numeric matrix", call. = FALSE)
            simFun <- param@simFun
            if (identical(simFun, cor)) {
                parms <- list()
                if (!is.null(param@dots$use))
                    parms$use <- param@dots$use
                else parms$use <- "pairwise.complete.obs"
                if (!is.null(param@dots$method))
                    parms$method <- param@dots$method
            } else parms <- param@dots
            res <- do.call(
                simFun, args = c(list(param@transform(object)), parms))
            if (!(is.matrix(res) && nrow(res) == ncol(res) &&
                  nrow(res) == nrow(object) && is.numeric(res)))
                stop("'simFun' did not return the expected results, ",
                     "i.e. a symmetric numeric matrix with ", nrow(object),
                     " columns and rows", call. = FALSE)
            do.call(param@groupFun,
                    args = c(list(res), threshold = param@threshold,
                             param@dots))
        })

## groupFeatures on a `SummarizedExperiment`.
