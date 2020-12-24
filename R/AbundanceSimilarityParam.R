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
#' `transform`). `simFun` defaults to `cor` but any function that
#' returns a (symmetric) numeric similarity matrix can be used. Parameter
#' `groupFun` allows to specify the function to group the features based on the
#' similarity function. It defaults to `groupSimilarityMatrix`. See
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
#' @param simFun `function` to be used to calculate the (pairwise)
#'     similarities. Defaults to `simFun = cor` in which case also
#'     `use = "pairwise.complete.obs"` is automatically passed along to the
#'     function.
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
#'     to the similarity calculation. Defaults to `transform = log2`. To use
#'     feature values *as is* use `transform = identity`.
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
             simFun = cor,
             groupFun = groupSimilarityMatrix,
             subset = integer(),
             transform = log2,
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
AbundanceSimilarityParam <- function(threshold = 0.9, simFun = cor,
                                     groupFun = groupSimilarityMatrix,
                                     subset = integer(), transform = log2,
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

## groupFeatures on a `matrix`.
## groupFeatures on a `SummarizedExperiment`.