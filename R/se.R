#' @title Quantified LC-MS preprocessing result test data
#'
#' @description
#'
#' The `se` variable is a [SummarizedExperiment::SummarizedExperiment()]
#' object representing the
#' results from a `xcms`-based pre-processing of an LC-MS untargeted
#' metabolomics data set. The raw data files are provided in the `faahKO`
#' package. The pre-processing of this data set is described in detail in
#' the *xcms* vignette of the `xcms` package. This object was created from the
#' `XCMSnExp` result object with the `quantify` method.
#'
#' @name se
#'
#' @examples
#'
#' ## Load the data
#' data(se)
#'
#' library(SummarizedExperiment)
#'
#' ## Access row (feature) data
#' rowData(se)
NULL
