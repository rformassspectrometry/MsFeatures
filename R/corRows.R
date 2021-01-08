#' @title Correlate rows of a numeric matrix
#'
#' @description
#'
#' `corRows` is a simple function to perform a pairwise correlation between
#' **rows** of a numeric matrix by calling [stats::cor()] on the transposed
#' input matrix `x`.
#'
#' @param x `numeric` `matrix`.
#'
#' @param method see information for parameter `method` in [stats::cor()].
#'
#' @param use see information for parameter `use` in [stats::cor()]. Defaults
#'     to `use = "pairwise.complete.obs"`.
#'
#' @param y not supported (ignored).
#'
#' @param ... additional parameters (ignored).
#'
#' @return `matrix` with correlation coefficients between rows in `x`.
#'
#' @author Johannes Rainer
#'
#' @export
#'
#' @examples
#'
#' ## Define a simple numeric matrix
#' x <- rbind(
#'     c(12, 34, 231, 234, 9, 5, 7),
#'     c(900, 900, 800, 10, 12, 9, 4),
#'     c(25, 70, 400, 409, 15, 8, 4),
#'     c(12, 13, 14, 15, 16, 17, 18),
#'     c(14, 36, 240, 239, 12, 7, 8)
#'     )
#'
#' corRows(x)
#'
#' corRows(x, method = "spearman")
corRows <- function(x, y = NULL, use = "pairwise.complete.obs",
                    method = c("pearson", "kendall", "spearman"), ...) {
    method <- match.arg(method)
    cor(t(x), use = use, method = method)
}
