% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corRows.R
\name{corRows}
\alias{corRows}
\title{Correlate rows of a numeric matrix}
\usage{
corRows(
  x,
  y = NULL,
  use = "pairwise.complete.obs",
  method = c("pearson", "kendall", "spearman"),
  ...
)
}
\arguments{
\item{x}{\code{numeric} \code{matrix}.}

\item{y}{not supported (ignored).}

\item{use}{see information for parameter \code{use} in \code{\link[stats:cor]{stats::cor()}}. Defaults
to \code{use = "pairwise.complete.obs"}.}

\item{method}{see information for parameter \code{method} in \code{\link[stats:cor]{stats::cor()}}.}

\item{...}{additional parameters (ignored).}
}
\value{
\code{matrix} with correlation coefficients between rows in \code{x}.
}
\description{
\code{corRows} is a simple function to perform a pairwise correlation between
\strong{rows} of a numeric matrix by calling \code{\link[stats:cor]{stats::cor()}} on the transposed
input matrix \code{x}.
}
\examples{

## Define a simple numeric matrix
x <- rbind(
    c(12, 34, 231, 234, 9, 5, 7),
    c(900, 900, 800, 10, 12, 9, 4),
    c(25, 70, 400, 409, 15, 8, 4),
    c(12, 13, 14, 15, 16, 17, 18),
    c(14, 36, 240, 239, 12, 7, 8)
    )

corRows(x)

corRows(x, method = "spearman")
}
\author{
Johannes Rainer
}
