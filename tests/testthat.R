library("testthat")
library("MsFeatures")

data("xdata")
## Update the path to the files for the local system
xcms::dirname(xdata) <- paste0(system.file("cdf/", package = "faahKO"),
                               rep(c("KO", "WT"), each = 4))

test_check("MsFeatures")
