#' @rdname featureGroups
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importMethodsFrom SummarizedExperiment rowData
#'
#' @export
setMethod("featureGroups", "SummarizedExperiment", function(object) {
    if (any(colnames(rowData(object)) == "feature_group"))
        as.character(rowData(object)$feature_group)
    else rep(NA_character_, nrow(rowData(object)))
})

#' @rdname featureGroups
#'
#' @importMethodsFrom SummarizedExperiment rowData<-
#'
#' @export
setReplaceMethod(
    "featureGroups", "SummarizedExperiment", function(object, value) {
        if (!(length(value) == 1 | length(value) == nrow(rowData(object))))
            stop("'value' has to be either of length 1 or equal to the number ",
                 "of rows in object")
        rowData(object)$feature_group <- as.character(value)
        object
    })
