#' @include hidden_aliases.R
NULL

#' @title General Feature Grouping Concept
#'
#' @name groupFeatures
#'
#' @rdname groupFeatures
#'
#' @description
#'
#' This documentation describes the general concepts of feature grouping, which
#' can be achieved by the different approaches described further below.
#'
#' The main function for the stepwise feature grouping is `groupFeatures`. The
#' selection of the actual grouping algorithm (along with the definition of its
#' parameters) is done by passing the respective *parameter* object, along with
#' the object containing the input data and optional additional arguments,
#' to the `groupFeatures` method.
#'
#' @section Single-step Feature Grouping:
#'
#' Each feature grouping algorithm can be applied individually as a single-step
#' approach, e.g. by grouping features only on a single feature property, such
#' as the retention time. Additional feature grouping approaches might also be
#' implemented that consider combination of different MS feature properties in
#' a single clustering process.
#'
#' @section Stepwise Feature Grouping Refinement:
#'
#' Stepwise feature grouping evaluates a single property of MS features (such
#' as their retention time or abundances) at a time to define the feature
#' groups. Each subsequent grouping step *builds* on the previous one by
#' eventually sub-grouping each feature group, if needed. Thus, feature groups
#' get refined in each step. As an example, grouping of features based on a
#' similar retention time would loosely group features from all compounds
#' eluting at about the same time from a e.g. liquid chromatography run. This
#' obviously would also group features representing ions from different
#' co-eluting compounds. Thus, calling `groupFeatures` on the previous feature
#' grouping result with a different parameter object would *refine* these
#' initial feature groups, splitting them based on another property of the
#' features (such as correlation of feature abundances across samples).
#'
#' The advantage of the stepwise approach is that results can be evaluated
#' after each grouping step and parameters adapted if needed. Also, it provides
#' flexibility by allowing to change the order of grouping approaches, or skip
#' individual steps if not suitable for the available data or the
#' experimental setup.
#'
#' The major disadvantage is that a wrong group assignment in one of the initial
#' steps can not be *corrected* for in later steps.
#'
#' @param object input data object on which (with which data) the feature
#'     grouping should be performed.
#'
#' @param param parameter object which type defines the selection of the
#'     grouping algorithm.
#'
#' @param ... additional arguments to be passed to the grouping algorithm.
#'
#' @author Johannes Rainer
#'
#' @seealso [featureGroups()] for the function to extract (defined) feature
#'     groups from a result object.
NULL

#' @rdname groupFeatures
#'
#' @exportMethod groupFeatures
setGeneric("groupFeatures", function(object, param, ...)
    standardGeneric("groupFeatures"))

#' @title Get or set feature group definitions from an object
#'
#' `featureGroups` and `featureGroups<-` allow to extract or set the feature
#' definitions from the input object. The implementations for
#' [SummarizedExperiment()] get or set the content of a column named
#' `"feature_group"` in the object's `rowData`.
#'
#'  This method should be implemented for all other object for which a
#' [groupFeatures()] method is defined.
#'
#' @param object the input object. In the `MsFeatures` package this method is
#'     implemented for `SummarizedExperiment`.
#'
#' @param value the new value for the *feature groups* variable.
#'
#' @param ... ignored.
#'
#' @return a `character` with the group assignment of the features. Has to have
#'     the same length as there are features in `object.`
#'
#' @author Johannes Rainer
#'
#' @rdname featureGroups
#'
#' @exportMethod featureGroups
#'
#' @exportMethod featureGroups<-
#'
#' @examples
#'
#' ## Load the test SummarizedExperiment
#' library(SummarizedExperiment)
#' data(se)
#'
#' ## No column "feature_group" present in the object, this NA is returned
#' featureGroups(se)
#'
#' ## Add a column "feature_group" to the `rowData` of the object
#' rowData(se)$feature_group <- seq_len(nrow(rowData(se)))
#'
#' featureGroups(se)
setGeneric("featureGroups", function(object, ...)
    standardGeneric("featureGroups"))

#' @rdname featureGroups
setGeneric("featureGroups<-", function(object, value)
    standardGeneric("featureGroups<-"))
