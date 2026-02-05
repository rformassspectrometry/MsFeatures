# General Feature Grouping Concept

This documentation describes the general concepts of feature grouping,
which can be achieved by the different approaches described further
below.

The main function for the stepwise feature grouping is `groupFeatures`.
The selection of the actual grouping algorithm (along with the
definition of its parameters) is done by passing the respective
*parameter* object, along with the object containing the input data and
optional additional arguments, to the `groupFeatures` method.

## Usage

``` r
groupFeatures(object, param, ...)
```

## Arguments

- object:

  input data object on which (with which data) the feature grouping
  should be performed.

- param:

  parameter object which type defines the selection of the grouping
  algorithm.

- ...:

  additional arguments to be passed to the grouping algorithm.

## Value

Depending on the implementation and the input object. Generally the
input object with grouping results added. See respective help pages for
more information.

## Single-step Feature Grouping

Each feature grouping algorithm can be applied individually as a
single-step approach, e.g. by grouping features only on a single feature
property, such as the retention time. Additional feature grouping
approaches might also be implemented that consider combination of
different MS feature properties in a single clustering process.

## Stepwise Feature Grouping Refinement

Stepwise feature grouping evaluates a single property of MS features
(such as their retention time or abundances) at a time to define the
feature groups. Each subsequent grouping step *builds* on the previous
one by eventually sub-grouping each feature group, if needed. Thus,
feature groups get refined in each step. As an example, grouping of
features based on a similar retention time would loosely group features
from all compounds eluting at about the same time from a e.g. liquid
chromatography run. This obviously would also group features
representing ions from different co-eluting compounds. Thus, calling
`groupFeatures` on the previous feature grouping result with a different
parameter object would *refine* these initial feature groups, splitting
them based on another property of the features (such as correlation of
feature abundances across samples).

The advantage of the stepwise approach is that results can be evaluated
after each grouping step and parameters adapted if needed. Also, it
provides flexibility by allowing to change the order of grouping
approaches, or skip individual steps if not suitable for the available
data or the experimental setup.

The major disadvantage is that a wrong group assignment in one of the
initial steps can not be *corrected* for in later steps.

## See also

[`featureGroups()`](https://rformassspectrometry.github.io/MsFeatures/reference/featureGroups.md)
for the function to extract (defined) feature groups from a result
object.

## Author

Johannes Rainer

## Examples

``` r

## For examples please refer to the help pages of the `SimilarRtimeParam` or
## `AbundanceSimilarityParam` objects.
NULL
#> NULL
```
