# Group features based on approximate retention times

Group features based on similar retention time. This method is supposed
to be used as an initial *crude* grouping of LC-MS features based on the
median retention time of all their chromatographic peaks. All features
with a difference in their retention time which is `<=` parameter
`diffRt` of the parameter object are grouped together.

If `object` is a
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):
if a column `"feature_group"` is found in
[`SummarizedExperiment::colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
feature groups defined in that column are further sub-grouped with this
method. See
[`groupFeatures()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.md)
for the general concept of this feature grouping. Also, it might be
required to specify the column in the object's `rowData` containing the
retention times with the `rtime` parameter (which defaults to
`rtime = "rtime"`.

Parameter `groupFun` allows to specify the function that should be used
for the actual grouping. Two possible choices are:

- `groupFun = groupClosest` (the default): this method creates groups of
  features with smallest differences in retention times between the
  individual group members. All differences between group members are
  `< diffRt` (in contrast to the other grouping functions listed below).
  See
  [`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)
  (which is used for the actual grouping on pairwise retention time
  differences) for more details.

- `groupFun = groupConsecutive`: the
  [`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md)
  function groups values together if their difference is smaller than
  `diffRt`. This function iterates over the sorted retention times
  starting the grouping from the lowest value. If the difference of a
  feature to more than one group is smaller `diffRt` it is assigned to
  the group to which its retention time is closest (most similar) to the
  mean retention time of that group. This leads to smaller group sizes.
  Be aware that with this grouping differences in retention times
  between individual features within a group could still be larger
  `diffRt`. See
  [`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md)
  for details and examples.

- `groupFun = MsCoreUtils::group`: this function consecutively groups
  elements together if their difference in retention time is smaller
  than `diffRt`. If two features are grouped into one group, also all
  other features with a retention time within the defined window to any
  of the two features are also included into the feature group. This
  grouping is recursively expanded which can lead, depending on
  `diffRt`, to very large feature groups spanning a large retention time
  window. See
  [`MsCoreUtils::group()`](https://rdrr.io/pkg/MsCoreUtils/man/group.html)
  for details.

Other grouping functions might be added in future. Alternatively it is
also possible to provide a custom function for the grouping operation.

## Usage

``` r
SimilarRtimeParam(diffRt = 1, groupFun = groupClosest)

# S4 method for class 'numeric,SimilarRtimeParam'
groupFeatures(object, param, ...)

# S4 method for class 'SummarizedExperiment,SimilarRtimeParam'
groupFeatures(object, param, rtime = "rtime", ...)
```

## Arguments

- diffRt:

  `numeric(1)` defining the retention time window within which features
  should be grouped. All features with a rtime difference smaller or
  equal than `diffRt` are grouped.

- groupFun:

  `function` that can be used to group values. Defaults to
  `groupFun = groupClosest`. See description for details and
  alternatives.

- object:

  input object that provides the retention times that should be grouped.
  The `MsFeatures` package defines a method for `object` being either a
  `numeric` or a `SummarizedExperiment`.

- param:

  `SimilarRtimeParam` object with the settings for the method.

- ...:

  additional parameters passed to the `groupFun` function.

- rtime:

  for `object` being a
  [`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):
  `character(1)` specifying the column in `rowData(object)` that
  contains the retention time values.

## Value

Depending on parameter `object`:

- for `object` being a `numeric`: returns a `factor` defining the
  feature groups.

- for `object` being `SummarizedExperiment`: returns the input object
  with the feature group definition added to a column `"feature_group"`
  in the result object's `rowData`.

## See also

[`groupFeatures()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.md)
for the general concept of feature grouping.

[`featureGroups()`](https://rformassspectrometry.github.io/MsFeatures/reference/featureGroups.md)
for the function to extract defined feature groups from a
`SummarizedExperiment`.

Other feature grouping methods:
[`groupFeatures-similar-abundance`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures-similar-abundance.md)

## Author

Johannes Rainer

## Examples

``` r

## Simple grouping of a numeric vector.
##
## Define a numeric vector that could represent retention times of features
x <- c(2, 3, 4, 5, 10, 11, 12, 14, 15)

## Group the values using a `group` function. This will create larger
## groups.
groupFeatures(x, param = SimilarRtimeParam(2, MsCoreUtils::group))
#> [1] 1 1 1 1 2 2 2 2 2
#> Levels: 1 2

## Group values using the default `groupClosest` function. This creates
## smaller groups in which all elements have a difference smaller than the
## defined `diffRt` with each other.
groupFeatures(x, param = SimilarRtimeParam(2, groupClosest))
#> [1] 1 1 1 4 2 2 2 3 3
#> Levels: 1 2 3 4

## Grouping on a SummarizedExperiment
##
## load the test SummarizedExperiment object
library(SummarizedExperiment)
data(se)

## No feature groups defined yet
featureGroups(se)
#>   [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [26] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [51] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#>  [76] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [101] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [126] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [151] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [176] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [201] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA

## Determine the column that contains retention times
rowData(se)
#> DataFrame with 225 rows and 11 columns
#>           mzmed     mzmin     mzmax     rtmed     rtmin     rtmax    npeaks
#>       <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
#> FT001     200.1     200.1     200.1   2901.63   2880.73   2922.53         2
#> FT002     205.0     205.0     205.0   2789.39   2782.30   2795.36         8
#> FT003     206.0     206.0     206.0   2788.73   2780.73   2792.86         7
#> FT004     207.1     207.1     207.1   2718.12   2713.21   2726.70         7
#> FT005     219.1     219.1     219.1   2518.82   2517.40   2520.81         3
#> ...         ...       ...       ...       ...       ...       ...       ...
#> FT221    591.30     591.3     591.3   3005.03   2992.87   3006.05         5
#> FT222    592.15     592.1     592.3   3022.11   2981.91   3107.59         6
#> FT223    594.20     594.2     594.2   3418.16   3359.10   3427.90         3
#> FT224    595.25     595.2     595.3   3010.15   2992.87   3013.77         6
#> FT225    596.20     596.2     596.2   2997.91   2992.87   3002.95         2
#>              KO        WT            peakidx  ms_level
#>       <numeric> <numeric>             <list> <integer>
#> FT001         2         0  287, 679,1577,...         1
#> FT002         4         4     47,272,542,...         1
#> FT003         3         4     32,259,663,...         1
#> FT004         4         3     19,249,525,...         1
#> FT005         1         2  639, 788,1376,...         1
#> ...         ...       ...                ...       ...
#> FT221         2         3    349,684,880,...         1
#> FT222         1         3     86,861,862,...         1
#> FT223         1         2  604, 985,1543,...         1
#> FT224         2         3     67,353,876,...         1
#> FT225         0         2  866,1447,1643,...         1

## Column "rtmed" contains the (median) retention time for each feature
## Group features that are eluting within 10 seconds
res <- groupFeatures(se, SimilarRtimeParam(10), rtime = "rtmed")

featureGroups(res)
#>   [1] "FG.021" "FG.001" "FG.001" "FG.057" "FG.062" "FG.019" "FG.062" "FG.075"
#>   [9] "FG.044" "FG.066" "FG.023" "FG.048" "FG.003" "FG.057" "FG.016" "FG.014"
#>  [17] "FG.048" "FG.066" "FG.030" "FG.004" "FG.015" "FG.011" "FG.069" "FG.033"
#>  [25] "FG.033" "FG.025" "FG.071" "FG.069" "FG.002" "FG.033" "FG.013" "FG.050"
#>  [33] "FG.023" "FG.038" "FG.003" "FG.043" "FG.012" "FG.012" "FG.013" "FG.004"
#>  [41] "FG.013" "FG.047" "FG.072" "FG.011" "FG.054" "FG.024" "FG.058" "FG.025"
#>  [49] "FG.032" "FG.012" "FG.003" "FG.001" "FG.005" "FG.005" "FG.006" "FG.004"
#>  [57] "FG.028" "FG.044" "FG.028" "FG.007" "FG.074" "FG.007" "FG.022" "FG.074"
#>  [65] "FG.007" "FG.025" "FG.025" "FG.045" "FG.052" "FG.072" "FG.048" "FG.020"
#>  [73] "FG.027" "FG.008" "FG.064" "FG.018" "FG.008" "FG.016" "FG.060" "FG.018"
#>  [81] "FG.076" "FG.047" "FG.018" "FG.011" "FG.053" "FG.037" "FG.013" "FG.002"
#>  [89] "FG.042" "FG.032" "FG.061" "FG.070" "FG.061" "FG.007" "FG.070" "FG.002"
#>  [97] "FG.008" "FG.046" "FG.065" "FG.008" "FG.022" "FG.022" "FG.033" "FG.061"
#> [105] "FG.077" "FG.047" "FG.014" "FG.023" "FG.064" "FG.045" "FG.021" "FG.078"
#> [113] "FG.007" "FG.042" "FG.079" "FG.059" "FG.040" "FG.025" "FG.027" "FG.040"
#> [121] "FG.027" "FG.018" "FG.009" "FG.065" "FG.045" "FG.020" "FG.018" "FG.021"
#> [129] "FG.024" "FG.052" "FG.015" "FG.037" "FG.038" "FG.052" "FG.043" "FG.036"
#> [137] "FG.067" "FG.035" "FG.067" "FG.053" "FG.053" "FG.032" "FG.073" "FG.071"
#> [145] "FG.073" "FG.037" "FG.054" "FG.056" "FG.028" "FG.056" "FG.006" "FG.053"
#> [153] "FG.014" "FG.060" "FG.034" "FG.049" "FG.044" "FG.049" "FG.041" "FG.026"
#> [161] "FG.059" "FG.036" "FG.008" "FG.006" "FG.008" "FG.019" "FG.035" "FG.016"
#> [169] "FG.015" "FG.015" "FG.006" "FG.041" "FG.024" "FG.045" "FG.055" "FG.029"
#> [177] "FG.030" "FG.015" "FG.059" "FG.068" "FG.068" "FG.055" "FG.055" "FG.031"
#> [185] "FG.031" "FG.031" "FG.020" "FG.080" "FG.025" "FG.081" "FG.029" "FG.049"
#> [193] "FG.029" "FG.082" "FG.040" "FG.010" "FG.010" "FG.068" "FG.046" "FG.040"
#> [201] "FG.039" "FG.009" "FG.058" "FG.063" "FG.083" "FG.063" "FG.026" "FG.017"
#> [209] "FG.065" "FG.009" "FG.009" "FG.017" "FG.017" "FG.034" "FG.084" "FG.011"
#> [217] "FG.072" "FG.011" "FG.026" "FG.051" "FG.051" "FG.019" "FG.034" "FG.039"
#> [225] "FG.050"

## Evaluating differences between retention times within each feature group
rts <- split(rowData(res)$rtmed, featureGroups(res))
lapply(rts, function(z) abs(diff(z)) <= 10)
#> $FG.001
#> [1] TRUE TRUE
#> 
#> $FG.002
#> [1] TRUE TRUE
#> 
#> $FG.003
#> [1] TRUE TRUE
#> 
#> $FG.004
#> [1] TRUE TRUE
#> 
#> $FG.005
#> [1] TRUE
#> 
#> $FG.006
#> [1] TRUE TRUE TRUE
#> 
#> $FG.007
#> [1] TRUE TRUE TRUE TRUE
#> 
#> $FG.008
#> [1] TRUE TRUE TRUE TRUE TRUE
#> 
#> $FG.009
#> [1] TRUE TRUE TRUE
#> 
#> $FG.010
#> [1] TRUE
#> 
#> $FG.011
#> [1] TRUE TRUE TRUE TRUE
#> 
#> $FG.012
#> [1] TRUE TRUE
#> 
#> $FG.013
#> [1] TRUE TRUE TRUE
#> 
#> $FG.014
#> [1] TRUE TRUE
#> 
#> $FG.015
#> [1] TRUE TRUE TRUE TRUE
#> 
#> $FG.016
#> [1] TRUE TRUE
#> 
#> $FG.017
#> [1] TRUE TRUE
#> 
#> $FG.018
#> [1] TRUE TRUE TRUE TRUE
#> 
#> $FG.019
#> [1] TRUE TRUE
#> 
#> $FG.020
#> [1] TRUE TRUE
#> 
#> $FG.021
#> [1] TRUE TRUE
#> 
#> $FG.022
#> [1] TRUE TRUE
#> 
#> $FG.023
#> [1] TRUE TRUE
#> 
#> $FG.024
#> [1] TRUE TRUE
#> 
#> $FG.025
#> [1] TRUE TRUE TRUE TRUE TRUE
#> 
#> $FG.026
#> [1] TRUE TRUE
#> 
#> $FG.027
#> [1] TRUE TRUE
#> 
#> $FG.028
#> [1] TRUE TRUE
#> 
#> $FG.029
#> [1] TRUE TRUE
#> 
#> $FG.030
#> [1] TRUE
#> 
#> $FG.031
#> [1] TRUE TRUE
#> 
#> $FG.032
#> [1] TRUE TRUE
#> 
#> $FG.033
#> [1] TRUE TRUE TRUE
#> 
#> $FG.034
#> [1] TRUE TRUE
#> 
#> $FG.035
#> [1] TRUE
#> 
#> $FG.036
#> [1] TRUE
#> 
#> $FG.037
#> [1] TRUE TRUE
#> 
#> $FG.038
#> [1] TRUE
#> 
#> $FG.039
#> [1] TRUE
#> 
#> $FG.040
#> [1] TRUE TRUE TRUE
#> 
#> $FG.041
#> [1] TRUE
#> 
#> $FG.042
#> [1] TRUE
#> 
#> $FG.043
#> [1] TRUE
#> 
#> $FG.044
#> [1] TRUE TRUE
#> 
#> $FG.045
#> [1] TRUE TRUE TRUE
#> 
#> $FG.046
#> [1] TRUE
#> 
#> $FG.047
#> [1] TRUE TRUE
#> 
#> $FG.048
#> [1] TRUE TRUE
#> 
#> $FG.049
#> [1] TRUE TRUE
#> 
#> $FG.050
#> [1] TRUE
#> 
#> $FG.051
#> [1] TRUE
#> 
#> $FG.052
#> [1] TRUE TRUE
#> 
#> $FG.053
#> [1] TRUE TRUE TRUE
#> 
#> $FG.054
#> [1] TRUE
#> 
#> $FG.055
#> [1] TRUE TRUE
#> 
#> $FG.056
#> [1] TRUE
#> 
#> $FG.057
#> [1] TRUE
#> 
#> $FG.058
#> [1] TRUE
#> 
#> $FG.059
#> [1] TRUE TRUE
#> 
#> $FG.060
#> [1] TRUE
#> 
#> $FG.061
#> [1] TRUE TRUE
#> 
#> $FG.062
#> [1] TRUE
#> 
#> $FG.063
#> [1] TRUE
#> 
#> $FG.064
#> [1] TRUE
#> 
#> $FG.065
#> [1] TRUE TRUE
#> 
#> $FG.066
#> [1] TRUE
#> 
#> $FG.067
#> [1] TRUE
#> 
#> $FG.068
#> [1] TRUE TRUE
#> 
#> $FG.069
#> [1] TRUE
#> 
#> $FG.070
#> [1] TRUE
#> 
#> $FG.071
#> [1] TRUE
#> 
#> $FG.072
#> [1] TRUE TRUE
#> 
#> $FG.073
#> [1] TRUE
#> 
#> $FG.074
#> [1] TRUE
#> 
#> $FG.075
#> logical(0)
#> 
#> $FG.076
#> logical(0)
#> 
#> $FG.077
#> logical(0)
#> 
#> $FG.078
#> logical(0)
#> 
#> $FG.079
#> logical(0)
#> 
#> $FG.080
#> logical(0)
#> 
#> $FG.081
#> logical(0)
#> 
#> $FG.082
#> logical(0)
#> 
#> $FG.083
#> logical(0)
#> 
#> $FG.084
#> logical(0)
#> 

## One feature group ("FG.053") has elements with a difference larger 10:
rts[["FG.053"]]
#> [1] 3473.784 3470.313 3472.164 3474.992
abs(diff(rts[["FG.053"]]))
#> [1] 3.470799 1.851390 2.828125

## But the difference between the **sorted** retention times is < 10:
abs(diff(sort(rts[["FG.053"]])))
#> [1] 1.851390 1.619410 1.208715

## Feature grouping with pre-defined feature groups: groupFeatures will
## sub-group the pre-defined feature groups, features with the feature group
## being `NA` are skipped. Below we perform the feature grouping only on
## features 40 to 70
fgs <- rep(NA, nrow(rowData(se)))
fgs[40:70] <- "FG"
featureGroups(se) <- fgs
res <- groupFeatures(se, SimilarRtimeParam(10), rtime = "rtmed")
featureGroups(res)
#>   [1] NA       NA       NA       NA       NA       NA       NA       NA      
#>   [9] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [17] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [25] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [33] NA       NA       NA       NA       NA       NA       NA       "FG.001"
#>  [41] "FG.009" "FG.010" "FG.007" "FG.001" "FG.011" "FG.012" "FG.006" "FG.004"
#>  [49] "FG.013" "FG.002" "FG.001" "FG.001" "FG.002" "FG.002" "FG.014" "FG.001"
#>  [57] "FG.005" "FG.015" "FG.005" "FG.003" "FG.008" "FG.003" "FG.016" "FG.008"
#>  [65] "FG.003" "FG.004" "FG.004" "FG.006" "FG.017" "FG.007" NA       NA      
#>  [73] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [81] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [89] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [97] NA       NA       NA       NA       NA       NA       NA       NA      
#> [105] NA       NA       NA       NA       NA       NA       NA       NA      
#> [113] NA       NA       NA       NA       NA       NA       NA       NA      
#> [121] NA       NA       NA       NA       NA       NA       NA       NA      
#> [129] NA       NA       NA       NA       NA       NA       NA       NA      
#> [137] NA       NA       NA       NA       NA       NA       NA       NA      
#> [145] NA       NA       NA       NA       NA       NA       NA       NA      
#> [153] NA       NA       NA       NA       NA       NA       NA       NA      
#> [161] NA       NA       NA       NA       NA       NA       NA       NA      
#> [169] NA       NA       NA       NA       NA       NA       NA       NA      
#> [177] NA       NA       NA       NA       NA       NA       NA       NA      
#> [185] NA       NA       NA       NA       NA       NA       NA       NA      
#> [193] NA       NA       NA       NA       NA       NA       NA       NA      
#> [201] NA       NA       NA       NA       NA       NA       NA       NA      
#> [209] NA       NA       NA       NA       NA       NA       NA       NA      
#> [217] NA       NA       NA       NA       NA       NA       NA       NA      
#> [225] NA      
```
