# Group features based on abundance similarities across samples

Group features based on similar abundances (i.e. *feature values*)
across samples. Parameter `subset` allows to define a sub set of samples
on which the similarity calculation should be performed. It might for
example be better to exclude QC samples from the analysis because
feature values are supposed to be constant in these samples.

The function first calculates a nxn similarity matrix with n being the
number of features and subsequently groups features for which the
similarity is higher than the user provided threshold. Parameter
`simFun` allows to specify the function to calculate the pairwise
similarities on the feature values (eventually transformed by the
function specified with parameter `transform`). `simFun` defaults to a
function that uses `cor` to calculate similarities between rows in
`object` but any function that calculates similarities between rows and
that returns a (symmetric) numeric similarity matrix can be used.

If `object` is a
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):
if a column `"feature_group"` is found in
[`SummarizedExperiment::colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
feature groups defined in that column are further sub-grouped with this
method. See
[`groupFeatures()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.md)
for the general concept of this feature grouping.

Parameter `groupFun` allows to specify the function to group the
features based on the similarity function. It defaults to
`groupSimilarityMatrix`. See
[`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)
for details.

Additional settings for the `groupFun` and `simFun` functions can be
passed to the **parameter object** with the `...` in the
`AbundanceSimilarityParam` constructor function. Other additional
parameters specific for the type of `object` can be passed *via* `...`
in the `groupFeatures` call.

## Usage

``` r
AbundanceSimilarityParam(
  threshold = 0.9,
  simFun = corRows,
  groupFun = groupSimilarityMatrix,
  subset = integer(),
  transform = identity,
  ...
)

# S4 method for class 'matrix,AbundanceSimilarityParam'
groupFeatures(object, param, ...)

# S4 method for class 'SummarizedExperiment,AbundanceSimilarityParam'
groupFeatures(object, param, i = 1L, ...)
```

## Arguments

- threshold:

  `numeric(1)` defining the (similarity) threshold to be used for the
  feature grouping. This parameter is passed to the `groupFun` function.

- simFun:

  `function` to be used to calculate (pairwise) similarities (between
  **rows**). Defaults to `simFun = corRows`. See description or
  [`corRows()`](https://rformassspectrometry.github.io/MsFeatures/reference/corRows.md)
  for more details.

- groupFun:

  `function` to group features based on the calculated similarity
  matrix. Defaults to `groupFun = groupSimilarityMatrix`. See
  [`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)
  for details.

- subset:

  `integer` or `logical` defining a subset of samples (at least 2) on
  which the similarity calculation should be performed. By default the
  calculation is performed on all samples.

- transform:

  `function` to be used to transform feature abundances prior to the
  similarity calculation. Defaults to `transform = identity`.
  Alternatively, values could e.g. transformed into log2 scale with
  `transform = log2`.

- ...:

  for `AbundanceSimilarityParam`: optional parameters to be passed along
  to `simFun` and `groupFun`. For `groupFeatures`: optional parameters
  for the extraction/definition of the feature values from `object`.

- object:

  object containing the feature abundances on which features should be
  grouped.

- param:

  `AbundanceSimilarityParam` defining the settings for the grouping
  based on feature values.

- i:

  for `object` being a
  [`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):
  `integer(1)` or `character(1)` specifying either the index or name of
  the the *assay* in `object` that contains the feature values that
  should be used. Use
  [`SummarizedExperiment::assayNames()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  on `object` to list all available assays.

## Value

for object being a `SummarizedExperiment`: a `SummarizedExperiment` with
the grouping results added to a column `"feature_group"` in the object's
`rowData`. For object being a `matrix`: an `integer` of length equal to
the number of rows with the group identifiers.

## See also

[`groupFeatures()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.md)
for the general concept of feature grouping.

[`featureGroups()`](https://rformassspectrometry.github.io/MsFeatures/reference/featureGroups.md)
for the function to extract defined feature groups from a
`SummarizedExperiment`.

Other feature grouping methods:
[`groupFeatures-similar-rtime`](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures-similar-rtime.md)

## Author

Johannes Rainer

## Examples

``` r

## Define a simple numeric matrix on which we want to group the rows
x <- rbind(
    c(12, 34, 231, 234, 9, 5, 7),
    c(900, 900, 800, 10, 12, 9, 4),
    c(25, 70, 400, 409, 15, 8, 4),
    c(12, 13, 14, 15, 16, 17, 18),
    c(14, 36, 240, 239, 12, 7, 8),
    c(100, 103, 80, 2, 3, 1, 1)
    )

## Group rows based on similarity calculated with Pearson's correlation
## on the actual data values (without transforming them).
res <- groupFeatures(x, AbundanceSimilarityParam())
res
#> [1] 1 2 1 3 1 2

## Use Spearman's rho to correlate rows of the log2 transformed x matrix
res <- groupFeatures(x, AbundanceSimilarityParam(method = "spearman",
    transform = log2))
res
#> [1] 2 1 2 3 2 1

## Perform the grouping on a SummarizedExperiment
library(SummarizedExperiment)
data(se)

## Group features based on log2 transformed feature values in the first
## assay of the SummarizedExperiment
res <- groupFeatures(se, param = AbundanceSimilarityParam(threshold = 0.7,
    transform = log2))

featureGroups(res)
#>   [1] "FG.001" "FG.017" "FG.001" "FG.002" "FG.041" "FG.003" "FG.037" "FG.021"
#>   [9] "FG.021" "FG.021" "FG.028" "FG.045" "FG.017" "FG.002" "FG.014" "FG.026"
#>  [17] "FG.048" "FG.027" "FG.008" "FG.007" "FG.008" "FG.017" "FG.002" "FG.033"
#>  [25] "FG.011" "FG.016" "FG.015" "FG.002" "FG.027" "FG.049" "FG.010" "FG.019"
#>  [33] "FG.006" "FG.010" "FG.006" "FG.018" "FG.004" "FG.004" "FG.026" "FG.013"
#>  [41] "FG.020" "FG.008" "FG.018" "FG.013" "FG.008" "FG.011" "FG.011" "FG.022"
#>  [49] "FG.047" "FG.004" "FG.006" "FG.003" "FG.004" "FG.004" "FG.043" "FG.016"
#>  [57] "FG.032" "FG.001" "FG.038" "FG.025" "FG.003" "FG.045" "FG.026" "FG.003"
#>  [65] "FG.050" "FG.041" "FG.015" "FG.022" "FG.008" "FG.002" "FG.011" "FG.013"
#>  [73] "FG.015" "FG.051" "FG.014" "FG.001" "FG.052" "FG.031" "FG.044" "FG.017"
#>  [81] "FG.047" "FG.019" "FG.017" "FG.025" "FG.046" "FG.035" "FG.053" "FG.012"
#>  [89] "FG.012" "FG.008" "FG.054" "FG.034" "FG.009" "FG.028" "FG.030" "FG.030"
#>  [97] "FG.002" "FG.010" "FG.031" "FG.028" "FG.039" "FG.024" "FG.055" "FG.044"
#> [105] "FG.004" "FG.006" "FG.016" "FG.029" "FG.044" "FG.035" "FG.006" "FG.020"
#> [113] "FG.024" "FG.008" "FG.027" "FG.004" "FG.005" "FG.004" "FG.014" "FG.010"
#> [121] "FG.004" "FG.004" "FG.010" "FG.014" "FG.042" "FG.039" "FG.056" "FG.029"
#> [129] "FG.036" "FG.008" "FG.004" "FG.005" "FG.022" "FG.010" "FG.003" "FG.020"
#> [137] "FG.020" "FG.004" "FG.003" "FG.042" "FG.010" "FG.023" "FG.025" "FG.006"
#> [145] "FG.057" "FG.010" "FG.022" "FG.014" "FG.058" "FG.023" "FG.020" "FG.046"
#> [153] "FG.037" "FG.021" "FG.008" "FG.010" "FG.023" "FG.013" "FG.013" "FG.024"
#> [161] "FG.032" "FG.010" "FG.005" "FG.020" "FG.010" "FG.020" "FG.022" "FG.013"
#> [169] "FG.010" "FG.010" "FG.020" "FG.040" "FG.023" "FG.043" "FG.034" "FG.038"
#> [177] "FG.007" "FG.007" "FG.025" "FG.022" "FG.022" "FG.022" "FG.017" "FG.010"
#> [185] "FG.010" "FG.010" "FG.010" "FG.029" "FG.010" "FG.059" "FG.024" "FG.030"
#> [193] "FG.010" "FG.019" "FG.020" "FG.020" "FG.026" "FG.022" "FG.010" "FG.020"
#> [201] "FG.020" "FG.040" "FG.013" "FG.027" "FG.060" "FG.027" "FG.036" "FG.020"
#> [209] "FG.013" "FG.022" "FG.022" "FG.005" "FG.010" "FG.019" "FG.029" "FG.017"
#> [217] "FG.002" "FG.033" "FG.009" "FG.020" "FG.020" "FG.018" "FG.009" "FG.020"
#> [225] "FG.020"

## Perform feature grouping only on a subset of rows/features:
featureGroups(res) <- NA_character_
featureGroups(res)[40:80] <- "FG"
res <- groupFeatures(res, AbundanceSimilarityParam(transform = log2))
featureGroups(res)
#>   [1] NA       NA       NA       NA       NA       NA       NA       NA      
#>   [9] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [17] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [25] NA       NA       NA       NA       NA       NA       NA       NA      
#>  [33] NA       NA       NA       NA       NA       NA       NA       "FG.001"
#>  [41] "FG.002" "FG.003" "FG.008" "FG.009" "FG.006" "FG.007" "FG.006" "FG.001"
#>  [49] "FG.010" "FG.004" "FG.005" "FG.005" "FG.004" "FG.004" "FG.011" "FG.012"
#>  [57] "FG.013" "FG.001" "FG.014" "FG.015" "FG.003" "FG.016" "FG.017" "FG.003"
#>  [65] "FG.018" "FG.019" "FG.020" "FG.021" "FG.007" "FG.002" "FG.022" "FG.023"
#>  [73] "FG.024" "FG.025" "FG.026" "FG.027" "FG.028" "FG.029" "FG.030" "FG.031"
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
