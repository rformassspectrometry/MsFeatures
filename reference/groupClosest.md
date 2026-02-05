# Group values with differences below threshold

Group values with a difference between them being smaller than a user
defined threshold. This function uses the
[`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)
or the
[`groupSimilarityMatrixTree()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrixTree.md)
function to create groups with smallest differences between its members.
Differences between **all** members of one group are below the user
defined threshold `maxDiff`. This is a more stringent grouping than what
[`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md)
performs leading thus to smaller groups (with smaller differences
between its members).

## Usage

``` r
groupClosest(x, maxDiff = 1, FUN = groupSimilarityMatrix, ...)
```

## Arguments

- x:

  `numeric` of values that should be grouped.

- maxDiff:

  `numeric(1)` defining the threshold for difference between values in
  `x` to be grouped into the same group.

- FUN:

  supported similarity calculation function. Can be either
  [groupSimilarityMatrix](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)
  (the default) or
  [groupSimilarityMatrixTree](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrixTree.md),
  the latter being faster and less memory demanding.

- ...:

  additional parameters passed to `FUN`.

## Value

`integer` with the group assignment (values grouped together have the
same return value).

## See also

Other grouping operations:
[`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md),
[`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md),
[`groupSimilarityMatrixTree()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrixTree.md)

## Author

Johannes Rainer, Johan Lassen

## Examples

``` r

x <- c(1.1, 1.9, 2.2)
groupClosest(x)
#> [1] 2 1 1
## Although the difference between the 1st and 2nd element would be smaller
## than the threshold, they are not grouped because the difference between
## the 2nd and 3rd element is even smaller. The first element is also not
## put into the same group, because it has a difference > diffRt to the 3rd
## element.

x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7,
    9, 9.5, 15)

groupClosest(x)
#>  [1] 2 2 2 5 5 1 1 3 3 3 4 4 7 6 6 8
```
