# Group rows of a diagonal matrix using a threshold

This function groups elements (rows or columns) of a diagonal matrix,
such as a pairwise correlation matrix or similarity matrix, with a value
`>= threshold`. This creates clusters of elements in which **all**
elements have a value `>= threshold` with **any** other element in that
cluster. On a correlation matrix (such as created with `cor`) it will
generate small clusters of highly correlated elements. Note however that
single elements in one cluster could also have a correlation
`>= threshold` to another element in another cluster. The average
similarity to its own cluster will however be higher to that of the
other.

## Usage

``` r
groupSimilarityMatrix(x, threshold = 0.9, full = TRUE, ...)
```

## Arguments

- x:

  symmetrix `numeric` `matrix`.

- threshold:

  `numeric(1)` above which rows in `x` should be grouped.

- full:

  `logical(1)` whether the full matrix should be considered, or just the
  upper triangular matrix (including the diagonal).

- ...:

  ignored.

## Value

`integer` same length than `nrow(x)`, grouped elements (rows) defined by
the same value.

## Details

The algorithm is defined as follows:

- all pairs of values in `x` which are `>= threshold` are identified and
  sorted decreasingly.

- starting with the pair with the highest correlation, groups are
  defined:

- if none of the two is in a group, both are put into the same new
  group.

- if one of the two is already in a group, the other is put into the
  same group if **all** correlations of it to that group are
  `>= threshold` (and are not `NA`).

- if both are already in the same group nothing is done.

- if both are in different groups: an element is put into the group of
  the other if a) all correlations of it to members of the other's group
  are not `NA` and `>= threshold` **and** b) the average correlation to
  the other group is larger than the average correlation to its own
  group.

This ensures that groups are defined in which all elements have a
correlation `>= threshold` with each other and the correlation between
members of the same group is maximized.

## See also

Other grouping operations:
[`groupClosest()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupClosest.md),
[`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md),
[`groupSimilarityMatrixTree()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrixTree.md)

## Author

Johannes Rainer

## Examples

``` r

x <- rbind(
    c(1, 0.9, 0.6, 0.8, 0.5),
    c(0.9, 1, 0.7, 0.92, 0.8),
    c(0.6, 0.7, 1, 0.91, 0.7),
    c(0.8, 0.92, 0.91, 1, 0.9),
    c(0.5, 0.8, 0.7, 0.9, 1)
    )

groupSimilarityMatrix(x, threshold = 0.9)
#> [1] 2 1 3 1 4

groupSimilarityMatrix(x, threshold = 0.1)
#> [1] 1 1 1 1 1

## Add also a correlation between 3 and 2
x[2, 3] <- 0.9
x[3, 2] <- 0.9
x
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  1.0 0.90 0.60 0.80  0.5
#> [2,]  0.9 1.00 0.90 0.92  0.8
#> [3,]  0.6 0.90 1.00 0.91  0.7
#> [4,]  0.8 0.92 0.91 1.00  0.9
#> [5,]  0.5 0.80 0.70 0.90  1.0
groupSimilarityMatrix(x, threshold = 0.9)
#> [1] 2 1 1 1 3

## Add a higher correlation between 4 and 5
x[4, 5] <- 0.99
x[5, 4] <- 0.99
x
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  1.0 0.90 0.60 0.80 0.50
#> [2,]  0.9 1.00 0.90 0.92 0.80
#> [3,]  0.6 0.90 1.00 0.91 0.70
#> [4,]  0.8 0.92 0.91 1.00 0.99
#> [5,]  0.5 0.80 0.70 0.99 1.00
groupSimilarityMatrix(x, threshold = 0.9)
#> [1] 2 2 3 1 1

## Increase correlation between 2 and 3
x[2, 3] <- 0.92
x[3, 2] <- 0.92
x
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  1.0 0.90 0.60 0.80 0.50
#> [2,]  0.9 1.00 0.92 0.92 0.80
#> [3,]  0.6 0.92 1.00 0.91 0.70
#> [4,]  0.8 0.92 0.91 1.00 0.99
#> [5,]  0.5 0.80 0.70 0.99 1.00
groupSimilarityMatrix(x, threshold = 0.9) ## Don't break previous cluster!
#> [1] 3 2 2 1 1
```
