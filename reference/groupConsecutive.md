# Grouping of sorted values into sets with smallest differences

`groupConsecutive` groups **sorted** values in `x` for which the
difference is smaller than `maxDiff`. As a result, the mean difference
between the groups will always be larger than `maxDiff`, but difference
between individual values within the same group (e.g. between the first
and last) can be larger `maxDiff`.

In detail, from the sorted `x`, the function starts from the smallest
value defining the first group as the one containing all values in `x`
with a difference to this first value which is `<= maxDiff`. The next
group is the defined based on the next larger value that is not part of
the first group and includes all values with a difference `<= maxDiff`
to this value. For values fulfilling this criteria but being already
part of a previous group, the differences to the mean value of the
current group and to the mean of previous groups are compared and values
are assigned to the group to which they have the smallest difference.

Example: values `1.1, 1.9, 2.2` should be grouped with a `maxDiff = 1`.
The first group is defined to include all values for which the
difference to the first value (`1.1`) is smaller `maxDiff`. Thus, the
first group is defined to contain values `1.1 and 1.9`. Then the next
group is defined based on the next larger value not part of any group,
`2.2`. This group contains values `1.9` and `2.2` with the value `1.9`
being already assigned to the first group. The difference between this
value `1.9` and the mean of the current group (`mean(c(1.9, 2.2)`) is
then compared to the difference of `1.9` to the mean value of the group
`1.9` is already part of (which is `mean(c(1.1, 1.9))`). Since the
difference to the second group is smaller, `1.9` is removed from the
first group and assigned to the second one.

## Usage

``` r
groupConsecutive(x, maxDiff = 1)
```

## Arguments

- x:

  `numeric` of values that should be grouped.

- maxDiff:

  `numeric(1)` defining the threshold for difference between values in
  `x` to be grouped into the same group.

## Value

`integer` with the group assignment (values grouped together have the
same return value).

## Note

The difference between consecutive (ordered) values within a defined
group is always `<= maxDiff`, but the difference between e.g. the first
and the last of the (ordered) values can be larger than `maxDiff`. See
[`groupClosest()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupClosest.md)
for a more stringent grouping function.

## See also

Other grouping operations:
[`groupClosest()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupClosest.md),
[`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md),
[`groupSimilarityMatrixTree()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrixTree.md)

## Author

Johannes Rainer

## Examples

``` r

## The example described above
x <- c(1.1, 1.9, 2.2)
groupConsecutive(x)
#> [1] 1 2 2

x <- c(1.1, 1.5, 1.7, 2.3, 2.7, 4.3, 4.4, 4.9, 5.2, 5.4, 5.8, 6, 7,
    9, 9.5, 15)

groupConsecutive(x)
#>  [1] 1 1 1 2 2 3 3 3 4 4 4 4 5 6 6 7
## value 5.2 was initially grouped with 4.3 (because their difference is
## smaller 1, but then re-grouped together with 5.4 because the difference
## between 5.4 (the next value outside the group of 4.3) and 5.2 is smaller
## than its difference to the mean value of the group for value 4.3

## Example for a case in which values are NOT grouped into the same group
## even if the difference between them is <= maxDiff
a <- c(4.9, 5.2, 5.4)
groupConsecutive(a, maxDiff = 0.3)
#> [1] 1 2 2
```
