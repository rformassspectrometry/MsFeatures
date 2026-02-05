# Identifies groups with intragroup pairwise distances less than maxDiff

Searches through a complete clustering tree for branches with internal
distances less than `maxDiff`.

## Usage

``` r
groupSimilarityMatrixTree(dists, maxDiff = 1)
```

## Arguments

- dists:

  a *distance* object such as returned by a call to
  [`dist()`](https://rdrr.io/r/stats/dist.html) with the pairwise
  retention time distances.

- maxDiff:

  `numeric(1)` with the maximum pairwise distance between members of a
  group.

## Value

`integers` representing groups of peaks based on retention times. Should
be of the same length as the feature definitions.

## Details

The search is performed top-down so that every node is evaluated for
qualifying as a group. If the node does not qualify the algorithm
iterates to the next node in the tree. If the node qualifies, all
members are added as a group to the *groups vector* and all following
child nodes are skipped.

## See also

Other grouping operations:
[`groupClosest()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupClosest.md),
[`groupConsecutive()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupConsecutive.md),
[`groupSimilarityMatrix()`](https://rformassspectrometry.github.io/MsFeatures/reference/groupSimilarityMatrix.md)

## Author

Johan Lassen
