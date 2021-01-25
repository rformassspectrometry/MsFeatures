# MsFeatures 0.0

## MsFeatures 0.0.4

- Add `AbundanceSimilarityParam` and functionality to group features based on
  similarity of abundances across samples.

## MsFeatures 0.0.3

- Add functionality to group features based on similar retention times.
- Rename *old* `groupClosest` function to `groupConsecutive`.
- Add new `groupClosest` function which uses `groupSimilarityMatrix`.

## MsFeatures 0.0.2

- Add `groupDiagonalMatrix` function, remove `groupByCorrelation`.

## MsFeatures 0.0.1

- Add `groupClosest` and `groupByCorrelation` functions.
