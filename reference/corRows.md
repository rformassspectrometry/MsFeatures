# Correlate rows of a numeric matrix

`corRows` is a simple function to perform a pairwise correlation between
**rows** of a numeric matrix by calling
[`stats::cor()`](https://rdrr.io/r/stats/cor.html) on the transposed
input matrix `x`.

## Usage

``` r
corRows(
  x,
  y = NULL,
  use = "pairwise.complete.obs",
  method = c("pearson", "kendall", "spearman"),
  ...
)
```

## Arguments

- x:

  `numeric` `matrix`.

- y:

  not supported (ignored).

- use:

  see information for parameter `use` in
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html). Defaults to
  `use = "pairwise.complete.obs"`.

- method:

  see information for parameter `method` in
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

- ...:

  additional parameters (ignored).

## Value

`matrix` with correlation coefficients between rows in `x`.

## Author

Johannes Rainer

## Examples

``` r

## Define a simple numeric matrix
x <- rbind(
    c(12, 34, 231, 234, 9, 5, 7),
    c(900, 900, 800, 10, 12, 9, 4),
    c(25, 70, 400, 409, 15, 8, 4),
    c(12, 13, 14, 15, 16, 17, 18),
    c(14, 36, 240, 239, 12, 7, 8)
    )

corRows(x)
#>             [,1]        [,2]       [,3]       [,4]       [,5]
#> [1,]  1.00000000  0.09610556  0.9995093 -0.2120207  0.9999276
#> [2,]  0.09610556  1.00000000  0.1152920 -0.8822763  0.1001638
#> [3,]  0.99950928  0.11529202  1.0000000 -0.2362897  0.9993837
#> [4,] -0.21202072 -0.88227627 -0.2362897  1.0000000 -0.2136497
#> [5,]  0.99992759  0.10016376  0.9993837 -0.2136497  1.0000000

corRows(x, method = "spearman")
#>            [,1]       [,2]       [,3]       [,4]       [,5]
#> [1,]  1.0000000  0.4865062  0.9642857 -0.6071429  0.9642857
#> [2,]  0.4865062  1.0000000  0.5225437 -0.9549937  0.5585812
#> [3,]  0.9642857  0.5225437  1.0000000 -0.6428571  0.9285714
#> [4,] -0.6071429 -0.9549937 -0.6428571  1.0000000 -0.6428571
#> [5,]  0.9642857  0.5585812  0.9285714 -0.6428571  1.0000000
```
