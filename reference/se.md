# Quantified LC-MS preprocessing result test data

The `se` variable is a
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object representing the results from a `xcms`-based pre-processing of an
LC-MS untargeted metabolomics data set. The raw data files are provided
in the `faahKO` package. The pre-processing of this data set is
described in detail in the *xcms* vignette of the `xcms` package. This
object was created from the `XCMSnExp` result object with the `quantify`
method.

## Examples

``` r

## Load the data
data(se)

library(SummarizedExperiment)

## Access row (feature) data
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
```
