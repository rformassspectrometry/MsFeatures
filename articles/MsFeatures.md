# Grouping Mass Spectrometry Features

**Package**:
*[MsFeatures](https://bioconductor.org/packages/3.23/MsFeatures)*  
**Authors**: Johannes Rainer \[aut, cre\] (ORCID:
<https://orcid.org/0000-0002-6977-7147>), Johan Lassen \[ctb\] (ORCID:
<https://orcid.org/0000-0002-0452-566X>)  
**Last modified:** 2026-02-05 12:00:04.669913  
**Compiled**: Thu Feb 5 12:06:32 2026

## Introduction

Electrospray ionization (ESI) is commonly used in mass spectrometry
(MS)-based metabolomics to generate ions from the compounds to enable
their detection by the MS instrument. Ionization can generate different
ions (adducts) of the same original compound which are then reported as
separate *MS features* with different mass-to-charge ratios (m/z). To
reduce data set complexity (and to aid subsequent annotation steps) it
is advisable to group features which supposedly represent signal from
the same original compound into a single entity.

The `MsFeatures` package provides key concepts and functions for this
feature grouping. Methods are implemented for base R objects as well as
for Bioconductor’s `SummarizedExperiment` class. See also the
description of the [general grouping
concept](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.html)
on the package webpage for more information. Additional grouping
methodology is expected to be implemented in other R packages for data
objects with additional LC-MS related information, such as the
`XCMSnExp` object in the `xcms` package. The implementation for the
`SummarizedExperiment` provided in this package can be used as a
reference for these additional methodology.

After definition of the feature groups, the
*[QFeatures](https://bioconductor.org/packages/3.23/QFeatures)* package
could be used to aggregate their abundances into a single signal.

## Installation

The package can be installed with the `BiocManager` package. To install
`BiocManager` use `install.packages("BiocManager")` and, after that,
`BiocManager::install("MsFeatures")` to install this package.

## Mass Spectrometry Feature Grouping

Features from the same originating compound inherit its characteristics
including its retention time (for LC or GC-MS experiments) and
abundance/intensity. For the latter it is expected that features from
the same compound have the same pattern of feature values/abundances
across samples.

The `MsFeatures` package defines the `groupFeatures` method to perform
MS feature grouping based on the provided input data and a parameter
object which selects and defines the feature grouping algorithm. This
algorithm is supposed to assign individual features to a (single)
feature group. Currently two feature grouping approaches are
implemented:

- `SimilarRtimeParam`: group features based on similar retention times.
- `AbundanceSimilarityParam`: group features based on similar feature
  values/abundances across samples.

Additional algorithms, e.g. by considering also differences in features’
m/z values matching expected ions/adducts or isotopes, may be
implemented in future in this or other packages.

In this document we demonstrate the feature grouping functionality on a
simple toy data set used also in the
*[xcms](https://bioconductor.org/packages/3.23/xcms)* package with the
raw data being provided in the `faahKO` data package. This data set
consists of samples from 4 mice with knock-out of the fatty acid amide
hydrolase (FAAH) and 4 wild type mice. Pre-processing of this data set
is described in detail in the *xcms* vignette of the `xcms` package.
Below we load all required packages and the result from this
pre-processing which is provided as a `SummarizedExperiment` within this
package and can be loaded with `data(se)`.

``` r

library(MsFeatures)
library(SummarizedExperiment)

data("se")
```

Before performing the feature grouping we inspect the result object.
Feature properties and definitions can be accessed with `rowData`, the
feature abundances with `assay`.

``` r

rowData(se)
```

    ## DataFrame with 225 rows and 11 columns
    ##           mzmed     mzmin     mzmax     rtmed     rtmin     rtmax    npeaks
    ##       <numeric> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
    ## FT001     200.1     200.1     200.1   2901.63   2880.73   2922.53         2
    ## FT002     205.0     205.0     205.0   2789.39   2782.30   2795.36         8
    ## FT003     206.0     206.0     206.0   2788.73   2780.73   2792.86         7
    ## FT004     207.1     207.1     207.1   2718.12   2713.21   2726.70         7
    ## FT005     219.1     219.1     219.1   2518.82   2517.40   2520.81         3
    ## ...         ...       ...       ...       ...       ...       ...       ...
    ## FT221    591.30     591.3     591.3   3005.03   2992.87   3006.05         5
    ## FT222    592.15     592.1     592.3   3022.11   2981.91   3107.59         6
    ## FT223    594.20     594.2     594.2   3418.16   3359.10   3427.90         3
    ## FT224    595.25     595.2     595.3   3010.15   2992.87   3013.77         6
    ## FT225    596.20     596.2     596.2   2997.91   2992.87   3002.95         2
    ##              KO        WT            peakidx  ms_level
    ##       <numeric> <numeric>             <list> <integer>
    ## FT001         2         0  287, 679,1577,...         1
    ## FT002         4         4     47,272,542,...         1
    ## FT003         3         4     32,259,663,...         1
    ## FT004         4         3     19,249,525,...         1
    ## FT005         1         2  639, 788,1376,...         1
    ## ...         ...       ...                ...       ...
    ## FT221         2         3    349,684,880,...         1
    ## FT222         1         3     86,861,862,...         1
    ## FT223         1         2  604, 985,1543,...         1
    ## FT224         2         3     67,353,876,...         1
    ## FT225         0         2  866,1447,1643,...         1

``` r

head(assay(se))
```

    ##        ko15.CDF   ko16.CDF   ko21.CDF  ko22.CDF  wt15.CDF  wt16.CDF  wt21.CDF
    ## FT001  159738.1  506848.88  113441.08  169955.6  216096.6  145509.7  230477.9
    ## FT002 1924712.0 1757150.96 1383416.72 1180288.2 2129885.1 1634342.0 1623589.2
    ## FT003  213659.3  289500.67  162897.19  178285.7  253825.6  241844.4  240606.0
    ## FT004  349011.5  451863.66  343897.76  208002.8  364609.8  360908.9  223322.5
    ## FT005  135978.5   25524.79   71530.84  107348.5  223951.8  134398.9  190203.8
    ## FT006  286221.4  289908.23  164008.97  149097.6  255697.7  311296.8  366441.5
    ##         wt22.CDF
    ## FT001  140551.30
    ## FT002 1354004.93
    ## FT003  185399.47
    ## FT004  221937.53
    ## FT005   84772.92
    ## FT006  271128.02

Columns `"mzmed"` and `"rtmed"` in the object’s `rowData` provide the
m/z and retention time which characterizes each feature. In total 225
features are available in the present data set, with many of them most
likely representing signal from different ions of the same compound. We
aim to identify these based on the following assumptions of the LC-MS
data:

- Features (ions) of the same compound should have similar retention
  time.
- The abundance of features (ions) of the same compound should have a
  similar pattern across samples, i.e. if a compound is highly
  concentrated in one sample and low in another, all ions from it also
  should follow the same pattern.

As detailed in the [general grouping
concept](https://rformassspectrometry.github.io/MsFeatures/reference/groupFeatures.html),
the feature grouping implemented in `MsFeatures` is by default intended
to be used as a stepwise approach in which each `groupFeatures` call
further sub-groups (and thus refines) previously defined feature groups.
This enables to either use a single algorithm for the feature grouping
or to build a feature grouping *pipeline* by combining different
algorithms. In our example we perform first a initial grouping of
features based on similar retention time and subsequently further refine
these feature groups by requiring also similarity of feature values
across samples.

Note that it would also be possible to perform the grouping only on a
subset of features instead of the full data set. An example is provided
in the last section of this vignette.

### Grouping of features by similar retention time

The most intuitive and simple way to group LC-MS features is based on
their retention times: ionization of the compounds happens after the LC
and thus all ions from the same compound should have the same retention
time. The plot below shows the retention times (and m/z) of all features
from the present data set.

``` r

plot(rowData(se)$rtmed, rowData(se)$mzmed,
     xlab = "retention time", ylab = "m/z", main = "features",
     col = "#00000060")
grid()
```

![Plot of retention times and m/z for all features in the data
set.](MsFeatures_files/figure-html/feature-rt-mz-plot-1.png)

Plot of retention times and m/z for all features in the data set.

As we can see there are several features with a similar retention time,
especially for lower retention times. By using `groupFeatures` with the
`SimilarRtimeParam` we can next group features if their difference in
retention time is below a certain threshold. This approach will however
not only group features representing ions from the same compound
together, but also features from different, but co-eluting compounds
(i.e. different compounds with the same retention time). Thus feature
groups defined by this algorithm should be further *refined* based on
another feature property to reduce false positives.

For the present example, we group features with a maximal difference in
retention time of 10 seconds into a feature group. We also have to
specify the column in the object’s `rowData` which contains the
retention times for the features.

``` r

se <- groupFeatures(se, param = SimilarRtimeParam(10), rtime = "rtmed")
```

The `groupFeatures` call on the `SummarizedExperiment` added the results
of the grouping into a new column called `"feature_group"` in the
object’s `rowData`. This column can also be directly accessed with the
`featureGroups` function. Below we print the number of features for each
feature grouped defined by the `SimilarRtimeParam` approach.

``` r

table(featureGroups(se))
```

    ## 
    ## FG.001 FG.002 FG.003 FG.004 FG.005 FG.006 FG.007 FG.008 FG.009 FG.010 FG.011 
    ##      3      3      3      3      2      4      5      6      4      2      5 
    ## FG.012 FG.013 FG.014 FG.015 FG.016 FG.017 FG.018 FG.019 FG.020 FG.021 FG.022 
    ##      3      4      3      5      3      3      5      3      3      3      3 
    ## FG.023 FG.024 FG.025 FG.026 FG.027 FG.028 FG.029 FG.030 FG.031 FG.032 FG.033 
    ##      3      3      6      3      3      3      3      2      3      3      4 
    ## FG.034 FG.035 FG.036 FG.037 FG.038 FG.039 FG.040 FG.041 FG.042 FG.043 FG.044 
    ##      3      2      2      3      2      2      4      2      2      2      3 
    ## FG.045 FG.046 FG.047 FG.048 FG.049 FG.050 FG.051 FG.052 FG.053 FG.054 FG.055 
    ##      4      2      3      3      3      2      2      3      4      2      3 
    ## FG.056 FG.057 FG.058 FG.059 FG.060 FG.061 FG.062 FG.063 FG.064 FG.065 FG.066 
    ##      2      2      2      3      2      3      2      2      2      3      2 
    ## FG.067 FG.068 FG.069 FG.070 FG.071 FG.072 FG.073 FG.074 FG.075 FG.076 FG.077 
    ##      2      3      2      2      2      3      2      2      1      1      1 
    ## FG.078 FG.079 FG.080 FG.081 FG.082 FG.083 FG.084 
    ##      1      1      1      1      1      1      1

We also calculate the mean retention time for all the feature groups and
order them increasingly.

``` r

split(rowData(se)$rtmed, featureGroups(se)) |>
vapply(FUN = mean, numeric(1)) |>
sort()
```

    ##   FG.005   FG.012   FG.062   FG.076   FG.069   FG.007   FG.018   FG.027 
    ## 2506.301 2511.821 2519.745 2595.503 2623.786 2684.671 2686.835 2688.495 
    ##   FG.064   FG.057   FG.066   FG.074   FG.003   FG.023   FG.004   FG.001 
    ## 2694.089 2718.799 2731.648 2751.041 2787.291 2788.022 2788.625 2788.949 
    ##   FG.011   FG.048   FG.014   FG.075   FG.084   FG.080   FG.060   FG.021 
    ## 2790.232 2793.065 2799.149 2832.560 2860.496 2871.735 2881.899 2899.457 
    ##   FG.063   FG.033   FG.071   FG.079   FG.050   FG.051   FG.039   FG.010 
    ## 2915.107 2924.922 2936.602 2953.627 2998.303 3005.543 3009.864 3011.159 
    ##   FG.040   FG.019   FG.006   FG.067   FG.036   FG.078   FG.037   FG.013 
    ## 3015.635 3022.568 3027.254 3038.254 3044.682 3057.242 3075.361 3079.912 
    ##   FG.017   FG.083   FG.044   FG.073   FG.070   FG.008   FG.068   FG.082 
    ## 3088.353 3114.721 3128.047 3142.710 3161.133 3170.219 3184.355 3203.298 
    ##   FG.029   FG.009   FG.055   FG.002   FG.042   FG.016   FG.035   FG.054 
    ## 3217.009 3226.662 3242.698 3258.465 3261.208 3263.141 3265.660 3284.004 
    ##   FG.043   FG.038   FG.052   FG.077   FG.056   FG.028   FG.061   FG.022 
    ## 3289.081 3295.654 3300.775 3311.799 3321.354 3324.717 3335.200 3353.939 
    ##   FG.031   FG.081   FG.041   FG.026   FG.015   FG.030   FG.047   FG.034 
    ## 3358.597 3372.908 3383.014 3395.552 3405.585 3407.192 3410.212 3417.383 
    ##   FG.024   FG.049   FG.059   FG.058   FG.045   FG.046   FG.065   FG.053 
    ## 3422.675 3428.323 3435.389 3442.423 3447.411 3457.582 3465.580 3472.813 
    ##   FG.032   FG.020   FG.025   FG.072 
    ## 3478.609 3483.910 3497.730 3510.664

Note that the differences in retention times between the feature groups
can be smaller than the used cut-off (10 seconds in our case). If we
were not happy with this feature grouping and would like to repeat it we
would need to drop the `"feature_group"` column in the object’s
`rowData` with `rowData(se)$feature_group <- NULL` and repeat the
feature grouping with different settings. This is required, because by
default `groupFeatures` will *refine* previous feature grouping results
but not overwrite them.

As stated above, this initial grouping on retention times put features
from the same, but also from different co-eluting compounds into the
same feature group. We thus next refine the feature groups requiring
also feature abundances across samples to be correlated.

### Grouping of features by abundance correlation across samples

Features representing ions of the same compound are expected to have
correlated feature values (intensities, abundances) across samples.
`groupFeatures` with `AbundanceSimilarityParam` allows to group features
with similar abundance patterns. This approach performs a pairwise
similarity calculation and puts features with a similarity
`>= threshold` into the same feature group. By calling this function on
the previous result object the initial feature groups will be refined,
by eventually splitting them based on the (missing) correlation of
feature abundances.

We below evaluate the correlation between individual features indicating
also the previously defined feature groups.

``` r

library(pheatmap)
fvals <- log2(assay(se))

cormat <- cor(t(fvals), use = "pairwise.complete.obs")
ann <- data.frame(fgroup = featureGroups(se))
rownames(ann) <- rownames(cormat)

res <- pheatmap(cormat, annotation_row = ann, cluster_rows = TRUE,
                cluster_cols = TRUE)
```

![Correlation of features based on their
abundances.](MsFeatures_files/figure-html/abundance-correlation-heatmap-1.png)

Correlation of features based on their abundances.

As expected, the clustering based on the feature abundances does not
perfectly match the retention time-based feature grouping. Many features
grouped based on retention time have a low, or even negative correlation
of feature abundances across samples hence most likely representing
features from different, but co-eluting compounds. On the other hand,
many features are highly correlated, but have a different retention time
and can thus also not represent signal from ions of the same compound.
Thus, each single approach has its drawbacks, but combination them can
reduce the number of wrongly grouped features.

We thus next perform the feature grouping with
`AbundanceSimilarityParam` on the result object to refine the retention
time-based feature groups. The approach can be further customized by
providing a function to calculate feature similarities with parameter
`simFun` (by default `cor` will be used to calculate similarities using
Pearson’s correlation). Parameter `transform` allows to specify a
function to transform feature abundances prior similarity calculation.
By default the feature values are taken *as-is*, but below we use
`transform = log2` to perform the calculations in log2 scale. With
`threshold = 0.7` we ensure that only features with a correlation
coefficient `>= 0.7` are assigned to the same feature group. Finally,
parameter `i` would allow to specify the assay in the
`SummarizedExperiment` that contains the feature abundances on which
similarities should be calculated. See the `AbundanceSimilarityParam`
help page for a full listing of the parameters and more details.

``` r

se <- groupFeatures(se, AbundanceSimilarityParam(threshold = 0.7,
                                                 transform = log2), i = 1)
table(featureGroups(se))
```

    ## 
    ## FG.001.001 FG.001.002 FG.002.001 FG.002.002 FG.002.003 FG.003.001 FG.003.002 
    ##          2          1          1          1          1          2          1 
    ## FG.004.001 FG.005.001 FG.006.001 FG.006.002 FG.007.001 FG.007.002 FG.007.003 
    ##          3          2          3          1          2          1          1 
    ## FG.007.004 FG.008.001 FG.008.002 FG.008.003 FG.008.004 FG.009.001 FG.009.002 
    ##          1          3          1          1          1          3          1 
    ## FG.010.001 FG.011.001 FG.011.002 FG.011.003 FG.012.001 FG.013.001 FG.013.002 
    ##          2          2          2          1          3          2          1 
    ## FG.013.003 FG.014.001 FG.014.002 FG.015.001 FG.015.002 FG.015.003 FG.015.004 
    ##          1          2          1          2          1          1          1 
    ## FG.016.001 FG.016.002 FG.017.001 FG.017.002 FG.018.001 FG.018.002 FG.018.003 
    ##          2          1          2          1          2          1          1 
    ## FG.018.004 FG.019.001 FG.019.002 FG.019.003 FG.020.001 FG.020.002 FG.021.001 
    ##          1          1          1          1          2          1          1 
    ## FG.021.002 FG.021.003 FG.022.001 FG.022.002 FG.022.003 FG.023.001 FG.023.002 
    ##          1          1          1          1          1          2          1 
    ## FG.024.001 FG.024.002 FG.024.003 FG.025.001 FG.025.002 FG.025.003 FG.025.004 
    ##          1          1          1          2          1          1          1 
    ## FG.025.005 FG.026.001 FG.026.002 FG.026.003 FG.027.001 FG.027.002 FG.027.003 
    ##          1          1          1          1          1          1          1 
    ## FG.028.001 FG.028.002 FG.029.001 FG.029.002 FG.030.001 FG.030.002 FG.031.001 
    ##          2          1          2          1          1          1          3 
    ## FG.032.001 FG.032.002 FG.032.003 FG.033.001 FG.033.002 FG.033.003 FG.033.004 
    ##          1          1          1          1          1          1          1 
    ## FG.034.001 FG.034.002 FG.034.003 FG.035.001 FG.035.002 FG.036.001 FG.036.002 
    ##          1          1          1          1          1          1          1 
    ## FG.037.001 FG.037.002 FG.038.001 FG.039.001 FG.040.001 FG.040.002 FG.041.001 
    ##          2          1          2          2          2          2          1 
    ## FG.041.002 FG.042.001 FG.042.002 FG.043.001 FG.044.001 FG.044.002 FG.045.001 
    ##          1          1          1          2          2          1          1 
    ## FG.045.002 FG.045.003 FG.045.004 FG.046.001 FG.047.001 FG.047.002 FG.047.003 
    ##          1          1          1          2          1          1          1 
    ## FG.048.001 FG.048.002 FG.048.003 FG.049.001 FG.049.002 FG.050.001 FG.050.002 
    ##          1          1          1          2          1          1          1 
    ## FG.051.001 FG.052.001 FG.052.002 FG.053.001 FG.053.002 FG.053.003 FG.054.001 
    ##          2          2          1          2          1          1          1 
    ## FG.054.002 FG.055.001 FG.055.002 FG.055.003 FG.056.001 FG.056.002 FG.057.001 
    ##          1          1          1          1          1          1          2 
    ## FG.058.001 FG.058.002 FG.059.001 FG.059.002 FG.059.003 FG.060.001 FG.060.002 
    ##          1          1          1          1          1          1          1 
    ## FG.061.001 FG.061.002 FG.062.001 FG.062.002 FG.063.001 FG.064.001 FG.064.002 
    ##          2          1          1          1          2          1          1 
    ## FG.065.001 FG.065.002 FG.066.001 FG.066.002 FG.067.001 FG.067.002 FG.068.001 
    ##          2          1          1          1          1          1          3 
    ## FG.069.001 FG.070.001 FG.070.002 FG.071.001 FG.071.002 FG.072.001 FG.072.002 
    ##          2          1          1          1          1          2          1 
    ## FG.073.001 FG.073.002 FG.074.001 FG.075.001 FG.076.001 FG.077.001 FG.078.001 
    ##          1          1          2          1          1          1          1 
    ## FG.079.001 FG.080.001 FG.081.001 FG.082.001 FG.083.001 FG.084.001 
    ##          1          1          1          1          1          1

Many of the larger retention time-based feature groups have been
splitted into two or more sub-groups based on the correlation of their
feature abundances. We evaluate this for one specific feature group
`"FG.003"` by plotting their pairwise correlation.

``` r

fts <- grep("FG.003", featureGroups(se))
pairs(t(fvals[fts, ]), gap = 0.1, main = "FG.003")
```

![Pairwise correlation plot for features initially grouped into the
feature group
FG.003.](MsFeatures_files/figure-html/abundance-correlation-fg003-1.png)

Pairwise correlation plot for features initially grouped into the
feature group FG.003.

A high correlation can be observed between *FT035* and *FT051* while
they are not correlated with feature *FT013*. We next evaluate the
feature grouping for another example.

``` r

fts <- grep("FG.008", featureGroups(se))
pairs(t(fvals[fts, ]), gap = 0.1, main = "FG.008")
```

![Pairwise correlation plot for features initially grouped into the
feature group
FG.008.](MsFeatures_files/figure-html/abundance-correlation-fg008-1.png)

Pairwise correlation plot for features initially grouped into the
feature group FG.008.

The results are less clear than for the previous example, still, some
features seem to be correlated with each other while others are not.
Generally, the abundance correlation approach in this data set suffers
from the low number of sample (8 in total). Also, the approach works
better for features with a high variance (biologically or technically)
across samples.

The table below lists the retention time, m/z and group assignment for
these features.

``` r

tmp <- as.data.frame(rowData(se)[fts, c("rtmed", "mzmed", "feature_group")])
tmp <- tmp[order(tmp$feature_group), ]
knitr::kable(tmp)
```

|       |    rtmed | mzmed | feature_group |
|:------|---------:|------:|:--------------|
| FT097 | 3173.240 | 385.1 | FG.008.001    |
| FT163 | 3170.171 | 502.1 | FG.008.001    |
| FT165 | 3170.171 | 503.1 | FG.008.001    |
| FT074 | 3167.106 | 361.1 | FG.008.002    |
| FT077 | 3170.675 | 362.1 | FG.008.003    |
| FT100 | 3169.951 | 386.1 | FG.008.004    |

The difference in m/z between features *FT163* and *FT165*, both being
assigned to the same feature group, is ~ 1 suggesting that one of the
two is in fact a (C13) isotope of the other feature.

### Performing feature grouping on a subset of features

Sometimes it might not be needed or required to perform the feature
grouping on the full data set but only on a subset of *interesting*
features (i.e. those with significant differences in feature abundances
between sample groups). This has also the advantage of a larger range of
feature values across samples which supports the abundance
similarity-based feature grouping.

Feature grouping on a subset of features can be performed by manually
assigning all features of interest to an initial feature group and
setting the feature group for all other features to `NA`. As an example
we perform below the feature grouping only features 30-60.

``` r

featureGroups(se) <- NA_character_
featureGroups(se)[30:60] <- "FG"

se <- groupFeatures(se, SimilarRtimeParam(10), rtime = "rtmed")
```

This did not *refine* this initial, manually specified feature group by
the retention time-based grouping. Features with `NA` value in their
feature group column are skipped. As a result we get the following
grouping:

``` r

featureGroups(se)
```

    ##   [1] NA       NA       NA       NA       NA       NA       NA       NA      
    ##   [9] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [17] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [25] NA       NA       NA       NA       NA       "FG.009" "FG.005" "FG.010"
    ##  [33] "FG.002" "FG.011" "FG.001" "FG.007" "FG.004" "FG.004" "FG.005" "FG.002"
    ##  [41] "FG.005" "FG.012" "FG.008" "FG.002" "FG.007" "FG.013" "FG.014" "FG.008"
    ##  [49] "FG.015" "FG.004" "FG.001" "FG.002" "FG.003" "FG.003" "FG.016" "FG.002"
    ##  [57] "FG.006" "FG.017" "FG.006" "FG.018" NA       NA       NA       NA      
    ##  [65] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [73] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [81] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [89] NA       NA       NA       NA       NA       NA       NA       NA      
    ##  [97] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [105] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [113] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [121] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [129] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [137] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [145] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [153] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [161] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [169] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [177] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [185] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [193] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [201] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [209] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [217] NA       NA       NA       NA       NA       NA       NA       NA      
    ## [225] NA

## Session information

    ## R Under development (unstable) (2026-02-01 r89366)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] pheatmap_1.0.13             SummarizedExperiment_1.41.0
    ##  [3] Biobase_2.71.0              GenomicRanges_1.63.1       
    ##  [5] Seqinfo_1.1.0               IRanges_2.45.0             
    ##  [7] S4Vectors_0.49.0            BiocGenerics_0.57.0        
    ##  [9] generics_0.1.4              MatrixGenerics_1.23.0      
    ## [11] matrixStats_1.5.0           MsFeatures_1.19.0          
    ## [13] BiocStyle_2.39.0           
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] sass_0.4.10         MsCoreUtils_1.23.2  SparseArray_1.11.10
    ##  [4] lattice_0.22-7      digest_0.6.39       evaluate_1.0.5     
    ##  [7] grid_4.6.0          RColorBrewer_1.1-3  bookdown_0.46      
    ## [10] fastmap_1.2.0       jsonlite_2.0.0      Matrix_1.7-4       
    ## [13] ProtGenerics_1.43.0 BiocManager_1.30.27 scales_1.4.0       
    ## [16] textshaping_1.0.4   jquerylib_0.1.4     abind_1.4-8        
    ## [19] cli_3.6.5           rlang_1.1.7         XVector_0.51.0     
    ## [22] cachem_1.1.0        DelayedArray_0.37.0 yaml_2.3.12        
    ## [25] otel_0.2.0          S4Arrays_1.11.1     tools_4.6.0        
    ## [28] R6_2.6.1            lifecycle_1.0.5     fs_1.6.6           
    ## [31] htmlwidgets_1.6.4   clue_0.3-66         MASS_7.3-65        
    ## [34] ragg_1.5.0          cluster_2.1.8.2     desc_1.4.3         
    ## [37] pkgdown_2.2.0.9000  bslib_0.10.0        gtable_0.3.6       
    ## [40] glue_1.8.0          systemfonts_1.3.1   xfun_0.56          
    ## [43] knitr_1.51          farver_2.1.2        htmltools_0.5.9    
    ## [46] rmarkdown_2.30      compiler_4.6.0
