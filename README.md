# `MsFeatures` - Functionality for Mass Spectrometry Features

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsFeatures/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsFeatures/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov](https://codecov.io/gh/rformassspectrometry/MsFeatures/branch/main/graph/badge.svg?token=zUtxxzOqMT)](https://codecov.io/gh/rformassspectrometry/MsFeatures)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)
[![years in bioc](http://bioconductor.org/shields/years-in-bioc/MsFeatures.svg)](https://bioconductor.org/packages/release/bioc/html/MsFeatures.html)
[![Ranking by downloads](http://bioconductor.org/shields/downloads/release/MsFeatures.svg)](https://bioconductor.org/packages/stats/bioc/MsFeatures/)
[![build release](http://bioconductor.org/shields/build/release/bioc/MsFeatures.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/MsFeatures/)
[![build devel](http://bioconductor.org/shields/build/devel/bioc/MsFeatures.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/MsFeatures/)

The `MsFeatures` package defines functionality for Mass Spectrometry (MS)
features. These features are characterized by their specific mass-to-charge
ration m/z and eventually a retention time rt (for LC-MS features).

This package defines thus basic concepts and functions for such MS features
which includes grouping of (LC-MS) features based on their properties, such as
retention time (coeluting compounds), or correlation of signals across samples.
Provided functions are for base data types and core MS data representations
(such as `Chromatogram` objects from the
[`MSnbase`](https://github.com/lgatto/MSnbase) package). Implementation of
feature grouping methods for more specific data objects, such as LC-MS
preprocessing results stored in `XCMSnExp` objects from the
[`xcms`](https://github.com/sneumann/xcms) package are implemented in the
respective packages.

See the package [homepage](https://rformassspectrometry.github.io/MsFeatures)
for more information.


# Installation

The package can be installed with

```r
install.packages("BiocManager")
BiocManager::install("MsFeatures")
```


# Contributions

Contributions are highly welcome and should follow the [contribution
guidelines](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html#contributions).
Also, please check the coding style guidelines in the [RforMassSpectrometry
vignette](https://rformassspectrometry.github.io/RforMassSpectrometry/articles/RforMassSpectrometry.html).
