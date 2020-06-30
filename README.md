# `MsFeatures` - Functionality for Mass Spectrometry Features

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsFeatures/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsFeatures/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/MsFeatures/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/MsFeatures?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

Functionality for Mass Spectrometry Features. An MS feature is characterized by
its specific m/z and retention time range and its abundance in a set of samples
is supposed to represent the concentration of the same ion species of a specific
compound in each of these samples.

This package builds upon and extends the
[`Features`](https://github.com/RforMassSpectrometry/Features) and the
[`xcms`](https://github.com/sneumann/xcms) packages. The package specifically
provides methods to perform a grouping of (LC-MS) features defined by the latter
package that are based on the properties of these (i.e. similar retention time,
correlation of signal across samples or similar chromatographic peak shape).

