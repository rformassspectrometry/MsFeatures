# `MsFeatures` - Functionality for Mass Spectrometry Features

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![R-CMD-check-bioc](https://github.com/RforMassSpectrometry/MsFeatures/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/RforMassSpectrometry/MsFeatures/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](http://codecov.io/github/rformassspectrometry/MsFeatures/coverage.svg?branch=master)](http://codecov.io/github/rformassspectrometry/MsFeatures?branch=master)
[![license](https://img.shields.io/badge/license-Artistic--2.0-brightgreen.svg)](https://opensource.org/licenses/Artistic-2.0)

Functionality for Mass Spectrometry Features. An MS feature is characterized by
its specific m/z and retention time range and its abundance in a set of samples
is supposed to represent the concentration of the same ion species of a specific
compound in each of these samples.

This package defines basic concepts and functions to work with MS features such
as grouping of features based on their properties. Provided functions are for
base data types and core MS data representations (such as `Chromatogram` objects
from the [`MSnbase`](https://github.com/lgatto/MSnbase) package). Implementation
of feature grouping methods for more specific data objects, such as LC-MS
preprocessing results stored in `XCMSnExp` objects from the
[`xcms`](https://github.com/sneumann/xcms) package are implemented in the
respective packages.
