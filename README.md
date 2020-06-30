# MsFeatures

Functionality for Mass Spectrometry Features. An MS feature is characterized by
its specific m/z and retention time range and its abundance in a set of samples
is supposed to represent the concentration of the same ion species of a specific
compound in each of these samples.

This package builds upon and extends the
`[Features](https://github.com/RforMassSpectrometry/Features)` and the
`[xcms](https://github.com/sneumann/xcms)` packages. The package specifically
provides methods to perform a grouping of (LC-MS) features defined by the latter
package that are based on the properties of these (i.e. similar retention time,
correlation of signal across samples or similar chromatographic peak shape).

