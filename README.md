# RNAtools

An R package for processing RNA-seq data using multiple methods. Provides
consistent interfaces for normalisation, testing, visualisation etc.

Multiple packages exist for differential expression testing using RNA-seq data.
This package aims to provide a consistent interface to multiple packages that
allows them to be applied simultaneously to a single data set with minimal
setup. In addition consistent, easily comparable visualisations can be produced
using _ggplot2_, allowing for easy modification by the user.

## Installation

Several of the packages used by RNAtools are only available from
[Bioconductor](www.bioconductor.org) so we need to install those first.

```r
source("http://bioconductor.org/biocLite.R")
biocLite(c("edgeR", "DESeq", "DESeq2", "limma", "genefilter", "HTSFilter"))
```

RNAtools can now easily be installed from Github using the _devtools_ package.

```r
install.packages("devtools")
library("devtools")
install_github("lazappi/RNAtools")
```

If you are using Windows you will be prompted to also install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html). While
this is recommended RNAtools should install fine without it.

## Guide

See the vignette for a walkthrough of using RNAtools.

```r
library("RNAtools")
browseVignettes("RNAtools")
```
