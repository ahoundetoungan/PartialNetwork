# An R Package for Estimating Peer Effects Using Partial Network Data

<!-- badges: start -->
  [![Lifecycle: stable](https://img.shields.io/badge/Lifecycle-Stable-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![R-CMD-check](https://github.com/ahoundetoungan/PartialNetwork/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/ahoundetoungan/PartialNetwork/actions/workflows/R-CMD-check.yml)

  [![R-universe](https://ahoundetoungan.r-universe.dev/badges/PartialNetwork)](https://ahoundetoungan.r-universe.dev/PartialNetwork)
  [![CRAN](https://www.r-pkg.org/badges/version/PartialNetwork)](https://CRAN.R-project.org/package=PartialNetwork)
  [![CRAN Downloads](https://img.shields.io/endpoint?url=https://ahoundetoungan.github.io/cranlogs/badges/PartialNetwork.json)](https://cran.r-project.org/package=PartialNetwork)

  [![Vignette](https://img.shields.io/badge/Vignette-blue.svg)](https://nbviewer.org/github/ahoundetoungan/PartialNetwork/blob/master/doc/PartialNetwork_vignette.pdf)

<!-- badges: end -->


The **PartialNetwork** package includes all functions for the replication of the results in Boucher and Houndetoungan (2023). The exact replication codes are located in the folder [**Results**](https://github.com/ahoundetoungan/PartialNetwork/tree/master/Results). Below, we also provide detailed examples on how to use the estimators described in the paper.

## Installation
### CRAN version
**PartialNetwork** can be directly installed from CRAN.
```R
install.packages("PartialNetwork")
```

### GitHub version
It may be possible that we updated the package without submitting the new version to CRAN. The latest version (*but not necessary stable*) of **PartialNetwork** can be installed from this GitHub repos.
```R
remotes::install_github("ahoundetoungan/PartialNetwork", build_vignettes = TRUE)
```
## How to use the `PartialNetwork` package?
See our [vignette in pdf](https://nbviewer.jupyter.org/github/ahoundetoungan/PartialNetwork/blob/master/doc/PartialNetwork_vignette.pdf).
