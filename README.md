# VBphenoR
Variational Bayes for Latent Patient Phenotypes in Electronic Health Records (EHR)

# Overview

`VBphenoR` is an R package for discovering latent patient phenotypes from realistically large EHR data using Bayesian statistics. 
In order to computationally support EHR data, we employ variational Bayes (VB). Currently, it supports latent class discovery
using VB Gaussian Mixture Model implemented with Coordinate-ascent Variational Inference (CAVI) and VB Logistic Regression for
biomarker levels shifted for the latent phenotype. Please note this package is still under development.

# Installation

Prior to analyzing your EHR data, the R package needs to be installed. The
easiest way to install `VBphenoR` is through CRAN:

``` r
install.packages("VBphenoR")
```

`VBphenoR` can also be downloaded and installed via GitHub. This is most useful for downloading a *specific* version of the package (which
can be found at <https://github.com/buckleybrian/VBphenoR/releases>):

``` r
devtools::install_github("buckleybrian/VBphenoR@vx.xx.x")

# or 

devtools::install_version("VBphenoR", version = "x.x.x", repos = "http://cran.us.r-project.org")
```

The latest developmental version of the package can be downloaded and
installed by running:

``` r
devtools::install_github("buckleybrian/VBphenoR", build_vignettes = TRUE, build_manual=TRUE)
```

After successful installation, the package must be loaded into the
working space:

``` r
library(VBphenoR)
```

# Usage

See the [vignettes](https://buckleybrian.github.io/VBphenoR/articles/VBpheno.html) for usage instructions and example.


# License

`VBphenoR` is available under the open source [MIT License](https://www.r-project.org/Licenses/MIT)

