
# mzAnnotation

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/jasenfinch/mzAnnotation/workflows/R-CMD-check/badge.svg)](https://github.com/jasenfinch/mzAnnotation/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/jasenfinch/mzAnnotation/devel.svg)](https://codecov.io/github/jasenfinch/mzAnnotation?branch=devel)
[![DOI](https://zenodo.org/badge/33083554.svg)](https://zenodo.org/badge/latestdoi/33083554)
<!-- badges: end -->

An R package containing tools for putative annotation of accurate m/z
from electrospray ionisation mass spectrometry data.

### Installation

Due to the dependency of `mzAnnotation` on the [Open Babel](http://openbabel.org/wiki/Main_Page), installation can vary depending on your OS.

#### Linux

The package can be installed from source on Linux by first installing the Open Babel library headers.
On Debian flavour distributions, these are available by installing the `libopenbabel-dev` package.

`mzAnnotation` can then be installed directly from this repository using:

``` r
devtools::install_github('aberHRML/mzAnnotation')
```

#### Windows

It is recommended to install the package using pre-compiled binaries using the following:

``` r
install.packages('mzAnnotation',repo = 'https://aberhrml.github.io/drat/')
```

#### macOS

`mzAnnotation` is not currently supported on macOS.

### Learn more

The package documentation can be browsed online at <https://aberHRML.github.io/mzAnnotation/>. 

If this is your first time using `mzAnnotation` see the [introduction vignette](https://aberHRML.github.io/binneR/articles/mzAnnotation.html) for information on how to get started.

If you believe you've found a bug in `mzAnnotation`, please file a bug (and, if
possible, a [reproducible example](https://reprex.tidyverse.org)) at
<https://github.com/aberHRML/mzAnnotation/issues>.
