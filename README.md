[![Build Status](https://travis-ci.org/jasenfinch/mzAnnotation.svg)](https://travis-ci.org/jasenfinch/mzAnnotation)[![Build status](https://ci.appveyor.com/api/projects/status/b9wgaej0u690ls20/branch/master?svg=true)](https://ci.appveyor.com/project/jasenfinch/mzannotation/branch/master)

### mzAnnotation

An R package for explanatory FIE-HRMS m/z annotation

#### Installation

To install run:
```R
source("https://bioconductor.org/biocLite.R")
biocLite("mzR")

library(devtools)

install_github('jasenfinch/OrbiFIEproc')

install_github('jasenfinch/mzAnnotation')
```