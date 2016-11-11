[![Build Status](https://travis-ci.org/jasenfinch/mzAnnotation.svg)](https://travis-ci.org/jasenfinch/mzAnnotation) [![Build status](https://ci.appveyor.com/api/projects/status/b9wgaej0u690ls20/branch/master?svg=true)](https://ci.appveyor.com/project/jasenfinch/mzannotation/branch/master)

### mzAnnotation

An R package for putative accurate m/z annotation

#### Installation

To install run:
```R
library(devtools)

install_github('jasenfinch/mzAnnotation')
```

#### Tools

Available tools include:

* Molecular formular generator
```r
  res <- mfGen(341.10894,
              max=c(C = 12,iC = 0,H = 22,iH = 0,N = 0,iN = 0,O = 11,iO = 0,F = 0 ,Na = 0,
                   Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0),
              min=c(C = 0,iC = 0,H = 0,iH = 0,N = 0,iN = 0,O = 0,iO = 0,F = 0 ,Na = 0,
                   Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0),
             tolerance=0.01,charge=-1,applygr=TRUE)
```

* Isotope distribution calculator
```r
res <- isoDistr('C4H5O5',chrg=-1)
```
* Putative ionisation product searches
```r
res <- getPIP(133.01378,'n',5)
```