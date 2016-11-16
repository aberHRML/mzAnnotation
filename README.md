[![Build Status](https://travis-ci.org/jasenfinch/mzAnnotation.svg)](https://travis-ci.org/jasenfinch/mzAnnotation) [![Build status](https://ci.appveyor.com/api/projects/status/b9wgaej0u690ls20/branch/master?svg=true)](https://ci.appveyor.com/project/jasenfinch/mzannotation/branch/master)

### mzAnnotation

An R package for putative annotation of accurate m/z

#### Installation

To install run:
```R
devtools::install_github('jasenfinch/mzAnnotation')
```

#### Tools

Available tools include:

* Adduct, isotope and biotransfromation relationship prediction
```r
res <- relationshipPredictor(mz = c(132.03023,133.01425,133.03359,168.00691),mode = 'n')
```

* Molecular formula generator
```r
res <- generateMF(mz = 341.10894,ppm = 5,charge = -1, 
                  applygr = TRUE, 
                  composition = c(C = 12,iC = 0,H = 22,iH = 0,
                  N = 0,iN = 0,O = 11,iO = 0,F = 0 ,Na = 0,
                  Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,
                  Br = 0,iBr = 0,K = 0,iK = 0))
```

* Isotope distribution calculator
```r
res <- isoDistr(mf = 'C4H5O5',chrg=-1)
```
* Putative ionisation product searches
```r
res <- PIPsearch(mz = 133.01378,mode = 'n',ppm = 5)
```