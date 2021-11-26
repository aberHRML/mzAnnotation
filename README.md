
# mzAnnotation

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R build
status](https://github.com/jasenfinch/mzAnnotation/workflows/R-CMD-check/badge.svg)](https://github.com/jasenfinch/mzAnnotation/actions)
[![Coverage
Status](https://img.shields.io/codecov/c/github/jasenfinch/mzAnnotation/devel.svg)](https://codecov.io/github/jasenfinch/mzAnnotation?branch=devel)
[![DOI](https://zenodo.org/badge/33083554.svg)](https://zenodo.org/badge/latestdoi/33083554)
<!-- badges: end -->

An R package containing tools for putative annotation of accurate m/z
from electrospray ionisation mass spectrometry data.

### Installation

To install run:

``` r
devtools::install_github('jasenfinch/mzAnnotation')
```

### Tools

Available tools include:

-   Adduct, isotope and biotransfromation relationship prediction

``` r
res <- relationshipCalculator(c(132.03023,168.00691))

res
#> # A tibble: 1 x 9
#>   `m/z1` `m/z2` Adduct1 Adduct2  Isotope1 Isotope2 Transformation1
#>    <dbl>  <dbl> <chr>   <chr>    <lgl>    <lgl>    <lgl>          
#> 1   132.   168. [M-H]1- [M+Cl]1- NA       NA       NA             
#> # … with 2 more variables: Transformation2 <lgl>, Error <dbl>
```

-   Molecular formula generation

``` r
res <- generateMF(342.11621,
                  element_max = c(C = 12,H = 22,N = 0,
                                  O = 11,P = 0,S = 0))
res
#> # A tibble: 1 x 3
#>   MF         Mass `PPM Error`
#>   <chr>     <dbl>       <dbl>
#> 1 C12H22O11  342.           0
```

-   Isotope distribution calculation

``` r
res <- isotopeDistribution(MF = 'C4H5O5',charge = -1)
res
#> # A tibble: 8 x 4
#>   Isotope      `m/z` `Relative Abundance` Probability
#>   <chr>        <dbl>                <dbl>       <dbl>
#> 1 <NA>          133.            1           0.945    
#> 2 13C 1         134.            0.0449      0.0424   
#> 3 18O 1         135.            0.0100      0.00947  
#> 4 17O 1         134.            0.00200     0.00189  
#> 5 13C 2         135.            0.000756    0.000714 
#> 6 2H 1          134.            0.000500    0.000472 
#> 7 13C 1; 18O 1  136.            0.000450    0.000425 
#> 8 13C 1; 17O 1  135.            0.0000900   0.0000850
```

-   Putative ionisation product searches

``` r
db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
res <- PIPsearch(db = db,
                 mz = 132.03023,
                 ppm = 5,
                 adduct = '[M-H]1-',
                 isotope = NA)
res
#> # A tibble: 1 x 12
#>      ID NAME   InChI       InChIKey   SMILES  MF    Accurate_Mass Isotope Adduct
#>   <int> <chr>  <chr>       <chr>      <chr>   <chr>         <dbl> <lgl>   <chr> 
#> 1     4 L-Asp… InChI=1S/C… CKLJMWTZI… C([C@@… C4H7…          133. NA      [M-H]…
#> # … with 3 more variables: Measured m/z <dbl>, Theoretical m/z <dbl>,
#> #   PPM Error <dbl>
```
