# mzAnnotation 2.0.0

* All chemical structure (SMILES, InChI) and chemoinformatics functionality based on OpenBabel has now be moved to the [cheminf](https://github.com/jasenfinch/cheminf) package. 

* The adduct, isotope and transformation rules functions `adduct()`, `isotope()` and `transformations()` have been renamed to `adduct_rules()`, `isotope_rules()` and `transformation_rules()` respectively.

* Added the `adduct_names()`, `isotope_names()` and `transformation_names()` functions to list the available adducts, isotopes and transformations in the default rules tables.

* Added `suitableElementRanges()` to calculate minimum and maximum elemental ranges for a specified molecular mass, suitable for molecular formula generation.

* Added `elementRatios()` for calculating the ratio between all elements present in specified molecular formulas.

* Added `CHOproportion()` for calculating the proportion of CHO present in specified molecular formulas.

* Added functions for calculating the rings plus double bonds equivalent (`rdbe()`), LEWIS valence test `lewis()` and `senior()` molecular graph test for molecular formulas.

* Added functions (`goldenRules()`,`elementCountCheck()`,`elementRatioCheck()`) for calculating the Golden Rules 2, 3, 5 and 6 from [Kind et al. 2007](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-105).

* Added `ipMF()` for generating molecular formulas for electrospray ionisation products.

* Added `transformationPossible()` to test if a transformation between two molecular formulas is possible.

* `count.elements()` is now re-exported from the [CHNOSZ](https://chnosz.net/) package.

* The `digits` option is now set to 7 on package load.

* The specified `adduct`, `isotope` and `transformation` arguments are now checked against the specified `adduct_rules_table`, `isotope_rules_table` and `transformation_rules_table` in `calcM()` and `calcMZ()`. An error is thrown if an are not present.

* Numerous function documentation improvements.

# mzAnnotation 1.7.5

* Fixed error in `relationshipCalculator()` for empty returns. 

# mzAnnotation 1.7.4

* Remove `modes` argument form `relationshipCalcualtor()`.

* Speed up`relationshipCalculator()` by ~x20.

# mzAnnotation 1.7.3

* Added a `NEWS.md` file to track changes to the package.

* Added a `pkgdown` site available at https://jasenfinch.github.io/mzAnnotation/.
