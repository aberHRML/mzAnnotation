#' calcM
#' @description calculate an m/z for a given M, adduct, isotope and transformation
#' @param M M for which to calculate an m/z
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adductTable adduct table containing available adduct rules. Defaults to \code{adducts()}.
#' @param isotopeTable isotope table containing available isotope rules. Defaults to \code{isotopes()}.
#' @param transformationTable transformations table containing available transformations rules. Defaults to \code{transformations()}.
#' @importFrom dplyr filter
#' @examples calcMZ(116.05182,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
#' @export

calcMZ <- function(M, adduct = '[M+H]1+', isotope = NA, transformation = NA, adductTable = adducts(), isotopeTable = isotopes(), transformationTable = transformations()){
  
  addRule <- filter(adductTable,Name == adduct)
  
  if (!is.na(transformation)) {
    transformRule <- filter(transformationTable, `MF Change` == transformation)
    M <- M + transformRule$Difference
  }
  
  if (!is.na(isotope)) {
    isoRule <- filter(isotopeTable,Isotope == isotope)
    mz <- (M * addRule$xM) + isoRule$`Mass Difference`
  } else {
    mz <- (M * addRule$xM)
  }
  
   mz <- mz / addRule$Charge + addRule$Add
  
  return(round(mz,5))
}