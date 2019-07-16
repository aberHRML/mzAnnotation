#' calcM
#' @description calculate M for a given m/z, adduct, isotope and transformation
#' @param mz m/z for which to calculate M
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adductTable adduct table containing available adduct rules. Defaults to \code{Adducts}.
#' @param isotopeTable isotope table containing available isotope rules. Defaults to \code{Isotopes}.
#' @param transformationTable transformations table containing available transformations rules. Defaults to \code{Transformations}.
#' @importFrom dplyr filter
#' @examples calcM(118.08626,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
#' @export

calcM <- function(mz, adduct = '[M+H]1+', isotope = NA, transformation = NA, adductTable = adducts(), isotopeTable = isotopes(), transformationTable = transformations()){
  
  addRule <- filter(adductTable,Name == adduct)
  
  M <- ((mz - addRule$Add) * addRule$Charge) 
  
  if (!is.na(isotope)) {
    isoRule <- filter(isotopeTable,Isotope == isotope)
    M <- (M - isoRule$`Mass Difference`)
  }
  
  M <- M / addRule$xM
  
  if (!is.na(transformation)) {
    transformRule <- filter(transformationTable, `MF Change` == transformation)
    M <- M - transformRule$Difference
  }
  return(round(M,5))
}
