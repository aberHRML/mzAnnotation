#' calcM
#' @description calculate M for a given m/z, adduct, isotope and transformation
#' @param mz m/z for which to calculate M
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adducts adduct table containing available adduct rules. Defaults to table returned by \code{availableAdducts()}.
#' @param isotopes isotope table containing available isotope rules. Defaults to table returned by \code{availableIsotopes()}.
#' @param transformations transformations table containing available transformations rules. Defaults to table returned by \code{availableTransformations()}.
#' @importFrom dplyr filter
#' @examples calcM(118.08626)
#' @export

calcM <- function(mz, adduct = '[M+H]1+', isotope = NA, transformation = NA, adducts = availableAdducts(), isotopes = availableIsotopes(), transformations = availableTransformations()){
  
  addRule <- filter(adducts,Name == adduct)
  
  M <- ((mz - addRule$Add) * addRule$Charge) / addRule$xM 
  
  if (!is.na(isotope)) {
    isoRule <- filter(isotopes,Isotope == isotope)
    M <- M - isoRule$Mass.Difference
  }
  
  if (!is.na(transformation)) {
    transformRule <- filter(transformations, MF.Change == transformation)
    M <- M - transformRule$Difference
  }
  return(round(M,5))
}
