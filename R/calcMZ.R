#' Calculate the m/z of a mass
#' @description calculate an m/z for a given M, adduct, isotope and transformation
#' @param M M for which to calculate an m/z
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adduct_rules_table adduct table containing available adduct rules. Defaults to \code{adducts()}.
#' @param isotope_rules_table isotope table containing available isotope rules. Defaults to \code{isotopes()}.
#' @param transformation_rules_table transformations table containing available transformations rules. Defaults to \code{transformations()}.
#' @importFrom dplyr filter
#' @examples calcMZ(116.05182,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
#' @export

calcMZ <- function(M, 
                   adduct = '[M+H]1+', 
                   isotope = NA, 
                   transformation = NA, 
                   adduct_rules_table = adduct_rules(), 
                   isotope_rules_table = isotope_rules(), 
                   transformation_rules_table = transformation_rules()){
  
  addRule <- filter(adduct_rules_table,Name == adduct)
  
  if (!is.na(transformation)) {
    transformRule <- filter(transformation_rules_table, 
                            `MF Change` == transformation)
    M <- M + transformRule$Difference
  }
  
  if (!is.na(isotope)) {
    isoRule <- filter(isotope_rules_table,
                      Isotope == isotope)
    mz <- (M * addRule$xM) + isoRule$`Mass Difference`
  } else {
    mz <- (M * addRule$xM)
  }
  
   mz <- mz / addRule$Charge + addRule$Add
  
  return(round(mz,5))
}