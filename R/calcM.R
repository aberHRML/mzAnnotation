#' Calculate the mass of an m/z
#' @description calculate M for a given m/z, adduct, isotope and transformation
#' @param mz m/z for which to calculate M
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adduct_rules_table tibble containing available adduct formation rules. Defaults to `adduct_rules()`.
#' @param isotope_rules_table tibble containing available isotopic rules. Defaults to `isotope_rules()`.
#' @param transformation_rules_table tibble containing available transformation rules. Defaults to `transformation_rules()`.
#' @importFrom dplyr filter
#' @examples calcM(118.08626,adduct = '[M+H]1+',isotope = '13C',transformation = 'M - [O] + NH2]')
#' @export

calcM <- function(mz, 
                  adduct = '[M+H]1+', 
                  isotope = NA, 
                  transformation = NA, 
                  adduct_rules_table = adduct_rules(), 
                  isotope_rules_table = isotope_rules(), 
                  transformation_rules_table = transformation_rules()){
  
  addRule <- filter(adduct_rules_table,Name == adduct)
  
  M <- ((mz - addRule$Add) * addRule$Charge) 
  
  if (!is.na(isotope)) {
    isoRule <- filter(isotope_rules_table,Isotope == isotope)
    M <- (M - isoRule$`Mass Difference`)
  }
  
  M <- M / addRule$xM
  
  if (!is.na(transformation)) {
    transformRule <- filter(transformation_rules_table, `MF Change` == transformation)
    M <- M - transformRule$Difference
  }
  return(round(M,5))
}
