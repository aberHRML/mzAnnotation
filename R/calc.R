#' Calculate the molecular mass of an m/z
#' @description Calculate the molecular mass for a given mass to charge ratio, adduct, isotope and transformation.
#' @param mz mass to charge ratio
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adduct_rules_table tibble containing available adduct formation rules. Defaults to `adduct_rules()`.
#' @param isotope_rules_table tibble containing available isotopic rules. Defaults to `isotope_rules()`.
#' @param transformation_rules_table tibble containing available transformation rules. Defaults to `transformation_rules()`.
#' @importFrom dplyr filter
#' @examples 
#' calcM(
#'   mz = 118.08626,
#'   adduct = '[M+H]1+',
#'   isotope = '13C',
#'   transformation = 'M - [O] + [NH2]')
#' @export

calcM <- function(mz, 
                  adduct = '[M+H]1+', 
                  isotope = NA, 
                  transformation = NA, 
                  adduct_rules_table = adduct_rules(), 
                  isotope_rules_table = isotope_rules(), 
                  transformation_rules_table = transformation_rules()){
  
  checkAdductTable(adduct_rules_table)
  checkIsotopeTable(isotope_rules_table)
  checkTransformationTable(transformation_rules_table)
  
  if (!adduct %in% adduct_rules_table$Name){
    stop('The specified `adduct` is not found in the adduct rules table',
         call. = FALSE)
  }
  
  if (!is.na(isotope) & !isotope %in% isotope_rules_table$Isotope){
    stop('The specified `isotope` is not found in the isotope rules table',
         call. = FALSE)
  }
  
  if (!is.na(transformation) & !transformation %in% transformation_rules_table$`MF Change`){
    stop('The specified `transformation` is not found in the transformation rules table',
         call. = FALSE)
  }
  
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

#' Calculate the m/z of a mass
#' @description Calculate the mass to charge ratio for a given molecular mass, adduct, isotope and transformation.
#' @param M molecular mass
#' @param adduct adduct to apply
#' @param isotope isotope to apply
#' @param transformation transformation to apply
#' @param adduct_rules_table adduct table containing available adduct rules. Defaults to \code{adducts()}.
#' @param isotope_rules_table isotope table containing available isotope rules. Defaults to \code{isotopes()}.
#' @param transformation_rules_table transformations table containing available transformations rules. Defaults to \code{transformations()}.
#' @importFrom dplyr filter
#' @examples
#' calcMZ(
#'   M = 116.05182,
#'   adduct = '[M+H]1+',
#'   isotope = '13C',
#'   transformation = 'M - [O] + [NH2]')
#' @export

calcMZ <- function(M, 
                   adduct = '[M+H]1+', 
                   isotope = NA, 
                   transformation = NA, 
                   adduct_rules_table = adduct_rules(), 
                   isotope_rules_table = isotope_rules(), 
                   transformation_rules_table = transformation_rules()){
  
  checkAdductTable(adduct_rules_table)
  checkIsotopeTable(isotope_rules_table)
  checkTransformationTable(transformation_rules_table)
  
  if (!adduct %in% adduct_rules_table$Name){
    stop('The specified `adduct` is not found in the adduct rules table',
         call. = FALSE)
  }
  
  if (!is.na(isotope) & !isotope %in% isotope_rules_table$Isotope){
    stop('The specified `isotope` is not found in the isotope rules table',
         call. = FALSE)
  }
  
  if (!is.na(transformation) & !transformation %in% transformation_rules_table$`MF Change`){
    stop('The specified `transformation` is not found in the transformation rules table',
         call. = FALSE)
  }
  
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

#' Calculate a molecular formula  accurate mass
#' @description Calculate the accurate mass of a specified molecular formula.
#' @param MF the molecular formula for which to calculate the accurate mass
#' @param charge the total charge of the specified molecular formula
#' @param element_table elemental mass information. Defaults to \code{\link{elements}()}.
#' @examples 
#' calcAccurateMass('C4H5O5',charge = 0)
#' @importFrom CHNOSZ count.elements
#' @importFrom dplyr mutate select
#' @export

calcAccurateMass <- function(MF,charge = 0, element_table = elements()) {
  e <- element_table$AtomicMass[element_table$Element == 'e']
  elementFreq <- tibble(Element = names(count.elements(MF)),Frequency = count.elements(MF))
  elementMasses <- element_table %>%
    filter(RelativeAbundance == 1) %>%
    select(Element,AtomicMass)
  
  suppressMessages(elementFreq <- left_join(elementFreq,elementMasses) %>%
                     mutate(M = Frequency * AtomicMass))
  
  M <- sum(elementFreq$M)
  
  if (charge > 0) {
    M <- M - e * abs(charge)
  } 
  if (charge < 0) {
    M <- M + e * abs(charge)
  }
  
  if (charge == 0) {
    charge <- 1
  }
  
  M <- M/abs(charge)
  
  M <- round(M,5)
  return(M)
}