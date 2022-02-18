#' Molecular formula generation
#' @description Molecular formula generation
#' @param mass accurate mass
#' @param ppm ppm tolerance
#' @param charge charge
#' @param element_ranges named list of element
#' @return A \code{tibble} containing the generated MFs, their theoretical mass and their PPM error.
#' @examples
#' generateMF(342.11621,
#'            element_ranges = list(C = c(0,12),
#'                                  H = c(0,22),
#'                                  O = c(0,11)))
#' @importFrom magrittr set_names
#' @export

generateMF <- function(mass, 
                       ppm = 1, 
                       charge = 0, 
                       element_ranges = suitableElementRanges(mass)){
  
  elements_symbols <- c('C','H','N','O','P','S')
  default_element_ranges <- map(seq_along(elements_symbols),
                                ~ rep(0,2)) %>% 
    set_names(elements_symbols)
  default_element_ranges[names(element_ranges)] <- element_ranges
  
  molecular_formulas <- generate(mass,ppm,charge,default_element_ranges) %>% 
    as_tibble() %>% 
    mutate(Mass = round(Mass,5),
           `PPM error` = ppmError(mass,Mass)) %>% 
    arrange(`PPM error`)
  
  return(molecular_formulas)
}

#' Ionisation product molecular formula generation
#' @description Generate molecular formulas for a given ionisation product accurate *m/z*.
#' @param mz Accurate *m/z*
#' @param adduct ionisation product adduct
#' @param isotope ionisation product isotope
#' @param ppm ppm error tolerance threshold
#' @param adduct_rules_table tibble containing available adduct formation rules. Defaults to `adduct_rules()`.
#' @param isotope_rules_table tibble containing available isotopic rules. Defaults to `isotope_rules()`.
#' @examples 
#' ipMF(118.08626,adduct = '[M+H]1+')
#' @importFrom dplyr arrange
#' @export

ipMF <- function(mz,
                 adduct = "[M+H]1+",
                 isotope = NA,
                 ppm = 5, 
                 adduct_rules_table = adduct_rules(),
                 isotope_rules_table = isotope_rules()){
  
  M <- calcM(mz,
             adduct,
             isotope,
             adduct_rules_table = adduct_rules_table,
             isotope_rules_table = isotope_rules_table)
  
  if (M < 100 & ppm < 10) {
    ppm <- 10
  } else {
    ppm <- ppm
  }
  
  ppm <- (ppm/10^6 * mz)/M * 10^6
  
  mfs <- generateMF(M,
                    ppm = ppm) 
  
  if (nrow(mfs) > 0){
    mfs <- mfs %>%
      rename(`Theoretical M` = Mass) %>% 
      rowwise() %>% 
      mutate(`Measured m/z` = mz,
             `Measured M` = M,
             `Theoretical m/z` = calcMZ(`Theoretical M`,
                                        adduct = adduct,
                                        isotope = isotope,
                                        adduct_rules_table = adduct_rules_table,
                                        isotope_rules_table = isotope_rules_table),
             Adduct = adduct,
             Isotope = isotope) %>% 
      select(MF,
             Adduct,
             Isotope,
             `Measured m/z`,
             `Measured M`,
             `Theoretical m/z`,
             `Theoretical M`,
             `PPM error`) %>% 
      ungroup()
    
    gr_scores <- mfs$MF %>% 
      goldenRules() %>% 
      goldenRulesScore()
    
    mfs <- mfs %>% 
      left_join(select(gr_scores,MF,
                       `Plausibility (%)`), 
                by = "MF") %>% 
      mutate(`PPM error` = abs(`PPM error`)) %>% 
      arrange(desc(`Plausibility (%)`),`PPM error`)
  } else {
    mfs <- tibble(MF = character(),
                  Adduct = character(),
                  Isotope = character(),
                  `Measured m/z` = numeric(),
                  `Measured M` = numeric(),
                  `Theoretical m/z` = numeric(),
                  `Theoretical M` = numeric(),
                  `PPM error` = numeric(),
                  `Plausibility (%)` = numeric())
  }
  
  return(mfs)
}
