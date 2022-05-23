#' Molecular formula generation
#' @description Molecular formula generation
#' @param mass accurate mass
#' @param ppm ppm tolerance
#' @param charge charge
#' @param element_ranges named list of element ranges
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

#' Molecular formula isotopic possiblity
#' @description Check if an isotope is possible for a vector of molecular formulas.
#' @param MF a character vector of molecular formulas
#' @param isotope the isotope to check
#' @param isotope_rules_table tibble containing available isotopic rules. Defaults to `isotope_rules()`.
#' @return A boolean vector specifying if the specified isotope is possible for the molecular formulas.
#' @examples isotopePossible(c('C12H22O11','H2O'))
#' @importFrom stringr str_remove_all
#' @export

isotopePossible <- function(MF,
                             isotope = '13C',
                             isotope_rules_table = isotope_rules()){
  
  checkIsotopeTable(isotope_rules_table)
  
  if (!is.na(isotope)){
    if (!(isotope %in% isotope_rules_table$Isotope)){
      stop('Specified isotope not found in isotopic rules table.',
           call. = FALSE)
    }
    
    element <- isotope %>% 
      str_remove_all('[0-9]')
    
    isotope_rule <- isotope_rules_table %>% 
      filter(Isotope == isotope) %>% 
      .$Rule %>% 
      parse_expr()
    
    element_frequencies <- MF %>% 
      elementFrequencies() 
    
    if (!(element %in% colnames(element_frequencies))){
      element_frequencies <- element_frequencies %>% 
        mutate(!!element := NA)
    }
    
    
    element_frequencies %>% 
      mutate(possible = !!isotope_rule %>% 
             replace_na(FALSE)
             ) %>% 
      select(MF,possible) %>% 
      .$possible 
  }
  else NA
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
  
  empty <- tibble(MF = character(),
                  Adduct = character(),
                  Isotope = character(),
                  `Measured m/z` = numeric(),
                  `Measured M` = numeric(),
                  `Theoretical m/z` = numeric(),
                  `Theoretical M` = numeric(),
                  `PPM error` = numeric(),
                  `Plausibility (%)` = numeric())
  
  M <- calcM(mz,
             adduct,
             isotope,
             adduct_rules_table = adduct_rules_table,
             isotope_rules_table = isotope_rules_table)
  
  ppm_M <- (ppm/10^6 * mz)/M * 10^6
  
  mfs <- generateMF(M,
                    ppm = ppm_M,
                    element_ranges = suitableElementRanges(M)) 
  
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
             Isotope = isotope,
             isotope_possible = isotopePossible(MF,
                                                Isotope,
                                                isotope_rules_table),
             `PPM error` = ppmError(`Measured m/z`,
                                    `Theoretical m/z`)
    ) %>%
      filter(isTRUE(isotope_possible) | 
               is.na(isotope_possible),
             `PPM error` < ppm) 
  }
  else return(empty)
  
  if (nrow(mfs) > 0){
    mfs <- mfs %>% 
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
  }
  else return(empty)
  
  return(mfs)
}

#' Molecular formula adduct transformation
#' @description adduct transform a molecular formula.
#' @param MF molecular formula to transform.
#' @param adduct adduct to use for transformation.
#' @param adduct_rules_table Adduct formation rules table to use. Defaults to `adduct_rules()`.
#' @examples 
#' adductTransformMF('C6H12O6','[M+H]1+')
#' @importFrom CHNOSZ count.elements
#' @export

adductTransformMF <- function(MF,adduct,adduct_rules_table = adduct_rules()){
  
  checkAdductTable(adduct_rules_table)
  
  tMF <- function(freq,trans,expres){
    if (!is.na(trans)) {
      expAt <- trans %>% count.elements()
      if (F %in% (names(expAt) %in% names(freq))) {
        tmp <- rep(0,length(which(!(names(expAt) %in% names(freq)))))
        names(tmp) <- names(expAt)[!(names(expAt) %in% names(freq))]
        freq <- c(freq,tmp)
        freq <- freq[order(names(freq))]
      }
      if (expres == '+') {
        freq[names(expAt)] <- freq[names(expAt)] + expAt
      }
      if (expres == '-') {
        freq[names(expAt)] <- freq[names(expAt)] - expAt
      }
    }
    return(freq)
  }
  
  adductRule <- adduct_rules_table %>% 
    filter(Name == adduct)
  
  freq <- MF %>% count.elements()
  
  freq <- freq %>% 
    {. * adductRule$xM} %>%
    tMF(adductRule$RemAt,'-') %>%
    tMF(adductRule$AddEx,'+') %>%
    tMF(adductRule$RemEx,'-')
  
  if (T %in% (freq < 0)) {
    freq <- ''
  } else {
    ch <- freq[names(freq) %in% c('C','H')]
    freq <- freq[!(names(freq) %in% c('C','H'))]
    freq <- c(ch,freq)
    freq <- freq[!(freq == 0)]
    freq[freq == 1] <- ''
    freq <- str_c(names(freq),freq) %>%
      str_c(collapse = '') 
  }
  
  return(freq)
  
}

#' Transform a molecular formula
#' @description transform a molecular formula
#' @param MF molecular formula to transform
#' @param transformation transformation to apply
#' @param transformation_rules_table transformations table containing available transformations rules. Defaults to `transformation_rules()`.
#' @details \code{NA} will be returned if \code{MF} cannot be transformed.
#' @examples 
#' transformMF('C4H5O5')
#' transformMF('C4H5N',transformation = 'M - [OH] + [NH2]')
#' @importFrom tidyr gather
#' @importFrom stringr str_c str_replace
#' @export

transformMF <- function(MF, 
                        transformation = 'M - [O] + [NH2]', 
                        transformation_rules_table = transformation_rules()){
  
  checkTransformationTable(transformation_rules_table)
  
  if (!is.na(transformation)) {
    elements <- c('C','H','O','N','P','S')
    
    transformation <- filter(transformation_rules_table,
                             `MF Change` == transformation) %>%
      select(C:S)
    
    MF <- count.elements(MF)
    
    if (length(which(!(elements %in% names(MF)))) > 0) {
      MF <- c(MF,rep(0,length(which(!(elements %in% names(MF))))))
      names(MF)[names(MF) == ''] <- elements[!(elements %in% names(MF))]
    }
    
    MF <- MF[order(names(MF))]
    MF <- MF + transformation
    MF <- gather(MF,'Element','Frequency') 
    
    if (T %in% (MF$Frequency < 0)) {
      MF <- NA
    } else {
      MF <- MF %>%
        filter(Frequency > 0)
      
      MF$Frequency[MF$Frequency == 1] <- ''
      
      MF <- str_c(MF$Element,MF$Frequency) %>%
        str_c(collapse = '') 
    }
  }
  return(MF)
}

#' Transformation check
#' @description Check if a transformation between two molecular formulas is possible
#' @param from molecular formula to be transformed
#' @param to molecular formula product of the transformation
#' @param transformation transformation to apply. As found in column `MF Change` of the table supplied to the arguement `transformation_rules_table`.
#' @param transformation_rules_table the transformation rules table. Defaults to the returned value of `transformation_rules()`. Alternative tables should be supplied in the same format.
#' @return TRUE/FALSE depending whether the transformation is possible
#' @export

transformationPossible <- function(from,
                                   to,
                                   transformation,
                                   transformation_rules_table = transformation_rules()){
  
  checkTransformation(
    transformation,
    transformation_rules_table
  )
  
  product_MF <- transformMF(from,
                            transformation,
                            transformation_rules_table)
  
  possible <- product_MF == to
  
  possible <- ifelse(is.na(possible),FALSE,possible)
  
  return(possible)
}