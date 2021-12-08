#' Calculate suitable elemental frequency ranges
#' @description Calculate elemental frequency ranges for a given mass which are suitable for molecular formula generation.
#' @param mass accurate mass
#' @return named numeric vector of element frequencies
#' @examples
#' suitableElementRanges(342.11621)
#' @export

suitableElementRanges <- function(mass){
  carb <- round(mass/12)
  Hs <- round(carb * 2)
  NO <- round(carb / 2)
  PS <- round(carb / 4)
  
  maxi <- list(C = c(0,carb),
               H = c(0,Hs),
               N = c(0,NO),
               O = c(0,NO),
               P = c(0,PS),
               S = c(0,PS))
  
  return(maxi)
}

#' Calculate the rings plus double bonds equivalent
#' @description Calculate the rings plus double bonds equivalent (RDBE) for molecular formulas.
#' @param element_frequencies table of element frequencies for a set of molecular formulas as returned by `elementFrequencies()`.
#' @param valences named list of element valences
#' @return A vector of RDBE values.
#' @examples 
#' element_frequencies <- elementFrequencies(c('C12H22O11','C12H22NO11'))
#' rdbe(element_frequencies)
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_split
#' @importFrom tidyr replace_na
#' @export

rdbe <- function(element_frequencies,
                 valences = list(C = 4,
                                 H = 1,
                                 N = 3,
                                 O = 2,
                                 P = 3,
                                 S = 4)){
  valences <- valences %>% 
    as_tibble() %>% 
    gather(element,valence) %>% 
    mutate(valence = valence - 2)
  
  element_frequencies <- element_frequencies %>% 
    gather(element,
           frequency,
           -MF) %>% 
    mutate(frequency = replace_na(frequency,0)) %>% 
    left_join(valences,by = 'element')
  
  element_frequencies %>% 
    group_by(MF) %>% 
    mutate(value = frequency * valence) %>%
    ungroup() %>% 
    group_split(MF) %>% 
    map_dbl(~{
      .x$value %>% 
        sum() %>% 
        {. + 2} %>% 
        {./2}
    })
}

#' LEWIS and SENIOR checks
#' @rdname lewis_senior
#' @description LEWIS molecular formula valence test and SENIOR test for the existence of molecular graphs. 
#' @param element_frequencies table of element frequencies for a set of molecular formulas as returned by `elementFrequencies()`.
#' @param valences named list of element valences
#' @return Boolean vector of check results for each molecular formula.
#' @examples 
#' element_frequencies <- elementFrequencies(c('C12H22O11','C12H22NO11'))
#' lewis(element_frequencies)
#' senior(element_frequencies)
#' @export

lewis <- function(element_frequencies,
                  valences = list(C = 4,
                                  H = 1,
                                  N = 3,
                                  O = 2,
                                  P = 3,
                                  S = 4)){
  tibble(RDBE = rdbe(element_frequencies,
                     valences = valences)) %>% 
    rowwise() %>% 
    mutate(remainder = RDBE %% 1,
           LEWIS =  ifelse(
             all(
               RDBE >= 0,
               remainder != 0.5
             ),
             TRUE,
             FALSE)) %>% 
    .$LEWIS
}

#' @rdname lewis_senior
#' @export

senior <- function(element_frequencies,
                   valences = list(C = 4,
                                   H = 1,
                                   N = 3,
                                   O = 2,
                                   P = 3,
                                   S = 4)){
  valences <- valences %>% 
    as_tibble() %>% 
    gather(element,valence)
  
  element_frequencies <- element_frequencies %>%  
    gather(element,frequency,-MF) %>% 
    mutate(frequency = replace_na(frequency,0)) %>% 
    left_join(valences,by = 'element') %>% 
    mutate(total_valence = frequency * valence) %>% 
    group_by(MF)
  
  mfs <- MF
  
  element_frequencies %>% 
    summarise(sum_valence = sum(total_valence)) %>% 
    select(MF,sum_valence) %>% 
    left_join(element_frequencies %>%
                filter((valence %% 2) != 0) %>% 
                summarise(odd_valence_total = sum(frequency)),
              by = 'MF') %>% 
    left_join(element_frequencies %>% 
                filter(frequency > 0) %>% 
                summarise(twice_maximum_valence = max(valence) * 2),
              by = 'MF') %>% 
    left_join(element_frequencies %>% 
                summarise(twice_atoms_minus_1 = sum(frequency) %>% 
                            {. * 2 - 1}),
              by = 'MF') %>% 
    rowwise() %>% 
    mutate(SENIOR = ifelse(
      all(
        (sum_valence %% 2) == 0 | (odd_valence_total %% 2) == 0,
        sum_valence >= twice_maximum_valence,
        sum_valence >= twice_atoms_minus_1
      ),
      TRUE,
      FALSE
    )) %>% 
    mutate(MF = factor(MF,levels = mfs)) %>% 
    arrange(MF) %>% 
    .$SENIOR
}


#' Molecular formula generation
#' @description Molecular formula generation
#' @param mass accurate mass
#' @param ppm ppm tolerance
#' @param charge charge
#' @param element_ranges named list of element
#' @param LEWIS return only molecular formulas that pass the LEWIS check
#' @param SENIOR return only molecular formulas that pass the SENiOR check
#' @author Jasen Finch
#' @importFrom rcdk generate.formula
#' @export
#' @return A \code{tibble} containing the generated MFs, their theoretical mass and their PPM error.
#' @examples
#' generateMF(342.11621,
#'            element_ranges = list(C = c(0,12),
#'                                  H = c(0,22),
#'                                  O = c(0,11)))

generateMF <- function(mass, 
                       ppm = 1, 
                       charge = 0, 
                       element_ranges = suitableElementRanges(mass),
                       LEWIS = TRUE,
                       SENIOR = TRUE){
  
  element_ranges <- element_ranges %>%
    names() %>% 
    map(~{
      c(.x,element_ranges[[.x]])
    })
  
  window <- ppmRange(mass,ppm) %>% 
    {.$upper - .$lower} %>% 
    {./2}
  
  molecular_formulas <- generate.formula(mass,
                                         window = window,
                                         elements = element_ranges,
                                         validation = FALSE,
                                         charge = charge)	 %>% 
    map(~{
      
      tibble(MF = .x@string,
             Mass = round(.x@mass,
                          5)) %>% 
          mutate(`PPM error` = ppmError(mass,Mass))
    }) %>% 
    bind_rows()
  
  mf_checks <- molecular_formulas$MF %>% 
    elementFrequencies() %>% 
    tibble(
      RDBE = rdbe(.),
      LEWIS = lewis(.),
      SENIOR = senior(.)
    ) %>% 
    select(MF,RDBE,LEWIS,SENIOR)
  
  molecular_formulas <- molecular_formulas %>% 
    left_join(mf_checks,
              by = 'MF')
  
  if (isTRUE(LEWIS)){
    molecular_formulas <- molecular_formulas %>% 
      filter(LEWIS == TRUE)
  }
  
  if (isTRUE(SENIOR)){
    molecular_formulas <- molecular_formulas %>% 
      filter(SENIOR == TRUE)
  }
  
  return(molecular_formulas)
}

#' Ionisation product molecular formula generation
#' @description Generate molecular formulas for a given ionisation product accurate *m/z*.
#' @param mz Accurate *m/z*
#' @param adduct ionisation product adduct
#' @param isotope ionisation product isotope
#' @param ppm ppm error tolerance threshold
#' @param LEWIS return only molecular formulas that pass the LEWIS check
#' @param SENIOR return only molecular formulas that pass the SENiOR check
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
                 LEWIS = TRUE,
                 SENIOR = TRUE,
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
             Isotope = isotope,
             Score = MFscore(MF)) %>% 
      select(MF,
             Adduct,
             Isotope,
             `Measured m/z`,
             `Measured M`,
             `Theoretical m/z`,
             `Theoretical M`,
             `PPM error`,
             Score) %>% 
      ungroup() %>% 
      arrange(Score)
  } else {
    mfs <- tibble(MF = character(),
                  Adduct = character(),
                  Isotope = character(),
                  `Measured m/z` = numeric(),
                  `Measured M` = numeric(),
                  `Theoretical M` = numeric(),
                  `PPM error` = numeric(),
                  Score = numeric())
  }
  
  return(mfs)
}

#' Molecular formula scoring
#' @description Molecular formula scoring based on elemental ratios
#' @param mf Molecular formula
#' @examples
#' MFscore('C12H22O11')
#' @importFrom CHNOSZ count.elements
#' @export

MFscore <- function(mf){
  elements <- c('C','H','N','O','P','S')
  
  eleFreq <- as.vector(count.elements(mf))
  ele <- names(count.elements(mf))
  
  if (length(which(!(elements %in% ele))) > 0) {
    eleFreq <- c(eleFreq,rep(0,length(which(!(elements %in% ele)))))
    ele <- c(ele,elements[!(elements %in% ele)])
  }
  
  colnames(eleFreq) <- NULL
  
  eleRatios <- c(
    `H/C` = if ('H' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'H']/eleFreq[ele == 'C']
    },
    `N/C` = if ('N' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'N']/eleFreq[ele == 'C']
    },
    `O/C` = if ('O' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'O']/eleFreq[ele == 'C']
    },
    `P/C` = if ('P' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'P']/eleFreq[ele == 'C']
    },
    `S/C` = if ('S' %in% ele & 'C' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'C']
    },
    `N/O` = if ('N' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'N']/eleFreq[ele == 'O']
    },
    `P/O` = if ('P' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'P']/eleFreq[ele == 'O']
    },
    `S/O` = if ('S' %in% ele & 'O' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'O']
    },
    `O/P` = if ('O' %in% ele & 'P' %in% ele) {
      eleFreq[ele == 'O']/eleFreq[ele == 'P']
    },
    `S/P` = if ('S' %in% ele & 'P' %in% ele) {
      eleFreq[ele == 'S']/eleFreq[ele == 'P']
    }
  )
  if (!is.null(eleRatios)) {
    if ('H/C' %in% names(eleRatios)) {
      eleRatios['H/C'] <- abs(eleRatios['H/C'] - 1.6)
    }
    if ('O/C' %in% names(eleRatios)) {
      eleRatios['O/C'] <- abs(eleRatios['O/C'] - 0.3)
    }
    if  (is.nan(eleRatios['N/O'])) {
      eleRatios['N/O'] <- 0
    }
    if (is.infinite(eleRatios['N/O'])) {
      eleRatios['N/O'] <- eleFreq[ele == 'N']
    }
    if (is.nan(eleRatios['P/O'])) {
      eleRatios['P/O'] <- 0
    }
    if (is.infinite(eleRatios['P/O'])) {
      eleRatios['P/O'] <- eleFreq[ele == 'P']
    }
    if (is.nan(eleRatios['S/O'])) {
      eleRatios['S/O'] <- 0
    }
    if (is.infinite(eleRatios['S/O'])) {
      eleRatios['S/O'] <- eleFreq[ele == 'S']
    }
    if (is.nan(eleRatios['O/P'])) {
      eleRatios['O/P'] <- 0
    }
    if ('O/P' %in% names(eleRatios) & eleRatios['O/P'] >= 3) {
      eleRatios['O/P'] <- 0
    }
    if (is.nan(eleRatios['S/P'])) {
      eleRatios['S/P'] <- 0
    }
    if (is.infinite(eleRatios['S/P'])) {
      eleRatios['S/P'] <- eleFreq[ele == 'S']
    }
    score <- sum(eleRatios) 
  } else {
    score <- NA
  }
  
  return(score)
}
