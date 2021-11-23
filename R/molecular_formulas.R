#' Calculate suitable maximum elemental frequencies
#' @description Calculate maximum elemental frequencies for a given mass which are suitable for molecular formula generation.
#' @param mass accurate mass
#' @return named numeric vector of element frequencies
#' @examples
#' suitableMaxElements(342.11621)
#' @export

suitableMaxElements <- function(mass){
  carb <- round(mass/12)
  Hs <- round(carb * 2)
  NO <- round(carb / 2)
  PS <- round(carb / 4)
  
  maxi <- c(C = carb,
            H = Hs,
            N = NO,
            O = NO,
            P = PS,
            S = PS)
  
  return(maxi)
}


#' Molecular formula generation
#' @description Molecular formula generation
#' @param mass accurate mass
#' @param ppm ppm tolerance
#' @param charge charge
#' @param validation \code{boolean}, apply validation rules
#' @param element_max numeric \code{vector} of maximum elemental composition
#' @param element_min numeric \code{vector} of minimum elemental composition
#' @details this uses the HR2 molecular formula generator available at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/supp/rhrmet.html}.
#' @author Jasen Finch
#' @importFrom CHNOSZ count.elements
#' @importFrom tibble as_tibble
#' @importFrom dplyr arrange
#' @export
#' @return A \code{tibble} containing the generated MFs, their theoretical mass and their PPM error.
#' @examples
#' res <- generateMF(342.11621,
#'                   element_max = c(C = 12,H = 22,N = 0,
#'                                 O = 11,P = 0,S = 0))

generateMF <- function(mass, 
                       ppm = 1, 
                       charge = 0, 
                       validation = TRUE, 
                       element_max = suitableMaxElements(mass), 
                       element_min = c(C = 0,H = 0,N = 0,O = 0,P = 0,S = 0)){

    comp_max <- c(C = 0,iC = 0,H = 0,iH = 0,N = 0,iN = 0,O = 0,iO = 0,F = 0 ,Na = 0,
             Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0)
    comp_min <- c(C = 0,iC = 0,H = 0,iH = 0,N = 0,iN = 0,O = 0,iO = 0,F = 0 ,Na = 0,
                  Si = 0,P = 0,S = 0,Cl = 0,iCl = 0,Br = 0,iBr = 0,K = 0,iK = 0)
    
    comp_max[names(element_max)] <- element_max
    comp_min[names(element_min)] <- element_min
    
    mmu <- ppm/10^6*mass*1000
    res <- HR2(mass,comp_max,comp_min,mmu,charge,validation)	
    
    if (length(res) > 0) {
      res <- res %>%
        map(~{
          tibble(MF = .[2],
                 Mass = .[4] %>%
                   as.numeric() %>%
                   round(5),
                 `PPM Error` = .[4] %>%
                   as.numeric() %>%
                   round(5) %>%
                   ppmError(measured = mass))
        }) %>%
        bind_rows() %>%
        arrange(`PPM Error`)
      
      res$MF <- sapply(res$MF,function(x){
        mf <- count.elements(x)
        if (1 %in% mf) {
          mf <- sapply(names(mf),function(y,freq){
            freq <- freq[y]
            if (freq == 1) {
              y
            } else {
              paste(y,freq,sep = '')
            }
          },freq = mf)
          mf <- paste(mf,collapse = '')
          return(mf)
        } else {
          return(x)
        }
      })
    } else {
      res <- tibble(MF = character(),Mass = numeric(),`PPM Error` = numeric())
    }
    
  return(res)
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
#' @export

ipMF <- function(mz,
                 adduct = "[M+H]1+",
                 isotope = NA,
                 ppm = 1, 
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
  
  if (M < 200) {
    gr <- FALSE
  } else {
    gr <- TRUE
  }
  
  mfs <- generateMF(M,
             ppm = ppm,
             validation = gr) 
  
  if (nrow(mfs) > 0){
    mfs <- mfs %>%
      rename(`Theoretical M` = Mass) %>% 
      rowwise() %>% 
      mutate(`Measured m/z` = mz,
             `Measured M` = M,
             Adduct = adduct,
             Isotope = isotope,
             Score = MFscore(MF)) %>% 
      select(MF,
             Adduct,
             Isotope,
             `Measured m/z`,
             `Measured M`,
             `Theoretical M`,
             `PPM Error`,
             Score) %>% 
      arrange(Score)
  } else {
    mfs <- tibble(MF = character(),
                  Adduct = character(),
                  Isotope = character(),
                  `Measured m/z` = numeric(),
                  `Measured M` = numeric(),
                  `Theoretical M` = numeric(),
                  `PPM Error` = numeric(),
                  Score = numeric())
  }
  
  return(mfs)
}

#' Molecular formula scoring
#' @description Molecular formula scoring based on elemental ratios
#' @param mf Molecular formula
#' @examples
#' MFscore('C12H22O11')
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
