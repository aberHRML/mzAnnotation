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
