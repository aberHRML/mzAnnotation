#' Calculate molecular formula accurate mass
#' @description calculate the accurate mass of a given molecular formula
#' @param MF molecular formula for which to calculate the accrate mass
#' @param charge charge of the given molecular formula
#' @param elementTable element information table. Defaults to \code{\link{elements}()}.
#' @examples 
#' calcAccurateMass('C4H5O5',charge = 0)
#' @importFrom CHNOSZ count.elements
#' @export

calcAccurateMass <- function(MF,charge = 0, elementTable = elements()) {
  e <- elementTable$AtomicMass[elementTable$Element == 'e']
  elementFreq <- tibble(Element = names(count.elements(MF)),Frequency = count.elements(MF))
  elementMasses <- elementTable %>%
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