#' calcAccurateMass
#' @description calculate the accurate mass of a given molecular formula
#' @param MF molecular formula for which to calculate the accrate mass
#' @param charge charge of the given molecular formula
#' @param elements element information table. Defaults to \code{\link{Elements}}.
#' @examples 
#' calcAccurateMass('C4H5O5',charge = 0)
#' @export

calcAccurateMass <- function(MF,charge = 0, elements = mzAnnotation::Elements) {
  e <- elements$AtomicMass[elements$Element == 'e']
  elementFreq <- tibble(Element = names(makeup(MF)),Frequency = makeup(MF))
  elementMasses <- elements %>%
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