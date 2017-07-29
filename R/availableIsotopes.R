#' availableIsotopes
#' @description returns isotopes available in mzAnnotation
#' @importFrom tibble as_tibble
#' @export

availableIsotopes <- function(){
  return(as_tibble(MZedDB$ISOTOPE_RULES))
}

