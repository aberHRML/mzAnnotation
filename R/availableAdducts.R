#' availableAdducts
#' @description returns adducts rules available in mzAnnotation
#' @importFrom tibble as_tibble
#' @export

availableAdducts <- function(){
  return(as_tibble(MZedDB$ADDUCT_FORMATION_RULES))
}

