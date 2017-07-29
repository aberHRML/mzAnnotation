#' availableTransfromations
#' @description returns transformations available in mzAnnotation
#' @importFrom tibble as_tibble
#' @export

availableTransformations <- function(){
  return(as_tibble(MZedDB$BIOTRANSFORMATION_RULES))
}

