#' ppmRange
#' @description Calculate the upper and lower ppm boundaries for a given m/z
#' @param mz the m/z for which to calculate the range
#' @param ppm the ppm 
#' @examples ppmRange(118.08626,5)
#' @export

ppmRange <- function(mz,ppm) {
  amu <- mz * ppm * 10^-6
  lower <- round(mz - amu,5)
  upper <- round(mz + amu,5)
  return(c(lower = lower, upper = upper))
}