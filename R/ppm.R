#' Calculate PPM error
#' @description Calculate ppmError between a measured and theoretical m/z
#' @param measured measured m/z
#' @param theoretical theoretical m/z
#' @examples ppmError(118.08626,118.08647)
#' @export

ppmError <- function(measured,theoretical) {
  error <- (measured - theoretical) / theoretical * 10^6
  return(error)
}
