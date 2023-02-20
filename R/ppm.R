#' Calculate ppm error
#' @description Calculate the parts per million error between a measured *m/z* and a theoretical *m/z*.
#' @param measured measured mass to charge ratio
#' @param theoretical theoretical mass to charge ratio
#' @return The parts per million error.
#' @examples ppmError(118.08626,118.08647)
#' @export

ppmError <- function(measured,theoretical) {
  error <- (measured - theoretical) / theoretical * 10^6
  return(error)
}
