.onLoad <- function(libname, pkgname) {
  if (getOption('digits') < 10) options(digits = 10)
  
  invisible()
}
