.onLoad <- function(libname,pkgname) {
  if(.Platform$OS == "windows")
    Sys.setenv(BABEL_DATADIR=system.file("openbabel_data",package=pkgname))
  
}