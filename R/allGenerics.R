#' @rdname filterMR
setGeneric("filterMR", function(db,lower,upper) {
  standardGeneric("filterMR")
})

#' @rdname filterIP
setGeneric("filterIP", function(db,rule) {
  standardGeneric("filterIP")
})

#' @rdname filterER
setGeneric("filterER", function(db,rule) {
  standardGeneric("filterER")
})

#' @rdname filterACCESSIONS
setGeneric("filterACCESSIONS", function(db,ids) {
  standardGeneric("filterACCESSIONS")
})

setGeneric("elementFrequencies", function(db) {
  standardGeneric("elementFrequencies")
})

#' @rdname filterMF
setGeneric('filterMF', function(db,mf){
  standardGeneric('filterMF')
})

#' @rdname getAccessions
setGeneric('getAccessions',function(db) {
  standardGeneric('getAccessions')
})

#' @rdname getDescriptors
setGeneric('getDescriptors',function(db) {
  standardGeneric('getDescriptors')
})

#' @rdname calcAdducts
setGeneric('calcAdducts',function(db,id,adductTable = adducts()) {
  standardGeneric('calcAdducts')
})

#' @rdname PIPsearch
setGeneric('PIPsearch',function(db,mz,ppm,adduct,isotope = NA, isotopeTable = isotopes(), adductTable = adducts()) {
  standardGeneric('PIPsearch')
})

