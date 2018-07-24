
setGeneric("filterMR", function(db,lower,upper) {
  standardGeneric("filterMR")
})

setGeneric("filterIP", function(db,rule) {
  standardGeneric("filterIP")
})

setGeneric("filterIR", function(db,rule) {
  standardGeneric("filterIR")
})

setGeneric("filterACCESSIONS", function(db,ids) {
  standardGeneric("filterACCESSIONS")
})

setGeneric("elementFrequencies", function(db) {
  standardGeneric("elementFrequencies")
})

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