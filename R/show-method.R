#' show-MetaboliteDatabase
#' @description show method for MetaboliteDatabase class.
#' @param object S4 object of class MetaboliteDatabase
#' @importFrom methods show
#' @export

setMethod('show',signature = 'MetaboliteDatabase',
          function(object){
            type <- object@type
            accessions <- object %>%
              getAccessions() %>%
              nrow()
            cat('\n',type,' MetaboliteDatabase object containing ',accessions,' accessions\n\n',sep = '')
          }
)