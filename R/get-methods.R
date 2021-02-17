#' Retrieve accessions
#' @rdname getAccessions
#' @description return accession data table from MetaboliteDatabase object
#' @param db MetaboliteDatabase
#' @export

setMethod('getAccessions',signature = 'MetaboliteDatabase',
          function(db){
            db@accessions[[1]]
          }
)

#' Retrieve descriptors
#' @rdname getDescriptors
#' @description return descriptor data table from MetaboliteDatabase object
#' @param db MetaboliteDatabase
#' @export

setMethod('getDescriptors',signature = 'MetaboliteDatabase',
          function(db){
            db@descriptors[[1]]
          }
)