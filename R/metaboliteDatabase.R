#' metaboliteDatabase
#' @description Build a metabolite database ready for use.
#' @examples 
#' db <- metaboliteDatabase(aminoAcids,descriptors(aminoAcids$SMILE))
#' @export

metaboliteDatabase <- function(accessions,descriptors,connection = NULL,type = 'local'){
  db <- new('MetaboliteDatabase')
  db@type <- type
  if (type == 'local'){
    db@accessions <- accessions
    db@descriptors <- descriptors
  }
  if (type == 'remote'){
    if (!is.null(connection)){
      
    } else {
      stop('No database connection specified')
    }
  }
  return(db)
}