#' metaboliteDatabase
#' @description Build a metabolite database ready for use.
#' @param accessions
#' @param descriptors
#' @param connection
#' @param type
#' @examples 
#' db <- metaboliteDatabase(aminoAcids,descriptors(aminoAcids$SMILE))
#' @importFrom dplyr tbl
#' @importFrom purrr map_chr
#' @importFrom methods new
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
      db@accessions <- tbl(connection,accessions)
      db@descriptors <- tbl(connection,descriptors)
    } else {
      stop('No database connection specified')
    }
  }
  return(db)
}