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
    db@accessions <- list(accessions)
    db@descriptors <- list(descriptors)
  }
  if (type == 'remote'){
    if (!is.null(connection)){
      db@accessions <- list(tbl(connection,accessions))
      db@descriptors <- list(tbl(connection,descriptors))
    } else {
      stop('No database connection specified')
    }
  }
  return(db)
}