#' Create a metabolite database
#' @description Build a metabolite database ready for use.
#' @param accessions tibble containing accession information. If \code{type = 'remote'} this should be the name of the table containing the accession information within the SQL database.
#' @param descriptors tibble containing descriptor information as returned by \code{descriptors()}. If \code{type = 'remote'} this should be the name of the table containing the descriptor information within the SQL database.
#' @param connection If \code{type = 'remote'} this should be a valid database connection as returned by \code{DBI::dbconnect()}.
#' @param type set to either \code{'local'} for in-memory databases or \code{remote} for SQL database connections.
#' @examples 
#' db <- metaboliteDB(amino_acids,descriptors(amino_acids$SMILES))
#' @importFrom dplyr tbl
#' @importFrom purrr map_chr
#' @importFrom methods new
#' @export

metaboliteDB <- function(accessions,descriptors,connection = NULL,type = 'local'){
  db <- new('MetaboliteDatabase')
  db@type <- type
  if (!(type %in% c('local','remote'))){
    stop('type not recognised!',call. = FALSE)
  }
  
  if (nrow(accessions) != nrow(descriptors)) {
   stop('Number of rows in accessions and descriptors do not match')  
  }
  
  if (type == 'local'){
    db@accessions <- list(accessions %>% as_tibble())
    db@descriptors <- descriptors %>% 
      as_tibble() %>%
      bind_cols(db %>%
                  getAccessions() %>%
                  select(ID)) %>%
      list(.)
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
