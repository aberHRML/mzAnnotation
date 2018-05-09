#' MetaboliteDatabase
#' @export

setClass('MetaboliteDatabase',
         slots = list(
          type = 'character',
          accessions = 'tbl_df',
          descriptors = 'tbl_df',
          connection = 'list'
         ))