# setOldClass('tbl_dbi')
# setClassUnion('metabolite_table',c('tbl_df','tbl_dbi'))

#' MetaboliteDatabase
#' @export

setClass('MetaboliteDatabase',
         slots = list(
          type = 'character',
          accessions = 'list',
          descriptors = 'list',
          connection = 'list'
         ))