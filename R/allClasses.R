setOldClass('tbl_dbi')
setClassUnion('metabolite_table',c('tbl_df','tbl_dbi'))

#' MetaboliteDatabase
#' @export

setClass('MetaboliteDatabase',
         slots = list(
          type = 'character',
          accessions = 'metabolite_table',
          descriptors = 'metabolite_table',
          connection = 'list'
         ))