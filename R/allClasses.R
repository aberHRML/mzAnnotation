setOldClass('tbl_dbi')
setOldClass('tbl_df')
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