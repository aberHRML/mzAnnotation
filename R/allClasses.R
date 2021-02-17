# setOldClass('tbl_dbi')
# setClassUnion('metabolite_table',c('tbl_df','tbl_dbi'))

#' Metabolite database class
#' @export

setClass('MetaboliteDatabase',
         slots = list(
          type = 'character',
          accessions = 'list',
          descriptors = 'list',
          connection = 'list'
         ))