

checkTableNames <- function(tab,type,correct_names){
  table_names <- colnames(tab)

  if (!all(correct_names %in% table_names)) {
    stop('Columns '
         ,paste0(columns$column[!(columns$column %in% nam)],collapse = ', '),
         ' not found in ',
         type,
         ' rules table!',
         call. = FALSE)
  }
  
  invisible()
}

checkAdductTable <- function(tab,
                             type = 'adduct',
                             correct_names = colnames(adduct_rules())){
  checkTableNames(tab,type,correct_names)
}

checkIsotopeTable <- function(tab
                              type = 'isotope',
                              correct_names = colnames(isotope_rules())){
  checkTableNames(tab,type,correct_names)
}

checkTransformationTable <- function(tab,
                                     type = 'transformation',
                                     correct_names = colnames(transformation_rules())){
  checkTableNames(tab,type,correct_names)
}

checkTransformation <- function(transformation,
                                transformation_rules_table = transformation_rules()){
  checkTransformationTable(transformation_rules_table)
  
  if (!(transformation %in% transformation_rules_table$`MF Change`)){
    stop('Specified transformation not found in transformation rules table.',
         call. = FALSE)
  }
}

checkEntries <- function(entries){
  necessary_names <- c('ID','NAME','SMILES')
  
  entry_names <- colnames(entries)
  
  presence <- necessary_names %in% entry_names
  
  if (FALSE %in% presence) {
    stop(paste0('The table of metabolite entries should contain the following column names: ',
          paste(necessary_names,collapse = ', '),
          '.'))
  }
  
  invisible()
}
