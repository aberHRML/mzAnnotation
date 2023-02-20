

checkTableNames <- function(tab,type,correct_names){
  table_names <- colnames(tab)

  if (!all(correct_names %in% table_names)) {
    stop('Columns '
         ,paste0(correct_names[!(correct_names %in% table_names)],collapse = ', '),
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

checkIsotopeTable <- function(tab,
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
