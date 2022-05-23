
checkAdductTable <- function(tab){
  result <- TRUE
  nam <- colnames(tab)
  columns <- tibble(column = c('Name',
                               'Charge', 
                               'xM',
                               'Add',
                               'Nelec',
                               'AddAt',
                               'RemAt', 
                               'AddEx', 
                               'RemEx', 
                               'Rule', 
                               'Default'),
                    type = c('character', 
                             'integer', 
                             'integer', 
                             'numeric',
                             'integer',
                             'character',
                             'character',
                             'character',
                             'character',
                             'character',
                             'numeric'),
               )
  if (F %in% (columns$column %in% nam)) {
    cat('Columns ',str_c(columns$column[!(columns$column %in% nam)],collapse = ', '),' not found!')
    result <- FALSE
  }
  
  return(result)
}

checkIsotopeTable <- function(tab){
  result <- TRUE
  return(result)
}

checkTransformationTable <- function(tab){
  result <- TRUE
  return(result)
}

checkAccessionTable <- function(tab){
  result <- TRUE
  return(result)
}

checkDescriptorsTable <- function(tab){
  result <- TRUE
  return(result)
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
