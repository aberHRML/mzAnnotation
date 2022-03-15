
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
