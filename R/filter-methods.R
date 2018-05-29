
filterMR <- function(db,lower,upper){
  desc <- db@descriptors[[1]]
  desc <- desc %>%
    filter(Accurate_Mass > lower & Accurate_Mass < upper)
  acc <- db@accessions[[1]]
  acc <- acc %>%
    filter(SMILE %in% desc$SMILE)
  db@descriptors <- list(desc)
  db@accessions <- list(acc)
  return(db)
}


filterIR <- function(db,rule){
  ef <- elementFrequencies(db)
  if (str_extract(rule,'[:alpha:]') %in% colnames(ef)){
    ef <- ef %>%
      filter(eval(parse(text = rule)))   
  } else {
    ef[0,]
  }
  db@descriptors[[1]] <- db@descriptors[[1]] %>%
    filter(ACCESSION_ID %in% ef$ACCESSION_ID)
  db@accessions[[1]] <- db@accessions[[1]] %>%
    filter(ACCESSION_ID %in% ef$ACCESSION_ID)
  return(db)
}


filterIP <- function(db,rule){
  desc <- db@descriptors[[1]]
  desc <- desc %>%
    filter(eval(parse(text = rule)))
  acc <- db@accessions[[1]]
  acc <- acc %>%
    filter(SMILE %in% desc$SMILE)
  db@descriptors <- list(desc)
  db@accessions <- list(acc)
  return(db)
}

filterACCESSIONS <- function(db,ids){
  db@accessions[[1]] <- db@accessions[[1]] %>%
    filter(ACCESSION_ID %in% ids)
  db@descriptors[[1]] <- db@descriptors[[1]] %>%
    filter(ACCESSION_ID %in% ids)
  return(db)
}

filterMF <- function(db,MF){
  db@accessions[[1]] <- db@accessions[[1]] %>%
    filter(ACCESSION_ID %in% MF)
  db@descriptors[[1]] <- db@descriptors[[1]] %>%
    filter(ACCESSION_ID %in% MF)
  return(db)
}