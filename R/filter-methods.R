
filterMR <- function(db,lower,upper){
  desc <- db@descriptors
  desc <- desc %>%
    filter(Accurate_Mass > lower & Accurate_Mass < upper)
  acc <- db@accessions
  acc <- acc %>%
    filter(SMILE %in% desc$SMILE)
  db@descriptors <- desc
  db@accessions <- acc
  return(db)
}


filterIR <- function(db,rule){
  ef <- elementFrequencies(db)
  ef <- ef %>%
    filter(eval(parse(text = rule)))
  db@descriptors <- db@descriptors %>%
    filter(ACCESSION_ID %in% ef$ACCESSION_ID)
  db@accessions <- db@accessions %>%
    filter(ACCESSION_ID %in% ef$ACCESSION_ID)
  return(db)
}


filterIP <- function(db,rule){
  desc <- db@descriptors
  desc <- desc %>%
    filter(eval(parse(text = rule)))
  acc <- db@accessions
  acc <- acc %>%
    filter(SMILE %in% desc$SMILE)
  db@descriptors <- desc
  db@accessions <- acc
  return(db)
}

filterACCESSIONS <- function(db,ids){
  db@accessions <- db@accessions %>%
    filter(ACCESSION_ID %in% ids)
  db@descriptors <- db@descriptors %>%
    filter(ACCESSION_ID %in% ids)
  return(db)
}

filterMF <- function(db,MF){
  db@accessions <- db@accessions %>%
    filter(ACCESSION_ID %in% MF)
  db@descriptors <- db@descriptors %>%
    filter(ACCESSION_ID %in% MF)
  return(db)
}