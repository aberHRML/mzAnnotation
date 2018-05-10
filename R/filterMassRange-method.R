
filterMassRange <- function(db,lower,upper){
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