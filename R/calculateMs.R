#' @importFrom tibble tibble rowid_to_column
#' @importFrom dplyr rowwise

calculateMs <- function(mz,add,iso,trans,adductTable = adducts(), isotopeTable = isotopes(), transformationTable = transformations()) {
  M <- map(trans,~{
    t <- .
    t <- map(iso,~{
      i <- .
      if (i == 'NA') {i <- NA}
      if (t == 'NA') {t <- NA}
      i <- tibble(`m/z` = mz, Adduct = add, Isotope = i, Transformation = t) %>%
        rowwise() %>%
        mutate(M = calcM(`m/z`,adduct = Adduct,isotope = Isotope,transformation = Transformation,adductTable = adductTable,isotopeTable = isotopeTable,transformationTable = transformationTable))
      return(i)
    })
    t <- bind_rows(t)
    return(t)
  })
  M <- bind_rows(M)
  M <- rowid_to_column(M, 'ID')
  return(M)
}