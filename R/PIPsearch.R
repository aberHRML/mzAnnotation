#' Putative Ionisation Product searching
#' @param mz the accurate m/z to search.
#' @param db object of class \code{MetaboliteDatabase}.
#' @param ppm the parts per million threshold to search.
#' @param adduct the adduct name to search.
#' @param isotope the isotope name to search. Defaults to NA for non-isotopic searches.
#' @param adducts adduct table containing available adduct rules. Defaults to table returned by \code{availableAdducts()}.
#' @param isotopes isotope table containing available isotope rules. Defaults to table returned by \code{availableIsotopes()}.
#' @export
#' @author  Jasen Finch
#' @importFrom dplyr bind_rows select filter
#' @importFrom magrittr %>%
#' @examples
#' res <- PIPsearch(132.03023,metaboliteDB(aminoAcids,descriptors(aminoAcids)),5,'[M-H]1-')

PIPsearch <- function(mz,db,ppm,adduct,isotope = NA, isotopes = mzAnnotation::Isotopes, adducts = mzAnnotation::Adducts){
  M <- calcM(mz,adduct = adduct,isotope = isotope,adducts = adducts,isotopes = isotopes)
  mr <- ppmRange(M,ppm)
  
  res <- db %>%
    filterMR(mr$lower,mr$upper)
  
  if (!is.na(isotope) & nrow(res@accessions[[1]]) > 0) {
    isoRule <- isotopes$Rule[isotopes$Isotope == isotope]
    res <- res %>%
      filterIR(isoRule)
  }
  
  addRule <- adducts$Rule[adducts$Name == adduct]
  
  res <- res %>%
    filterIP(addRule)
  
  res <- res %>%
  {left_join(.@accessions[[1]],.@descriptors[[1]],by = c("ACCESSION_ID", "SMILE"))} %>%
    select(ACCESSION_ID:Accurate_Mass) %>%
    mutate(Isotope = isotope,
           Adduct = adduct,
           `Measured m/z` = mz,
           `Theoretical m/z` = calcMZ(Accurate_Mass,adduct),
           `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)
    ) 
  return(res)
}