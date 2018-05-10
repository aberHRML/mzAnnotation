#' Putative Ionisation Product searching
#' @param mz the accurate m/z to search
#' @param ppm the parts per million threshold to search
#' @param mode either 'p' or 'n' specifiying the acquisition mode. Can be overrided using \code{add} to do mixed mode searches.
#' @param add a character \code{vector} specifying the adducts to search. If \code{NA}, all adducts for that acquisition mode will be used.
#' @param iso a character \code{vector} specifying the isotopes to search. If \code{NA} isotopes will not be searched.
#' @param adducts adduct table containing available adduct rules. Defaults to table returned by \code{availableAdducts()}.
#' @param isotopes isotope table containing available isotope rules. Defaults to table returned by \code{availableIsotopes()}.
#' @param DB object of class \code{MetaboliteDatabase}.
#' @export
#' @author  Jasen Finch
#' @importFrom dplyr bind_rows
#' @examples
#' res <- PIPsearch(132.03023,metaboliteDatabase(aminoAcids,descriptors(aminoAcids)),5,'[M-H]1-')

PIPsearch <- function(mz,db,ppm,adduct,isotope = NA,isotopes = mzAnnotation::Isotopes,adducts = mzAnnotation::Adducts){
  M <- calcM(mz,adduct = adduct,isotope = isotope,adducts = adducts,isotopes = isotopes)
  mr <- ppmRange(M,ppm)
  
  res <- db %>%
    filterMR(mr$lower,mr$upper)
  
  if (!is.na(isotope) & nrow(res@accessions) > 0){
    isoRule <- isotopes$Rule[isotopes$Isotope == isotope]
    res <- res %>%
      filterIR(isoRule)
  }
  
  addRule <- adducts$Rule[adducts$Name == adduct]
  
  res <- res %>%
    filterIP(addRule)
  
  res <- res %>%
  {left_join(.@accessions,.@descriptors,by = c("ACCESSION_ID", "SMILE"))} %>%
    select(ACCESSION_ID:Accurate_Mass) %>%
    mutate(Isotope = isotope,
           Adduct = adduct,
           `Measured m/z` = mz,
           `Theoretical m/z` = calcMZ(Accurate_Mass,adduct),
           `PPM Error` = ppmError(`Measured m/z`,`Theoretical m/z`)
    ) 
  return(res)
}