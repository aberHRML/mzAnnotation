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
#' #res <- PIPsearch(133.01378,5)

PIPsearch <- function(mz, DB, ppm = 5, mode = 'n', add = NA, iso = NA, adducts = mzAnnotation::Adducts, isotopes = mzAnnotation::Isotopes){
  
  if (class(DB) != 'MetaboliteDatabase'){
    stop('DB needs to be of class "MetaboliteDatabase"!')
  }
  
  if (T %in% is.na(add)) {
	  if (mode == "p") {
	    adductList <- adducts$Name[adducts$Nelec < 0] 
	  }
	  if (mode == "n") {
	    adductList <- adducts$Name[adducts$Nelec > 0]
	  }
    if (mode == "ne") {
      adductList <- c("M") 
    }
  } else {
    adductList <- add
  }
  if (!(T %in% is.na(iso))) {
    iso <- c(NA,iso)
  }
	
	res <- map(iso,~{
	  i <- .
	  map(adductList,~{
	    queryPIP(mz,ppm,.,i,adducts,isotopes,DB)
	  })
	}) %>%
	  map(bind_rows) %>%
	  bind_rows()

	return(res)
}
