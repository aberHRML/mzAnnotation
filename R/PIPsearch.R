#' Putative Ionisation Product searching
#' @param mz the accurate m/z to search
#' @param ppm the parts per million threshold to search
#' @param mode either 'p' or 'n' specifiying the acquisition mode. Can be overrided using \code{add} to do mixed mode searches.
#' @param add a character \code{vector} specifying the adducts to search. If \code{NA}, all adducts for that acquisition mode will be used.
#' @param iso a character \code{vector} specifying the isotopes to search. If \code{NA} isotopes will not be searched.
#' @param adducts adduct table containing available adduct rules. Defaults to table returned by \code{availableAdducts()}.
#' @param isotopes isotope table containing available isotope rules. Defaults to table returned by \code{availableIsotopes()}.
#' @param DB Database to use. Defaults to \code{MZedDB}
#' @details The underlying database is that of MZedDB (\url{http://maltese.dbs.aber.ac.uk:8888/hrmet/index.html}). 
#' A list of available adducts can be found at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/search/disprules.php}. 
#' Isotopic adducts have also been added and include [2M+K41]1+, [M+K41]1+, [M+K41-2H]1- and [M+Cl37]1-. 
#' Available isotopes include 13C, 13C2, 13C3, 13C4, 18O, 37Cl, 41K and 34S.
#' @export
#' @author  Jasen Finch
#' @importFrom dplyr bind_rows
#' @examples
#' res <- PIPsearch(133.01378,5)

PIPsearch <-
function(mz, ppm = 5, mode = 'n', add = NA, iso = NA, adducts = mzAnnotation::Adducts, isotopes = mzAnnotation::Isotopes, DB = mzAnnotation::MZedDB){
  
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
	  lapply(adductList,function(a){
	    queryPIP(mz,ppm,a,i,adducts,isotopes,DB)
	  })
	}) %>%
	  map(bind_rows) %>%
	  bind_rows() %>%
	  select(ID,Name,MF,`Accurate Mass`,Smile = `Smile 1`,Adduct:`PPM Error`) %>%
	  mutate(`PPM Error` = round(`PPM Error`,5)) 
	

	return(res)
}
