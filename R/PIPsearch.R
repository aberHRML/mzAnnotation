#' Putative Ionisation Product searching
#' @param mz the accurate m/z to search
#' @param mode either 'p' or 'n' specifiying the acquisition mode
#' @param ppm the parts per million threshold to search
#' @param add a character \code{vector} specifying the adducts to search. If \code{NULL}, all adducts for that acquisition mode will be used.
#' @param iso a character \code{vector} specifying the isotopes to search. If \code{NULL} isotopes will not be searched.
#' @details The underlying database is that of MZedDB (\url{http://maltese.dbs.aber.ac.uk:8888/hrmet/index.html}). 
#' A list of available adducts can be found at \url{http://maltese.dbs.aber.ac.uk:8888/hrmet/search/disprules.php}. 
#' Isotopic adducts have also been added and include [2M+K41]1+, [M+K41]1+, [M+K41-2H]1- and [M+Cl37]1-. 
#' Available isotopes include C13, 2C13, 3C13, 4C13, O18, Cl37, K41 and S34.
#' @export
#' @author  Jasen Finch
#' @importFrom plyr ldply
#' @examples
#' res <- PIPsearch(133.01378,'n',5)

PIPsearch <-
function(mz,mode,ppm,add = NULL, iso = NULL){
  adducts <- MZedDB$ADDUCT_FORMATION_RULES
  if (is.null(add)) {
	  if (mode == "p" & is.null(add)) {
	    adducts <- adducts$Name[adducts$Nelec < 0] 
	  }
	  if (mode == "n" & is.null(add)) {
	    adducts <- adducts$Name[adducts$Nelec > 0]
	  }
    if (mode == "ne") {
      adducts <- c("M") 
    }
  } else {
    adducts <- add
  }
  if (!is.null(iso)) {
    iso <- c('',iso)
  } else {
    iso <- ''
  }
	
	## Loop over adducts irrespective of isotopes
	res <- lapply(iso,function(iso,adducts,mz,ppm,bio){
	  if (iso == '') {
	    iso <- NULL
	  }
	  lapply(adducts,queryPIP,mz = mz,ppm = ppm,iso = iso)
	  },adducts = adducts,mz = mz,ppm = ppm)
	res <- lapply(res,ldply,stringsAsFactors = F)
	res <- ldply(res,stringsAsFactors = F)
	names(res) <- c("ID","Name","MF","Accurate Mass","Smiles","Adduct",'Isotope',"Adduct m/z","PPM Error")

	return(res)
}
