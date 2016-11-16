#' Putative Ionisation Product searching
#' @param mz the accurate m/z to search
#' @param mode either 'p' or 'n' specifiying the acquisition mode
#' @param ppm the parts per million threshold to search
#' @param add a character \code{vector} specifying the adducts to search. If \code{NULL} a default selection will be used.
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
  if (is.null(add)) {
	  if (mode == "p" & is.null(add)) {
				  adducts <- c("[M+2Na]2+","[M+H+K]2+","[M+H+Na]2+","[M+H+NH4]2+","[M+2H]2+","[M+2H-H2O-NH3]2+","[2M+K]1+","[2M+K41]1+","[2M+Na]1+","[2M+NH4]1+",
				         "[2M+H]1+","[M+2K-H]1+","[M+2Na-H]1+","[M+K]1+","[M+K41]1+","[M+K41]1+","[M+Na]1+","[M+H]1+","[M1+.]1+","[M+H-NH3]1+","[M+H-H2O]1+",
				         "[M+H-FA]1+")
	  }
	  if (mode == "n" & is.null(add)) {
		  	adducts <- c("[2M+Na-2H]1-","[2M-H]1-","[M+K-2H]1-","[M+K41-2H]1-","[M+Cl]1-",'[M+Cl37]1-',"[M+Na-2H]1-","[M1-.]1-","[M-H]1-","[M-2H]2-")
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
