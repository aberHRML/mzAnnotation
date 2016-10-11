
getPIP <-
function(mz,mode,ppm,add=NULL,iso=NULL,bio=NULL){
  if(!is.null(add)){
	  if (mode=="p" & is.null(add)){
				  adducts <- c("[M+2Na]2+","[M+H+K]2+","[M+H+Na]2+","[M+H+NH4]2+","[M+2H]2+","[M+2H-H2O-NH3]2+","[2M+K]1+","[2M+Na]1+","[2M+NH4]1+",
				         "[2M+H]1+","[M+2K-H]1+","[M+2Na-H]1+","[M+K]1+","[M+Na]1+","[M+Na]1+","[M+H]1+","[M1+.]1+","[M+H-NH3]1+","[M+H-H2O]1+",
				         "[M+H-FA]1+")
	  }
	  if (mode=="n" & is.null(add)){
		  	adducts <- c("[2M+Na-2H]1-","[2M-H]1-","[M+K-2H]1-","[M+Cl]1-","[M+Na-2H]1-","[M1-.]1-","[M-H]1-","[M-2H]2-")
	  }
    if (mode=="ne"){
      adducts <- c("M") 
    }
  } else {
    adducts <- add
  }
  if(!is.null(iso)){
    iso <- c('',iso)
  }
	data("MZedDB")
	
	## Loop over adducts irrespective of isotopes
	res <- lapply(iso,function(iso,adducts,mz,ppm,MZedDB,bio){
	  if(iso==''){
	    iso <- NULL
	  }
	  lapply(adducts,queryPIP,mz=mz,ppm=ppm,MZedDB=MZedDB,iso=iso,bio=bio)
	  },adducts=adducts,mz=mz,ppm=ppm,MZedDB=MZedDB,bio=bio)
	res <- lapply(res,ldply)
	res <- ldply(res)
	#names(res) <- c("ID","Name","MF","Accurate Mass","Smiles","Adduct",'Isotope','Biotransformation',"Adduct m/z","PPM Error")

	return(res)
}
