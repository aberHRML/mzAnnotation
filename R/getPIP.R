
getPIP <-
function(mz,mode,ppm,add=NULL,iso=F,bio=F){
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
  mass_low <- mz - (mz * ppm * 10^-6)
	mass_high <- mz + (mz * ppm * 10^-6)
	data("MZedDB")
	res <- lapply(adducts,queryPIP,mass_low=mass_low,mass_high=mass_high,mz=mz,MZedDB=MZedDB)
	res <- lapply(res,function(x){names(x) <- c("ID","Name","MF","Accurate Mass","Smiles","Adduct","Adduct m/z","PPM Error");return(x)})
	res <- ldply(res,data.frame)
	return(res)
}
