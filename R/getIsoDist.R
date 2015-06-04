
  
getIsoDist <- function(){
	mf_res <- apply(explan.mz,getMF,ppm,charge=char,latmin=latmin,latmax=latmax,applygr=applygr)
	mf_res <- filterMF(mf_res)
	iso <- grep("i",mf_res[,"Clean MF"],fixed=TRUE)
	if(length(iso)>0){
		mf <- mf[-iso]
	}
	if (length(mf) > 0){
		## get isotopic distributions
		iso_ab <- lapply(mf,isoDistr,chrg=char)
		iso_ab <- filterIso(iso_ab)
		names(iso_ab) <- mf
	}
}