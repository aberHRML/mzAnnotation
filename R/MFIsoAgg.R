#' Molecular Formula and Isotope Distribution aggregation

MFIsoAgg <-
function(mz,mo,ppm,latmin=c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),latmax=c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),...){
  ## get molecular formulas
  if (mo=="p"){
		char <- 1
	}
	if (mo=="n"){
		char <- -1
	}	
	if (mz<180){
		applygr <- FALSE
	} else {
		applygr <- TRUE
	}
	mf_res <- getMF(mz,ppm,charge=char,latmin=latmin,latmax=latmax,applygr=applygr)
	mf_res <- filterMF(mf_res)
	if(nrow(mf_res)==0){
			mf_res <- "No MF hits"
	} else {
		iso <- grep("i",mf_res[,"Clean MF"],fixed=TRUE)
		if(length(iso)>0){
				mf <- mf[-iso]
			}
			if (length(mf) > 0){
				## get isotopic distributions
				iso_ab <- lapply(mf,iso.distr,chrg=char)
				iso_ab <- filter.iso(iso_ab)
				names(iso_ab) <- mf
				 ## calculate measured relative isotope intensities for C13
		}
	}
	## get mzeddb hits for accurate mass
	cat("\n","Querying MZedDB..."); flush.console()
	mz.hits <- mzeddb_get(as.numeric(m[1,2]),mo,ppm,MZedDB.2)	
	res <- list(ACCURATE_MASS=m[1,2],CORRELATIONS=cor,MF=mf_res,ISO=iso_ab,M_ISO=all.max,SUM_ISO=sum.max,MZEDDB_HITS=mz.hits)
	cat("\n","...done","\n",sep=""); flush.console()
	return(MFres)
}
