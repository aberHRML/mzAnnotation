MF_iso_agg <-
function(m,acc_mat,mo,ppm,cor,MZedDB.2,lmin=c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),lmax=c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),...){	
	cat("\n",as.character(m[1,1]),sep=""); flush.console()
	if (mo=="p"){
		char=1
	}
	if (mo=="n"){
		char=-1
	}	
	## get molecular formulas
	if (m[1,2]<180){
		appgr <- FALSE
	} else {
		appgr <- TRUE
	}
	mf_res <- getmf_3(m[1,2],ppm,charge=char,latmin=lmin,latmax=lmax,adduct=1,applygr=appgr)
	if(mf_res=="No MF hits"){
		iso_ab <- "No theoretical isotope distributions"
		all.max <- "No measured isotope distributions"
		sum.max <- "No composition prediction"
		cat("\n","No MF's generated"); flush.console()
	} else {
		mf_res <- filter.MF(mf_res)
		if(nrow(mf_res)==0){
			iso_ab <- "No theoretical isotope distributions"
			all.max <- "No measured isotope distributions"
			sum.max <- "No composition prediction"
			mf_res <- "No MF hits"
			cat("\n","No MF's generated"); flush.console()
		} else {
			mf <- as.character(mf_res[,3])
			if (length(mf)==1){
				cat("\n",length(mf), "MF generated...",sep=" "); flush.console()
			} else {
			cat("\n",length(mf), "MF's generated...",sep=" "); flush.console()
			}
			iso <- grep("i",mf_res[,2],fixed=TRUE)
			if(length(iso)>0){
				mf <- mf[-iso]
			}
			if (length(mf) > 0){
				## get isotopic abundances
				iso_ab <- lapply(mf,iso.distr,chrg=char)
				iso_ab <- filter.iso(iso_ab)
				if (length(iso_ab)==1){
					cat("\n","theoretical isotope distributions calculated for",length(iso_ab),"MF...",sep=" "); flush.console()
				} else {
					cat("\n","theoretical isotope distributions calculated for",length(iso_ab),"MF's...",sep=" "); flush.console()
				}
				names(iso_ab) <- mf
				 ## calculate measured relative isotope intensities for C13
				all.max <- meas_iso_abun(m[1,2],acc_mat,mo)
				if (length(all.max)>1){
					all.max.1 <- na.omit(all.max)
					all.max.1 <- all.max.1[all.max.1$No_Element>0,]
					all.max.1 <- all.max.1[is.finite(all.max.1$No_Element),]
					if(sum(all.max.1$No_Element)>0){
						sum.max <- summarySE(data=all.max.1,measurevar="No_Element")
						names(sum.max)[3] <- "Mean"
						cat("\n","measured isotope distributions calculated..."); flush.console()
			 	 } else {
			 		 	all.max <- "No measured isotope distributions"
					  sum.max <- "No composition prediction"
			 	 }
				} else {
					all.max <- "No measured isotope distributions"
					sum.max <- "No composition prediction"
				}
			} else {
				iso_ab <- "No theoretical isotope distributions"
				all.max <- "No measured isotope distributions"
				sum.max <- "No composition prediction"
			}
		}
	}
	## get mzeddb hits for accurate mass
	cat("\n","Querying MZedDB..."); flush.console()
	mz.hits <- mzeddb_get(as.numeric(m[1,2]),mo,ppm,MZedDB.2)	
	res <- list(ACCURATE_MASS=m[1,2],CORRELATIONS=cor,MF=mf_res,ISO=iso_ab,M_ISO=all.max,SUM_ISO=sum.max,MZEDDB_HITS=mz.hits)
	cat("\n","...done","\n",sep=""); flush.console()
	return(res)
}
