#' Molecular formula generator

getMF <-
function(explan.mz,ppm,mode,MFfilter=T,lat=list(c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0))){
	if (mode=="p"){
		charge <- 1
	}
	if (mode=="n"){
		charge <- -1
	}	
	if (explan.mz[1]<180){
		applygr <- FALSE
	} else {
		applygr <- TRUE
	}
  prec <- ppm/10^6*explan.mz[1]*1000
	resp <- mfGen(explan.mz[1],lat[[2]],lat[[1]],prec,charge,applygr)	
	if(length(resp)>0){
	  resp <- matrix(unlist(resp), nrow = length(resp),byrow=T)
	  resp[,7] <- sapply(resp[,4],function(x,mass){x <- as.numeric(x);x <- (x-mass)/mass*10^6; return(round(x,5))},mass=explan.mz[[1]])
	  resp <- cbind(rep(explan.mz[[1]],nrow(resp)),resp)
    resp <- resp[,-c(6:7)]
    resp <- matrix(resp,ncol=6)
	} else {
	  resp <- matrix(ncol=6,nrow=0)
	}
	colnames(resp) <- c("Measured m/z", "Clean MF", "MF", "RDB", "m/z","PPM Error")
	if (MFfilter==T){
	  resp <- filterMF(resp)
	}
	return(data.frame(resp))
}
