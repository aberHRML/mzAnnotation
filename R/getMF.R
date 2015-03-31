#' Molecular formula generator

getMF <-
function(mass,ppm=1,charge=1,applygr=TRUE,latmin=c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),latmax=c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0)){
	prec <- ppm/10^6*mass*1000
	resp <- mfGen(mass,latmax,latmin,prec,charge,applygr)	
	if(length(resp)>0){
	  resp <- matrix(unlist(resp), nrow = length(resp),byrow=T)
	  resp[,7] <- sapply(resp[,4],function(x,mass){x <- as.numeric(x);x <- (x-mass)/mass*10^6; return(round(x,5))},mass=mass)
	  resp <- cbind(rep(mass,nrow(resp)),resp)
    resp <- resp[,-(6:7)]
	  colnames(resp) <- c("Measured m/z", "Clean MF", "MF", "RDB", "m/z","PPM Error")
  } else {
		resp <- "No MF hits"
	}
	return(resp)
}
