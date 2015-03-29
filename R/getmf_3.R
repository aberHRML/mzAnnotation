getmf_3 <-
function(mass,ppm,charge=1,adduct=1,applygr=TRUE,latmin=c(C=0,iC=0,H=0,iH=0,N=0,iN=0,O=0,iO=0,F=0,Na=0,Si=0,P=0,S=0,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0),latmax=c(C=41,iC=0,H=72,iH=0,N=15,iN=0,O=30,iO=0,F=0,Na=0,Si=0,P=2,S=2,Cl=0,iCl=0,Br=0,iBr=0,K=0,iK=0)){
	prec <- ppm/10^6*mass*1000
	resp=mf_gen(mass,prec,charge,applygr,mini= latmin,maxi = latmax)
	if(length(resp)>0){
		resp=cbind(rep(mass,nrow(resp)),resp)
		resp <- resp[,c(1,2,3,5,6,7)]
		if (class(resp)=="character"){
			resp <- matrix(resp,nrow=1)
		}
		resp[,6] <- sapply(seq(1,nrow(resp)),function(x,meamass,mass){meamass <- as.numeric(meamass);mass <- as.numeric(mass);y <- sqrt(((meamass[x]-mass[x])/mass[x]*10^6)^2);return(y)},resp[,1],resp[,4])
		resp <- data.frame(as.numeric(resp[,1]),resp[,2],resp[,3],as.numeric(resp[,4]),as.numeric(resp[,5]),as.numeric(resp[,6]))
		resp[,6] <- round(resp[,6],digits=2)
		colnames(resp)<-c("Measured Mass","Original MF","Clean MF","Mass","No.Isotopes","PPM Error")
		resp <- resp[order(resp[,6]),]
	} else {
		resp <- "No MF hits"
	}
	return(resp)
}
