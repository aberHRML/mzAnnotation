getmf_2 <-
function(mass,ppm,charge=1,adduct=1,applygr=TRUE,hr="./myhr.exe",latms=NULL){
	if(is.null(latms))
	 latms="-C 0-41 -H 0-72 -1 0-2 -K 0-2 -A 0-2 -P 0-4 -S 0-4 -N 0-10 -O 0-15 -J 0-2"
	
	if(!applygr)
	 latms=paste(latms,"-g")

	prec <- ppm/10^6*mass*1000
	cmd=paste(hr,"-t",prec,latms," -s -c",charge,"-a",adduct,"-m",mass)
	resp=system(cmd,intern=T)
	if(length(resp)>0){
		if(class(resp)=="character"){
			resp <- matrix(data=resp,nrow=1)
		}
		resp<-do.call("rbind",lapply(resp,function(x) strsplit(x,"\t")[[1]]))
		lso=sort(abs(as.numeric(resp[,3])-mass),ind=T)$ix
		resp=resp[lso,c(1,2,4,6,7),drop=F]
		resp=cbind(rep(mass,nrow(resp)),resp)
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
