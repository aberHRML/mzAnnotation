mZedDB <-
function(rand_mz,ppm,Path){    # wrapper function for mZedDB searches of feature lists in rand_mz at specified ppm
	for(x in 1:length(ppm)) {
		for (i in 1:length(rand_mz)){
			list.mass <- rand_mz[[i]]
			savefile <- paste(Path,paste(names(rand_mz)[i],"_MS_",ppm[x],".csv",sep=""),sep="/")
		
			if (names(rand_mz)[i]=="pos"){
				add=c("[M]+","[M+H]1+","[M+NH4]1+","[M+Na]1+","[M+K]1+","[M+H-NH3]1+","[M+H-FA]1+","[M+H-H2O]1+","[2M+H]1+","[2M+NH4]1+",
				"[2M+Na]1+","[2M+K]1+","[M+2H]2+","[M+2Na]2+","[M+2K-H]1+","[M+2Na-H]1+","[M+H+Na]2+","[M+H+K]2+","[M+H+NH4]2+")
			} 
			if (names(rand_mz)[i]=="neg"){
			add=c("[M-H]1-","[M+Na-2H]1-","[M+Cl]1-","[M+K-2H]1-","[M]1-","[2M-H]1-","[M-2H]2-")
			}
		
  	  resMS=getmany(host="maltese.dbs.aber.ac.uk",accmas=list.mass,add=add,massacc=ppm[x],printlong=TRUE)  ##Be sure to change mode when needed
  	  write.table(resMS,file=savefile,col.names = TRUE,sep=",")
  		if (length(resMS[,1]) > 0 ) {
  		   umz=as.factor(unique(resMS[,1]))
   		 	 xy<-NULL
   	  for(i in c(1:length(umz)))  {
   	     xy[i]<-which(list.mass==umz[i])
   	   }
   	 list.mass=list.mass[-xy]
     }
	 }
  }
}
