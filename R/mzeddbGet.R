mzeddb_get <-
function(acc.mass,mode,ppm,MZedDB){
	if (mode=="p"){
				add=c("[M1+.]1+","[M+H]1+","[M+NH4]1+","[M+Na]1+","[M+K]1+","[M+H-FA]1+","[M+H-H2O]1+","[2M+H]1+","[2M+NH4]1+",
				"[2M+Na]1+","[2M+K]1+","[M+2H]2+","[M+2Na]2+","[M+2K-H]1+","[M+2Na-H]1+","[M+H+Na]2+","[M+H+K]2+","[M+H+NH4]2+")
			} 
	if (mode=="n"){
			add=c("[M-H]1-","[M+Na-2H]1-","[M+Cl]1-","[M+K-2H]1-","[M1-.]1-","[2M-H]1-","[M-2H]2-")
			}
	res <- query_db(acc.mass,add,ppm,MZedDB)
	if (length(res)>0){
		if(class(res)=="character"){
			res <- matrix(data=res,nrow=1)
		}
		res <- cbind(res[,8],res[,4],res[,3],res[,2],res[,7])
		mz <- rep(acc.mass,nrow(res))					
		ppm.err <- round((as.numeric(res[,1])-mz)/mz*10^6,4)
		res <- cbind(res,ppm.err)
		res <- res[order(res[,6]),]
		if(class(res)=="character"){
			res <- matrix(data=res,nrow=1)
		}
		colnames(res)<-c('m/z','Accurate Mass','Molecular Formula','Name','Adduct','PPM Error')
		patterns <- c("ic acid","keto","<i>","</i>","(R)-","(S)-","-L-","L-","(+)-",
                  "cis,","cis-","trans-","&beta;-","-D-",
                  "D-","D.","&gamma;-","-n-","N-","&alpha;-","&alpha;","&beta;")
		replacements <- c("ate","oxo","","","","","","","","","","","","","","","","","","","","")
		for(i in 1:length(patterns)){
			res[,4] <- sapply(res[,4],gsub,pattern=patterns[i],replacement=replacements[i],fixed=TRUE)
		}
		res[,4] <- tolower(res[,4])
		res <- res[!duplicated(res),]
	} else {
		res <- "No MZedDB hits"
		cat("\n","No MZedDB hits",sep=" ")
	}
	return(res)
}
