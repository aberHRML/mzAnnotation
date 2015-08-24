
getIsoDist <- function(mf_res,mode){
  if(nrow(mf_res)>0){
    if (mode=="p"){
	  	charge <- 1
	  }
	  if (mode=="n"){
	  	charge <- -1
	  }	
	  iso <- grepl("i",mf_res[,"Clean.MF"],fixed=TRUE)
	  if(length(which(iso==T))>0){
		  mf_res <- mf_res[-which(iso==T),]
	  }
	  if(nrow(mf_res)>0){
		  iso_ab <- lapply(as.character(mf_res[,"MF"]),isoDistr,chrg=charge)
		  iso_ab <- filterIso(iso_ab)
		  names(iso_ab) <- mf_res[,"MF"]
		  return(iso_ab)
	  }
  }
}