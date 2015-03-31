#' Filter generated MFs
#' @description Filter MFs that have more than one adduct or more than one isotope.

filterMF <-
function(MFs){
  mf <- as.character(MFs[,"Clean MF"])
	mf <- sapply(mf,function(x){y <- strsplit(x,"[.]");return(y)})
  mf <- lapply(mf,function(x){y <- sapply(x,gsub,pattern="iCl",replacement="ICL");return(y)})
  mf <- lapply(mf,function(x){y <- sapply(x,gsub,pattern="iO",replacement="Io");return(y)})
  mf <- lapply(mf,function(x){y <- sapply(x,gsub,pattern="iK",replacement="Ik");return(y)})
	mf <- unlist(lapply(mf,function(x){y <- !(TRUE %in% (grepl("iC",x) & grepl("Ik",x) || 
                                                grepl("iC",x) & grepl("Io",x) || 
                                                grepl("iC",x) & grepl("ICL",x) || 
                                                grepl("Io",x) & grepl("ICL",x) || 
                                                grepl("Io",x) & grepl("Ik",x) || 
                                                grepl("K",x) & grepl("Na",x) || 
                                                grepl("Na",x) & grepl("Ik",x) || 
                                                grepl("K",x) & grepl("Cl",x) || 
                                                grepl("ICL",x) & grepl("K",x) || 
                                                grepl("K",x) & grepl("Ik",x) || 
                                                grepl("Cl",x) & grepl("ICL",x)));return(y)}))
  filtered.mfs <- MFs[which(mf==T),]
	return (filtered.mfs)
}
