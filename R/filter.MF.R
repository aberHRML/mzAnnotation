filter.MF <-
function(MFs){
	filtered.mfs <- data.frame(matrix(nrow=0,ncol=6))
	for(i in 1:nrow(MFs)){
		mf <- as.character(MFs[i,2])
		x <- strsplit(mf,"[.]")
		# change iK and iCl to avoid confusion with K, Cl and iC when greping
		if ((grepl("iCl",x))){
			x <- list(c(sapply(x,gsub,pattern="iCl",replacement="ICL")))
		}
		if (grepl("iK",x)){
			x <- list(c(sapply(x,gsub,pattern="iK",replacement="Ik")))
		}
		
		if(!(grepl("iC",x) & grepl("K",x) || grepl("iC",x) & grepl("Na",x) || grepl("iC",x) & grepl("Ik",x) || grepl("iC",x) & grepl("Cl",x) || grepl("iC",x) & grepl("ICL",x) || grepl("K",x) & grepl("Na",x) || grepl("Na",x) & grepl("Ik",x) || grepl("K",x) & grepl("Cl",x) || grepl("ICL",x) & grepl("K",x) || grepl("K",x) & grepl("Ik",x) || grepl("Cl",x) & grepl("ICL",x))){
			filtered.mfs <- rbind(filtered.mfs,MFs[i,])
		}
	}
	names(filtered.mfs) <- names(MFs)
	return (filtered.mfs)
}
