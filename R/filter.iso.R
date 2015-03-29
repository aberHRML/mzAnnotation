filter.iso <-
function(isos){
	for (i in 1:length(isos)){
		curr.iso <- isos[[i]]
		iso.mat <- matrix(nrow=0,ncol=4)
		for(x in 1:nrow(curr.iso)){
			iso <- curr.iso[x,4]
			iso <- unlist(strsplit(iso," "))
			if (length(iso) < 3){
				iso.mat <- rbind(iso.mat,curr.iso[x,])
			}
		}
		isos[[i]] <- iso.mat
	}
 return(isos)
}
