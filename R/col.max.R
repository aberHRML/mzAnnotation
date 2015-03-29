col.max <-
function(x.mass,name=NULL){
	if(is.null(dim(x.mass))){
		nam <- rep(name,length(x.mass))	
		x.max <- cbind(nam,x.mass)
	} else {
		x.max <- matrix(0,ncol=2,nrow=nrow(x.mass))
		for (i in 1:nrow(x.mass)){
			x <- x.mass[i,]
			ma <- max(x)
			if (ma==0){
				nam <- names(x)[1]	
			} else {
				nam <- names(x)[x %in% ma]
			}
			if (length(nam)>1){
				nam <- nam[1]
			} 
			a <- c(nam,ma)
			x.max[i,] <- a
		}
	}
	if(ncol(x.max)==1){
		colnames(x.max) <- c("intensity")
	} else {
		colnames(x.max) <- c("mz","intensity")
	}	
	return(data.frame(x.max))
}
