cor.lists <-
function(x) {
	col.1 = NULL
	for (i in 1:ncol(x)){
		col <- x[,i]
		col <- col[col>0]
		col <- sort(col,decreasing=T)
		col <- data.frame(names(col),col)
		col.1[i] <- list(col)
	}
	len <- NULL
	for (i in 1:length(col.1)){
		len[i] <- nrow(col.1[[i]])
	}
	mat <- matrix(0,ncol=length(col.1)*2,nrow=max(len))
	seq.1 <- seq(1,length(col.1)*2,2)
	nam <- NULL
	for (i in 1:length(col.1)){
		col.2 <- col.1[[i]]
		col.2[,1] <- as.character(col.2[,1])
		nam <- c(nam ,colnames(x)[i],"r")
		mat[1:nrow(col.2),seq.1[i]] <- col.2[,1]
		mat[1:nrow(col.2),seq.1[i]+1] <- col.2[,2]
	}
	colnames(mat) <- nam
	return(mat)
}
