acc.mass.final <-
function(acc_mat){		# Make final accurate mass matrix from acc_mat using the most intense mass
	final.mat.1 <- NULL
	for (x in 1:length(acc_mat)){
		mat <- acc_mat[[x]]
		final.mat <- matrix(0,nrow=length(rownames(mat)),ncol=2)
		for (i in 1:nrow(mat)){
			s <- seq(1,ncol(mat),2)
			row <- mat[i,]
			mass <- row[s]
			int <- row[s+1]
			m <- row[s[grep(max(int),int)]:(s[grep(max(int),int)]+1)]
			final.mat[i,] <- m
		}
	rownames(final.mat) <- rownames(mat)
	colnames(final.mat) <- c("mz","intensity")
	final.mat <- final.mat[is.finite(final.mat[,1]),]
	final.mat <- final.mat[order(final.mat[,1]),]
	final.mat.1[x] <- list(final.mat)
	}
	return(final.mat.1)
}
