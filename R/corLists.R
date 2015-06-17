#' Aggregate significant correlations into lists for each variable
#'
#' @param x 

corLists <-
function(x) {
	col.1 = NULL
	for (i in 1:ncol(x)){
		col <- x[,i]
		col <- col[col>0]
		col <- sort(col,decreasing=T)
		col <- data.frame(names(col),col)
		colnames(col) <- c(colnames(x)[i],"r")
		col.1[i] <- list(col)
	}
	names(col.1) <- colnames(x)
	return(col.1)
}
