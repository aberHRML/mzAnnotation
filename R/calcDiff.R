#' Identify structural relationships based on m/z differences
#'
#' @param mat 
#' @param mo 

calcDiff <-
function(mat,mo){
  data("adduct_rel")
  if(nrow(mat)>0){
    data('adduct_rel')
  	par <- colnames(mat)[1]
	  par <- as.numeric(gsub(mo,"",par))
	  if(mo=="n"){
      add <- adduct_rel[which(adduct_rel$Mode=="n"|adduct_rel$Mode=="b"),]
	  }
	  if(mo=="p"){
	  	add <- adduct_rel[which(adduct_rel$Mode=="p"|adduct_rel$Mode=="b"),]
	  }
	  cors <- mat[,1]
	  cors <- as.numeric(gsub(mo,"",cors))
	  rel <- round(sapply(add[,"Equation"],function(x,M){y <- eval(parse(text=x));return(y)},M=par),2)
	  cors.1 <- round(cors-par,2)
	  diff <- round(rel-par,2)
	  mat.1 <- cbind(mat,cors.1)
	  colnames(mat.1)[3] <- "Difference"
	  posi.add <- sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=diff,ad=add[,"MF"])
	  add.itself <- sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=-diff,ad=add[,"MF"])
	  if (class(posi.add)=="list"){
	    posi.add <- sapply(posi.add,paste,collapse=", ")
	  }
	  if (class(add.itself)=="list"){
	    add.itself <- sapply(add.itself,paste,collapse=", ")
	  }
	  mat.1 <- cbind(mat.1,posi.add,add.itself)
	  rownames(mat.1) <- NULL
  } else {
    mat.1 <- matrix(nrow=0,ncol=5)
    colnames(mat.1) <- c(colnames(mat),"Difference","posi.add","add.itself") 
  }
	return(data.frame(mat.1))
}
