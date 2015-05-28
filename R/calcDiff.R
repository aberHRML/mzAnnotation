#' Identify structural relationships based on m/z differences

calcDiff <-
function(mat,mo){
  data(adduct_rel)
	mat <- na.omit(mat)
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
  cors.1 <- round(cors-par,2)
	rel <- round(sapply(add[,"Equation"],function(x,M){y <- eval(parse(text=x));return(y)},M=par),2)
  names(rel) <- add[,"MF"]
	mat.1 <- cbind(mat,cors.1)
	colnames(mat.1)[3] <- "Difference"
	posi.add <-sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=diff,ad=add)
	add.itself <- sapply(mat.1[,3],function(x,dif,ad){matc <- x==dif;if (TRUE %in% matc) {add <- ad[matc]} else {add <- NA};return(add)},dif=-diff,ad=add)
	mat.1 <- cbind(mat.1,posi.add,add.itself)
	rownames(mat.1) <- NULL
	return(data.frame(mat.1))
}
